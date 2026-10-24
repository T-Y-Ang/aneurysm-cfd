
from pathlib import Path
from collections import Counter, defaultdict

import gmsh
import numpy as np
import trimesh

# -------------------------------------------------------------
# BUILD A CLOSED PATIENT-SPECIFIC VESSEL SURFACE
# -------------------------------------------------------------
#
# Input:
#   mesh/patient_wall_classified_remesh_test.msh
#
# Output (only if validation passes):
#   geometry/patient_aneurysm_closed_remeshed.stl
#
# Units: millimetres
# -------------------------------------------------------------

INPUT = Path("mesh/patient_wall_classified_remesh_test.msh")
OUTPUT = Path("geometry/patient_aneurysm_closed_remeshed.stl")


def polygon_area(points):
    """Signed area of a polygon in two dimensions."""
    x = points[:, 0]
    y = points[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def cross2d(a, b):
    """Two-dimensional cross product."""
    return a[0] * b[1] - a[1] * b[0]


def point_in_triangle(p, a, b, c, tolerance=1e-12):
    """Check whether a point lies inside or on a CCW triangle."""
    return (
        cross2d(b - a, p - a) >= -tolerance
        and cross2d(c - b, p - b) >= -tolerance
        and cross2d(a - c, p - c) >= -tolerance
    )


def triangulate_polygon(points_2d):
    """
    Triangulate a simple planar polygon by ear clipping.

    Returns triangles indexed into the original ordered
    boundary vertices. No new vertices are introduced.
    """
    count = len(points_2d)
    remaining = list(range(count))

    # Ear clipping expects counterclockwise ordering.
    if polygon_area(points_2d) < 0:
        remaining.reverse()

    faces = []
    tolerance = 1e-12

    while len(remaining) > 3:
        ear_found = False

        for j in range(len(remaining)):
            i0 = remaining[j - 1]
            i1 = remaining[j]
            i2 = remaining[(j + 1) % len(remaining)]

            a = points_2d[i0]
            b = points_2d[i1]
            c = points_2d[i2]

            # An ear must be convex and nondegenerate.
            if cross2d(b - a, c - a) <= tolerance:
                continue

            other_vertices = [
                k for k in remaining
                if k not in (i0, i1, i2)
            ]

            if any(
                point_in_triangle(points_2d[k], a, b, c)
                for k in other_vertices
            ):
                continue

            faces.append([i0, i1, i2])
            remaining.pop(j)
            ear_found = True
            break

        if not ear_found:
            raise RuntimeError(
                "Ear clipping failed: check polygon geometry."
            )

    faces.append(remaining)
    return np.asarray(faces, dtype=int)


def get_boundary_loops(triangles):
    """Extract ordered loops from open triangle edges."""
    edge_counts = Counter()

    for a, b, c in triangles:
        for u, v in ((a, b), (b, c), (c, a)):
            edge_counts[tuple(sorted((int(u), int(v))))] += 1

    adjacency = defaultdict(set)

    for (a, b), count in edge_counts.items():
        if count == 1:
            adjacency[a].add(b)
            adjacency[b].add(a)

    if any(len(neighbours) != 2 for neighbours in adjacency.values()):
        raise RuntimeError("Non-manifold opening boundary detected.")

    loops = []
    unvisited = set(adjacency)

    while unvisited:
        start = min(unvisited)
        ordered = [start]
        previous = None
        current = start

        while True:
            candidates = sorted(
                v for v in adjacency[current]
                if v != previous
            )

            next_vertex = candidates[0]

            if next_vertex == start:
                break

            if next_vertex in ordered:
                raise RuntimeError("Repeated vertex in boundary loop.")

            ordered.append(next_vertex)
            previous, current = current, next_vertex

        unvisited.difference_update(ordered)
        loops.append(ordered)

    return sorted(loops, key=len, reverse=True)


def main():
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    gmsh.initialize()

    try:
        gmsh.open(str(INPUT))

        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        node_tags = np.asarray(node_tags, dtype=np.int64)
        vertices = np.asarray(coordinates).reshape(-1, 3)

        tag_to_index = {
            int(tag): i for i, tag in enumerate(node_tags)
        }

        _, triangle_nodes = gmsh.model.mesh.getElementsByType(2)
        wall_tags = triangle_nodes.reshape(-1, 3)

    finally:
        gmsh.finalize()

    wall_faces = np.array([
        [tag_to_index[int(tag)] for tag in face]
        for face in wall_tags
    ], dtype=int)

    loops = get_boundary_loops(wall_faces)

    if len(loops) != 3:
        raise RuntimeError(
            f"Expected 3 openings; found {len(loops)}."
        )

    print("\nRECONSTRUCTING THREE VESSEL CAPS")

    cap_faces = []

    for number, loop in enumerate(loops):
        points = vertices[loop]

        # Find the best-fit opening plane using SVD.
        centre = points.mean(axis=0)
        _, _, vh = np.linalg.svd(
            points - centre, full_matrices=False
        )

        plane_coordinates = (points - centre) @ vh[:2].T

        deviation = np.max(
            np.abs((points - centre) @ vh[2])
        )

        if deviation > 1e-4:
            raise RuntimeError(
                f"Opening {number} is insufficiently planar."
            )

        local_faces = triangulate_polygon(plane_coordinates)

        # Check that all original opening edges are retained.
        cap_edge_counts = Counter()

        for a, b, c in local_faces:
            for u, v in ((a, b), (b, c), (c, a)):
                cap_edge_counts[tuple(sorted((int(u), int(v))))] += 1

        actual_boundary = {
            edge for edge, count in cap_edge_counts.items()
            if count == 1
        }

        expected_boundary = {
            tuple(sorted((i, (i + 1) % len(loop))))
            for i in range(len(loop))
        }

        if actual_boundary != expected_boundary:
            raise RuntimeError(
                f"Opening {number}: cap boundary mismatch."
            )

        global_faces = np.asarray(loop)[local_faces]
        cap_faces.append(global_faces)

        triangle_points = vertices[global_faces]

        areas = 0.5 * np.linalg.norm(
            np.cross(
                triangle_points[:, 1] - triangle_points[:, 0],
                triangle_points[:, 2] - triangle_points[:, 0]
            ),
            axis=1
        )

        if np.any(areas <= 1e-12):
            raise RuntimeError(
                f"Opening {number}: degenerate triangles."
            )

        print(
            f"Opening {number}: "
            f"{len(loop)} vertices, "
            f"{len(global_faces)} triangles, "
            f"area={areas.sum():.6f} mm²"
        )

    combined_faces = np.vstack([wall_faces, *cap_faces])

    mesh = trimesh.Trimesh(
        vertices=vertices,
        faces=combined_faces,
        process=False
    )

    print("\nCLOSED SURFACE VALIDATION")
    print("Watertight:", mesh.is_watertight)
    print("Consistent winding:", mesh.is_winding_consistent)
    print("Enclosed volume (mm³):", mesh.volume)
    print("Surface area (mm²):", mesh.area)

    # Trimesh can correct inconsistent face orientation.
    if not mesh.is_winding_consistent:
        trimesh.repair.fix_winding(mesh)

    # Ensure outward-facing orientation for positive volume.
    if mesh.is_watertight and mesh.volume < 0:
        mesh.invert()

    if not mesh.is_watertight:
        raise RuntimeError(
            "Closed surface is not watertight. STL not exported."
        )

    if not mesh.is_winding_consistent:
        raise RuntimeError(
            "Surface winding is inconsistent. STL not exported."
        )

    if mesh.volume <= 0:
        raise RuntimeError(
            "Invalid enclosed volume. STL not exported."
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(str(OUTPUT))

    print("\nSUCCESS")
    print("Saved:", OUTPUT)


if __name__ == "__main__":
    main()
