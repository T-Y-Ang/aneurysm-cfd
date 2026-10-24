
from pathlib import Path
from collections import Counter, defaultdict

import gmsh
import numpy as np

# -------------------------------------------------------------
# BUILD PATIENT-SPECIFIC ANEURYSM VOLUME MESH
# -------------------------------------------------------------
#
# Starting geometry:
#   mesh/patient_wall_classified_remesh_test.msh
#
# Outputs:
#   mesh/patient_aneurysm_volume.msh
#
# Units: millimetres (conversion to metres will occur
# when importing into OpenFOAM).
#
# Physical groups:
#   vessel_wall
#   opening_0
#   opening_1
#   opening_2
#   fluid
# -------------------------------------------------------------

INPUT = Path("mesh/patient_wall_classified_remesh_test.msh")
OUTPUT = Path("mesh/patient_aneurysm_volume.msh")

# Gmsh curve entities belonging to each opening.
OPENING_CURVES = {
    "opening_0": [36, 66],
    "opening_1": [67],
    "opening_2": [51, 52],
}


def get_ordered_boundary_loops(triangles):
    """Find ordered boundary loops from surface triangles."""

    edge_counts = Counter()

    for a, b, c in triangles:
        for u, v in ((a, b), (b, c), (c, a)):
            edge_counts[tuple(sorted((int(u), int(v))))] += 1

    adjacency = defaultdict(set)

    for (a, b), count in edge_counts.items():
        if count == 1:
            adjacency[a].add(b)
            adjacency[b].add(a)

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
                raise RuntimeError("Invalid boundary loop.")

            ordered.append(next_vertex)
            previous, current = current, next_vertex

        unvisited.difference_update(ordered)
        loops.append(ordered)

    return sorted(loops, key=len, reverse=True)


def cross2d(a, b):
    return a[0] * b[1] - a[1] * b[0]


def polygon_area(points):
    x = points[:, 0]
    y = points[:, 1]

    return 0.5 * np.sum(
        x * np.roll(y, -1) - np.roll(x, -1) * y
    )


def point_in_triangle(p, a, b, c, tolerance=1e-12):
    return (
        cross2d(b - a, p - a) >= -tolerance
        and cross2d(c - b, p - b) >= -tolerance
        and cross2d(a - c, p - c) >= -tolerance
    )


def triangulate_polygon(points_2d):
    """Triangulate a planar polygon without adding vertices."""

    remaining = list(range(len(points_2d)))

    if polygon_area(points_2d) < 0:
        remaining.reverse()

    faces = []

    while len(remaining) > 3:
        found = False

        for j in range(len(remaining)):
            i0 = remaining[j - 1]
            i1 = remaining[j]
            i2 = remaining[(j + 1) % len(remaining)]

            a, b, c = points_2d[[i0, i1, i2]]

            if cross2d(b - a, c - a) <= 1e-12:
                continue

            others = [
                k for k in remaining
                if k not in (i0, i1, i2)
            ]

            if any(
                point_in_triangle(points_2d[k], a, b, c)
                for k in others
            ):
                continue

            faces.append([i0, i1, i2])
            remaining.pop(j)
            found = True
            break

        if not found:
            raise RuntimeError("Cap triangulation failed.")

    faces.append(remaining)

    return np.asarray(faces, dtype=int)


def main():

    gmsh.initialize()

    try:
        gmsh.open(str(INPUT))

        # -----------------------------------------------------
        # 1. READ EXISTING WALL TRIANGLES
        # -----------------------------------------------------

        wall_surfaces = [
            tag for dim, tag in gmsh.model.getEntities(2)
        ]

        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()

        xyz = {
            int(tag): point
            for tag, point in zip(
                node_tags,
                np.asarray(coordinates).reshape(-1, 3)
            )
        }

        _, triangle_nodes = gmsh.model.mesh.getElementsByType(2)
        wall_triangles = triangle_nodes.reshape(-1, 3)

        loops = get_ordered_boundary_loops(wall_triangles)

        if len(loops) != 3:
            raise RuntimeError("Expected three vessel openings.")

        # -----------------------------------------------------
        # 2. ADD THREE CAP SURFACES
        # -----------------------------------------------------

        cap_surfaces = {}

        for index, loop in enumerate(loops):

            points = np.array([xyz[tag] for tag in loop])
            centre = points.mean(axis=0)

            _, _, vh = np.linalg.svd(
                points - centre,
                full_matrices=False
            )

            points_2d = (points - centre) @ vh[:2].T

            local_faces = triangulate_polygon(points_2d)

            global_faces = np.asarray(loop)[local_faces]

            # -----------------------------------------------------
            # CREATE CAP SURFACE WITH EXPLICIT BOUNDARY CURVES
            # -----------------------------------------------------
            #
            # The cap must reference the existing Gmsh curves forming its perimeter.
            #
            # This establishes topological connectivity between the cap and the classified vessel wall.
            # -----------------------------------------------------
            opening_names = [
                 "opening_0",
                 "opening_1",
                 "opening_2"
            ]

            opening_name = opening_names[index]

            boundary_curves = OPENING_CURVES[opening_name]

            surface_tag = gmsh.model.addDiscreteEntity(
                2,
                -1,
                boundary_curves
            )

            # Add cap triangles using the existing node tags.
            gmsh.model.mesh.addElementsByType(
                surface_tag,
                2,
                [],
                global_faces.flatten().tolist()
            )

            cap_surfaces[index] = surface_tag

            print(
                f"Cap {index}: "
                f"{len(global_faces)} triangles, "
                f"surface entity {surface_tag}"
            )

        # -----------------------------------------------------
        # 3. DEFINE THE ENCLOSED FLUID VOLUME
        # -----------------------------------------------------

        all_surfaces = wall_surfaces + list(
            cap_surfaces.values()
        )

        volume_tag = gmsh.model.addDiscreteEntity(
            3,
            -1,
            all_surfaces
        )

        # -----------------------------------------------------
        # 4. DEFINE PHYSICAL GROUPS
        # -----------------------------------------------------

        group = gmsh.model.addPhysicalGroup(
            2, wall_surfaces
        )
        gmsh.model.setPhysicalName(
            2, group, "vessel_wall"
        )

        for index, surface_tag in cap_surfaces.items():
            group = gmsh.model.addPhysicalGroup(
                2, [surface_tag]
            )

            gmsh.model.setPhysicalName(
                2, group, f"opening_{index}"
            )

        group = gmsh.model.addPhysicalGroup(
            3, [volume_tag]
        )
        gmsh.model.setPhysicalName(
            3, group, "fluid"
        )

        # -----------------------------------------------------
        # 5. GENERATE TETRAHEDRA
        # -----------------------------------------------------

        gmsh.option.setNumber(
            "Mesh.MeshSizeMin", 0.10
        )
        gmsh.option.setNumber(
            "Mesh.MeshSizeMax", 0.15
        )

        gmsh.option.setNumber(
            "Mesh.MeshSizeFromPoints", 0
        )
        gmsh.option.setNumber(
            "Mesh.MeshSizeFromCurvature", 0
        )
        gmsh.option.setNumber(
            "Mesh.MeshSizeExtendFromBoundary", 0
        )

        print("\nGenerating tetrahedral volume mesh...")

        gmsh.model.mesh.generate(3)

        gmsh.write(str(OUTPUT))

        _, tetra_nodes = gmsh.model.mesh.getElementsByType(4)

        print("\nVOLUME MESH RESULT")
        print("Tetrahedra:", len(tetra_nodes) // 4)
        print("Saved:", OUTPUT)

    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
