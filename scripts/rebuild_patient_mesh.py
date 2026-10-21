
# -------------------------------------------------------------
# PATIENT-SPECIFIC ANEURYSM MESH RECONSTRUCTION
# -------------------------------------------------------------
#
# Purpose:
# Rebuild Gmsh geometric topology from the patient-specific
# vessel wall and three separately exported boundary caps.
#
# The eventual objective is a tetrahedral fluid-domain mesh
# with four distinct boundary patches for OpenFOAM.
#
# At this stage, we only load and inspect the surfaces.

from pathlib import Path
import gmsh
import numpy as np
from collections import Counter

# Locate the repository independently of the working directory.
project_dir = Path(__file__).resolve().parent.parent

surface_files = [
    "patient_vessel_wall.stl",
    "patient_cap_0.stl",
    "patient_cap_1.stl",
    "patient_cap_2.stl",
]

gmsh.initialize()
gmsh.model.add("patient_aneurysm_rebuilt")

# Import the vessel wall and each cap as separate surfaces.
for filename in surface_files:
    gmsh.merge(str(project_dir / "geometry" / filename))

# Merge coincident vertices along the cap perimeters.
gmsh.model.mesh.removeDuplicateNodes()

# Inspect the imported surfaces before reconstruction.
surfaces = gmsh.model.getEntities(2)

print("Imported surfaces:", surfaces)

for _, surface_tag in surfaces:
    element_types, element_tags, _ = (
        gmsh.model.mesh.getElements(2, surface_tag)
    )

    triangle_count = sum(
        len(tags) for tags in element_tags
    )

    print(
        f"Surface {surface_tag}: {triangle_count} triangles"
    )

# -------------------------------------------------------------
# EXTRACT CAP PERIMETER EDGES
# -------------------------------------------------------------
#
# Each cap consists of triangular elements.
#
# Interior edges are shared by two triangles.
# Perimeter edges belong to only one triangle.
#
# Store the perimeter edges for later reconstruction
# of geometric boundary curves.

cap_perimeters = {}

for cap_tag in [2, 3, 4]:

    # Retrieve the triangular elements belonging to this cap.
    element_types, element_tags, node_groups = (
        gmsh.model.mesh.getElements(2, cap_tag)
    )

    triangles = np.asarray(
        node_groups[0]
    ).reshape(-1, 3)

    edge_counts = Counter()

    # Count how often each undirected edge occurs.
    for a, b, c in triangles:

        for u, v in [(a, b), (b, c), (c, a)]:

            edge = tuple(sorted((int(u), int(v))))
            edge_counts[edge] += 1

    # Edges occurring once form the cap perimeter.
    perimeter_edges = [
        edge
        for edge, count in edge_counts.items()
        if count == 1
    ]

    cap_perimeters[cap_tag] = perimeter_edges

    print(
        f"Cap {cap_tag}: "
        f"{len(perimeter_edges)} perimeter edges"
    )

# -------------------------------------------------------------
# ORDER CAP PERIMETER VERTICES
# -------------------------------------------------------------
#
# Convert unordered perimeter edges into continuous loops.
# Each vertex should have exactly two neighbours.

cap_loops = {}

for cap_tag, perimeter_edges in cap_perimeters.items():

    # Build a lookup table of neighbouring vertices.
    adjacency = {}

    for u, v in perimeter_edges:
        adjacency.setdefault(u, []).append(v)
        adjacency.setdefault(v, []).append(u)

    # Verify that the perimeter has no branches or endpoints.
    if any(len(neighbours) != 2 for neighbours in adjacency.values()):
        raise ValueError(
            f"Cap {cap_tag} has an invalid perimeter."
        )

    # Start at an arbitrary vertex and follow the loop.
    start = min(adjacency)
    loop = [start]

    previous = None
    current = start

    while True:

        neighbours = adjacency[current]

        # Continue to the neighbour we did not just visit.
        next_vertex = (
            neighbours[0]
            if neighbours[0] != previous
            else neighbours[1]
        )

        if next_vertex == start:
            break

        if next_vertex in loop:
            raise ValueError(
                f"Cap {cap_tag} contains a repeated vertex."
            )

        loop.append(next_vertex)

        previous = current
        current = next_vertex

    # Confirm that the traversal included every perimeter vertex.
    if len(loop) != len(adjacency):
        raise ValueError(
            f"Cap {cap_tag} contains disconnected loops."
        )

    cap_loops[cap_tag] = loop

    print(
        f"Cap {cap_tag}: ordered loop with "
        f"{len(loop)} vertices"
    )

# -------------------------------------------------------------
# CHECK CAP ORIENTATION AGAINST VESSEL WALL
# -------------------------------------------------------------
#
# For a consistently oriented surface, adjacent triangles
# traverse their shared edge in opposite directions.
#
# Compare the orientation of each cap perimeter edge with
# the corresponding vessel-wall edge.

def get_directed_edges(surface_tag):

    _, _, node_groups = gmsh.model.mesh.getElements(
        2, surface_tag
    )

    triangles = np.asarray(
        node_groups[0]
    ).reshape(-1, 3)

    directed_edges = {}

    for a, b, c in triangles:

        for u, v in [(a, b), (b, c), (c, a)]:

            u, v = int(u), int(v)

            edge = tuple(sorted((u, v)))

            directed_edges.setdefault(edge, []).append(
                (u, v)
            )

    return directed_edges


wall_edges = get_directed_edges(1)

for cap_tag, perimeter_edges in cap_perimeters.items():

    cap_edges = get_directed_edges(cap_tag)

    same_direction = 0
    opposite_direction = 0

    for edge in perimeter_edges:

        wall_direction = wall_edges[edge][0]
        cap_direction = cap_edges[edge][0]

        if wall_direction == cap_direction:
            same_direction += 1
        else:
            opposite_direction += 1

    print(
        f"Cap {cap_tag}: "
        f"{same_direction} same-direction edges, "
        f"{opposite_direction} opposite-direction edges"
    )

# -------------------------------------------------------------
# CORRECT CAP TRIANGLE ORIENTATION
# -------------------------------------------------------------
#
# All three caps have triangle winding reversed relative
# to the vessel wall.
#
# Reverse the orientation of their existing mesh elements.
# The original STL files remain unchanged.

for cap_tag in [2, 3, 4]:

    gmsh.model.mesh.reverse([
        (2, cap_tag)
    ])

    print(
        f"Reversed triangle orientation for cap {cap_tag}"
    )

# -------------------------------------------------------------
# VERIFY CORRECTED CAP ORIENTATION
# -------------------------------------------------------------

wall_edges = get_directed_edges(1)

for cap_tag, perimeter_edges in cap_perimeters.items():

    cap_edges = get_directed_edges(cap_tag)

    opposite_direction = sum(
        wall_edges[edge][0] != cap_edges[edge][0]
        for edge in perimeter_edges
    )

    print(
        f"Cap {cap_tag}: "
        f"{opposite_direction}/{len(perimeter_edges)} "
        f"oppositely directed boundary edges"
    )

# -------------------------------------------------------------
# CHECK SIGNED VOLUME OF THE CLOSED SURFACE
# -------------------------------------------------------------
#
# Calculate the signed enclosed volume using oriented
# surface triangles.
#
# Positive volume indicates outward orientation, assuming
# the mesh is consistently oriented and non-self-intersecting.
#
# The STL coordinates are in millimetres, so the resulting
# volume is in cubic millimetres.

node_tags, node_coordinates, _ = gmsh.model.mesh.getNodes()

coordinates = np.asarray(node_coordinates).reshape(-1, 3)

node_lookup = {
    int(tag): coordinates[i]
    for i, tag in enumerate(node_tags)
}

signed_volume = 0.0

for surface_tag in [1, 2, 3, 4]:

    _, _, node_groups = gmsh.model.mesh.getElements(
        2, surface_tag
    )

    triangles = np.asarray(
        node_groups[0]
    ).reshape(-1, 3)

    for a, b, c in triangles:

        p0 = node_lookup[int(a)]
        p1 = node_lookup[int(b)]
        p2 = node_lookup[int(c)]

        signed_volume += (
            np.dot(p0, np.cross(p1, p2)) / 6.0
        )

print(f"Signed enclosed volume: {signed_volume:.6f} mm³")

# -------------------------------------------------------------
# INSPECT CAP BOUNDARY COORDINATES
# -------------------------------------------------------------
#
# Calculate each cap perimeter's centroid and bounding-box
# dimensions in millimetres.
#
# This helps verify that the ordered loops represent
# plausible vascular openings before rebuilding topology.

for cap_tag, loop in cap_loops.items():

    points = np.array([
        node_lookup[int(node_tag)]
        for node_tag in loop
    ])

    centroid = points.mean(axis=0)

    dimensions = points.max(axis=0) - points.min(axis=0)

    print(f"\nCap {cap_tag}:")

    print(
        "  Centroid (mm):",
        np.round(centroid, 3)
    )

    print(
        "  Bounding-box dimensions (mm):",
        np.round(dimensions, 3)
    )

# -------------------------------------------------------------
# CALCULATE CAP CROSS-SECTIONAL AREAS
# -------------------------------------------------------------
#
# Calculate the area of each triangular cap by summing
# the areas of its constituent triangles.
#
# Triangle area = 0.5 * |(p1 - p0) x (p2 - p0)|
#
# Coordinates are in millimetres, so areas are in mm².

cap_areas = {}

for cap_tag in [2, 3, 4]:

    _, _, node_groups = gmsh.model.mesh.getElements(
        2, cap_tag
    )

    triangles = np.asarray(
        node_groups[0]
    ).reshape(-1, 3)

    total_area = 0.0

    for a, b, c in triangles:

        p0 = node_lookup[int(a)]
        p1 = node_lookup[int(b)]
        p2 = node_lookup[int(c)]

        triangle_area = 0.5 * np.linalg.norm(
            np.cross(p1 - p0, p2 - p0)
        )

        total_area += triangle_area

    cap_areas[cap_tag] = total_area

    print(
        f"Cap {cap_tag}: "
        f"cross-sectional area = {total_area:.4f} mm²"
    )

# -------------------------------------------------------------
# PRESERVE CORRECTED SURFACE MESH DATA
# -------------------------------------------------------------
#
# Store the existing node coordinates and triangle connectivity
# before rebuilding Gmsh's geometric entities.
#
# The cap triangles have already had their orientation corrected.

saved_nodes = {
    int(tag): coordinates[i].copy()
    for i, tag in enumerate(node_tags)
}

saved_triangles = {}

for surface_tag in [1, 2, 3, 4]:

    element_types, element_tags, node_groups = (
        gmsh.model.mesh.getElements(2, surface_tag)
    )

    # Type 2 represents a three-node triangular element.
    triangle_index = list(element_types).index(2)

    saved_triangles[surface_tag] = np.asarray(
        node_groups[triangle_index],
        dtype=np.int64
    ).reshape(-1, 3).copy()

print("\nSaved mesh data for reconstruction:")
print("Saved nodes:", len(saved_nodes))

for surface_tag, triangles in saved_triangles.items():
    print(
        f"Surface {surface_tag}: "
        f"{len(triangles)} saved triangles"
    )

gmsh.finalize()
