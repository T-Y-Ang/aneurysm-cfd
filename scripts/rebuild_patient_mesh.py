
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

# -------------------------------------------------------------
# BEGIN DISCRETE TOPOLOGY RECONSTRUCTION
# -------------------------------------------------------------
#
# We have already saved the original nodes, corrected triangles, # and three ordered cap boundary loops.
#
# Clear the current Gmsh model so that we can rebuild its # geometric entities with explicit shared boundaries.
#
# The saved Python dictionaries remain available in memory.

gmsh.clear()

gmsh.model.add("patient_aneurysm_topology")

print("\nNew Gmsh model created for topology reconstruction.")
print("Existing entities:", gmsh.model.getEntities())
print("Saved boundary loops:", list(cap_loops.keys()))

# -------------------------------------------------------------
# CREATE DISCRETE ENDPOINTS FOR THE CAP BOUNDARY CURVES
# -------------------------------------------------------------
#
# Split each closed cap perimeter at two approximately opposite # positions. These nodes become the endpoints of two curves.
#
# We retain the original node identifiers and coordinates.

cap_endpoints = {}
point_tags = {}
next_point_tag = 1

for cap_tag, loop in cap_loops.items():

    first_node = int(loop[0])
    opposite_node = int(loop[len(loop) // 2])

    cap_endpoints[cap_tag] = (
        first_node,
        opposite_node
    )

    for node_tag in [first_node, opposite_node]:

        point_tag = next_point_tag
        next_point_tag += 1

        gmsh.model.addDiscreteEntity(0, point_tag)

        gmsh.model.mesh.addNodes(
            0,
            point_tag,
            [node_tag],
            saved_nodes[node_tag].tolist()
        )

        point_tags[node_tag] = point_tag

    print(
        f"Cap {cap_tag}: endpoints "
        f"{first_node}, {opposite_node}"
    )

print(
    "Created geometric points:",
    gmsh.model.getEntities(0)
)

# -------------------------------------------------------------
# CREATE DISCRETE CURVES ALONG THE CAP PERIMETERS
# -------------------------------------------------------------
#
# Each closed perimeter is divided into two arcs.
# The arcs share the two geometric endpoints created above.
#
# Interior arc nodes are classified on their respective curves.
# Existing node coordinates and identifiers are preserved.

cap_curves = {}
next_curve_tag = 1

for cap_tag, loop in cap_loops.items():

    midpoint = len(loop) // 2

    # Two arcs connecting the same endpoints.
    arcs = [
        loop[:midpoint + 1],
        loop[midpoint:] + [loop[0]]
    ]

    cap_curves[cap_tag] = []

    for arc in arcs:

        curve_tag = next_curve_tag
        next_curve_tag += 1

        # Create a curve bounded by its two endpoint points.
        start_point = point_tags[int(arc[0])]
        end_point = point_tags[int(arc[-1])]

        gmsh.model.addDiscreteEntity(
            1,
            curve_tag,
            [start_point, end_point]
        )

        # Add interior nodes only.
        # Endpoint nodes already belong to geometric points.
        interior_nodes = [int(n) for n in arc[1:-1]]

        interior_coordinates = [
            coordinate
            for node in interior_nodes
            for coordinate in saved_nodes[node]
        ]

        if interior_nodes:
            gmsh.model.mesh.addNodes(
                1,
                curve_tag,
                interior_nodes,
                interior_coordinates
            )

        # Create a line element for every perimeter edge.
        line_nodes = []

        for a, b in zip(arc[:-1], arc[1:]):
            line_nodes.extend([int(a), int(b)])

        gmsh.model.mesh.addElementsByType(
            curve_tag,
            1,       # Two-node line elements
            [],      # Automatically assign element tags
            line_nodes
        )

        cap_curves[cap_tag].append(curve_tag)

    print(
        f"Cap {cap_tag}: curves {cap_curves[cap_tag]}"
    )

print(
    "Created geometric curves:",
    gmsh.model.getEntities(1)
)

# -------------------------------------------------------------
# RECONSTRUCT DISCRETE SURFACES
# -------------------------------------------------------------
#
# Create four surfaces with explicit shared curve boundaries:
#   1 = vessel wall
#   2 = cap_0
#   3 = cap_1
#   4 = cap_2
#
# The vessel wall shares all six boundary curves with the caps.

surface_boundaries = {
    1: [1, 2, 3, 4, 5, 6],
    2: cap_curves[2],
    3: cap_curves[3],
    4: cap_curves[4],
}

for surface_tag, boundaries in surface_boundaries.items():

    gmsh.model.addDiscreteEntity(
        2,
        surface_tag,
        boundaries
    )

    print(
        f"Created surface {surface_tag} "
        f"with boundary curves {boundaries}"
    )

print(
    "Reconstructed surfaces:",
    gmsh.model.getEntities(2)
)

# -------------------------------------------------------------
# RESTORE ORIGINAL SURFACE NODES AND TRIANGLES
# -------------------------------------------------------------
#
# Nodes on points and curves have already been registered.
# Only the remaining interior surface nodes need to be added.
#
# Preserve the original triangle connectivity and orientation.

assigned_nodes = set(point_tags.keys())

for loop in cap_loops.values():
    assigned_nodes.update(int(node) for node in loop)

for surface_tag, triangles in saved_triangles.items():

    # Find all node IDs referenced by this surface.
    surface_nodes = set(
        int(node) for node in triangles.flatten()
    )

    # Add only nodes not already assigned to another entity.
    interior_nodes = sorted(surface_nodes - assigned_nodes)

    if interior_nodes:

        interior_coordinates = [
            coordinate
            for node in interior_nodes
            for coordinate in saved_nodes[node]
        ]

        gmsh.model.mesh.addNodes(
            2,
            surface_tag,
            interior_nodes,
            interior_coordinates
        )

    assigned_nodes.update(interior_nodes)

    # Restore the original three-node triangles.
    gmsh.model.mesh.addElementsByType(
        surface_tag,
        2,
        [],
        triangles.flatten().tolist()
    )

    print(
        f"Surface {surface_tag}: "
        f"{len(interior_nodes)} interior nodes, "
        f"{len(triangles)} triangles restored"
    )

print(
    "Total uniquely assigned nodes:",
    len(assigned_nodes)
)

# -------------------------------------------------------------
# CREATE THE ENCLOSED FLUID VOLUME
# -------------------------------------------------------------
#
# The four reconstructed surfaces form the boundary of the # patient-specific aneurysm fluid domain.
#
# Surface 1 = vessel wall
# Surface 2 = cap_0
# Surface 3 = cap_1
# Surface 4 = cap_2

###fluid_volume_tag = gmsh.model.addDiscreteEntity(
###    3,
###    1,
###    [1, 2, 3, 4]
###)

# -------------------------------------------------------------
# CREATE A GEOMETRIC FLUID VOLUME
# -------------------------------------------------------------
#
# Test whether Gmsh's built-in geometry kernel can construct
# a meshable volume from the reconstructed surface entities.

surface_loop_tag = gmsh.model.geo.addSurfaceLoop(
    [1, 2, 3, 4]
)

fluid_volume_tag = gmsh.model.geo.addVolume(
    [surface_loop_tag]
)

gmsh.model.geo.synchronize()

print(
    "Geometric volume entities:",
    gmsh.model.getEntities(3)
)

print("\nCreated fluid volume:", fluid_volume_tag)

print(
    "Volume entities:",
    gmsh.model.getEntities(3)
)

print(
    "Volume boundary surfaces:",
    gmsh.model.getBoundary(
        [(3, fluid_volume_tag)],
        oriented=False
    )
)

# -------------------------------------------------------------
# GENERATE PATIENT-SPECIFIC TETRAHEDRAL MESH
# -------------------------------------------------------------
#
# The original surface triangles have already been restored.
# Gmsh will now attempt to fill the enclosed fluid volume with tetrahedral elements.

print("\nStarting tetrahedral mesh generation...")

gmsh.model.mesh.generate(3)

# -------------------------------------------------------------
# OPTIMISE THE 3D TETRAHEDRAL MESH
# -------------------------------------------------------------
#
# Apply Gmsh's Netgen optimisation algorithm to improve tetrahedral element quality after volume mesh generation.
#
# This aims to reduce highly distorted or nearly degenerate tetrahedra while retaining the existing surface geometry.

gmsh.model.mesh.optimize("Netgen")

# Element type 4 represents a four-node tetrahedron.
tetrahedron_tags, tetrahedron_nodes = (
    gmsh.model.mesh.getElementsByType(4)
)

tetrahedron_count = len(tetrahedron_tags)

print("\nTetrahedral meshing result:")
print("Number of tetrahedra:", tetrahedron_count)

if tetrahedron_count == 0:
    raise RuntimeError(
        "No tetrahedra generated. "
        "The reconstructed topology requires further investigation."
    )

print("SUCCESS: A non-empty tetrahedral mesh was generated.")

# -------------------------------------------------------------
# DEFINE PHYSICAL GROUPS FOR CFD BOUNDARIES
# -------------------------------------------------------------
#
# Physical groups identify the vessel wall, the three
# arterial openings, and the internal fluid volume.
#
# We have not yet assigned inlet/outlet roles to the caps.

physical_surfaces = {
    1: "vessel_wall",
    2: "cap_0",
    3: "cap_1",
    4: "cap_2",
}

for surface_tag, name in physical_surfaces.items():

    group_tag = gmsh.model.addPhysicalGroup(
        2,
        [surface_tag]
    )

    gmsh.model.setPhysicalName(
        2,
        group_tag,
        name
    )

    print(f"Created surface group: {name}")


# Identify the enclosed fluid region.
fluid_group_tag = gmsh.model.addPhysicalGroup(
    3,
    [fluid_volume_tag]
)

gmsh.model.setPhysicalName(
    3,
    fluid_group_tag,
    "fluid"
)

print("Created volume group: fluid")

# -------------------------------------------------------------
# SAVE THE INITIAL PATIENT-SPECIFIC VOLUME MESH
# -------------------------------------------------------------
#
# This is our first successful tetrahedral mesh.
# It is saved as an intermediate checkpoint because mesh quality and CFD boundary conditions are not yet validated.

# Create the output directory if it does not already exist.
mesh_directory = Path("mesh")
mesh_directory.mkdir(parents=True, exist_ok=True)

# Define the output filename.
# Save the optimised mesh separately so the original remains available for comparison and reproducibility.
mesh_path = mesh_directory / "patient_aneurysm_optimized.msh"

# Use Gmsh MSH version 2.2 for OpenFOAM compatibility.
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# Save the current mesh.
gmsh.write(str(mesh_path))

print("\nSaved initial tetrahedral mesh:")
print(mesh_path.resolve())
print("File size:", mesh_path.stat().st_size, "bytes")

gmsh.finalize()
