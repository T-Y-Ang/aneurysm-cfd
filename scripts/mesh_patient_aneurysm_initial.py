# =============================================================
# PATIENT-SPECIFIC INTRACRANIAL ANEURYSM — VOLUME MESHING
# =============================================================
#
# Dataset: OpenNeuro ds005096
# Case: sub-038, October 2018
#
# Objective:
# Generate a tetrahedral CFD mesh from the segmented
# aneurysm and parent-vessel surface.
#
# This script will be developed incrementally, with
# geometry and mesh validation at each stage.
# =============================================================

from pathlib import Path
import gmsh
import numpy as np

# Resolve paths relative to the repository.
project_dir = Path(__file__).resolve().parent.parent

geometry_file = (
    project_dir
    / "geometry"
    / "patient_aneurysm_closed_corrected.stl"
)

# -------------------------------------------------------------
# IMPORT INDIVIDUAL CFD BOUNDARY SURFACES
# -------------------------------------------------------------
#
# Import the vessel wall and three caps separately.
# This allows us to inspect whether Gmsh retains the
# intended four surface entities.
#
# The cap numbers are provisional. Anatomical inlet/outlet
# assignments have not yet been verified.

surface_files = [
    "patient_vessel_wall.stl",
    "patient_cap_0.stl",
    "patient_cap_1.stl",
    "patient_cap_2.stl",
]

gmsh.initialize()
gmsh.model.add("patient_aneurysm")

for filename in surface_files:
    gmsh.merge(str(project_dir / "geometry" / filename))

# Inspect the imported surface entities.
surfaces = gmsh.model.getEntities(2)

print("Imported surface entities:", surfaces)

for dimension, surface_tag in surfaces:
    element_types, element_tags, _ = gmsh.model.mesh.getElements(
        dimension, surface_tag
    )

    triangle_count = sum(len(tags) for tags in element_tags)

    print(f"Surface {surface_tag}: {triangle_count} triangles")

# -------------------------------------------------------------
# CHECK NODE CONNECTIVITY BETWEEN IMPORTED SURFACES
# -------------------------------------------------------------
#
# Compare the imported node count with our previously
# validated closed surface containing 3,909 vertices.
#
# A larger count may indicate duplicated nodes along
# the boundaries between the vessel wall and caps.

node_tags, _, _ = gmsh.model.mesh.getNodes()

print("Total imported nodes:", len(node_tags))
print("Expected unique nodes:", 3909)

# -------------------------------------------------------------
# MERGE DUPLICATE NODES BETWEEN CFD BOUNDARY SURFACES
# -------------------------------------------------------------
#
# The four STL files were imported independently.
# Merge coincident mesh nodes so adjacent surfaces share
# the same vertices along their boundaries.

gmsh.model.mesh.removeDuplicateNodes()

# Verify the resulting number of nodes.
node_tags, _, _ = gmsh.model.mesh.getNodes()

print("Nodes after merging:", len(node_tags))

# -------------------------------------------------------------
# VERIFY CLOSED SURFACE CONNECTIVITY
# -------------------------------------------------------------
#
# Each edge in a closed triangular surface should belong
# to exactly two triangles.
#
# Count how many edges occur once, twice or more than twice.

from collections import Counter

edge_counts = Counter()

for _, surface_tag in gmsh.model.getEntities(2):
    element_types, _, node_groups = gmsh.model.mesh.getElements(
        2, surface_tag
    )

    for element_type, nodes in zip(element_types, node_groups):
        if element_type != 2:
            continue

        # Convert node IDs into a NumPy array of triangles.
        triangles = np.asarray(nodes).reshape(-1, 3)

        for a, b, c in triangles:
            edge_counts.update([
                tuple(sorted((int(a), int(b)))),
                tuple(sorted((int(b), int(c)))),
                tuple(sorted((int(c), int(a))))
            ])

print("Open edges:", sum(n == 1 for n in edge_counts.values()))
print("Manifold edges:", sum(n == 2 for n in edge_counts.values()))
print("Non-manifold edges:", sum(n > 2 for n in edge_counts.values()))

# -------------------------------------------------------------
# CREATE A VOLUME BOUNDED BY THE FOUR SURFACES
# -------------------------------------------------------------
#
# The vessel wall and three caps form the boundary of
# the intended fluid domain.
#
# First create the geometric representation needed by
# Gmsh for meshing the imported discrete surfaces.

gmsh.model.mesh.createGeometry()

# -------------------------------------------------------------
# INSPECT GEOMETRIC CURVES
# -------------------------------------------------------------
#
# Check whether Gmsh has created the shared boundary
# curves between the vessel wall and three caps.

curves = gmsh.model.getEntities(1)

print("Geometric curves:", curves)
print("Number of geometric curves:", len(curves))

print(
    "Geometric surfaces:",
    gmsh.model.getEntities(2)
)

# -------------------------------------------------------------
# DEFINE THE PATIENT-SPECIFIC FLUID VOLUME
# -------------------------------------------------------------
#
# The four surface entities enclose the intended blood-flow
# domain. Create a surface loop and the corresponding volume.

surface_tags = [tag for _, tag in gmsh.model.getEntities(2)]

surface_loop = gmsh.model.geo.addSurfaceLoop(surface_tags)

fluid_volume = gmsh.model.geo.addVolume([surface_loop])

gmsh.model.geo.synchronize()

print("Surface loop:", surface_loop)
print("Fluid volume:", fluid_volume)
print("Volume entities:", gmsh.model.getEntities(3))

# -------------------------------------------------------------
# DEFINE PHYSICAL GROUPS FOR CFD
# -------------------------------------------------------------
#
# Preserve the identities of the vessel wall, three caps,
# and enclosed fluid volume in the exported Gmsh mesh.
#
# Cap numbers remain provisional until anatomical flow
# direction has been verified.

boundary_groups = {
    "vessel_wall": 1,
    "cap_0": 2,
    "cap_1": 3,
    "cap_2": 4,
}

for name, surface_tag in boundary_groups.items():
    group_tag = gmsh.model.addPhysicalGroup(
        2, [surface_tag]
    )

    gmsh.model.setPhysicalName(2, group_tag, name)

# Assign the three-dimensional fluid region.
fluid_group = gmsh.model.addPhysicalGroup(
    3, [fluid_volume]
)

gmsh.model.setPhysicalName(
    3, fluid_group, "fluid"
)

# Verify the registered physical groups.
for dimension, group_tag in gmsh.model.getPhysicalGroups():
    print(
        "Physical group:",
        dimension,
        group_tag,
        gmsh.model.getPhysicalName(dimension, group_tag)
    )

# -------------------------------------------------------------
# PRESERVE THE IMPORTED SURFACE MESH
# -------------------------------------------------------------
#
# The STL surfaces already contain a closed, conforming
# triangular mesh.
#
# Attempt to prevent Gmsh from remeshing these surfaces
# while generating tetrahedra in the enclosed volume.

gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)

# -------------------------------------------------------------
# GENERATE THE FIRST PATIENT-SPECIFIC TETRAHEDRAL MESH
# -------------------------------------------------------------
#
# Generate a 3D mesh of the enclosed aneurysm fluid domain.
# This is an initial meshing test, not the final CFD mesh.

gmsh.model.mesh.generate(3)

# Report the number of tetrahedral elements.
element_types, element_tags, _ = gmsh.model.mesh.getElements(3)

tetrahedra = sum(
    len(tags)
    for tags in element_tags
)

print("Generated 3D elements:", tetrahedra)
print("3D element types:", element_types)

gmsh.finalize()