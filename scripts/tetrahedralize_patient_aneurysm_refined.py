
from pathlib import Path

import numpy as np
import pyvista as pv
import tetgen

# -------------------------------------------------------------
# PATIENT-SPECIFIC ANEURYSM TETRAHEDRAL MESH
# -------------------------------------------------------------
#
# Input:
#   Validated, watertight closed aneurysm surface (STL)
#
# Method:
#   TetGen constrained tetrahedralization
#
# Output:
#   Volume mesh in VTK UnstructuredGrid format (.vtu)
#
# Coordinate units:
#   Millimetres
#
# This is an initial volume-meshing attempt. Mesh quality
# and CFD suitability must be assessed separately.
# -------------------------------------------------------------

INPUT = Path(
    "geometry/patient_aneurysm_closed_remeshed.stl"
)

# Save the refined candidate separately from the baseline mesh.
# This preserves the original CFD mesh for comparison.
OUTPUT = Path(
    "mesh/patient_aneurysm_tetgen_refined.vtu"
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)


# -------------------------------------------------------------
# 1. LOAD THE CLOSED SURFACE
# -------------------------------------------------------------

surface = pv.read(str(INPUT))

# Ensure that all surface cells are triangles.
surface = surface.triangulate()

print("\nINPUT SURFACE")
print("Points:", surface.n_points)
print("Triangles:", surface.n_cells)

# Check for open or non-manifold edges.
boundary_edges = surface.extract_feature_edges(
    boundary_edges=True,
    non_manifold_edges=True,
    feature_edges=False,
    manifold_edges=False
)

print("Boundary/non-manifold edges:", boundary_edges.n_cells)

if boundary_edges.n_cells != 0:
    raise RuntimeError(
        "Surface is not watertight or has non-manifold edges."
    )


# -------------------------------------------------------------
# 2. GENERATE TETRAHEDRA
# -------------------------------------------------------------

print("\nGenerating tetrahedral volume mesh...")

tet = tetgen.TetGen(surface)

# TetGen generates tetrahedra inside the closed surface.
# Quality constraints are deliberately moderate initially.
# Generate the tetrahedral mesh.
# TetGen versions can differ in their returned values,
# so we access the generated mesh through the TetGen object.

# -------------------------------------------------------------
# REFINED TETRAHEDRAL MESH
# -------------------------------------------------------------
#
# Retain the baseline tetrahedral shape-quality settings.
#
# Introduce a maximum element volume of 0.002 mm³ to produce a finer candidate mesh.
#
# This is an experimental refinement target.
# Mesh quality and CFD suitability must be checked separately.
# -------------------------------------------------------------

# Request global tetrahedral volume refinement.
# The maximum volume is expressed in mm³ because the input vascular surface uses millimetres.
tet.tetrahedralize(
    order=1,
    minratio=1.2,
    mindihedral=15,
    fixedvolume=True,
    maxvolume=0.002
)

# Retrieve the generated node coordinates and
# four-node tetrahedral connectivity.
nodes = tet.node
elements = tet.elem


# -------------------------------------------------------------
# 3. VALIDATE BASIC MESH PROPERTIES
# -------------------------------------------------------------

nodes = np.asarray(nodes)
elements = np.asarray(elements)

print("\nTETRAHEDRAL MESH REPORT")
print("Nodes:", len(nodes))
print("Tetrahedra:", len(elements))

if len(elements) == 0:
    raise RuntimeError("TetGen generated zero tetrahedra.")

if elements.ndim != 2 or elements.shape[1] != 4:
    raise RuntimeError(
        "Expected four-node linear tetrahedra."
    )

# Calculate signed tetrahedral volumes.
a = nodes[elements[:, 0]]
b = nodes[elements[:, 1]]
c = nodes[elements[:, 2]]
d = nodes[elements[:, 3]]

signed_volumes = (
    np.einsum(
        "ij,ij->i",
        np.cross(b - a, c - a),
        d - a
    ) / 6.0
)

absolute_volumes = np.abs(signed_volumes)

print(
    "Minimum tetrahedron volume (mm³):",
    absolute_volumes.min()
)
print(
    "Total tetrahedral volume (mm³):",
    absolute_volumes.sum()
)
print(
    "Non-positive tetrahedron volumes:",
    np.count_nonzero(absolute_volumes <= 1e-12)
)

if np.any(absolute_volumes <= 1e-12):
    raise RuntimeError(
        "Mesh contains degenerate tetrahedra."
    )


# -------------------------------------------------------------
# 4. SAVE VOLUME MESH
# -------------------------------------------------------------

# TetGen exposes its volume mesh as a PyVista grid.
volume_mesh = tet.grid

volume_mesh.save(str(OUTPUT))


# -------------------------------------------------------------
# 5. VALIDATE THE EXTERIOR BOUNDARY OF THE VOLUME MESH
# -------------------------------------------------------------
#
# Extract the exterior triangular faces of the tetrahedral
# mesh and compare them against the original input surface.
#
# This checks boundary preservation before OpenFOAM conversion.
# -------------------------------------------------------------

exterior = volume_mesh.extract_surface().triangulate()

print("\nEXTERIOR BOUNDARY VALIDATION")
print("Exterior triangles:", exterior.n_cells)
print("Original triangles:", surface.n_cells)

# Check for open or non-manifold edges in the extracted boundary.
exterior_problem_edges = exterior.extract_feature_edges(
    boundary_edges=True,
    non_manifold_edges=True,
    feature_edges=False,
    manifold_edges=False
)

print(
    "Boundary/non-manifold edges:",
    exterior_problem_edges.n_cells
)

# Compare surface areas.
original_area = surface.area
exterior_area = exterior.area

relative_area_difference = (
    abs(exterior_area - original_area) / original_area
)

print("Original surface area (mm²):", original_area)
print("Exterior surface area (mm²):", exterior_area)
print(
    "Relative area difference (%):",
    100 * relative_area_difference
)

if exterior_problem_edges.n_cells != 0:
    raise RuntimeError(
        "Extracted tetrahedral boundary is not closed."
    )


print("\nSUCCESS")
print("Saved:", OUTPUT)
