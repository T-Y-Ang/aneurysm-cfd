
import numpy as np
import pyvista as pv

# -------------------------------------------------------------
# VALIDATE TETGEN EXTERIOR BOUNDARY
# -------------------------------------------------------------
#
# Check whether the three known vessel opening planes
# remain represented in the tetrahedral mesh.
#
# This script reads the existing volume mesh.
# It does not regenerate or modify the mesh.
# -------------------------------------------------------------

INPUT = "mesh/patient_aneurysm_tetgen.vtu"

mesh = pv.read(INPUT)

# Extract the triangular exterior surface.
exterior = mesh.extract_surface(
    algorithm="dataset_surface"
).triangulate()

points = exterior.points

# Each triangle is represented by three vertex indices.
faces = exterior.faces.reshape(-1, 4)[:, 1:]

# Obtain the coordinates of each triangle's vertices.
triangles = points[faces]

# -------------------------------------------------------------
# IDENTIFY TRIANGLES LYING ON THE THREE OPENING PLANES
# -------------------------------------------------------------

# Coordinates are in millimetres.
tolerance = 1e-3

# Opening 0: x = 27 mm.
opening_0 = np.all(
    np.abs(triangles[:, :, 0] - 27.0) < tolerance,
    axis=1
)

# Opening 2: x = 42 mm.
opening_2 = np.all(
    np.abs(triangles[:, :, 0] - 42.0) < tolerance,
    axis=1
)

# Opening 1: oblique plane.
#
# Previously identified centroid:
# approximately (35.8347, -10.3963, -26.5043)
#
# Normal approximately (0.7071, -0.7071, 0).
#
# Therefore x - y is approximately constant.
oblique_constant = 35.8347 - (-10.3963)

opening_1 = np.all(
    np.abs(
        triangles[:, :, 0]
        - triangles[:, :, 1]
        - oblique_constant
    ) < tolerance,
    axis=1
)

# -------------------------------------------------------------
# CALCULATE TRIANGLE AREAS
# -------------------------------------------------------------

edge_1 = triangles[:, 1] - triangles[:, 0]
edge_2 = triangles[:, 2] - triangles[:, 0]

areas = 0.5 * np.linalg.norm(
    np.cross(edge_1, edge_2),
    axis=1
)

print("\nTETGEN OPENING VALIDATION")

for name, mask in [
    ("opening_0", opening_0),
    ("opening_1", opening_1),
    ("opening_2", opening_2),
]:
    print(f"\n{name}")
    print("Triangles:", np.count_nonzero(mask))
    print("Surface area (mm²):", areas[mask].sum())

# -------------------------------------------------------------
# CHECK FOR OVERLAPPING OPENING CLASSIFICATIONS
# -------------------------------------------------------------

overlap = (
    opening_0.astype(int)
    + opening_1.astype(int)
    + opening_2.astype(int)
)

print("\nTriangles assigned to multiple openings:",
      np.count_nonzero(overlap > 1))

print(
    "Remaining vessel-wall triangles:",
    np.count_nonzero(overlap == 0)
)
