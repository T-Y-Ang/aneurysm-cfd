
from pathlib import Path

import numpy as np
import pyvista as pv
import meshio

# -------------------------------------------------------------
# EXPORT TETGEN VOLUME MESH TO GMSH FORMAT
# -------------------------------------------------------------
#
# Preserve four distinct boundary groups:
#   vessel_wall
#   opening_0
#   opening_1
#   opening_2
#
# Coordinates remain in millimetres.
# Conversion to metres occurs after OpenFOAM import.
# -------------------------------------------------------------

# Read the refined TetGen mesh, preserving the baseline files.
INPUT = Path("mesh/patient_aneurysm_tetgen_refined.vtu")

# Export to a separate Gmsh file for OpenFOAM conversion.
OUTPUT = Path("mesh/patient_aneurysm_tetgen_refined.msh")

volume = pv.read(str(INPUT))

# -------------------------------------------------------------
# 1. EXTRACT TETRAHEDRAL ELEMENTS
# -------------------------------------------------------------

tetrahedra = volume.cells.reshape(-1, 5)[:, 1:]

if len(tetrahedra) != volume.n_cells:
    raise RuntimeError("Unexpected volume-cell structure.")

points = np.asarray(volume.points)

# -------------------------------------------------------------
# 2. EXTRACT EXTERIOR TRIANGULAR FACES
# -------------------------------------------------------------

exterior = volume.extract_surface(
    algorithm="dataset_surface"
).triangulate()

surface_faces = exterior.faces.reshape(-1, 4)[:, 1:]

# Map exterior-surface point indices to original volume indices.
original_point_ids = exterior.point_data["vtkOriginalPointIds"]

triangles = original_point_ids[surface_faces]

triangle_xyz = points[triangles]

# -------------------------------------------------------------
# 3. CLASSIFY THE THREE OPENINGS
# -------------------------------------------------------------

tolerance = 1e-3

opening_0 = np.all(
    np.abs(triangle_xyz[:, :, 0] - 27.0) < tolerance,
    axis=1
)

opening_2 = np.all(
    np.abs(triangle_xyz[:, :, 0] - 42.0) < tolerance,
    axis=1
)

oblique_constant = 35.8347 - (-10.3963)

opening_1 = np.all(
    np.abs(
        triangle_xyz[:, :, 0]
        - triangle_xyz[:, :, 1]
        - oblique_constant
    ) < tolerance,
    axis=1
)

assigned = (
    opening_0.astype(int)
    + opening_1.astype(int)
    + opening_2.astype(int)
)

if np.any(assigned > 1):
    raise RuntimeError("Overlapping boundary classifications.")

vessel_wall = assigned == 0

# -------------------------------------------------------------
# 4. DEFINE PHYSICAL GROUPS
# -------------------------------------------------------------

groups = {
    "vessel_wall": (1, vessel_wall),
    "opening_0": (2, opening_0),
    "opening_1": (3, opening_1),
    "opening_2": (4, opening_2),
}

triangle_physical = np.zeros(len(triangles), dtype=np.int32)

for name, (physical_id, mask) in groups.items():
    triangle_physical[mask] = physical_id

    print(
        f"{name}: {np.count_nonzero(mask)} triangles"
    )

if np.any(triangle_physical == 0):
    raise RuntimeError("Unassigned exterior triangles.")

# -------------------------------------------------------------
# 5. WRITE GMSH MESH
# -------------------------------------------------------------

# Gmsh physical ID 5 represents the fluid volume.
fluid_physical = np.full(
    len(tetrahedra), 5, dtype=np.int32
)

mesh = meshio.Mesh(
    points=points,
    cells=[
        ("triangle", triangles),
        ("tetra", tetrahedra),
    ],
    cell_data={
        "gmsh:physical": [
            triangle_physical,
            fluid_physical,
        ],
        "gmsh:geometrical": [
            triangle_physical,
            fluid_physical,
        ],
    },
    field_data={
        "vessel_wall": np.array([1, 2]),
        "opening_0": np.array([2, 2]),
        "opening_1": np.array([3, 2]),
        "opening_2": np.array([4, 2]),
        "fluid": np.array([5, 3]),
    },
)

meshio.write(
    str(OUTPUT),
    mesh,
    file_format="gmsh22",
    binary=False,
)

print("\nEXPORT COMPLETE")
print("Tetrahedra:", len(tetrahedra))
print("Exterior triangles:", len(triangles))
print("Saved:", OUTPUT)
