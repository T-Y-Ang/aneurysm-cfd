
# -------------------------------------------------------------
# ANEURYSM NECK REFINEMENT — GEOMETRY INSPECTION
# -------------------------------------------------------------
#
# Purpose:
# Load our existing aneurysm geometry and identify the
# edges that might be suitable for geometric smoothing.
#
# We will inspect the geometry before applying a fillet.

from pathlib import Path
import gmsh


# -------------------------------------------------------------
# 1. DEFINE INPUT FILE
# -------------------------------------------------------------

geometry_directory = Path(__file__).resolve().parent

input_file = geometry_directory / "sidewall_aneurysm.step"


# -------------------------------------------------------------
# 2. INITIALISE GMSH
# -------------------------------------------------------------

gmsh.initialize()
gmsh.model.add("aneurysm_neck_refinement")


# -------------------------------------------------------------
# 3. IMPORT EXISTING GEOMETRY
# -------------------------------------------------------------

gmsh.model.occ.importShapes(str(input_file))

gmsh.model.occ.synchronize()


# -------------------------------------------------------------
# 4. IDENTIFY VOLUMES AND EDGES
# -------------------------------------------------------------
#
# Dimension 3 = volumes
# Dimension 1 = curves or edges
#
# Not every edge represents the aneurysm neck.
# We will inspect their locations before selecting one.

volumes = gmsh.model.getEntities(dim=3)
edges = gmsh.model.getEntities(dim=1)

print("Number of 3D volumes:", len(volumes))
print("Number of geometric edges:", len(edges))

print("\nGeometric edges:")

for dim, tag in edges:
    # Bounding box coordinates:
    # xmin, ymin, zmin, xmax, ymax, zmax
    bbox = gmsh.model.getBoundingBox(dim, tag)

    # Convert metres to millimetres for easier interpretation.
    bbox_mm = tuple(round(value * 1000, 2) for value in bbox)

    print(f"Edge {tag}: {bbox_mm} mm")



# -------------------------------------------------------------
# 5. TEST ANEURYSM NECK FILLET
# -------------------------------------------------------------
#
# Edge 4 was identified as the likely intersection between
# the spherical aneurysm and cylindrical parent artery.
#
# We will attempt a small fillet to smooth this junction.
#
# The radius is expressed in metres.

fillet_radius = 0.0003  # 0.3 mm

try:
    refined_volumes = gmsh.model.occ.fillet(
        [volumes[0][1]],  # Tag of the existing 3D volume
        [4],              # Candidate aneurysm neck edge
        [fillet_radius]
    )

    gmsh.model.occ.synchronize()

    print("\nFillet operation completed.")
    print("Resulting volume entities:", refined_volumes)

    # Export the refined geometry.
    output_file = geometry_directory / "sidewall_aneurysm_refined.step"

    gmsh.write(str(output_file))

    print("Refined geometry saved to:", output_file)

except Exception as error:
    print("\nFillet operation failed.")
    print("Error:", error)


# -------------------------------------------------------------
# 6. CLOSE GMSH
# -------------------------------------------------------------

gmsh.finalize()

