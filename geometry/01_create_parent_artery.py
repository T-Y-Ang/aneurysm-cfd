
# -------------------------------------------------------------
# SYNTHETIC PARENT ARTERY GEOMETRY
# -------------------------------------------------------------
#
# Purpose:
# Create a simple cylindrical fluid domain representing
# a straight intracranial artery.
#
# This is the first geometric component of our aneurysm
# computational fluid dynamics (CFD) project.
#
# All dimensions are expressed in metres (SI units).

from pathlib import Path
import gmsh


# -------------------------------------------------------------
# 1. DEFINE GEOMETRIC PARAMETERS
# -------------------------------------------------------------

artery_diameter = 0.005   # 5 mm
artery_radius = artery_diameter / 2
artery_length = 0.030     # 30 mm


# -------------------------------------------------------------
# 2. INITIALISE GMSH
# -------------------------------------------------------------

gmsh.initialize()
gmsh.model.add("parent_artery")


# -------------------------------------------------------------
# 3. CREATE CYLINDRICAL FLUID DOMAIN
# -------------------------------------------------------------
#
# The cylinder extends along the x-axis.
#
# Its inlet is located at x = 0.
# Its outlet is located at x = 0.030 m.
#
# The cylinder represents the blood-filled region,
# not the arterial wall itself.

cylinder = gmsh.model.occ.addCylinder(
    0, 0, 0,              # Starting coordinates (x, y, z)
    artery_length, 0, 0,  # Direction and length
    artery_radius         # Radius
)

gmsh.model.occ.synchronize()


# -------------------------------------------------------------
# 4. SAVE GEOMETRY
# -------------------------------------------------------------
#
# STEP format stores the CAD geometry.
# We will use this later for meshing and visualisation.

output_directory = Path(__file__).resolve().parent
output_file = output_directory / "parent_artery.step"

gmsh.write(str(output_file))


# -------------------------------------------------------------
# 5. VERIFY GEOMETRY
# -------------------------------------------------------------

volumes = gmsh.model.getEntities(dim=3)

print("Number of 3D volumes:", len(volumes))
print("Volume entities:", volumes)
print("Geometry saved to:", output_file)


# -------------------------------------------------------------
# 6. CLOSE GMSH
# -------------------------------------------------------------

gmsh.finalize()
