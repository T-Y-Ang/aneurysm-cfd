
# -------------------------------------------------------------
# SYNTHETIC SIDEWALL INTRACRANIAL ANEURYSM
# -------------------------------------------------------------
#
# Purpose:
# Construct a simplified three-dimensional aneurysm geometry
# for computational fluid dynamics (CFD) simulation.
#
# The model will consist of:
# 1. A straight cylindrical parent artery
# 2. A spherical aneurysm sac projecting from its side
#
# Both structures will eventually be combined into a single
# continuous fluid domain.
#
# All dimensions are expressed in metres (SI units).

from pathlib import Path
import gmsh


# -------------------------------------------------------------
# 1. DEFINE PARENT ARTERY PARAMETERS
# -------------------------------------------------------------

artery_diameter = 0.005       # 5 mm
artery_radius = artery_diameter / 2
artery_length = 0.030         # 30 mm


# -------------------------------------------------------------
# 2. DEFINE ANEURYSM PARAMETERS
# -------------------------------------------------------------

aneurysm_diameter = 0.006     # 6 mm
aneurysm_radius = aneurysm_diameter / 2

# Position the aneurysm halfway along the parent artery.
aneurysm_center_x = artery_length / 2

# The aneurysm projects in the positive Y direction.
#
# Its centre is placed 4.5 mm above the artery centreline.
# This creates an overlap between the sphere and cylinder,
# allowing them to be joined into one fluid domain.

aneurysm_center_y = 0.0045

# Both structures are centred at z = 0.
aneurysm_center_z = 0.0


# -------------------------------------------------------------
# 3. DISPLAY PARAMETERS FOR VERIFICATION
# -------------------------------------------------------------

print("Parent artery diameter:", artery_diameter * 1000, "mm")
print("Parent artery length:", artery_length * 1000, "mm")

print("Aneurysm diameter:", aneurysm_diameter * 1000, "mm")
print(
    "Aneurysm centre:",
    (
        aneurysm_center_x,
        aneurysm_center_y,
        aneurysm_center_z
    )
)


# -------------------------------------------------------------
# 4. INITIALISE GMSH
# -------------------------------------------------------------
#
# Create a new Gmsh model for our synthetic aneurysm.

gmsh.initialize()
gmsh.model.add("sidewall_aneurysm")


# -------------------------------------------------------------
# 5. CREATE THE PARENT ARTERY
# -------------------------------------------------------------
#
# A cylinder extending along the X-axis.
# This represents the blood-filled arterial lumen.

artery = gmsh.model.occ.addCylinder(
    0, 0, 0,
    artery_length, 0, 0,
    artery_radius
)


# -------------------------------------------------------------
# 6. CREATE THE ANEURYSM SAC
# -------------------------------------------------------------
#
# A sphere positioned above the parent artery.
# Part of the sphere overlaps the cylinder.

aneurysm = gmsh.model.occ.addSphere(
    aneurysm_center_x,
    aneurysm_center_y,
    aneurysm_center_z,
    aneurysm_radius
)


# -------------------------------------------------------------
# 7. SYNCHRONISE THE GEOMETRY
# -------------------------------------------------------------
#
# Make the OpenCASCADE geometry available to the Gmsh model.

gmsh.model.occ.synchronize()


# -------------------------------------------------------------
# 8. VERIFY BOTH VOLUMES
# -------------------------------------------------------------

volumes = gmsh.model.getEntities(dim=3)

print("Number of 3D volumes:", len(volumes))
print("Volume entities:", volumes)



# -------------------------------------------------------------
# 9. COMBINE ARTERY AND ANEURYSM
# -------------------------------------------------------------
#
# Perform a Boolean union (fuse).
#
# This combines the overlapping cylindrical artery and
# spherical aneurysm into one continuous fluid domain.
#
# removeObject=True and removeTool=True remove the
# original solids after the union.

fused_volumes, _ = gmsh.model.occ.fuse(
    [(3, artery)],
    [(3, aneurysm)],
    removeObject=True,
    removeTool=True
)

# Synchronise the updated geometry.
gmsh.model.occ.synchronize()


# -------------------------------------------------------------
# 10. VERIFY THE FUSED GEOMETRY
# -------------------------------------------------------------

final_volumes = gmsh.model.getEntities(dim=3)

print("Number of volumes after fusion:", len(final_volumes))
print("Fused volume entities:", final_volumes)


# -------------------------------------------------------------
# 11. EXPORT THE GEOMETRY
# -------------------------------------------------------------
#
# Save the combined aneurysm geometry as a STEP file.

output_directory = Path(__file__).resolve().parent
output_file = output_directory / "sidewall_aneurysm.step"

gmsh.write(str(output_file))

print("Geometry saved to:", output_file)


# -------------------------------------------------------------
# 12. CLOSE GMSH
# -------------------------------------------------------------

gmsh.finalize()

