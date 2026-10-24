
import math
import gmsh

# -------------------------------------------------------------
# STEP 103: CLASSIFY THE CLOSED ANEURYSM SURFACE
# -------------------------------------------------------------
#
# Separate the imported STL into discrete surface patches using a 40-degree feature-angle threshold.
#
# This does not modify the original STL.
# -------------------------------------------------------------

INPUT = "geometry/patient_aneurysm_closed_remeshed.stl"
OUTPUT = "mesh/patient_closed_classified.msh"

gmsh.initialize()

try:
    gmsh.model.add("patient_closed_classified")

    # Load the validated, watertight STL.
    gmsh.merge(INPUT)

    # Identify surface patches separated by sharp features.
    gmsh.model.mesh.classifySurfaces(
        math.radians(40),
        True,
        True,
        math.pi
    )

    # Construct parametrizations for the classified patches.
    gmsh.model.mesh.createGeometry()

    # Report the resulting surfaces.
    surfaces = gmsh.model.getEntities(2)

    print("\nCLASSIFIED SURFACE REPORT")
    print("Number of surface patches:", len(surfaces))

    for dimension, tag in surfaces:
        _, triangle_nodes = gmsh.model.mesh.getElementsByType(
            2, tag
        )

        print(
            f"Surface {tag}: "
            f"{len(triangle_nodes) // 3} triangles"
        )

    # Preserve the classified mesh for subsequent steps.
    gmsh.write(OUTPUT)

    print("\nSaved:", OUTPUT)

finally:
    gmsh.finalize()
