
# -------------------------------------------------------------
# ANEURYSM CFD — 3D MESH GENERATION
# -------------------------------------------------------------
#
# Purpose:
# Import our refined synthetic aneurysm geometry into Gmsh.
#
# We will subsequently generate a tetrahedral mesh for
# computational fluid dynamics (CFD) simulation.
#
# All geometric dimensions are expressed in metres.


from pathlib import Path
import gmsh


# -------------------------------------------------------------
# 1. DEFINE PROJECT PATHS
# -------------------------------------------------------------

# Locate the repository root relative to this script.
project_directory = Path(__file__).resolve().parents[1]

# Path to our refined aneurysm STEP geometry.
geometry_file = (
    project_directory
    / "geometry"
    / "sidewall_aneurysm_refined.step"
)


# -------------------------------------------------------------
# 2. INITIALISE GMSH
# -------------------------------------------------------------

gmsh.initialize()
gmsh.model.add("aneurysm_mesh")


# -------------------------------------------------------------
# 3. IMPORT THE REFINED GEOMETRY
# -------------------------------------------------------------

gmsh.model.occ.importShapes(str(geometry_file))

# Synchronise OpenCASCADE geometry with the Gmsh model.
gmsh.model.occ.synchronize()


# -------------------------------------------------------------
# 4. VERIFY THE IMPORTED FLUID DOMAIN
# -------------------------------------------------------------

volumes = gmsh.model.getEntities(dim=3)

print("Imported geometry:", geometry_file.name)
print("Number of 3D volumes:", len(volumes))
print("Volume entities:", volumes)

# Stop if the geometry does not contain exactly one volume.
assert len(volumes) == 1, (
    "Expected exactly one connected fluid volume."
)


# -------------------------------------------------------------
# 4A. INSPECT FLUID BOUNDARY SURFACES
# -------------------------------------------------------------
#
# Each surface forms part of the boundary of our fluid volume.
#
# We will identify:
# - Inlet: circular opening at x = 0 m
# - Outlet: circular opening at x = 0.030 m
# - Wall: remaining vessel and aneurysm surfaces
#
# Bounding boxes help us identify these surfaces
# without relying on their Gmsh entity numbers.

surfaces = gmsh.model.getBoundary(
    volumes,
    oriented=False,
    recursive=False
)

print("\nBoundary surface information:")

for dimension, surface_tag in surfaces:

    bounding_box = gmsh.model.getBoundingBox(
        dimension,
        surface_tag
    )

    print(
        f"Surface {surface_tag}: "
        f"x = {bounding_box[0]:.6f} to "
        f"{bounding_box[3]:.6f} m"
    )


# -------------------------------------------------------------
# 4B. IDENTIFY AND LABEL PHYSICAL BOUNDARIES
# -------------------------------------------------------------
#
# Physical groups tell OpenFOAM which mesh faces belong
# to the inlet, outlet and vessel walls.
#
# The parent artery extends from x = 0 to x = 0.030 m.
# We classify surfaces using their X-coordinate ranges.

inlet_surfaces = []
outlet_surfaces = []
wall_surfaces = []

# Numerical tolerance for comparing coordinates (metres).
tolerance = 1e-6

artery_length = 0.030

for dimension, surface_tag in surfaces:

    xmin, ymin, zmin, xmax, ymax, zmax = (
        gmsh.model.getBoundingBox(dimension, surface_tag)
    )

    # An inlet surface lies entirely at x = 0.
    if abs(xmin) < tolerance and abs(xmax) < tolerance:
        inlet_surfaces.append(surface_tag)

    # An outlet surface lies entirely at x = 0.030 m.
    elif (
        abs(xmin - artery_length) < tolerance
        and abs(xmax - artery_length) < tolerance
    ):
        outlet_surfaces.append(surface_tag)

    # All remaining boundary surfaces are vessel walls.
    else:
        wall_surfaces.append(surface_tag)


# Confirm that exactly one inlet and outlet were identified.
assert len(inlet_surfaces) == 1, (
    f"Expected one inlet, found {inlet_surfaces}"
)

assert len(outlet_surfaces) == 1, (
    f"Expected one outlet, found {outlet_surfaces}"
)

assert len(wall_surfaces) > 0, (
    "No vessel-wall surfaces were identified."
)


# -------------------------------------------------------------
# 4C. CREATE GMSH PHYSICAL GROUPS
# -------------------------------------------------------------

inlet_group = gmsh.model.addPhysicalGroup(
    2, inlet_surfaces
)
gmsh.model.setPhysicalName(2, inlet_group, "inlet")

outlet_group = gmsh.model.addPhysicalGroup(
    2, outlet_surfaces
)
gmsh.model.setPhysicalName(2, outlet_group, "outlet")

wall_group = gmsh.model.addPhysicalGroup(
    2, wall_surfaces
)
gmsh.model.setPhysicalName(2, wall_group, "walls")

fluid_group = gmsh.model.addPhysicalGroup(
    3, [volumes[0][1]]
)
gmsh.model.setPhysicalName(3, fluid_group, "fluid")


# -------------------------------------------------------------
# 4D. VERIFY PHYSICAL GROUP ASSIGNMENTS
# -------------------------------------------------------------

print("\nPhysical boundary groups:")
print("Inlet surfaces:", inlet_surfaces)
print("Outlet surfaces:", outlet_surfaces)
print("Wall surfaces:", wall_surfaces)
print("Fluid volume:", volumes[0][1])


# -------------------------------------------------------------
# 5. DEFINE INITIAL MESH SIZE
# -------------------------------------------------------------
#
# Start with a coarse mesh to test whether the geometry
# can be discretised successfully.
#
# 0.001 m = 1 mm
#
# This is a preliminary mesh, not the final CFD mesh.

mesh_size = 0.001

gmsh.option.setNumber(
    "Mesh.MeshSizeMin",
    mesh_size
)

gmsh.option.setNumber(
    "Mesh.MeshSizeMax",
    mesh_size
)


# -------------------------------------------------------------
# 6. GENERATE THREE-DIMENSIONAL MESH
# -------------------------------------------------------------
#
# Dimension 3 tells Gmsh to generate a volume mesh.
#
# Gmsh will first mesh the boundary surfaces and then
# generate tetrahedral elements inside the fluid domain.

print("\nGenerating 3D tetrahedral mesh...")

gmsh.model.mesh.generate(3)


# -------------------------------------------------------------
# 7. REPORT MESH STATISTICS
# -------------------------------------------------------------

node_tags, _, _ = gmsh.model.mesh.getNodes()

element_types, element_tags, _ = (
    gmsh.model.mesh.getElements(dim=3)
)

number_of_elements = sum(
    len(tags) for tags in element_tags
)

print("\nMesh generation completed.")
print("Number of nodes:", len(node_tags))
print("Number of 3D elements:", number_of_elements)
print("3D element types:", element_types)


# -------------------------------------------------------------
# 8. SAVE THE MESH
# -------------------------------------------------------------
#
# Save in Gmsh's native MSH format.
#
# Later, we will establish physical boundary groups
# for the inlet, outlet and arterial wall.

output_file = (
    project_directory
    / "mesh"
    / "aneurysm_coarse.msh"
)


# -------------------------------------------------------------
# EXPORT MESH IN OPENFOAM-COMPATIBLE FORMAT
# -------------------------------------------------------------

# Use Gmsh MSH 2.2 ASCII format for OpenFOAM compatibility.
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)

# Save the mesh with its physical boundary groups.
gmsh.write(str(output_file))

# Confirm the output file location.
print("Mesh saved to:", output_file)



# -------------------------------------------------------------
# 9. CLOSE GMSH
# -------------------------------------------------------------

gmsh.finalize()

