# Technical Methodology — Intracranial Aneurysm CFD

## 1. Project Scope

This project investigates computational fluid dynamics in a simplified, synthetic sidewall intracranial aneurysm.

The initial objective is to establish a reproducible workflow for geometry construction, computational meshing, numerical simulation and haemodynamic visualisation.

The synthetic geometry is not derived from patient imaging.

## 2. Geometry Construction

**Software:** Python and Gmsh 4.15.2

| Parameter | Value |
|---|---|
| Parent artery diameter | 5 mm |
| Parent artery length | 30 mm |
| Aneurysm sac diameter | 6 mm |
| Aneurysm type | Synthetic sidewall saccular aneurysm |
| Neck refinement | 0.3 mm geometric fillet |
| Geometry format | STEP |

The parent artery was constructed as a cylindrical volume, and the aneurysm sac was generated from a sphere intersecting the arterial wall.

Both structures were combined into a single continuous three-dimensional fluid domain.

A geometric fillet was applied at the aneurysm neck to smooth the transition between the sac and parent artery.

The resulting geometry was exported as a STEP file for subsequent mesh generation.

## 3. Computational Mesh Generation

**Software:** Gmsh 4.15.2

The refined geometry was imported into Gmsh and discretised into an unstructured tetrahedral mesh.

The initial mesh was generated using a nominal characteristic element size of 1 mm.

### Initial Mesh Statistics

| Parameter | Result |
|---|---|
| Mesh nodes | 1,079 |
| Tetrahedral volume cells | 3,882 |
| Computational fluid regions | 1 |
| Inlet boundary faces | 64 |
| Outlet boundary faces | 64 |
| Vessel-wall boundary faces | 1,476 |
| Mesh format for OpenFOAM conversion | Gmsh MSH 2.2 ASCII |

### Boundary Classification

Three surface boundary groups were defined:

- `inlet` — Proximal opening of the parent artery
- `outlet` — Distal opening of the parent artery
- `walls` — Parent artery and aneurysm surfaces

A single three-dimensional physical group, `fluid`, was defined for the computational volume.

## 4. OpenFOAM Mesh Conversion and Validation

**Software:** OpenFOAM Foundation 14

The Gmsh mesh was converted into OpenFOAM's native mesh representation using `gmshToFoam`.

Mesh validation was performed using `checkMesh`.

### Mesh Quality Results

| Quality metric | Result |
|---|---|
| Connected regions | 1 |
| Maximum aspect ratio | 7.97396 |
| Maximum non-orthogonality | 69.0638° |
| Mean non-orthogonality | 22.6217° |
| Maximum skewness | 1.03837 |
| Minimum cell volume | 1.14032 × 10⁻¹¹ m³ |
| Maximum cell volume | 4.91712 × 10⁻¹⁰ m³ |
| Total fluid volume | 6.83253 × 10⁻⁷ m³ |
| Mesh validation | Passed — `Mesh OK` |

The mesh passed OpenFOAM's standard geometric and topological checks.

These results establish that the initial computational mesh is structurally suitable for preliminary numerical experiments. They do not establish mesh-independent haemodynamic accuracy.

## 5. Planned Flow Simulation

The initial simulation will use a simplified incompressible, Newtonian, laminar blood-flow model with rigid vessel walls.

Provisional material properties:

- Blood density: 1,060 kg/m³
- Dynamic viscosity: 0.0035 Pa·s

The inlet velocity, outlet pressure conditions, numerical solver settings and convergence criteria will be documented once configured.

## 6. Numerical Verification and Limitations

The preliminary mesh is relatively coarse and does not include dedicated near-wall boundary-layer refinement.

Further work will include:

- Mesh refinement near the aneurysm neck and vessel walls
- Mesh-convergence assessment
- Solver convergence monitoring
- Assessment of mass conservation
- Evaluation of numerical sensitivity to boundary conditions

Wall shear stress and other haemodynamic measurements will not be interpreted quantitatively until appropriate numerical verification has been undertaken.

The model is an educational engineering demonstration and is not intended for patient-specific rupture-risk prediction or clinical decision-making.
