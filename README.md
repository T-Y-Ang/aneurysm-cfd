# Intracranial Aneurysm Haemodynamics — Computational Fluid Dynamics

**A computational investigation of blood-flow characteristics in a synthetic intracranial saccular aneurysm using OpenFOAM.**

**Project status:** In development

## Clinical Background

Intracranial aneurysms are pathological dilatations of cerebral arteries that may rupture, resulting in subarachnoid haemorrhage.

Although aneurysm size, location, morphology and patient-related factors inform clinical risk assessment, aneurysm haemodynamics may provide additional insights into the biological mechanisms underlying aneurysm development, progression and rupture.

Computational fluid dynamics (CFD) enables numerical investigation of blood-flow behaviour within vascular geometries, including velocity distributions, pressure gradients, flow recirculation and wall shear stress.

## Research Question

How do changes in intracranial aneurysm geometry influence local haemodynamic characteristics, particularly intra-aneurysmal flow patterns and wall shear stress?

## Project Objectives

1. Construct a simplified three-dimensional model of an intracranial sidewall saccular aneurysm.
2. Generate and validate a computational mesh.
3. Simulate incompressible blood flow using OpenFOAM.
4. Quantify velocity, pressure and wall shear stress distributions.
5. Visualise flow patterns and regions of intra-aneurysmal recirculation.
6. Compare haemodynamic characteristics between different synthetic aneurysm geometries.
7. Evaluate numerical reliability and document modelling assumptions and limitations.

## Patient-Derived Aneurysm CFD Demonstration (Planned)

In addition to the synthetic aneurysm model, this project aims to demonstrate the application of computational fluid dynamics to a single intracranial aneurysm identified on publicly available, de-identified neurovascular imaging.

### Objectives

1. Identify a suitable open-access CTA, MRA or 3D rotational angiography dataset containing an intracranial aneurysm.
2. Select one imaging case with sufficient anatomical resolution for vascular reconstruction.
3. Segment the aneurysm and parent vessel from the imaging data.
4. Generate a three-dimensional vascular surface model and computational mesh.
5. Apply simplified physiological flow conditions and perform CFD simulation using OpenFOAM.
6. Visualise velocity distributions, intra-aneurysmal flow patterns, pressure and wall shear stress.
7. Compare the practical challenges of image-derived modelling against the synthetic aneurysm simulation.

### Intended Outputs

- Representative angiographic images showing the aneurysm
- Three-dimensional reconstruction of the aneurysm and parent vessel
- Computational mesh visualisation
- Velocity streamlines and flow recirculation patterns
- Wall shear stress distribution over the aneurysm surface

### Data and Methodological Considerations

The imaging dataset will be selected based on accessibility, anatomical suitability, image quality and permissions for research use and publication.

The demonstration will use one de-identified case and will not attempt to predict patient-specific rupture risk.

Where patient-specific physiological measurements are unavailable, flow boundary conditions will be based on explicitly documented assumptions. The resulting haemodynamic measurements will therefore represent an exploratory simulation rather than a clinically validated patient-specific assessment.

**Status:** Planned extension, following completion of the synthetic aneurysm CFD model.

## From Aneurysm Geometry to Computational Mesh

### 1. Constructing a Synthetic Intracranial Aneurysm

A simplified three-dimensional model of a sidewall saccular intracranial aneurysm was constructed using Python and Gmsh.

The model consists of a straight parent artery measuring 5 mm in diameter and 30 mm in length, with a 6 mm diameter aneurysm sac arising from its sidewall.

The artery and aneurysm were combined into a continuous fluid domain. The aneurysm neck was subsequently refined using a geometric fillet to create a smoother transition between the parent artery and aneurysm sac.

![Comparison of the initial and refined synthetic aneurysm geometry](outputs/figures/aneurysm_neck_refinement_comparison.png)

*Figure 1. Comparison of the initial synthetic sidewall aneurysm and the model following geometric refinement of the aneurysm neck. The geometry is illustrative and does not represent patient-specific anatomy.*

### 2. Generating the Computational Mesh

To simulate blood flow, the continuous three-dimensional aneurysm geometry must be divided into many smaller computational cells. This process is called **meshing**.

Gmsh was used to generate an initial tetrahedral mesh, comprising 3,882 three-dimensional cells and 1,079 nodes.


![Surface mesh of the synthetic intracranial aneurysm](outputs/figures/aneurysm_coarse_surface_mesh.png)

*Figure 2. Preliminary computational mesh of the synthetic sidewall intracranial aneurysm. The triangular elements visible on the surface form the boundaries of the three-dimensional tetrahedral mesh. The aneurysm sac and parent artery are represented as a continuous fluid domain.*


The mesh was assigned three boundary regions:

- **Inlet:** The opening through which blood enters the parent artery.
- **Outlet:** The opening through which blood leaves the parent artery.
- **Vessel walls:** The surfaces of the parent artery and aneurysm sac.

These boundaries allow different physiological conditions to be specified during subsequent CFD simulations.

### 3. Validating the Computational Mesh

The mesh was exported from Gmsh and successfully imported into OpenFOAM 14.

OpenFOAM's `checkMesh` utility was used to assess the mesh's geometric and topological integrity.

| Mesh characteristic | Result |
|---|---|
| Three-dimensional cells | 3,882 tetrahedra |
| Mesh nodes | 1,079 |
| Connected fluid regions | 1 |
| Boundary regions | Inlet, outlet, walls |
| Maximum cell aspect ratio | 7.97 |
| Maximum non-orthogonality | 69.06° |
| Maximum skewness | 1.038 |
| OpenFOAM mesh validation | **Passed (`Mesh OK`)** |

The initial mesh passed OpenFOAM's standard quality checks, establishing a suitable starting point for preliminary numerical experiments.

However, passing these checks does not demonstrate that the mesh is sufficiently refined for accurate haemodynamic measurements. Further refinement, particularly near the aneurysm neck and vessel walls, and a mesh-convergence study will be required before interpreting quantitative wall shear stress results.

**Current milestone:** Synthetic geometry construction, preliminary tetrahedral meshing and initial OpenFOAM mesh validation completed. Blood-flow simulation remains in development.

---


## Computational Methodology

### Geometry

The initial model will consist of a straight cylindrical parent artery with a synthetic sidewall saccular aneurysm.

Preliminary geometric parameters:

- Parent artery diameter: 5 mm
- Aneurysm sac diameter: approximately 6 mm
- Aneurysm neck: to be defined during geometry construction

These dimensions are illustrative and do not represent patient-specific anatomy.

### Fluid Dynamics

Blood flow will initially be modelled using the incompressible Navier–Stokes equations.

Initial assumptions include:

- Newtonian fluid behaviour
- Incompressible flow
- Rigid vessel walls
- No-slip wall boundary conditions
- Laminar flow
- Simplified inlet and outlet boundary conditions

These assumptions may be refined in subsequent experiments.

### Numerical Simulation

OpenFOAM 14 will be used to solve the governing fluid-flow equations using the finite-volume method.

Mesh quality will be assessed before simulation. Numerical reliability will be investigated through mesh-refinement comparisons and monitoring of solver convergence and conservation properties.

### Haemodynamic Analysis

Planned outputs include:

- Velocity magnitude and vector fields
- Flow streamlines and recirculation patterns
- Pressure distribution
- Wall shear stress magnitude and spatial distribution
- Comparative haemodynamic measurements between aneurysm geometries


## Technical Documentation

For readers interested in the engineering methodology, detailed technical documentation is available:

- [Technical Methodology](docs/methodology.md) — Geometry construction, mesh generation, boundary classification, OpenFOAM mesh validation and numerical limitations.

The documentation will be expanded as the project progresses through blood-flow simulation, mesh refinement and haemodynamic analysis.


## Project Structure

- `geometry/` — Parametric geometry generation
- `mesh/` — Mesh generation scripts and configuration
- `cases/` — Reproducible OpenFOAM case files
- `scripts/` — Simulation and post-processing utilities
- `notebooks/` — Exploratory numerical analysis
- `outputs/figures/` — Selected visualisations
- `docs/` — Technical documentation and limitations

Large datasets, generated meshes and intermediate solver outputs will not be stored in the repository.

## Software

- OpenFOAM 14
- Python
- Gmsh 4.15.2
- ParaView (planned)
- Ubuntu 26.04 LTS under WSL2


## Current Progress

### Environment and Initial Validation
- [x] Install and configure OpenFOAM 14
- [x] Run a standard cavity-flow validation tutorial
- [x] Verify tutorial mesh quality and solver execution

### Synthetic Aneurysm Model
- [x] Construct a three-dimensional synthetic sidewall aneurysm
- [x] Refine the aneurysm neck geometry
- [x] Generate a preliminary tetrahedral computational mesh using Gmsh
- [x] Define inlet, outlet and vessel-wall boundary regions
- [x] Import the mesh into OpenFOAM 14
- [x] Validate the preliminary mesh using `checkMesh`
- [ ] Configure physiological flow conditions
- [ ] Run the initial aneurysm CFD simulation
- [ ] Visualise velocity, pressure and flow streamlines
- [ ] Refine the computational mesh and assess mesh convergence
- [ ] Analyse wall shear stress distributions
- [ ] Compare haemodynamic characteristics between aneurysm geometries

### Patient-Derived Aneurysm Demonstration
- [ ] Identify a suitable publicly available neurovascular imaging dataset
- [ ] Reconstruct an image-derived aneurysm geometry
- [ ] Generate and validate an image-derived computational mesh
- [ ] Perform an exploratory image-derived CFD simulation
- [ ] Generate comparative visualisations for the portfolio

### Documentation
- [x] Document synthetic geometry construction and initial mesh validation
- [ ] Create an annotated computational mesh figure
- [ ] Document simulation parameters and numerical verification
- [ ] Complete project results, limitations and interpretation



## Limitations

This project is intended as an educational and exploratory computational investigation.

The initial geometry is synthetic, and the model will use simplified assumptions regarding vascular anatomy, vessel-wall mechanics and blood rheology.

Wall shear stress and other CFD-derived measurements cannot independently establish an individual aneurysm's rupture risk.

Results should not be interpreted as clinically validated predictions or used to guide patient management.

## Future Directions

Potential extensions include:

- Pulsatile physiological flow conditions
- More advanced haemodynamic indices
- Fluid–structure interaction and vessel-wall mechanics
- Investigation of aneurysm morphology and haemodynamic relationships

---

*Part of a neurosurgery-focused portfolio exploring computational modelling, machine learning and biomedical data science.*
