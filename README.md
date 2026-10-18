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
- Gmsh (planned)
- ParaView (planned)
- Ubuntu 26.04 LTS under WSL2

## Current Progress

- [x] Install and configure OpenFOAM 14
- [x] Run a standard cavity-flow validation tutorial
- [x] Verify tutorial mesh quality and solver execution
- [ ] Generate synthetic aneurysm geometry
- [ ] Generate aneurysm computational mesh
- [ ] Configure physiological flow conditions
- [ ] Run aneurysm CFD simulation
- [ ] Analyse haemodynamic outputs
- [ ] Compare aneurysm geometries
- [ ] Complete numerical verification and documentation
- [ ] Identify a suitable public intracranial aneurysm imaging dataset
- [ ] Reconstruct one image-derived aneurysm geometry
- [ ] Perform a single image-derived CFD demonstration
- [ ] Generate comparative visualisations for the portfolio


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
