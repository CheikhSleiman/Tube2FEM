# Tube2FEM - Tubular Objects FEM Simulations (v1.0.0)

<p align="center">
  <img src="figs/FigureLiverPV_repoCroped.jpg" alt="PDF Preview" width="80%">
  <br>
  <em>Image-based computational pipeline outlining the complete workflow for CFD and advection–diffusion–reaction simulations, starting from segmented vascular geometries.</em>
</p>


**Tube2FEM** is a general-purpose highly-automated pipeline for flow related processes in (embedded) tubular objects. 
Addressing a critical gap in computational fluid dynamics (CFD) and simulation sciences, it facilitates the transition from raw three-dimensional imaging, graph networks, or Computer Aided Design (CAD) models of tubular objects to refined, simulation-ready meshes.
The pipeline leverages a range of open-source software and libraries, notably [GIBBON](https://www.gibboncode.org), [FEniCS](https://fenicsproject.org/), and [Paraview](https://www.paraview.org/), to provide
flexibility and broad applicability across different simulation scenarios, ranging from biomedical to industrial applications.

Features
--------

- Fully automated pipeline: from geometry pre-processing to simulation postprocessing
- Mesh generation from 3D images, graphs, or surfaces
- Centerline extraction
- 3D CFD, 3D-1D coupling, and advection-diffusion-reaction simulations
- Automated post-processing and visualization scripts
- Modular structure for custom problem setups


Interoperability
----------------
**Tube2FEM** uses multiple Python packages and MATLAB libraries (see Figure below)

<p align="center">
  <img src="figs/InteroperabilityDarkMode.jpg" alt="PDF Preview" width="80%">
  <br>
  <em>Tube2FEM packages/libraries interoperability.</em>
</p>

Installation
------------
Tube2FEM runs across both Windows and WSL environments. Please follow the installation guidlines [here](https://github.com/CheikhSleiman/Tube2FEM/blob/main/INSTALLATION_GUIDELINES.md)


Repository Structure
--------------------
```
Tube2FEM
├── CaseStudies
│   ├── Input
│   │   ├──── Graph
│   │	└──── segmentedCT
│   │
│   ├── Mesh
│   │   ├──── surfaceMesh
│   │	└──── volumeMesh 
│   │
│   ├── Centreline
│   │
│   ├── FiniteElement
│   │   ├──── CFD
│   │	├──── 3D-1D
│   │   └──── AdvectionDiffusionReaction
│   │
│   ├── Postprocessing
│   │   ├──── ParaView
│   │   ├──── Animations
│   │   ├──── Figures
│   │   └──── simOutput
│   │
│   └── Main.m
│
├── src
│   ├── FiniteElement
│   ├── boundaryConditions
│   ├── postprocessingParaView
│   ├── skeletonisation
│   ├── surfaceMeshProcessing
│   └── volumetricMesh
│
├── tests
├── figs
├── CONTRIBUTING.md
├── CODE_OF_CONDUCT.md
├── INSTALLING_GUIDLINES.md
├── README.md
├── LICENSE.md
└── .gitignore
```

Authors
-------
Tube2FEM is developed by:

  * Hani Cheikh Sleiman
  * Shiyu Wang

Licence
-------
This program is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License [GPL](https://github.com/CheikhSleiman/Tube2FEM/blob/main/LICENSE.md) as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.


Citing
------
If you wish to use Tube2FEM for journal publications, please refer to the [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.06.22.600203v1)

Contact
-------
Please report bugs and other issues through the issue tracker at:
[https://github.com/CheikhSleiman/Tube2FEM/issues](https://github.com/CheikhSleiman/Tube2FEM/issues)


