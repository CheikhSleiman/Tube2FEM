# Tube2FEM - Tubular Objects FEM Simulations

Tube2FEM is a general-purpose highly-automated pipeline for flow related processes in (embedded) tubular objects. 
Addressing a critical gap in computational fluid dynamics (CFD) and simulation sciences, it facilitates the transition from raw three-dimensional imaging, graph networks, or Computer Aided Design (CAD) models of tubular objects to refined, simulation-ready meshes.
The pipeline leverages a range of open-source software and libraries, notably [GIBBON](https://www.gibboncode.org), [FEniCS](https://fenicsproject.org/), and [Paraview](https://www.paraview.org/), to provide
flexibility and broad applicability across different simulation scenarios, ranging from biomedical to industrial applications.



Dependencies & Interoperability
-------------------------------
Tube2FEM uses multiple Python packages and MATLAB libraries (see Figure below)

![PDF Preview](figs/InteroperabilityDarkMode.jpg)
<div align="center">Tube2FEM packages/libraries interoperability.</div>

For Windows users, Ubuntu 18.04 LTS WSL has to be installed, mainly to run the *FEniCS* simulations
  * [Matlab >R2020a](https://www.mathworks.com/) --> Win
  * [Gibbon](https://www.gibboncode.org) --> Win
  * [FEniCS 2019.1.0](https://fenicsproject.org/download/archive/) --> WSL
  * [Paraview](https://www.paraview.org/) --> Win
  * [Trimesh](https://trimesh.org/) --> Win python
  * [PyVista](https://pyvista.org/) --> Win python 
  * [VesselVio](https://github.com/JacobBumgarner/VesselVio) --> Win
  * [TreeSkel](https://github.com/ashkanpakzad/TreeSkel) --> WSL
  * [meshio](https://github.com/nschloe/meshio) --> WSL or Win python
  * [tifffile](https://pypi.org/project/tifffile) --> WSL or Win python

Repository Structure
--------------------
```
Tube2FEM
├── Problem
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
├── README.md (this file)
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
This program is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License [GPL](https://github.com/CheikhSleiman/Tube2FEM/blob/main/LICENSE.txt) as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.


Citing
------
If you wish to use Tube2FEM for journal publications, please refer to the [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.06.22.600203v1)



