# Step-by-Step Installation Guide (Windows + WSL)
-----------------------------------------------
Tube2FEM runs across both Windows and WSL environments. 

## System Compatibility

We recommend:

* **Windows 10/11** with **WSL (Ubuntu 18.04 LTS)**
* **MATLAB R2020a or later**
* **Python 3.8–3.10**


## Step 1 — Install Required Software

### 1.1 On Windows (Host)

Install the following **on Windows**:

| Tool                                               | Link                        | Notes                              |
| -------------------------------------------------- | --------------------------- | ---------------------------------- |
| [MATLAB ≥ R2020a](https://www.mathworks.com/)      | Commercial license required | Geometry preprocessing and meshing |
| [GIBBON Toolbox](https://www.gibboncode.org)       | GitHub project              | MATLAB open-source toolbox         |
| [ParaView](https://www.paraview.org/)              | Official site               | Post-processing and visualization  |
| [Python (≥3.8)](https://www.python.org/downloads/) | Download installer          | Needed for Python utilities        |

Install required Python packages on Windows:

- Install [Trimesh](https://trimesh.org/),[PyVista](https://pyvista.org/) and [tifffile](https://pypi.org/project/tifffile) libraries. 
Using Windows Command Prompt:

```bash
pip install trimesh pyvista tifffile
```

- Install [VesselVio](https://github.com/JacobBumgarner/VesselVio):
Tube2FEM requires running VesselVio from terminal (using the single-line executable VesselVio.py file) 
Follow the Windows build instructions provided by the developper [here](https://jacobbumgarner.github.io/VesselVio/Build.html).


---

### 1.2 On WSL (Ubuntu 18.04 LTS)

Open your WSL terminal and install the following:

Install legacy **FEniCS 2019.1.0** using:

```bash
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
```
install [meshio](https://github.com/nschloe/meshio):

```bash
pip install meshio[all]
```

If you prefer to use a different Skeletonisation algorithm then VesselVio, we also recommend [TreeSkel](https://github.com/ashkanpakzad/TreeSkel).
Clone and install :

```bash
git clone https://github.com/ashkanpakzad/TreeSkel.git
cd TreeSkel
pip3 install .
```

---

## Step 2 — Share Files Across Windows and WSL

Clone or Download the **Tube2FEM** repository to a directory accessible from both Windows and WSL:

```bash
C:\Users\YourName\Documents\Tube2FEM
```
In WSL, this path will be accessible as:
```bash
/mnt/c/Users/YourName/Documents/Tube2FEM
```

This allows MATLAB (Windows) and FEniCS (WSL) to access the same files without needing manual transfer.


Create a "Problem" file in the Tube2FEM directory and create the required subfolders (e.g., FiniteElement, Input, Mesh) and a main.m file for your example. Finally, use the modular structure of Tube2FEM (presented in the Case Study main.m files) to customize your problem.


---
