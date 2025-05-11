"""
    Generates a 3D-1D mixed-dimension simulation example

    Parameters:
    -----------
    input_file : mesh file path 'tetraFileName'
        a tetrahedral mesh file (.xdmf)

    Returns:
    --------
    .vtk file  
        Simulation output in .vtk format

    """
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from dolfin import *


mesh = Mesh()
with XDMFFile(tetraFileName) as infile:
    infile.read(mesh)

mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile(triFileName) as infile:
    infile.read(mvc, "name_to_read")
surf_markers = cpp.mesh.MeshFunctionSizet(mesh, mvc)

V = FunctionSpace(mesh, 'P', 1)
u = Function(V)


# Read the content of the BC_values.txt
file_path = BCvalues

with open(file_path, 'r') as file:
    content = file.read()

# Split the content into lines
BC_values = content.split('\\n')
# Convert strings to float
BC_values = [float(value) for value in BC_values]
#print(BC_values)


# Read the content of the BC_points.txt
file_path = BCpoints

with open(file_path, 'r') as file:
    content = file.read()

# Split the content into lines
BC_points = content.split('\\n')
# Convert strings to tuple
BC_points = [tuple(map(float, s.split())) for s in BC_points]
#print(BC_points)


bcs = []
for i, (bc_pt, bc_val) in enumerate(zip(BC_points, BC_values)):
    bc_i = CompiledSubDomain(f'near(x[0], {bc_pt[0]}, tol) && near(x[1], {bc_pt[1]}, tol)',
                             tol=1E-10)
    bcs.append(DirichletBC(V, bc_val, bc_i, method="pointwise"))
    #print(f'{round(i / len(BC_points) * 100, 2):6.2f}%%')



# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution to file in VTK format
vtkfile = File(outputFileName)
vtkfile << u

