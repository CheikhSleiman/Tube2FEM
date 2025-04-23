from fenics import *
from mshr import *
import numpy as np
from dolfin import *

T = 1.0            # final time
num_steps = 100   # number of time steps
dt = T / num_steps # time step size
mu = 1.7984e-5         # dynamic viscosity
rho = 1.225       # density

# Read mesh

mesh = Mesh()
with XDMFFile(tetraFileName) as infile:
    infile.read(mesh)

mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile(triFileName) as infile:
    infile.read(mvc, "name_to_read")
surf_markers = dolfin.cpp.mesh.MeshFunctionSizet(mesh, mvc)

# dx = Measure('dx', domain = mesh, subdomain_data= vol_markers)
ds = Measure('ds', domain = mesh, subdomain_data = surf_markers)
n = FacetNormal(mesh)


# Define function spaces
V = VectorFunctionSpace(mesh, 'CG', 2, dim=3)
Q = FunctionSpace(mesh, 'CG', 1)






