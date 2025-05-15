tetraFileName="../../../../../../../../../../../../../../../../../mnt/c/Users/homeuser/Documents/GitHub/Tube2FEM/tests/Mesh/Tetra.xdmf"
triFileName="../../../../../../../../../../../../../../../../../mnt/c/Users/homeuser/Documents/GitHub/Tube2FEM/tests/Mesh/Tri.xdmf"
fileName='../../../../../../../../../../../../../../../../../mnt/c/Users/homeuser/Documents/GitHub/Tube2FEM/tests/simOutput_3D/velocity.pvd'
fileName1='../../../../../../../../../../../../../../../../../mnt/c/Users/homeuser/Documents/GitHub/Tube2FEM/tests/simOutput_3D/pressure.pvd'
fileName2='../../../../../../../../../../../../../../../../../mnt/c/Users/homeuser/Documents/GitHub/Tube2FEM/tests/simOutput_3D/wss.pvd'
timeSeriesFileName='../../../../../../../../../../../../../../../../../mnt/c/Users/homeuser/Documents/GitHub/Tube2FEM/tests/simOutput_3D/velocity_series'

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







u_in = Expression('-0.025*t', degree = 1, t=0)

p_out = Expression('0*t', degree = 1, t=0)

bcp_outflow2 = DirichletBC(Q, p_out, surf_markers, 2)

bcu_walls     = DirichletBC(V, Constant((0, 0, 0)),surf_markers, 1)
bcu_inflow1    = DirichletBC(V.sub(2), u_in , surf_markers, 3)
bcu = [bcu_walls,bcu_inflow1]
bcp =[bcp_outflow2]
# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx  + inner(sigma(U, p_n), epsilon(v))*dx + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Time-stepping
t = 0
count = 0



#fileName='output/velocity.pvd'
file = File(fileName)
#fileName1 = 'output/vpressure.pvd'
file1 = File(fileName1)
#fileName2 = 'output/wss.pvd'
file2 = File(fileName2)

# Create time series (for use in reaction_system.py)
timeseries_u = TimeSeries(timeSeriesFileName)

n_normal = FacetNormal(mesh)
for n in range(num_steps):

    # Update current time
    t += dt
    u_in.t = t
    p_out.t= t
    print(t)
    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, 'cg', 'sor')
    # Save nodal values to file
    timeseries_u.store(u_.vector(), t)
    
    
    V_boundary = VectorFunctionSpace(mesh, 'CG', 1)
    u_wss = TrialFunction(V_boundary)
    v_wss = TestFunction(V_boundary)
    wss = Function(V_boundary,name="wss")
    T1 = -mu * rho * dot((grad(u_) + grad(u_).T), n_normal)
    Tn = dot(T1, n_normal)
    Tt = T1 - Tn*n_normal
    LHS_wss = assemble(inner(u_wss, v_wss) * ds, keep_diagonal=True)
    RHS_wss = assemble(inner(Tt, v_wss) * ds)
    LHS_wss.ident_zeros()
    solve(LHS_wss, wss.vector(), RHS_wss,'gmres', 'default')


    



    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)
    if (n%10==0):
       file << (u_,t)
       file1 << (p_,t)
       file2 << (wss,t)
       
       
    
    

 

