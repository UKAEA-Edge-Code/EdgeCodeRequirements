# simple_advect_CG_fenicsx.py
# attempt at time-dependent solver of 2D advection
# based upon:
# https://jsdokken.com/dolfinx-tutorial/chapter2/diffusion_code.html
# generates h5 / xdmf files containing the history of the advection

# TODO:
# higher-order time-stepper
# periodic BC (and L2 error evaluation at end of calculation)
# DG implementation# add alternate function to advect, only Gaussian implemented so far

import ufl
import numpy as np

from petsc4py import PETSc

from dolfinx import fem, mesh, io, plot
from mpi4py import MPI
from dolfinx.fem.petsc import assemble_vector, assemble_matrix, create_vector, apply_lifting, set_bc

T = 40.0
timeres = 400
t = 0.0
dt = T/timeres

meshres = 64
nx = meshres
ny = int(meshres/4)
domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([0,0]), np.array([40,10])],[nx, ny], mesh.CellType.quadrilateral)
V = fem.functionspace(domain, ("Lagrange", 3))

# it's necessary to interpolate to lower order to output, see https://fenicsproject.discourse.group/t/interpolating-functions-on-higher-order-elements/12658
V_output = fem.functionspace(domain, ("Lagrange",1))
n_output = fem.Function(V_output)

n = fem.Function(V)
n0 = fem.Function(V)
n.name = "density"
n0.name = "density_existing"  # NOTE calling it "density existing" works for order-1 but not order-2 (!)

driftvel = ufl.as_vector([1.0, 0.0])

# Gaussian initial data
s = 2.0  # width parameter
#n.interpolate(lambda x: np.exp(-x[0]))
n.interpolate(lambda x: np.exp(-((x[0]-20.0)**2+(x[1]-5.0)**2)/(s*s)))
n0.interpolate(lambda x: np.exp(-((x[0]-20.0)**2+(x[1]-5.0)**2)/(s*s)))

# Create boundary condition
fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(
    domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool))
bc = fem.dirichletbc(PETSc.ScalarType(0), fem.locate_dofs_topological(V, fdim, boundary_facets), V)

# output
xdmf = io.XDMFFile(domain.comm, "simple_advect_CG_fenicsx.xdmf", "w")
xdmf.write_mesh(domain)
n_output.interpolate(n)
xdmf.write_function(n_output, t)

# variational problem and solver
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
a = u*v*ufl.dx - dt*ufl.dot(driftvel, ufl.grad(u))*v*ufl.dx
L = n0*v*ufl.dx
bilinear_form  = fem.form(a)
linear_form = fem.form(L)

A = assemble_matrix(bilinear_form)
A.assemble()
b = create_vector(linear_form)

# form problem
solver = PETSc.KSP().create(domain.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)

# time-loop + output
for i in range(timeres):
    t += dt

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    assemble_vector(b, linear_form)

    # Apply Dirichlet boundary condition to the vector
    #apply_lifting(b, [bilinear_form], [[bc]])
    #b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    #set_bc(b, [bc])

    # Solve linear problem
    solver.solve(b, n.vector)
    n.x.scatter_forward()

    # Update solution at previous time step (n)
    n0.x.array[:] = n.x.array

    # Write solution to file
    n_output.interpolate(n)
    xdmf.write_function(n_output, t)
xdmf.close()
