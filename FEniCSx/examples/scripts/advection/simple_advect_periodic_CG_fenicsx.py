# simple_advect_CG_fenicsx.py
# attempt at time-dependent solver of 2D advection
# based upon:
# https://jsdokken.com/dolfinx-tutorial/chapter2/diffusion_code.html
# generates h5 / xdmf files containing the history of the advection

# periodic BC adapted from:
# https://github.com/jorgensd/dolfinx_mpc/blob/main/python/demos/demo_periodic_geometrical.py

# TODO:
# higher-order time-stepper
# DG implementation# add alternate function to advect, only Gaussian implemented so far

import ufl
import numpy as np

from petsc4py import PETSc

from dolfinx import fem, mesh, io, plot, default_scalar_type
from mpi4py import MPI
from dolfinx.fem.petsc import assemble_vector, assemble_matrix, create_vector, apply_lifting, set_bc

import dolfinx_mpc.utils
from dolfinx_mpc import LinearProblem, MultiPointConstraint

from dolfinx.mesh import locate_entities_boundary


T = 40.0
timeres = 400
t = 0.0
dt = T/timeres

meshres = 64
nx = meshres
ny = int(meshres/4)
domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([0,0]), np.array([40,10])],[nx, ny], mesh.CellType.quadrilateral)
V = fem.functionspace(domain, ("Lagrange", 1))

# it's necessary to interpolate to lower order to output, see https://fenicsproject.discourse.group/t/interpolating-functions-on-higher-order-elements/12658
V_output = fem.functionspace(domain, ("Lagrange",1))
n_output = fem.Function(V_output)

n = fem.Function(V)
n0 = fem.Function(V)
n.name = "density"
n0.name = "density_existing"  # NOTE calling it "density existing" works for order-1 but not order-2 (!)

driftvel = ufl.as_vector([-1.0, 0.0])

# Gaussian initial data
s = 2.0  # width parameter
n.interpolate(lambda x: np.exp(-((x[0]-20.0)**2+(x[1]-5.0)**2)/(s*s)))
n0.interpolate(lambda x: np.exp(-((x[0]-20.0)**2+(x[1]-5.0)**2)/(s*s)))

# Create boundary condition
tol = 250 * np.finfo(default_scalar_type).resolution

def dirichletboundary(x):
    return np.logical_or(np.isclose(x[1], 0, atol=tol), np.isclose(x[1], 10.0, atol=tol))

# Create Dirichlet boundary condition
facets = locate_entities_boundary(domain, 1, dirichletboundary)
topological_dofs = fem.locate_dofs_topological(V, 1, facets)
zero = np.array(0, dtype=default_scalar_type)
bc = fem.dirichletbc(zero, topological_dofs, V)
bcs = [bc]


# periodic boundary geometrical constraint
def periodic_boundary(x):
    return np.isclose(x[0], 0.0, atol=tol)


def periodic_relation(x):
    out_x = np.zeros_like(x)
    out_x[0] = 40.0 + x[0]
    out_x[1] = x[1]
    out_x[2] = x[2]
    return out_x

mpc = MultiPointConstraint(V)
mpc.create_periodic_constraint_geometrical(V, periodic_boundary, periodic_relation, bcs)
mpc.finalize()


# output
xdmf = io.XDMFFile(domain.comm, "simple_advect_periodic_CG_fenicsx.xdmf", "w")
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
problem = LinearProblem(a, L, mpc, bcs = bcs)

solver = problem.solver
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
    apply_lifting(b, [bilinear_form], [[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, [bc])

    # Solve linear problem
    n = problem.solve()
    n.x.scatter_forward()

    # Update solution at previous time step (n)
    n0.x.array[:] = n.x.array

    # Write solution to file
    n_output.interpolate(n)
    xdmf.write_function(n_output, t)
xdmf.close()
