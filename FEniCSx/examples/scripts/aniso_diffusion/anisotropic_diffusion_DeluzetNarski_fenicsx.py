# anisotropic_diffusion_DeluzetNarski_fenicsx.py
# port of aniso_diffusion_DeluzetNarski.py Firedrake script
# various bits taken from FEniCSx examples
# not written by an expert, just messed with it until it worked

# What this does: tests straight- and curved-fieldlines cases of anisotropic diffusion
# after "A two field iterated AP method for highly anisotropic elliptic equations"
# by Deluzet and Narski.
# This shows the issues explored in the paper:
# 1. Null space problem for eps -> 0.
# 2. "Locking" problem for order-1 (comes if basis functions have no cpt in the direction of the B-field): solution amplitude comes out much smaller than should be.
# NOTE modified to use actual MMS solution from the paper
# NOTE can choose Cartesian mesh or use gmsh .msh unstructured triangle mesh
# NOTE triangle mesh version tested to highest res used in the paper, h = 0.00078125 / 3.8M elements, fine on laptop

# maybe not all these imports are necessary - ?
from mpi4py import MPI
from petsc4py import PETSc

import dolfinx.fem as fem
from dolfinx_mpc import LinearProblem, MultiPointConstraint
import numpy as np
import scipy.sparse.linalg
from dolfinx import default_scalar_type
from dolfinx.common import Timer, TimingType, list_timings
from dolfinx import io
from dolfinx.io import XDMFFile, gmshio
from dolfinx import mesh
from dolfinx.mesh import create_unit_square, locate_entities_boundary
from dolfinx.fem.petsc import assemble_vector, assemble_matrix, create_vector
from ufl import *
import ufl

import dolfinx_mpc.utils
from dolfinx_mpc import LinearProblem, MultiPointConstraint

from math import *

# MODEL PARAMETERS
alpha = 0.0  # parameter for varying B-field away from x-direction; 0 is x-direction only
eps = 0.1  # anisotropy parameter, 1.0 is isotropic, intended to be decreased to small positive values
m = 0  # number of sinusoidal wiggles in magnetic fieldlines
# END OF MODEL PARAMETERS

# load gmsh msh and use
#mesh0, cell_markers, facet_markers = gmshio.read_from_msh("aniso_diffusion_DeluzetNarski.msh", MPI.COMM_WORLD, gdim=2)

# ... or use the following to create a quad mesh
NX = 400
NY = 400
mesh0 = create_unit_square(MPI.COMM_WORLD, NX, NY, mesh.CellType.quadrilateral)

V = fem.functionspace(mesh0, ("Lagrange", 1))
tol = 250 * np.finfo(default_scalar_type).resolution

def dirichletboundary(x):
   return np.logical_or(np.isclose(x[1], 0.0, atol=tol), np.isclose(x[1], 1.0, atol=tol))

# create Dirichlet boundary condition
facets = locate_entities_boundary(mesh0, 1, dirichletboundary)
topological_dofs = fem.locate_dofs_topological(V, 1, facets)
zero = np.array(0, dtype=default_scalar_type)
bc = fem.dirichletbc(zero, topological_dofs, V)
bcs = [bc]

u = TrialFunction(V)
v = TestFunction(V)

x = SpatialCoordinate(mesh0)

## magnetic field lines direction unit vector
norm_denom = ufl.sqrt((alpha*(2*x[1]-1)*(ufl.cos(m*pi*x[0]))+pi)**2+(pi*alpha*m*(x[1]**2-x[1])*(ufl.sin(m*pi*x[0])))**2)
bhat = as_vector([(alpha*(2*x[1]-1)*ufl.cos(m*pi*x[0])+pi)/norm_denom, (pi*alpha*m*(x[1]**2-x[1])*ufl.sin(m*pi*x[0]))/norm_denom])

k_par = 1.0
k_per = eps  # THIS IS EPSILON IN ANISO DIFFUSION TENSOR

flux = k_par * bhat * dot(bhat, grad(u)) + k_per * (grad(u) - bhat * dot(bhat, grad(u)))

a = inner(flux, grad(v))*dx
f = fem.Function(V)

## got to work out f from MMS given the solution in "solution" above
## algebra is a total PITA so let's use a high-order function space to do it numerically ...
## can easily be tested for straight fieldlines case (MMS is easy algebra in that case)
V_S = fem.functionspace(mesh0, ("Lagrange", 10))
solution_S = fem.Function(V_S)
solution_S.interpolate(lambda x: np.sin(pi*x[1]+alpha*(x[1]**2-x[1])*np.cos(m*pi*x[0]))+eps*np.cos(2*pi*x[0])*np.sin(pi*x[1]))
mms_expr = fem.Expression(-div(k_par * bhat * dot(bhat, grad(solution_S)) + k_per * (grad(solution_S) - bhat * dot(bhat, grad(solution_S)))), V.element.interpolation_points())
f.interpolate(mms_expr)

L = inner(f,v)*dx

bilinear_form = fem.form(a)
A = assemble_matrix(bilinear_form, bcs=bcs)
A.assemble()
linear_form = fem.form(L)
b = create_vector(linear_form)

mpc = MultiPointConstraint(V)
mpc.finalize()
problem = LinearProblem(a, L, mpc, bcs=bcs)

solver = problem.solver

solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)

T = fem.Function(V)
T.name = "density"

nh = problem.solve()

# output has to be first order so need ot interpolate (problem is output will look crap on a low res high order calc even though it's correct ...)
# really need to interpolate onto a finer mesh in that case, but I got segfaults when I tried ...
V1 = fem.functionspace(mesh0, ("Lagrange",1))
n_out = fem.Function(V1)
n_out.interpolate(nh)

xdmf = io.XDMFFile(mesh0.comm, "anisotropic_diffusion_DeluzetNarski_fenicsx.xdmf", "w")
xdmf.write_mesh(mesh0)
xdmf.write_function(n_out)
