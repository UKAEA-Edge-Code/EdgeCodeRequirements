"""
Based on ../../../../Firedrake/examples/scripts/aniso_diffusion/aniso_diffusion_DeluzetNarski.py
"""

from dune.alugrid import aluSimplexGrid as leafGridView
from dune.fem import assemble, integrate
from dune.fem.function import gridFunction
from dune.fem.space import lagrange
from dune.grid import reader
from dune.ufl import Constant, DirichletBC
from scipy.sparse.linalg import spsolve as solver
import numpy as np
from ufl import (
    as_vector,
    conditional,
    cos,
    dx,
    div,
    dot,
    grad,
    inner,
    pi,
    sin,
    SpatialCoordinate,
    sqrt,
    TestFunction,
    TrialFunction,
)
import os.path

# ============================== MODEL PARAMETERS =============================
# parameter for varying B-field away from x-direction; 10 gives reasonable bent fieldlines, 0 is x-direction only
alpha = Constant(0.0)
# anisotropy parameter, 1.0 is isotropic, intended to be decreased to small positive values
eps = Constant(1.0e-1)
m = 0
# Note if changing mesh need to make sure indices in bcB and bcT are correct - see comment.
# =============================================================================

# Read .msh file from Firedrake eg
repo_root = os.path.abspath(os.path.dirname(__file__) + "../../../../")
msh_path = os.path.join(
    repo_root,
    "Firedrake/examples/scripts/aniso_diffusion/aniso_diffusion_DeluzetNarski.msh",
)
domain = (reader.gmsh, msh_path)
gridView = leafGridView(domain, dimgrid=2)

V = lagrange(gridView, order=1)

u = TrialFunction(V)
v = TestFunction(V)

x, y = SpatialCoordinate(V)

solution_expr = sin(pi * y + alpha * (y**2 - y) * cos(m * pi * x)) + eps * cos(
    2 * pi * x
) * sin(pi * y)

solution = gridFunction(V, solution_expr)

# Fix zero Dirichlet BC on top and bottom boundaries
bcs = DirichletBC(V, 0, abs(y - 0.5) + 1e-6 >= 0.5)
# optional Dirichlet BC for lam=0.0 case: removes null space problem
# bcLR = DirichletBC(
#     V, solution, (1, 2)
# )

# magnetic field lines direction unit vector
norm_denom = sqrt(
    (alpha * (2 * y - 1) * cos(m * pi * x) + pi) ** 2
    + (pi * alpha * m * (y**2 - y) * sin(m * pi * x)) ** 2
)
bhat = as_vector(
    [
        (alpha * (2 * y - 1) * cos(m * pi * x) + pi) / norm_denom,
        (pi * alpha * m * (y**2 - y) * sin(m * pi * x)) / norm_denom,
    ]
)

k_par = Constant(1.0)
k_per = eps  # THIS IS EPSILON IN ANISO DIFFUSION TENSOR

flux = k_par * bhat * dot(bhat, grad(u)) + k_per * (grad(u) - bhat * dot(bhat, grad(u)))

a = inner(flux, grad(v)) * dx

# Work out f numerically from MMS given the solution in "solution_expr" above
V_S = lagrange(gridView, order=10)
solution_S = V_S.interpolate(
    solution_expr,
    name="solution_s",
)
flux_S = k_par * bhat * dot(bhat, grad(solution_S)) + k_per * (
    grad(solution_S) - bhat * dot(bhat, grad(solution_S))
)
f = V.interpolate(-div(flux_S), name="f")

L = inner(f, v) * dx

mat, rhs = assemble([a == L, bcs])
T = V.interpolate(0, name="T")
T.as_numpy[:] = solver(mat.as_numpy, rhs.as_numpy)


# # output magnetic field lines
# VB = VectorFunctionSpace(mesh, "CG", 1)
# BField = Function(VB)
# BField.interpolate(bhat)
# File("aniso_diffusion_DeluzetNarski_BField.pvd").write(BField)

# Compute and print L2 error
T_exact = sin(pi * y + alpha * (y**2 - y) * cos(m * pi * x)) + eps * cos(
    2 * pi * x
) * sin(pi * y)
L2_error = np.sqrt(
    integrate(gridFunction(dot(T - T_exact, T - T_exact), name="l2error"))
)
print(f"L2 error norm {L2_error:.2e}")

# Plot the solution
T.plot()
Tdiff = V.interpolate(T - T_exact)
solution_S.plot()

# T.rename("numerical solution")
# solution.rename("analytic solution")

# File("aniso_diffusion_DeluzetNarski.pvd").write(T, solution)
