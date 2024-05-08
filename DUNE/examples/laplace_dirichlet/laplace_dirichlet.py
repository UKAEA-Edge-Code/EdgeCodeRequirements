import numpy as np

from dune.grid import structuredGrid

gridView = structuredGrid([0, 0], [1, 1], [20, 20])

from dune.fem import assemble, integrate
from dune.fem.space import lagrange
from dune.ufl import DirichletBC

from scipy.sparse.linalg import spsolve as solver

# Set pyplot color map
from matplotlib import pyplot as plt

plt.set_cmap("jet")


space = lagrange(gridView, order=1)
u_h = space.interpolate(0, name="u_h")

from ufl import (
    FacetNormal,
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    conditional,
    dx,
    grad,
    inner,
    dot,
    div,
    ds,
)

# Set up functions
x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)

exact = 1 / 2 * (x[0] ** 2 + x[1] ** 2) - 1 / 3 * (x[0] ** 3 - x[1] ** 3) + 1
a = dot(grad(u), grad(v)) * dx
f = -div(grad(exact))
g_N = grad(exact)
n = FacetNormal(space)
l = f * v * dx + dot(g_N, n) * conditional(x[0] >= 1e-8, 1, 0) * v * ds

# Set up BCs and assemble system matrix, rhs
dbc = DirichletBC(space, exact, x[0] <= 1e-8)
mat, rhs = assemble([a == l, dbc])

# Solve (using scipy, can also use petsc for linear solves) and calc errors
u_h.as_numpy[:] = solver(mat.as_numpy, rhs.as_numpy)
u_h.plot()
e_h = u_h - exact
squaredErrors = integrate([e_h**2, inner(grad(e_h), grad(e_h))])
print("L^2 and H^1 errors:", [np.sqrt(e) for e in squaredErrors])

# Re-solve for a different RHS (system matrix, mat, doesn't change)
l = conditional(dot(x, x) < 0.1, 1, 0) * v * dx
dbc = DirichletBC(space, x[1] * (1 - x[1]), x[0] <= 1e-8)
rhs = assemble([l, dbc])
u_h.as_numpy[:] = solver(mat.as_numpy, rhs.as_numpy)
u_h.plot()
