from dune.fem import assemble, integrate, lagrange
from dune.grid import structuredGrid

import numpy as np
from scipy.sparse.linalg import spsolve as solver

from ufl import (
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    dx,
    grad,
    inner,
    cos,
    pi,
)

# Grid
gridView = structuredGrid([0, 0], [1, 1], [20, 20])

# Function space
space = lagrange(gridView, order=1)
u_h = space.interpolate(0, name="u_h")

# Define the problem
x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)

f = (8 * pi**2 + 1) * cos(2 * pi * x[0]) * cos(2 * pi * x[1])
a = (inner(grad(u), grad(v)) + u * v) * dx
l = f * v * dx

# Assemble lhs, rhs
mat, rhs = assemble(a == l)

# Solve
A = mat.as_numpy
b = rhs.as_numpy
y = u_h.as_numpy
y[:] = solver(A, b)

# Plot the solution
u_h.plot()

# Compute and print errors
exact = cos(2 * pi * x[0]) * cos(2 * pi * x[1])
e_h = u_h - exact
squaredErrors = integrate([e_h**2, inner(grad(e_h), grad(e_h))])
print("L^2 and H^1 errors:", [np.sqrt(e) for e in squaredErrors])

print("average:", integrate(exact, gridView=gridView))
# since the integrand is scalar, the following is equivalent:
print("average:", assemble(exact * dx, gridView=gridView))
