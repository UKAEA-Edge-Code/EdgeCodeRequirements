import numpy as np
from dune.fem import globalRefine, integrate
from dune.fem.function import gridFunction
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.space import lagrange as solutionSpace
from dune.grid import structuredGrid
from dune.ufl import Constant
from ufl import (
    FacetNormal,
    ds,
    dx,
    inner,
    dot,
    div,
    grad,
    sqrt,
    exp,
    SpatialCoordinate,
    TrialFunction,
    TestFunction,
)

gridView = structuredGrid([0, 0], [1, 1], [8, 8])
space = solutionSpace(gridView, order=2)
u_sln = space.interpolate(0, name="u_sln")
u_sln_prev = u_sln.copy(name="previous")

# Compute ICs
x = SpatialCoordinate(space)
u_init = 1 / 2 * (x[0] ** 2 + x[1] ** 2) - 1 / 3 * (x[0] ** 3 - x[1] ** 3) + 1

# Set up (time-dependent) a = a(u,v)
dt = Constant(0, name="dt")  # time step; 0 for now
t = Constant(0, name="t")  # current time
u = TrialFunction(space)
v = TestFunction(space)

# Helper functions to make LHS more concise
abs_du = lambda u: sqrt(inner(grad(u), grad(u)))
K = lambda u: 2 / (1 + sqrt(1 + 4 * abs_du(u)))

# LHS
lhs = (
    dot((u - u_sln_prev) / dt, v)
    + 0.5 * dot(K(u) * grad(u), grad(v))
    + 0.5 * dot(K(u_sln_prev) * grad(u_sln_prev), grad(v))
) * dx

# (Time-dependent) exact solution
exact = lambda t: exp(-2 * t) * (u_init - 1) + 1

# Helper functions to make LHS more concise
#   Time-dependent forcing term, chosen to correspond to the exact solution above
f = lambda t: -2 * exp(-2 * t) * (u_init - 1) - div(K(exact(t)) * grad(exact(t)))
#   Boundary flux
g = lambda s: K(exact(s)) * grad(exact(s))
n = FacetNormal(space)

# RHS
rhs = 0.5 * (f(t) + f(t + dt)) * v * dx + 0.5 * dot(g(t) + g(t + dt), n) * v * ds

# Set up non-linear solver (uses Newton method)
scheme = solutionScheme(lhs == rhs, solver="cg")

# # Default quadrature order is 2*space.order for elements, 2*space.order+1 for surface integrals; can increase it if necessary
# scheme.setQuadratureOrders( 2*space.order, 2*space.order+1 )

# Compute the exact solution at the end, then set up functions to work out the error as the simulation evolves
t_final = 0.25
u_exact_final = exact(t_final)
l2error = gridFunction(
    dot(u_sln - u_exact_final, u_sln - u_exact_final), name="l2error"
)
h1error = gridFunction(
    dot(grad(u_sln - u_exact_final), grad(u_sln - u_exact_final)), name="h1error"
)

# Set timestep
dt.value = 0.01
time = 0
# Set ICs
u_sln.interpolate(u_init)
# Time integration loop
while time < (t_final - 1e-6):
    t.value = time
    u_sln_prev.assign(u_sln)
    scheme.solve(target=u_sln)
    time += dt.value

# Report errors at the end
errors = [np.sqrt(e) for e in integrate([l2error, h1error])]
print(f"Result for grid size {gridView.size(0)}:")
print("\t | u_sln - u | =", "{:0.5e}".format(errors[0]))
print("\t | grad(uh - u) | =", "{:0.5e}".format(errors[1]))

# Plot final state
u_sln.plot()
# Write final state to VTU
gridView.writeVTK(
    "forchheimer", pointdata={"u": u_sln, "l2error": l2error, "h1error": h1error}
)


# Refine grid, reset ICs and repeat
globalRefine(1, gridView.hierarchicalGrid)
dt.value /= 2
u_sln.interpolate(u_init)
time = 0
while time < (t_final - 1e-6):
    t.value = time
    u_sln_prev.assign(u_sln)
    scheme.solve(target=u_sln)
    time += dt.value

# Report errors at the end
errorsFine = [np.sqrt(e) for e in integrate([l2error, h1error])]
print(f"Result for grid size {gridView.size(0)}:")
print("\t | u_sln - u | =", "{:0.5e}".format(errorsFine[0]))
print("\t | grad(uh - u) | =", "{:0.5e}".format(errorsFine[1]))

# Plot final state
u_sln.plot()

# Print experimental order of convergence
# Should be cubic for L2 error and quadratic for the H1 error
eocs = [
    round(np.log(fine / coarse) / np.log(0.5), 2)
    for fine, coarse in zip(errorsFine, errors)
]
print("EOCs:", eocs)
