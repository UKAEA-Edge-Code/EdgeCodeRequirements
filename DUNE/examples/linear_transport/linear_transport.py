# import numpy, math, sys
from matplotlib import pyplot
from dune.grid import structuredGrid as leafGridView
from dune.fem.space import dglegendre as dgSpace

# from dune.fem.scheme import galerkin as solutionScheme
# from dune.ufl import Constant

from ufl import (
    FacetNormal,
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    triangle,
    FacetNormal,
    dx,
    ds,
    grad,
    grad,
    dot,
    conditional,
    as_vector,
    jump,
    dS,
)

# remove later, shows same error: from dune.fem.space import dgonb as dgSpace
N = 200
gridView = leafGridView([0], [1], [N], overlap=1)
space = dgSpace(gridView, order=1)

u = TrialFunction(space)
v = TestFunction(space)
n = FacetNormal(space)
x = SpatialCoordinate(space)

# transport direction and upwind flux
speed = 0.5
b = as_vector([speed])
hatb = (dot(b, n) + abs(dot(b, n))) / 2.0
# boundary values (for left boundary)
dD = conditional(x[0] < 0.01, 1, 0)
g = 1

# forms
advInternal = dot(-b * u, grad(v)) * dx
advSkeleton = jump(hatb * u) * jump(v) * dS
advBnd = (hatb * u + (dot(b, n) - hatb) * g) * v * dD * ds
form = advInternal + advSkeleton + advBnd


from dune.fem.operator import molGalerkin as molOperator

op = molOperator(form)

uh = space.interpolate(conditional(x[0] < 0.25, 1, 0), name="solution")
uh.plot()
w = uh.copy()
un = uh.copy()
cfl = 0.1
tau = cfl / (speed * N)

t = 0
while t < 1:
    un.assign(uh)

    op(uh, w)
    uh.axpy(
        -tau, w
    )  # with numpy backend equivalent to 'uh.as_numpy[:] -= tau*w.as_numpy[:]'
    op(uh, w)
    uh.axpy(
        -tau, w
    )  # with numpy backend equivalent to 'uh.as_numpy[:] -= tau*w.as_numpy[:]'

    # with numpy backend the following is equivalent to
    # 'uh.as_numpy[:] = 0.5*(uh.as_numpy[:] + un.as_numpy[:])'
    uh *= 0.5
    uh.axpy(0.5, un)
    t += tau
uh.plot()
