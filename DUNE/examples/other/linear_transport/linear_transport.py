from dune.grid import structuredGrid as leafGridView
from dune.fem.operator import molGalerkin as molOperator
from dune.fem.space import dglegendre as dgSpace
from ufl import (
    FacetNormal,
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    triangle,
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
gridView = leafGridView([0], [1], [N], periodic=[1], overlap=1)
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


op = molOperator(form)

# Set ICs - step function from 0 <= x <= 0.25
uh = space.interpolate(conditional(x[0] < 0.25, 1, 0), name="solution")

# # Plot ICs
# uh.plot()

w = uh.copy()
un = uh.copy()
# Set timestep using target cfl val
cfl = 0.1
dt = cfl / (speed * N)
t_final = 3.0

nsteps = int(t_final / dt)
nchks = int(5 * t_final)
chk_freq = nsteps / nchks

t = 0
step = 0
while t < t_final:
    un.assign(uh)

    op(uh, w)
    uh.axpy(
        -dt, w
    )  # with numpy backend equivalent to 'uh.as_numpy[:] -= tau*w.as_numpy[:]'
    op(uh, w)
    uh.axpy(
        -dt, w
    )  # with numpy backend equivalent to 'uh.as_numpy[:] -= tau*w.as_numpy[:]'

    # with numpy backend the following is equivalent to
    # 'uh.as_numpy[:] = 0.5*(uh.as_numpy[:] + un.as_numpy[:])'
    uh *= 0.5
    uh.axpy(0.5, un)
    if step % chk_freq == 0:
        print(f"t = {t}")
        uh.plot()

    t += dt
    step += 1
