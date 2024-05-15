from dune.grid import structuredGrid as leafGridView
from dune.fem.space import dglegendre as dgSpace
from dune.fem.scheme import galerkin as solutionScheme
from dune.ufl import Constant, DirichletBC
from ufl import (
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    FacetNormal,
    dx,
    ds,
    exp,
    grad,
    grad,
    dot,
    conditional,
    as_vector,
    jump,
    dS,
)

# User params
Nx = 64
Ny = 16
xmin = 0
ymin = 0
xmax = 40
ymax = 10
adv_vel_x = 1.0
target_cfl = 0.05
tfinal = 1.0
dt_plot = tfinal / 20

# Derived
delta_x = (xmax - xmin) / Nx
centre_x = (xmax - xmin) / 2
centre_y = (ymax - ymin) / 2

# overlap=1 is needed for parallel computations
gridView = leafGridView([xmin, ymin], [xmax, ymax], [Nx, Ny], overlap=0)
is_parallel = gridView.comm.size > 1
space = dgSpace(gridView, order=2)

u = TrialFunction(space)
v = TestFunction(space)
norms = FacetNormal(space)
x, y = SpatialCoordinate(space)

# Set advection velocity and compute upwind flux
adv_vel = as_vector([adv_vel_x, 0])
hatb = (dot(adv_vel, norms) + abs(dot(adv_vel, norms))) / 2.0

# Mask that picks out x boundaries
x_bdy_mask = conditional((1 + x) * (1 - x) < 1e-10, 1, 0)
# Boundary values
x_bdy_vals = 0.0

advInternal = dot(-adv_vel * u, grad(v)) * dx
advInteriorFacets = jump(hatb * u) * jump(v) * dS
advExteriorFacets = (
    (hatb * u + (dot(adv_vel, norms) - hatb) * x_bdy_vals) * v * x_bdy_mask * ds
)
rhs = advInternal + advInteriorFacets + advExteriorFacets

# Gaussian ICs with width s
s = 2.0
ICs = exp(
    -((x - centre_x) * (x - centre_x) + (y - centre_y) * (y - centre_y)) / (s * s)
)

# Weak form of the time derivative
dt = Constant(0, name="dt")
sln_cur = space.interpolate(ICs, name="u_h")
sln_prev = space.interpolate(0.0, name="previous")
time_deriv = dot((u - sln_prev) / dt, v) * dx

scheme = solutionScheme(
    time_deriv == rhs,
    solver="gmres",
    parameters={
        "newton.linear.preconditioning.method": "jacobi",
        "newton.verbose": True,
    },
)

# Time integration loop
time = 0
dt.value = target_cfl * delta_x / adv_vel_x
print(f"Start time int loop. dt = {dt.value}")
t_plot_next = dt_plot
while time < (tfinal - 1e-6):
    sln_prev.assign(sln_cur)
    print(f"t = {time}")
    scheme.solve(target=sln_cur)
    if time >= t_plot_next:
        sln_cur.plot()
        t_plot_next += dt_plot
    time += dt.value

print("Finished")

# Output/plot final state
if is_parallel:
    gridView.writeVTK("u_h", pointdata=[sln_cur])
else:
    sln_cur.plot()
