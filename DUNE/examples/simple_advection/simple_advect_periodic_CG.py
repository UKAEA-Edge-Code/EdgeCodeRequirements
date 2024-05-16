from dune.common import comm
from dune.fem import integrate
from dune.fem.function import gridFunction
from dune.fem.scheme import galerkin
from dune.fem.space import lagrange
from dune.grid import structuredGrid
from dune.grid import reader
from dune.alugrid import aluSimplexGrid
from functools import partial
import numpy as np
import os.path
from dune.ufl import Constant, DirichletBC
import time
from ufl import (
    as_vector,
    conditional,
    div,
    dot,
    dx,
    exp,
    le,
    SpatialCoordinate,
    TestFunction,
    TrialFunction,
)

# Store path of the directory containing this script for use in IO functions
this_dir = os.path.abspath(os.path.dirname(__file__))

# override print to work in parallel
print = partial(print, flush=True) if comm.rank == 0 else lambda *args, **kwargs: None

# ================================= Options ===================================
meshres_x = 64
meshres_y = 16
xmin = 0.0
ymin = 0.0
xmax = 40.0
ymax = 10.0
xperiodic = True
mesh_from_file = False
ICs_type = "gaussian"
# ICs_type = "frankenstein"
t_final = 40.0
output_dt = 1.0

target_cfl = 0.01
delta_x = (xmax - xmin) / meshres_x
adv_vel_x = 1.0

solver_type = "gmres"
solver_settings = {
    "newton.tolerance": 1e-10,
    "newton.verbose": False,
    "newton.linear.tolerance": 1e-12,
    "newton.linear.preconditioning.method": "jacobi",
    "newton.linear.verbose": False,
}
# =============================================================================


def main():
    start = time.time()

    if mesh_from_file:
        # Read periodic mesh from file
        domain = (
            reader.dgf,
            os.path.join(os.path.abspath(os.path.dirname(__file__)), "rectangle.dgf"),
        )
        gridView = aluSimplexGrid(domain, dimgrid=2)
    else:
        # Generate Cartesian structured periodic mesh
        gridView = structuredGrid(
            [xmin, ymin],
            [xmax, ymax],
            [meshres_x, meshres_y],
            periodic=[xperiodic, False],
        )

    # CG function space
    V = lagrange(gridView, order=3)

    x, y = SpatialCoordinate(V)

    # Fix zero Dirichlet BC on top and bottom boundaries
    bcs = DirichletBC(V, 0, abs(y - ymax / 2) + 1e-6 >= ymax / 2)

    t = Constant(0.0, name="t")
    # Set dt
    dt = Constant(0.0, name="dt")
    dt.value = target_cfl * delta_x / adv_vel_x
    print(f"dt = {dt.value}")

    x, y = SpatialCoordinate(V)
    n_tri = TrialFunction(V)
    n_test = TestFunction(V)
    n = V.interpolate(0, name="n")
    n_prev = V.interpolate(0, name="n_prev")

    # Fixed drift velocity: 1 space unit / time unit in the x direction
    driftvel = as_vector([adv_vel_x, 0.0])

    # Set ICs and compute corresponding solution
    centre_x = xmax / 2
    centre_y = ymax / 2

    if ICs_type == "gaussian":
        # Gaussian ICs centered on domain midpoint
        s = 2.0  # width
        n.interpolate(
            exp(
                -((x - centre_x) * (x - centre_x) + (y - centre_y) * (y - centre_y))
                / (s * s)
            )
        )
        n_exact = exp(
            -((x - centre_x) * (x - centre_x) + (y - centre_y) * (y - centre_y))
            / (s * s)
        )
    elif ICs_type == "frankenstein":
        # 'Frankenstein's monster' ICs centered on domain midpoint
        s = 2.0
        a = 1.0
        n.interpolate(
            conditional(
                le(
                    (x - centre_x) * (x - centre_x) + (y - centre_y) * (y - centre_y),
                    s * s,
                ),
                exp(-a * a / (s * s - (x - centre_x) ** 2 - (y - centre_y) ** 2)),
                0,
            )
        )
        n_exact = conditional(
            le(
                (x - centre_x) * (x - centre_x) + (y - centre_y) * (y - centre_y), s * s
            ),
            exp(-a * a / (s * s - (x - centre_x) ** 2 - (y - centre_y) ** 2)),
            0,
        )
    else:
        raise ValueError(f"{ICs_type} isn't a valid ICs type")

    n_prev = n.copy()

    lhs = (n_tri - n_prev) / dt * n_test * dx
    rhs = n_test * div(n_tri * driftvel) * dx

    scheme = galerkin([lhs == rhs, bcs], solver=solver_type, parameters=solver_settings)

    t_next_output = output_dt
    output_num = 1
    while t.value <= t_final:
        n_prev.assign(n)
        scheme.solve(target=n)
        if t.value >= t_next_output:
            print(f"t = {t.value:.1f}")
            write_vtk(gridView, n, output_num)
            t_next_output += output_dt
            output_num += 1
        t.value += dt.value
    # Final output
    write_vtk(gridView, n, output_num)

    end = time.time()
    print(f"Finished in {end-start:.3f} s")

    L2_error = np.sqrt(
        integrate(gridFunction(dot(n - n_exact, n - n_exact), name="l2error"))
    )
    print(f"L2 error norm {L2_error:.2e}")


def write_vtk(gridView, n, output_num, dir=this_dir, prefix=ICs_type, verbose=True):
    fpath = os.path.join(dir, "output", f"{prefix}_{output_num}")
    gridView.writeVTK(fpath, pointdata={"n": n})
    if verbose:
        print(f"Wrote {fpath}.vtu")


main()
