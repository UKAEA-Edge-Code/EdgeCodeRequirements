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
from matplotlib import pyplot as plt
import numpy as np
from ufl import (
    as_vector,
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

# Store path of the directory containing this script for use in IO functions
this_dir = os.path.abspath(os.path.dirname(__file__))


def aniso_diff(params={}):
    # parameter for varying B-field away from x-direction; 10 gives reasonable bent fieldlines, 0 is x-direction only
    alpha = Constant(params.pop("alpha", 0.0))
    # anisotropy parameter, 1.0 is isotropic, intended to be decreased to small positive values
    eps = Constant(params.pop("eps", 1.0e-15))
    m = Constant(params.pop("m", 0))

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

    # Fix zero Dirichlet BC on top and bottom boundaries
    bcs = DirichletBC(V, 0, abs(y - 0.5) + 1e-6 >= 0.5)
    # optional Dirichlet BC for lam=0.0 case: removes null space problem
    # bcLR = DirichletBC(
    #     V, solution_expr, abs(x - 0.5) + 1e-6 >= 0.5
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

    flux = k_par * bhat * dot(bhat, grad(u)) + k_per * (
        grad(u) - bhat * dot(bhat, grad(u))
    )

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

    # Compute and print L2 error
    T_exact = sin(pi * y + alpha * (y**2 - y) * cos(m * pi * x)) + eps * cos(
        2 * pi * x
    ) * sin(pi * y)
    L2_error = np.sqrt(
        integrate(gridFunction(dot(T - T_exact, T - T_exact), name="l2error"))
    )
    print(f"L2 error norm {L2_error:.2e}")

    # Plot the solution and save to png
    label = f"a{int(alpha.value):d}_m{int(m.value):d}"
    fig = plt.figure()
    T.plot(fig)
    fig.savefig(os.path.join(this_dir, f"output/aniso_diff{label}.png"))


# Run 3 cases
aniso_diff(params=dict(alpha=0.0))
aniso_diff(params=dict(alpha=2.0, m=1.0))
aniso_diff(params=dict(alpha=2.0, m=10.0))
