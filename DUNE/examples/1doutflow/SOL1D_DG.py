# Adapted from ET implementation; his original preamble included below
#
# Attempt at time-dependent solver of 1D SOL equations, with upwind flux
# based on:
# https://www.firedrakeproject.org/demos/DG_advection.py.html
# irksome implementation follows
# https://www.firedrakeproject.org/Irksome/demos/demo_cahnhilliard.py.html
# and
# https://www.firedrakeproject.org/Irksome/demos/demo_monodomain_FHN.py.html


import math
import ufl
from dune.grid import structuredGrid as leafGridView
from dune.fem import assemble
from dune.fem.space import dglegendre as dgSpace
from dune.ufl import Constant, DirichletBC


def check_dict_empty(d, key_desc, is_fatal=False):
    prefix = "ERROR: Invalid" if is_fatal else "WARNING: Ignoring"
    [print(f"***{prefix} {key_desc} '{k}'***\n") for k in d.keys()]
    if is_fatal:
        raise RuntimeError("check_dict_empty failed with is_fatal=True")


def run_1D_SOL(model_params={}, settings={}, **invalid_kws):
    check_dict_empty(invalid_kws, "keyword")

    # Set defaults for model params that weren't supplied
    nstar = Constant(model_params.pop("nstar", 1.0))
    Temp = Constant(model_params.pop("Temp", 1.0))
    check_dict_empty(model_params, "model parameter")
    # Fix Mach number of 1
    mach_num = Constant(1.0)

    # Set defaults for settings that weren't supplied
    # ic_settings = settings.pop("ICs",)
    meshres = settings.pop("meshres", 200)
    nsteps = settings.pop("nsteps", 250)
    tfinal = settings.pop("tfinal", 10.0)
    check_dict_empty(model_params, "setting")

    # Mesh
    gridView = leafGridView([-1], [1], [meshres], overlap=1)
    space = dgSpace(gridView, order=1)

    # Function spaces
    V1 = ufl.FunctionSpace(space, "DG", 1)
    # IMPORTANT: velocity space needs to be continuous
    V2 = ufl.VectorFunctionSpace(space, "CG", 1)
    V = V1 * V2

    # Set timestep
    t = Constant(0.0)
    dt = Constant(tfinal / nsteps)

    # parameters for irksome
    butcher_tableau = irk.GaussLegendre(2)

    (x,) = ufl.SpatialCoordinate(space)
    nu = ufl.Function(V)
    n, u = ufl.split(nu)
    v1, v2 = ufl.TestFunctions(V)

    # Density ICs - Gaussian blob
    mach_num = Constant(model_params.pop("mach_num", 1.0))
    width = 0.1
    nu.sub(0).interpolate(
        (
            1.0
            + 0.2
            * (1 / ufl.sqrt(2 * math.pi * width**2))
            * ufl.exp(-(x**2) / (2 * width**2))
        )
    )
    # Velocity ICs v=x
    nu.sub(1).interpolate(ufl.as_vector([0.0 + 1.0 * x]))

    # Density source - const nstar everywhere
    nstarFunc = ufl.Function(V)
    nstarFunc.sub(0).interpolate(nstar + 0.0 * x)

    zero = Constant(0.0)
    dbc = DirichletBC(space, zero, ufl.abs(x[0]) > 1 - 1e-8)

    # Construct a quantity that =u.n for faces with u.n +ve, 0 otherwise
    norm = ufl.FacetNormal(space)
    u_n = 0.5 * (ufl.dot(u, norm) + abs(ufl.dot(u, norm)))

    # outflow BC imposed weakly in here
    F = (
        -((dt * v1) * ufl.dx + (n * ufl.dot(irk.Dt(u), v2)) * ufl.dx)
        + (n * ufl.dot(u, ufl.grad(v1)) + v1 * nstar) * ufl.dx
        + (
            nstar * ufl.dot(u, v2)
            + n * u[0] * ufl.grad(ufl.dot(u, v2))[0]
            + n * u[0] * ufl.dot(u, ufl.grad(v2[0]))
            + Temp * n * ufl.grad(v2[0])[0]
        )
        * ufl.dx
        - (v1("+") - v1("-")) * (u_n("+") * n("+") - u_n("-") * n("-")) * ufl.dS
        + (u("+")[0] * v2("+")[0] - u("-")[0] * v2("-")[0])
        * (n("+") * u_n("+") - n("-") * u_n("-"))
        * ufl.dS
        - ufl.conditional(ufl.dot(u, norm) > 0, v1 * ufl.dot(u, norm) * n, 0.0) * ufl.ds
    )
    # Solver params taken from Cahn-Hilliard example cited above
    solver_params = {
        "snes_monitor": None,
        "snes_max_it": 100,
        "snes_linesearch_type": "l2",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "mat_type": "aij",
        "pc_factor_mat_solver_type": "mumps",
    }

    # Dirichlet BCs are needed for boundary velocity; only works for mach number = 1 for now
    u_bc_left = ufl.DirichletBC(
        V.sub(1), ufl.as_vector([-mach_num * ufl.sqrt(Temp)]), 1
    )
    u_bc_right = ufl.DirichletBC(
        V.sub(1), ufl.as_vector([mach_num * ufl.sqrt(Temp)]), 2
    )

    # stepper = irk.TimeStepper(
    #     F,
    #     butcher_tableau,
    #     t,
    #     dt,
    #     nu,
    #     solver_parameters=solver_params,
    #     bcs=[u_bc_left, u_bc_right],
    # )

    # # Set up output
    # outdir = os.path.dirname(__file__)
    # outfile = fd.output.VTKFile(os.path.join(outdir, "SOL_1D_DG_upwind.pvd"))

    # Calc solution vectors
    n_sln = nstar / ufl.sqrt(Temp) * (Constant(1) + ufl.sqrt(1 - x * x))
    u_sln = ufl.sqrt(Temp) * (Constant(1) - ufl.sqrt(1 - x * x)) / x

    # # Run
    # while float(t) < float(tfinal):
    #     if (float(t) + float(dt)) >= tfinal:
    #         dt.assign(tfinal - float(t))

    #     print(f"t = {t}")

    #     outfile.write(nu.sub(0), nu.sub(1))
    #     stepper.advance()

    #     t.assign(float(t) + float(dt))

    #     # Print fractional errors
    #     print(f"    n_err = {ufl.norm(n-n_sln)/ufl.norm(n_sln)}")
    #     print(f"    u_err = {ufl.norm(u[0]-u_sln)/ufl.norm(u_sln)}")

    print("done.\n")


run_1D_SOL()
