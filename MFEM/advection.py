"""Adapted from MFEM example 9

https://github.com/mfem/PyMFEM/blob/master/examples/ex9.py
"""
import mfem.ser as mfem
import numpy as np
from functools import cache


# Inflow boundary condition (zero for the problems considered in this example)
class inflow_coeff(mfem.PyCoefficient):
    def EvalValue(self, x):
        return 0


class AdvectionEvolution(mfem.PyTimeDependentOperator):
    def __init__(
        self, element_space, velocity, inflow
    ):  # mass_matrix, advection_term, boundary_term):
        mass_matrix = mfem.BilinearForm(element_space)
        mass_matrix.AddDomainIntegrator(mfem.MassIntegrator())
        advection_term = mfem.BilinearForm(element_space)
        # Confusingly named, as this operator is for ADvection, rather than CONvection
        advection_term.AddDomainIntegrator(mfem.ConvectionIntegrator(velocity, -1.0))
        advection_term.AddInteriorFaceIntegrator(
            mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5))
        )
        # Don't think this is actually needed here, as the boundaries
        # are all periodic anyway
        advection_term.AddBdrFaceIntegrator(
            mfem.TransposeIntegrator(mfem.DGTraceIntegrator(velocity, 1.0, -0.5))
        )

        boundary_term = mfem.LinearForm(element_space)
        boundary_term.AddBdrFaceIntegrator(
            mfem.BoundaryFlowIntegrator(inflow, velocity, -1.0, -0.5)
        )

        mass_matrix.Assemble()
        mass_matrix.Finalize()
        skip_zeros = 0
        advection_term.Assemble(skip_zeros)
        advection_term.Finalize(skip_zeros)
        boundary_term.Assemble()

        # Even those these aren't used, we need to keep a reference to
        # them. The SparseMatrix objects are just references to
        # components of the BilinearTerm objects, so if the latter
        # gets garbage collected then using the former will cause
        # segfaults.
        self._advection_bilinear_form = advection_term
        self._mass_bilinear_form = mass_matrix

        self.element_space = element_space
        self.advection_term = advection_term.SpMat()
        self.mass_matrix = mass_matrix.SpMat()
        self.boundary_term = boundary_term
        self.tmp = mfem.Vector(mass_matrix.Size())
        self.mass_preconditioner = mfem.DSmoother()
        self.mass_solver = mfem.CGSolver()
        self.mass_solver.SetPreconditioner(self.mass_preconditioner)
        self.mass_solver.SetOperator(self.mass_matrix)
        self.mass_solver.iterative_mode = False
        self.mass_solver.SetRelTol(1e-9)
        self.mass_solver.SetAbsTol(0.0)
        self.mass_solver.SetMaxIter(100)
        self.mass_solver.SetPrintLevel(0)
        self.implicit_preconditioner = mfem.BlockILU(element_space.GetFE(0).GetDof())
        mfem.PyTimeDependentOperator.__init__(self, self.mass_matrix.Size())

    def Mult(self, x, y):
        # Return the answer in y
        self.advection_term.Mult(x, self.tmp)
        self.tmp += self.boundary_term
        self.mass_solver.Mult(self.tmp, y)

    @cache
    def implicit_mass_solver(self, dt):
        matrix = mfem.SparseMatrix(self.advection_term)
        matrix *= -dt
        matrix += self.mass_matrix
        solver = mfem.GMRESSolver()
        solver.iterative_mode = False
        solver.SetRelTol(1e-9)
        solver.SetAbsTol(0.0)
        solver.SetMaxIter(100)
        solver.SetPrintLevel(0)
        solver.SetPreconditioner(self.implicit_preconditioner)
        solver.SetOperator(matrix)
        return solver

    def ImplicitSolve(self, dt, x, k):
        solver = self.implicit_mass_solver(dt)
        self.advection_term.Mult(x, self.tmp)
        self.tmp += self.boundary_term
        solver.Mult(self.tmp, k)


class GaussianInitial(mfem.PyCoefficient):
    def __init__(self, xcentre, ycentre, sigma, magnitude=1.0):
        super().__init__()
        self.xcentre = xcentre
        self.ycentre = ycentre
        self.sigma_sq = sigma**2
        self.magnitude = magnitude

    def EvalValue(self, x):
        dist_sq = (x[0] - self.xcentre) ** 2 + (x[1] - self.ycentre) ** 2
        return self.magnitude * np.exp(-dist_sq / self.sigma_sq)


class FrankensteinInitial(mfem.PyCoefficient):
    def __init__(self, xcentre, ycentre, radius, magnitude=1.0, a=1.0):
        super().__init__()
        self.asq = a * a
        self.xcentre = xcentre
        self.ycentre = ycentre
        self.rad_sq = radius**2
        self.magnitude = magnitude * np.exp(self.asq / self.rad_sq)

    def EvalValue(self, x):
        dist_sq = (x[0] - self.xcentre) ** 2 + (x[1] - self.ycentre) ** 2
        return (
            np.exp(-self.asq / self.rad_sq - dist_sq) if self.rad_sq > dist_sq else 0.0
        )


class AnalyticSolution(mfem.PyCoefficient):
    def __init__(self, velocity, initial, sx, sy):
        super().__init__()
        self.velocity = velocity.GetVec()
        self.initial = initial
        self.sx = sx
        self.sy = sy
        self.t = 0

    def set_t(self, t):
        self.t = t

    def EvalValue(self, x):
        xp = (x[0] - self.t * self.velocity[0]) % self.sx
        yp = (x[1] - self.t * self.velocity[1]) % self.sy
        return self.initial.EvalValue(mfem.Vector([xp, yp]))


def run(
    order=3,
    t_final=40,
    nx=64,
    ny=16,
    sx=40.0,
    sy=10.0,
    dt=0.05,
    velocity=1.0,
    angle=0.0,
    radius=0.2,
    magnitude=1.0,
    frankenstein_param=False,
    vis_steps=5,
    visualization=False,
    device="cpu",
    paraview=False,
):
    device = mfem.Device(device)
    device.Print()

    # 2. Read the mesh from the given mesh file. We can handle geometrically
    #    periodic meshes in this code.

    nonperiodic_mesh = mfem.Mesh(nx, ny, "QUADRILATERAL", False, sx, sy)
    mesh = mfem.Mesh.MakePeriodic(
        nonperiodic_mesh,
        nonperiodic_mesh.CreatePeriodicVertexMapping(
            (mfem.Vector([sx, 0.0]), mfem.Vector([0.0, sy]))
        ),
    )
    dim = mesh.Dimension()

    # 3. Define the ODE solver used for time integration.
    ode_solver = mfem.SDIRK23Solver(2)

    # 5. Define the discontinuous DG finite element space of the given
    #    polynomial order on the refined mesh.
    fec = mfem.DG_FECollection(order, dim, mfem.BasisType.GaussLobatto)
    fes = mfem.FiniteElementSpace(mesh, fec)

    print("Number of unknowns: " + str(fes.GetVSize()))

    # 6. Set up and assemble the bilinear and linear forms corresponding to the
    #    DG discretization. The DGTraceIntegrator involves integrals over mesh
    #    interior faces.

    velocity = mfem.VectorConstantCoefficient(
        [velocity * np.cos(angle), velocity * np.sin(angle)]
    )
    inflow = inflow_coeff()
    u0 = (
        FrankensteinInitial(sx / 2, sy / 2, radius, magnitude, frankenstein_param)
        if frankenstein_param != 0.0
        else GaussianInitial(sx / 2, sy / 2, radius, magnitude)
    )

    # 7. Define the initial conditions, save the corresponding grid function to
    #    a file
    u = mfem.GridFunction(fes)
    u.ProjectCoefficient(u0)
    mesh.Print("advection.mesh", 8)
    u.Save("advection-init.gf", 8)

    expected = AnalyticSolution(velocity, u0, sx, sy)

    if paraview:
        pd = mfem.ParaViewDataCollection("advection", mesh)
        pd.SetPrefixPath("ParaView")
        pd.RegisterField("solution", u)
        pd.SetLevelsOfDetail(order)
        pd.SetDataFormat(mfem.VTKFormat_BINARY)
        pd.SetHighOrderOutput(True)
        pd.SetCycle(0)
        pd.SetTime(0.0)
        pd.Save()

    if visualization:
        sol_sock = mfem.socketstream("localhost", 19916)
        sol_sock.precision(8)
        sol_sock.send_solution(mesh, u)

    adv = AdvectionEvolution(fes, velocity, inflow)

    ode_solver.Init(adv)
    t = 0.0
    ti = 0

    expected.set_t(t)
    l2 = u.ComputeL2Error(expected)
    print(f"Initial L2 error: {l2}")

    while True:
        if t > t_final - dt / 2:
            break
        t, dt = ode_solver.Step(u, t, dt)
        ti = ti + 1

        if ti % vis_steps == 0:
            expected.set_t(t)
            l2 = u.ComputeL2Error(expected)
            print(
                "time step: " + str(ti) + ", time: " + str(np.round(t, 3)),
                ", L2 error: ",
                str(l2),
            )
            if paraview:
                pd.SetCycle(ti)
                pd.SetTime(t)
                pd.Save()
            if visualization:
                sol_sock.send_solution(mesh, u)

    expected.set_t(t)
    l2 = u.ComputeL2Error(expected)
    print(f"Final L2 error: {l2}")

    u.Save("advection-final.gf", 8)


if __name__ == "__main__":
    from mfem.common.arg_parser import ArgParser

    parser = ArgParser(description="Simple advection app")
    parser.add_argument(
        "-vel",
        "--velocity",
        default=1.0,
        action="store",
        type=float,
        help="The magnitude of the advection velocity.",
    )
    parser.add_argument(
        "-a",
        "--angle",
        default=0,
        action="store",
        type=float,
        help="The angle of the velocity field relative to the x-axis, in degrees.",
    )
    parser.add_argument(
        "-rad",
        "--radius",
        default=1.0,
        action="store",
        type=float,
        help="Radius of the perturbation being advected.",
    )
    parser.add_argument(
        "-mag",
        "--magnitude",
        default=1.0,
        action="store",
        type=float,
        help="Magnitude of the perturbation being advected.",
    )
    parser.add_argument(
        "-fr",
        "--frankenstein-param",
        default=0,
        type=float,
        help="The 'a' parameter in the 'Frankenstein's monster' type function for the "
        "perturbation. If 0 then a guassian will be used instead.",
    )
    parser.add_argument(
        "-vis",
        "--visualization",
        action="store_true",
        help="Enable GLVis visualization",
    )
    parser.add_argument(
        "-o",
        "--order",
        action="store",
        default=3,
        type=int,
        help="Order (degree) of the finite elements.",
    )
    parser.add_argument(
        "-d",
        "--device",
        default="cpu",
        type=str,
        help="Device configuration string, see Device::Configure().",
    )
    parser.add_argument(
        "-tf", "--t-final", default=40, type=float, help="Final time; start time is 0."
    )
    parser.add_argument(
        "-dt", "--time-step", default=0.05, type=float, help="Time step."
    )
    parser.add_argument(
        "-nx",
        "--x-resolution",
        default=64,
        type=int,
        help="Number of elements in the x-direction.",
    )
    parser.add_argument(
        "-ny",
        "--y-resolution",
        default=16,
        type=int,
        help="Number of elements in the y-direction.",
    )
    parser.add_argument(
        "-sx",
        "--x-size",
        default=40.0,
        type=float,
        help="The length of the domain in the x-direction.",
    )
    parser.add_argument(
        "-sy",
        "--y-size",
        default=10.0,
        type=float,
        help="The length of the domain in the y-direction.",
    )
    parser.add_argument(
        "-vs",
        "--visualization-steps",
        default=5,
        type=int,
        help="Visualize every n-th timestep.",
    )
    parser.add_argument(
        "-paraview",
        "--paraview-datafiles",
        default=False,
        action="store_true",
        help="Save data files for ParaView (paraview.org) visualization.",
    )

    args = parser.parse_args()
    parser.print_options(args)

    run(
        nx=args.x_resolution,
        ny=args.y_resolution,
        sx=args.x_size,
        sy=args.y_size,
        order=args.order,
        visualization=args.visualization,
        angle=args.angle * np.pi / 180,
        velocity=args.velocity,
        radius=args.radius,
        magnitude=args.magnitude,
        frankenstein_param=args.frankenstein_param,
        device=args.device,
        t_final=args.t_final,
        dt=args.time_step,
        vis_steps=args.visualization_steps,
        paraview=args.paraview_datafiles,
    )
