# Simple advection test with fixed mesh.
#
# Run with (command line):
#    > julia simple_advect_periodic_DG.jl
# OR (in julia REPL):
#    julia> include("simple_advect_periodic_DG.jl")
# To plot result (and mesh):
#    julia> using Plots;
#    julia> pd = PlotData2D(sol)
#    julia> plot(pd, seriescolor=:heat)
#    julia> plot!(getmesh(pd))

using OrdinaryDiffEq, Trixi

# Linear advection equation comes with Trixi; just specify the advection velocity.
v_0 = (1.0, 0.0)
equation = LinearScalarAdvectionEquation2D(v_0)

Lx=40.
Ly=10.
cells_per_dimension = (64,16)
# Set up a 40x10 structured periodic mesh on rectangle with fixed size.
# Put (0,0) at the centre of the rectangle for convenience.
coordinates_min = (-Lx/2., -Ly/2.) # lower left corner of the rectangle
coordinates_max = ( Lx/2.,  Ly/2.) # upper right corner of the rectangle
mesh_static = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max; periodicity=true)

# Set up initial conditions.
# Here they are expressed as the exact solution so that a callback can be used to check the error.

# Translate coordinates for a periodic domain
function x_trans_periodic_2d(x, domain_length = SVector(Lx, Ly), center = SVector(0, 0))
    x_normalized = x .- center
    x_shifted = x_normalized .% domain_length
    x_offset = ((x_shifted .< -0.5 * domain_length) -
                (x_shifted .> 0.5 * domain_length)) .* domain_length
    return center + x_shifted + x_offset
end
function initial_condition_gauss(x, t, equation::LinearScalarAdvectionEquation2D)
    x_trans = x_trans_periodic_2d(x - equation.advection_velocity * t)
    s = 2.0
    scalar = exp(-(x_trans[1]^2 + x_trans[2]^2)/(s*s))
    return SVector(scalar)
end
function initial_condition_monster(x, t, equation::LinearScalarAdvectionEquation2D)
    x_trans = x_trans_periodic_2d(x - equation.advection_velocity * t)
    s = 2.0
    a = 0.1
    scalar = max(0,exp( -( a^2 / (s^2 - x_trans[1]^2 + x_trans[2]^2) ) ))
    return SVector(scalar)
end
initial_condition = initial_condition_gauss

# Select solver
solver = DGSEM(polydeg=3)

# Construct semidiscretization and ode
semi = SemidiscretizationHyperbolic(mesh_static, equation, initial_condition, solver)
tspan = (0.0, 40.0)
ode = semidiscretize(semi, tspan)

# Construct some callbacks
analysis_callback = AnalysisCallback(semi, interval=100)
save_solution = SaveSolutionCallback(dt = 10,
                                     save_initial_solution = true,
                                     save_final_solution = true)
callbacks = CallbackSet(analysis_callback,save_solution);

# Solve with `solve` from OrdinaryDiffEq
# (nonstiff explicit Runge-Kutta, Bogacki-Shampine 3/2 method, "ode23()" in MATLAB)
sol = solve(ode, BS3(), save_everystep = false, callback = callbacks);

# With initial_condition_gauss, s=2, t=40 (polydeg=3), v = (1,0):
#L2 error:    2.28576555e-04
#Linf error:  2.13218567e-03
