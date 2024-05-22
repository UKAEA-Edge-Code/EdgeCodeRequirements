using OrdinaryDiffEq
using Trixi: DGSEM, TreeMesh, SemidiscretizationHyperbolic, semidiscretize, SVector, SummaryCallback
using Trixi: SaveRestartCallback, SaveSolutionCallback, CallbackSet

include("DN_source.jl")             # Symbolic source term callable function.
include("nonconservative.jl")       # Dummy hyperbolic equations with non-conservative/nonlinear terms.
include("anisotropic_diffusion.jl") # Actual parabolic equations to solve

#
# DN solution parameters
#
aa = 2.0
mm = 10.0
eeps = 1.0e-6

# kappa
diffusivity_para() = 1.0
diffusivity_perp() = eeps

#
# Equations
#
equations           = DummyEquations2D()
equations_parabolic = AnisotropicDiffusion2D(diffusivity_para(),diffusivity_perp(),equations)

#
# Source terms (keep the initial Source term as 4th variable for simplicity)
#
@inline function source_terms_DN(u, x, t, equations)
    S_eps = u[4]
	return SVector(S_eps,0.0,0.0,0.0)
end

#
# Various initial conditions
#
pi_ = float(pi)
initial_condition_zero = (x, t, equations) -> SVector(0.0,0.0,0.0,0.0)
@inline function initial_condition_DN(x, t, equations)
  bx = (aa*(2*x[2]-1)*cos(mm*pi_*x[1])+ pi_)
  by = (pi_*aa*mm*(x[2]^2 - x[2])*sin(mm*pi_*x[1]))
  rho = sin(pi_*x[2]+aa*(x[2]^2 - x[2])*cos(mm*pi_*x[1]))+eeps*cos(2*pi_*x[1])*sin(pi_*x[2]);
  S_eps = SourceFunctionDN(x[1],x[2],aa,mm,eeps)
  return SVector(rho,bx,by,S_eps)
end
@inline function initial_condition_sinpix_sinpiy(x, t, equations)
  bx = (aa*(2*x[2]-1)*cos(mm*pi_*x[1])+ pi_)
  by = (pi_*aa*mm*(x[2]^2 - x[2])*sin(mm*pi_*x[1]))
  rho = sin(pi_*x[1])*sin(pi_*x[2]);
  S_eps = SourceFunctionDN(x[1],x[2],aa,mm,eeps)
  return SVector(rho,bx,by,S_eps)
end
@inline function initial_condition_DN_rhozero(x, t, equations)
  bx = (aa*(2*x[2]-1)*cos(mm*pi_*x[1])+ pi_)
  by = (pi_*aa*mm*(x[2]^2 - x[2])*sin(mm*pi_*x[1]))
  S_eps = SourceFunctionDN(x[1],x[2],aa,mm,eeps)
  return SVector(0.0,bx,by,S_eps) # rho = 0
end
@inline function initial_condition_DN_dumzero(x, t, equations)
  rho = sin(pi_*x[2]+aa*(x[2]^2 - x[2])*cos(mm*pi_*x[1]))+eeps*cos(2*pi_*x[1])*sin(pi_*x[2]);
  return SVector(rho,0.0,0.0,0.0) # dummy variables zero
end



#
# solver/mesh
#
dg = DGSEM(polydeg = 3)
coordinates_min = (0.0, 0.0)
coordinates_max = (1.0, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 100_000,
                periodicity = false) # squares

#
# Actual initial/boundary conditions
#
initial_condition    = initial_condition_DN
boundary_condition   = BoundaryConditionDirichlet(initial_condition_zero)
boundary_condition_2 = BoundaryConditionDirichlet(initial_condition_DN_dumzero) # additional BC for +/- x boundary

# Set condition at each boundary
boundary_conditions = (; x_neg = boundary_condition_2,
                         y_pos = boundary_condition,
                         y_neg = boundary_condition,
                         x_pos = boundary_condition_2)

boundary_conditions_parabolic = (; x_neg = boundary_condition_2,
                                   y_pos = boundary_condition,
                                   y_neg = boundary_condition,
                                   x_pos = boundary_condition_2)

#
# Semidiscretization (hyperbolic+parabolic)
#
semi = SemidiscretizationHyperbolicParabolic(mesh, (equations, equations_parabolic),
                                             initial_condition, dg;
                                             boundary_conditions = (boundary_conditions,
                                                                    boundary_conditions_parabolic),source_terms=source_terms_DN)
#
# Define set of ODEs and callbacks
#
tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()
alive_callback = AliveCallback(alive_interval = 300)
analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval, uEltype = real(dg))
callbacks = CallbackSet(summary_callback, alive_callback, analysis_callback)

###############################################################################
# run the simulation

time_int_tol = 1e-8
sol = solve(ode, RDPK3SpFSAL49(); abstol = time_int_tol, reltol = time_int_tol,
            ode_default_options()..., callback = callbacks)
summary_callback()
