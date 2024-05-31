using OrdinaryDiffEq
using Trixi
include("DN_source.jl")
include("anisotropic_hyperbolic_diffusion_2d_sourceterm.jl")

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
equations = AnisoHyperDiffusEquations2D(Kpara=diffusivity_para(),Kperp=diffusivity_perp())

#
# Source
#
@inline function source_terms_DN(u, x, t, equations::AnisoHyperDiffusEquations2D)
    (; inv_Tr) = equations
    phi, q1, q2, _, _, S_eps = u
    du2 = -inv_Tr * q1
    du3 = -inv_Tr * q2
    return SVector(S_eps, du2, du3, 0, 0, 0)
end

#
# ICs (also used for BCs)
#
pi_ = float(pi)
@inline function initial_condition_DN(x, t, equations::AnisoHyperDiffusEquations2D)
    if iszero(t)
        phi = sin(pi_*x[2]+aa*(x[2]^2 - x[2])*cos(mm*pi_*x[1]))+eeps*cos(2*pi_*x[1])*sin(pi_*x[2]);
        q1 = du_dx(x[1],x[2],aa,mm,eeps)
        q2 = du_dy(x[1],x[2],aa,mm,eeps)
        bx = (aa*(2*x[2]-1)*cos(mm*pi_*x[1])+ pi_)
        by = (pi_*aa*mm*(x[2]^2 - x[2])*sin(mm*pi_*x[1]))
        S_eps = SourceFunctionDN(x[1],x[2],aa,mm,eeps)
    else
        phi = sin(pi_*x[2]+aa*(x[2]^2 - x[2])*cos(mm*pi_*x[1]))+eeps*cos(2*pi_*x[1])*sin(pi_*x[2]);
        q1 = du_dx(x[1],x[2],aa,mm,eeps)
        q2 = du_dy(x[1],x[2],aa,mm,eeps)
        bx = (aa*(2*x[2]-1)*cos(mm*pi_*x[1])+ pi_)
        by = (pi_*aa*mm*(x[2]^2 - x[2])*sin(mm*pi_*x[1]))
        S_eps = SourceFunctionDN(x[1],x[2],aa,mm,eeps)
    end
    return SVector(phi, q1, q2, bx, by, S_eps)
end

@inline function initial_condition_DN_solzero(x, t, equations::AnisoHyperDiffusEquations2D)
        T = eltype(x)
        bx = (aa*(2*x[2]-1)*cos(mm*pi_*x[1])+ pi_)
        by = (pi_*aa*mm*(x[2]^2 - x[2])*sin(mm*pi_*x[1]))
        S_eps = SourceFunctionDN(x[1],x[2],aa,mm,eeps)
    return SVector(zero(T), zero(T), zero(T), bx, by, S_eps)
end

#
# BCs
#
@inline function boundary_condition_DN(u_inner, orientation, direction, x, t,
                                       surface_flux_function,
                                       equations::AnisoHyperDiffusEquations2D)
    # Known solution at boundary (including magnetic field and source term)
    u_boundary = initial_condition_DN(x, one(t), equations)

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end

@inline function boundary_condition_DN_upper_lower(u_inner, orientation, direction, x, t,
                                                   surface_flux_function,
                                                   equations::AnisoHyperDiffusEquations2D)
    # u, u_x, u_y = 0 on boundary (magnetic field and source terms as prescribed)
    u_boundary = initial_condition_DN_solzero(x, one(t), equations)

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end


###############################################################################
# semidiscretization of the hyperbolic diffusion equations


initial_condition = initial_condition_DN

boundary_condition   = boundary_condition_DN
boundary_condition_2 = boundary_condition_DN_upper_lower
# Set condition at each boundary
boundary_conditions = (  x_neg = boundary_condition,
                         y_pos = boundary_condition_2,
                         y_neg = boundary_condition_2,
                         x_pos = boundary_condition)

solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

x_shift = 0.0 # shift to change angle of field at boundary
coordinates_min = (0.0 + x_shift, 0.0)
coordinates_max = (1.0 + x_shift, 1.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 30_000,
                periodicity = (false,false))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_conditions,
                                    source_terms = source_terms_DN)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0)
ode = semidiscretize(semi, tspan);

summary_callback = SummaryCallback()

resid_tol = 5.0e-12
steady_state_callback = SteadyStateCallback(abstol = resid_tol, reltol = 0.0)

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback, steady_state_callback,
                        analysis_callback, alive_callback,
                        save_solution,
                        stepsize_callback)
#callbacks = CallbackSet(summary_callback)
###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep = false, callback = callbacks);
summary_callback() # print the timer summary
