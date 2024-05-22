#
# Advection and non-conservative terms (dummy equations)
# Parabolic equations will be built on top of these
#
using Trixi
using Trixi: AbstractEquations, Gradient, Divergence
import Trixi: varnames, default_analysis_integrals, have_nonconservative_terms

#
# variable coefficient equations
#
struct DummyEquations2D{} <:
       AbstractEquations{2, 4} # 2D, 4 variables
end

# retain bx, by and a source function as variables
varnames(::typeof(cons2cons), ::DummyEquations2D) = ("rho", "bx", "by", "sourcefun")
default_analysis_integrals(::DummyEquations2D) = ()

#
# There are different flavours of the following functions depending on context,
# eps. what type of mesh is being used.
#

#
# Fluxes
#
# Zero conservative flux
@inline function flux(u, orientation::Integer, equations::DummyEquations2D)
    return zero(u)
end
@inline function flux(u, normal::AbstractVector, equations::DummyEquations2D)
    return zero(u)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, ::DummyEquations2D)
    return 0.0
end
@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector, ::DummyEquations2D)
    return 0.0
end


# ENABLE IF nonlinear terms added
#have_nonconservative_terms(::DummyEquations2D) = Trixi.True()

# This "nonconservative numerical flux" implements the nonconservative terms.
# In general, nonconservative terms can be written in the form
#   g(u) ∂ₓ h(u)
# Thus, a discrete difference approximation of this nonconservative term needs
# - `u mine`:  the value of `u` at the current position (for g(u))
# - `u_other`: the values of `u` in a neighborhood of the current position (for ∂ₓ h(u))
function flux_nonconservative(u_mine, u_other, orientation,
                              equations::DummyEquations2D)
    _, bx, by, _             = u_mine
    scalar, _, _, _          = u_other

    return SVector(zero(scalar), zero(scalar), zero(scalar))
end

#
# Boundary conditions (same condition for different cases)
#
@inline function boundary_condition_fixed_value(u_inner, normal_direction::AbstractVector,
                                                x, t, surface_flux_function,
                                                equations::DummyEquations2D)
    return SVector(zero(eltype(u_inner)),0.0,0.0,0.0)
end

@inline function boundary_condition_fixed_value(u_inner, orientation, direction,
		                                x, t, surface_flux_function,
                                                equations::DummyEquations2D)
    # get the appropriate normal vector from the orientation
    if orientation == 1
        normal_direction = SVector(1, 0)
    else # orientation == 2
        normal_direction = SVector(0, 1)
    end

    return boundary_condition_fixed_value(u_inner, normal_direction, direction,
                                        x, t, surface_flux_function, equations)
end

@inline function boundary_condition_fixed_value(u_inner, normal_direction::AbstractVector, direction,
		                                x, t, surface_flux_function,
                                                equations::DummyEquations2D)
    # flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
    # to be inward pointing on the -x and -y sides due to the orientation convention used by StructuredMesh
    if isodd(direction)
        boundary_flux = -boundary_condition_fixed_value(u_inner, -normal_direction,
                                                        x, t, surface_flux_function,
                                                        equations)
    else
        boundary_flux = boundary_condition_fixed_value(u_inner, normal_direction,
                                                       x, t, surface_flux_function,
                                                       equations)
    end

    return boundary_flux
end

#
# Initial conditions
#
function initial_condition_constant(x, t, equations::DummyEquations2D)
    rho = 0.0
    bx  = 0.0
    by  = 0.0
    return SVector(rho,bx,by,0.0)
end

#
# Source terms
#
@inline function source_terms_constant(u, x, t, equations::DummyEquations2D)
    Srho = 0.0
    Sbx  = 0.0
    Sby  = 0.0
    return SVector(Srho,Sbx,Sby,0.0)
end

@inline function (boundary_condition::BoundaryConditionDirichlet)(flux_inner,
                                                                  u_inner,
                                                                  normal::AbstractVector,
                                                                  x, t,
                                                                  operator_type::Gradient,
                                                                  equations::DummyEquations2D)
    # BCs are usually specified as conservative variables so we convert them to primitive variables
    #  because the gradients are assumed to be with respect to the primitive variables
    u_boundary = boundary_condition.boundary_value_function(x, t, equations)

    return cons2prim(u_boundary, equations)
end

@inline function (boundary_condition::BoundaryConditionDirichlet)(flux_inner,
                                                                  u_inner,
                                                                  normal::AbstractVector,
                                                                  x, t,
                                                                  operator_type::Divergence,
                                                                  equations::DummyEquations2D)
    # for Dirichlet boundary conditions, we do not impose any conditions on the viscous fluxes
    return flux_inner
end
