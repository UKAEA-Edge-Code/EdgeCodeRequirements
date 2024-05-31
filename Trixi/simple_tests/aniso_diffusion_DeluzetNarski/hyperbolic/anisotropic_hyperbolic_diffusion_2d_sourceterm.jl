@doc raw"""
    AnisoHyperDiffusEquations2D
The linear field-aligned hyperbolic diffusion equations in two space dimensions.
A description of this system can be found in Sec. 2.5 of the book "I Do Like CFD, Too: Vol 1".
The book is freely available at [http://www.cfdbooks.com/](http://www.cfdbooks.com/) and further analysis can be found in
the paper by Nishikawa [DOI: 10.1016/j.jcp.2007.07.029](https://doi.org/10.1016/j.jcp.2007.07.029)
"""

using Trixi
using Trixi: AbstractEquations
import Trixi: varnames, cons2prim, cons2cons, cons2entropy, residual_steady_state, default_analysis_integrals, flux, max_abs_speeds, max_abs_speed_naive, have_nonconservative_terms

abstract type AbstractAnisoHyperDiffusEquations{NDIMS, NVARS} <:
              AbstractEquations{NDIMS, NVARS} end
struct AnisoHyperDiffusEquations2D{RealT <: Real} <:
       AbstractAnisoHyperDiffusEquations{2, 6}
    Lr::RealT         # reference length scale
    inv_Tr::RealT     # inverse of the reference time scale
    Kpara::RealT # diffusion para
    Kperp::RealT # diffusion perp
end

function AnisoHyperDiffusEquations2D(; Kpara = 1.0, Kperp=0.1, Lr = inv(2pi))
    Tr = Lr^2 / max(Kpara, Kperp)
    AnisoHyperDiffusEquations2D(promote(Lr, inv(Tr), Kpara, Kperp)...)
end

varnames(::typeof(cons2cons),    ::AnisoHyperDiffusEquations2D) = ("phi", "q1", "q2", "bx", "by", "source")
varnames(::typeof(cons2prim),    ::AnisoHyperDiffusEquations2D) = ("phi", "q1", "q2", "bx", "by", "source")
varnames(::typeof(cons2entropy), ::AnisoHyperDiffusEquations2D) = ("phi", "q1", "q2", "bx", "by", "source")

function default_analysis_errors(::AnisoHyperDiffusEquations2D)
    (:l2_error, :linf_error, :residual)
end

@inline function residual_steady_state(du, ::AnisoHyperDiffusEquations2D)
    abs(du[1])
end

# Calculate 1D flux in for a single point
@inline function flux(u, orientation::Integer,
                      equations::AnisoHyperDiffusEquations2D)
    phi, q1, q2, bx, by, _ = u
    (; inv_Tr, Kpara, Kperp) = equations
    bnorm = sqrt(bx*bx+by*by)
    bxh   = bx/bnorm
    byh   = by/bnorm
    dotbbG =  bxh*q1 + byh*q2
    zero_ = zero(phi)
    if orientation == 1
        f1 = - Kpara * bxh * dotbbG - Kperp * (q1 - bxh*dotbbG)
        f2 = -phi * inv_Tr
        f3 = zero_
    else
        f1 = - Kpara * byh * dotbbG - Kperp * (q2 - byh*dotbbG)
        f2 = zero_
        f3 = -phi * inv_Tr
    end

    return SVector(f1, f2, f3, zero_, zero_, zero_)
end

# Note, this directional vector is not normalized
@inline function flux(u, normal_direction::AbstractVector,
                      equations::AnisoHyperDiffusEquations2D)
    phi, q1, q2, bx, by, _ = u
    (; inv_Tr, Kpara, Kperp) = equations
    bnorm = sqrt(bx*bx+by*by)
    bxh   = bx/bnorm
    byh   = by/bnorm
    dotbbG =  bxh*q1 + byh*q2
    fq1 = - Kpara * bxh * dotbbG - Kperp * (q1 - bxh*dotbbG)
    fq2 = - Kpara * byh * dotbbG - Kperp * (q2 - byh*dotbbG)
    f1 =  normal_direction[1] * fq1 + normal_direction[2] * fq2 # projection of diffusive flux onto normal direction
    f2 = -phi * inv_Tr * normal_direction[1]
    f3 = -phi * inv_Tr * normal_direction[2]
    zero_ = zero(phi)
    return SVector(f1, f2, f3, zero_, zero_, zero_)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::AnisoHyperDiffusEquations2D)
    sqrt(equations.Kpara * equations.inv_Tr)
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::AnisoHyperDiffusEquations2D)
    sqrt(equations.Kpara * equations.inv_Tr) * norm(normal_direction)
end

@inline function flux_godunov(u_ll, u_rr, orientation::Integer,
                              equations::AnisoHyperDiffusEquations2D)
    # Obtain left and right fluxes
    phi_ll, q1_ll, q2_ll = u_ll
    phi_rr, q1_rr, q2_rr = u_rr
    f_ll = flux(u_ll, orientation, equations)
    f_rr = flux(u_rr, orientation, equations)

    # this is an optimized version of the application of the upwind dissipation matrix:
    #   dissipation = 0.5*R_n*|Λ|*inv(R_n)[[u]]
    λ_max = sqrt(equations.Kpara * equations.inv_Tr)
    f1 = 1 / 2 * (f_ll[1] + f_rr[1]) - 1 / 2 * λ_max * (phi_rr - phi_ll)
    if orientation == 1 # x-direction
        f2 = 1 / 2 * (f_ll[2] + f_rr[2]) - 1 / 2 * λ_max * (q1_rr - q1_ll)
        f3 = 1 / 2 * (f_ll[3] + f_rr[3])
    else # y-direction
        f2 = 1 / 2 * (f_ll[2] + f_rr[2])
        f3 = 1 / 2 * (f_ll[3] + f_rr[3]) - 1 / 2 * λ_max * (q2_rr - q2_ll)
    end

    return SVector(f1, f2, f3, 0, 0)
end

@inline function flux_godunov(u_ll, u_rr, normal_direction::AbstractVector,
                              equations::AnisoHyperDiffusEquations2D)
    # Obtain left and right fluxes
    phi_ll, q1_ll, q2_ll = u_ll
    phi_rr, q1_rr, q2_rr = u_rr
    f_ll = flux(u_ll, normal_direction, equations)
    f_rr = flux(u_rr, normal_direction, equations)

    # this is an optimized version of the application of the upwind dissipation matrix:
    #   dissipation = 0.5*R_n*|Λ|*inv(R_n)[[u]]
    λ_max = sqrt(equations.Kpara * equations.inv_Tr)
    f1 = 1 / 2 * (f_ll[1] + f_rr[1]) -
         1 / 2 * λ_max * (phi_rr - phi_ll) *
         sqrt(normal_direction[1]^2 + normal_direction[2]^2)
    f2 = 1 / 2 * (f_ll[2] + f_rr[2]) -
         1 / 2 * λ_max * (q1_rr - q1_ll) * normal_direction[1]
    f3 = 1 / 2 * (f_ll[3] + f_rr[3]) -
         1 / 2 * λ_max * (q2_rr - q2_ll) * normal_direction[2]

    return SVector(f1, f2, f3, 0, 0)
end

@inline have_constant_speed(::AnisoHyperDiffusEquations2D) = True()

@inline function max_abs_speeds(eq::AnisoHyperDiffusEquations2D)
    λ = sqrt(eq.Kpara * eq.inv_Tr)
    return λ, λ
end

@inline function max_abs_speeds(u, eq::AnisoHyperDiffusEquations2D)
    phi, _, _, _, _, _ = u
    λ = sqrt(eq.Kpara * eq.inv_Tr)
    return λ, λ
end

# Convert conservative variables to primitive
@inline cons2prim(u, equations::AnisoHyperDiffusEquations2D) = u

# Convert conservative variables to entropy found in I Do Like CFD, Too, Vol. 1
@inline function cons2entropy(u, equations::AnisoHyperDiffusEquations2D)
    phi, q1, q2, bx, by, Su = u
    w1 = phi
    w2 = equations.Lr^2 * q1
    w3 = equations.Lr^2 * q2

    return SVector(w1, w2, w3, bx, by, Su)
end

# Calculate entropy for a conservative state `u` (here: same as total energy)
@inline function entropy(u, equations::AnisoHyperDiffusEquations2D)
    energy_total(u, equations)
end

# Calculate total energy for a conservative state `u`
@inline function energy_total(u, equations::AnisoHyperDiffusEquations2D)
    # energy function as found in equations (2.5.12) in the book "I Do Like CFD, Vol. 1"
    phi, q1, q2, _, _, _ = u
    return 0.5 * (phi^2 + equations.Lr^2 * (q1^2 + q2^2))
end
