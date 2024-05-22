@doc raw"""
    AnisotropicDiffusion2D(diffusivity, equations)

`AnisotropicDiffusion2D` represents an anisotropic diffusion term ``\nabla \cdot (\kappa\nabla u))``
with diffusivity ``\kappa`` aligned along field lines.
"""
struct AnisotropicDiffusion2D{E, N, Kpara, Kperp} <: AbstractAnisotropicDiffusion{2, N}
    Kpara::Kpara
    Kperp::Kperp
    equations_hyperbolic::E
end

function AnisotropicDiffusion2D(Kpara,Kperp, equations_hyperbolic)
    AnisotropicDiffusion2D{typeof(equations_hyperbolic), nvariables(equations_hyperbolic),
			   typeof(Kpara),typeof(Kperp)}(Kpara, Kperp, equations_hyperbolic)
end

function varnames(variable_mapping, equations_parabolic::AnisotropicDiffusion2D)
    varnames(variable_mapping, equations_parabolic.equations_hyperbolic)
end

function flux(u, gradients, orientation::Integer, equations::AnisotropicDiffusion2D)
    (; Kpara, Kperp) = equations
    dudx, dudy = gradients
    _, bx, by, _ = u
    bnorm = sqrt(bx*bx+by*by)
    bxh   = bx/bnorm
    byh   = by/bnorm
    dotbbgradu = bxh*dudx[1] + byh*dudy[1]
    if orientation == 1
        fx    = Kpara * bxh * dotbbgradu + Kperp * (dudx[1] - bxh*dotbbgradu)
        return SVector(fx,0.0,0.0,0.0)
    else # if orientation == 2
        fy    = Kpara * byh * dotbbgradu + Kperp * (dudy[1] - byh*dotbbgradu)
        return SVector(fy,0.0,0.0,0.0)
    end
end
