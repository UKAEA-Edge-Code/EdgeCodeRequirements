#
# Compute the DN source term symbolically and convert to callable function
#
using Symbolics

MyPi = float(pi) # just want to make sure it gets converted correctly

@variables x y α m ϵ

function b_hat(x)
  let π = oftype(MyPi, π)
    [(α*(2*x[2]-1)*cos(m*π*x[1]) + π) /
        sqrt( (α*(2*x[2]-1)*cos(m*π*x[1])+ π)^2 + (π*α*m*(x[2]^2 - x[2])*sin(m*π*x[1]))^2 ),
      π*α*m*(x[2]^2 - x[2])*sin(m*π*x[1]) /
        sqrt( (α*(2*x[2]-1)*cos(m*π*x[1])+ π)^2 + (π*α*m*(x[2]^2 - x[2])*sin(m*π*x[1]))^2 )
    ]
  end
end

function u_eps(x)
  let π = oftype(MyPi, π)
    sin(π*x[2] + α*(x[2]^2 - x[2])*cos(m*π*x[1])) + ϵ*cos(2*π*x[1])*sin(π*x[2])
  end
end

bhat = b_hat([x,y]) # symbolic
ueps = u_eps([x,y]) # symbolic

dudx = Symbolics.derivative(ueps,x)
dudy = Symbolics.derivative(ueps,y)
dotbbgradu = bhat[1]*dudx + bhat[2]*dudy
fx = bhat[1]*dotbbgradu + ϵ * (dudx - bhat[1]*dotbbgradu)
fy = bhat[2]*dotbbgradu + ϵ * (dudy - bhat[2]*dotbbgradu)
SourceTerm = Symbolics.derivative(fx,x) + Symbolics.derivative(fy,y) # \nabla \cdot (\kappa grad(u^{\epsilon})), symbolic

#Params = Dict(m => 10.0, α => 2.0, ϵ => 1e-6)
#S_0 = substitute(SourceTerm,Params)

SourceFunctionDN=build_function(SourceTerm,x,y,α,m,ϵ,expression=Val{false}) # callable function
du_dx=build_function(dudx,x,y,α,m,ϵ,expression=Val{false})
du_dy=build_function(dudy,x,y,α,m,ϵ,expression=Val{false})
#SourceFunctionDN=build_function(S_0,x,y,expression=Val{false}) # callable function
