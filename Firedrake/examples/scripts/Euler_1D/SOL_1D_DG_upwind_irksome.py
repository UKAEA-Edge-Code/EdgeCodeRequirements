# SOL_1D_DG_upwind_irksome.py
# Attempt at time-dependent solver of 1D SOL equations, with upwind flux
# based on:
# https://www.firedrakeproject.org/demos/DG_advection.py.html
# irksome implementation follows
# https://www.firedrakeproject.org/Irksome/demos/demo_cahnhilliard.py.html
# and
# https://www.firedrakeproject.org/Irksome/demos/demo_monodomain_FHN.py.html

from firedrake import *
import math
from irksome import Dt, GaussLegendre, MeshConstant, TimeStepper

meshres = 200
mesh = IntervalMesh(meshres, -1, 1)

# velocity space V2 is intended to be CG(N), density space V1 needs to be DG(N-1) for N = 1,2,3, ... (element order)
# this is the discretization used in eq.(2) of https://arxiv.org/abs/1605.00551 by Natale, Shipton, and Cotter

V1 = FunctionSpace(mesh, "DG", 0)
V2 = VectorFunctionSpace(mesh, "CG", 1)
V = V1*V2

# time parameters
T = 4.0  # duration
timeres = 100  # number of timesteps
t = Constant(0.0)
dt = Constant(T/timeres)

# parameters for irksome
butcher_tableau = GaussLegendre(2)

# model parameters
nstar = Constant(1.0)
#diffusion = 1.0e0
Temp = 1.0

x = SpatialCoordinate(mesh)
nu = Function(V)
n, u = split(nu)
#n, u = nu.split()  # doesn't work with Irksome
v1, v2 = TestFunctions(V)

# sonic outflow equilibrium init data
#nu.sub(0).interpolate((nstar/sqrt(Temp))*(1+sqrt(1-x[0]*x[0])))
#nu.sub(1).interpolate((sqrt(Temp)/(x[0]))*(1-sqrt(1-x[0]*x[0])))

# Gaussian blob init data
width = 0.1
nu.sub(0).interpolate((1.0+0.2*(1/sqrt(2*math.pi*width**2))*exp(-x[0]**2/(2*width**2))))
nu.sub(1).interpolate(as_vector([0.0+1.0*x[0]]))

# source function
nstarFunc = Function(V)
nstarFunc.sub(0).interpolate(nstar + 0.0*x[0])

# TRIALCODE check init data
#File("SOL_1D_DG_upwind_irksome_init.pvd").write(nu.sub(0), nu.sub(1))
#quit()

norm = FacetNormal(mesh)
u_n = 0.5*(dot(u,norm)+abs(dot(u,norm)))

# outflow BC imposed weakly in here, penalty term in ds (note conditional not really needed here as both are outflow)

F = -((Dt(n)*v1)*dx + (n*dot(Dt(u),v2))*dx) \
   + (n*dot(u, grad(v1))+v1*nstar)*dx \
   + (nstar*dot(u,v2)+n*u[0]*grad(dot(u,v2))[0]+n*u[0]*dot(u, grad(v2[0]))+Temp*n*grad(v2[0])[0])*dx \
   - (v1('+') - v1('-'))*(u_n('+')*n('+') - u_n('-')*n('-'))*dS \
   + (u('+')[0]*v2('+')[0]-u('-')[0]*v2('-')[0])*(n('+')*u_n('+')-n('-')*u_n('-'))*dS \
   - conditional(dot(u, norm) > 0, v1*dot(u, norm)*n, 0.0)*ds \

# weirdly, penultimate term does not seem to be needed?


# params taken from Cahn-Hilliard example cited above
params = {'snes_monitor': None, 'snes_max_it': 100,
          'snes_linesearch_type': 'l2',
          'ksp_type': 'preonly',
          'pc_type': 'lu', 'mat_type': 'aij',
          'pc_factor_mat_solver_type': 'mumps'}

# Dirichlet BCs are needed for boundary velocity

bc_test1 = DirichletBC(V.sub(1),as_vector([-1.0]),1)
bc_test2 = DirichletBC(V.sub(1),as_vector([1.0]),2)

stepper = TimeStepper(F, butcher_tableau, t, dt, nu, solver_parameters=params, bcs=[bc_test1, bc_test2])

outfile = File("SOL_1D_DG_upwind_irksome.pvd")

nu.sub(0).rename("n")
nu.sub(1).rename("u")

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))
    outfile.write(nu.sub(0), nu.sub(1))
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))

print("done.")
print("\n")

# export final state if desired
#ns, us = nu.split()
#ns.rename("density")
#us.rename("velocity")
#File("SOL_1D_DG_upwind_irksome_final.pvd").write(ns, us)
