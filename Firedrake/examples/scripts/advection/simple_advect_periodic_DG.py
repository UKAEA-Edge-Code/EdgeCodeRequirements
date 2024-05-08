# simple_advect_periodic_DG.py
# Attempt at time-dependent solver of 2D advection

from firedrake import *
import math
from irksome import Dt, MeshConstant, TimeStepper, RadauIIA, GaussLegendre, LobattoIIIC
import time
import numpy

meshres = 64
mesh = PeriodicRectangleMesh(meshres, meshres/4, 40, 10, quadrilateral=True)  # note non-periodic for debug
V = FunctionSpace(mesh, "DG", 3)
V4 = VectorFunctionSpace(mesh, "CG", 3)  # for advection velocity

T = 40.0
timeres = 400

t = Constant(0.0)
dt = Constant(T/timeres)

# parameters for irksome
#butcher_tableau = LobattoIIIC(2)
#butcher_tableau = GaussLegendre(2)
butcher_tableau = RadauIIA(2)

# model parameters
v_0 = as_vector([1.0,0.0])  # constant advection velocity vector

x, y = SpatialCoordinate(mesh)
n = Function(V)
v = TestFunction(V)

# Gaussian init data
s = 2.0  # Gaussian width of initial data
n.interpolate(exp(-((x-20)*(x-20)+(y-5)*(y-5))/(s*s)))
n_exact = exp(-((x-20)*(x-20)+(y-5)*(y-5))/(s*s))

# Frankenstein's monster initial data
#s = 2.0
#a = 1.0
#n.interpolate(conditional(le((x-20)*(x-20)+(y-5)*(y-5), s*s), exp(-a*a/(s*s-(x-20)**2-(y-5)**2)), 0))
#n_exact = conditional(le((x-20)*(x-20)+(y-5)*(y-5), s*s), exp(-a*a/(s*s-(x-20)**2-(y-5)**2)), 0)


driftvel = Function(V4)
driftvel.interpolate(as_vector([1.0, 0.0]))

norms = FacetNormal(mesh)
driftvel_n = 0.5*(dot(driftvel, norms)+abs(dot(driftvel, norms)))

# TRIALCODE check init data
#File("simple_advect_periodic_DG_init.pvd").write(n)
#quit()

h = 2.0/meshres

F = Dt(n)*v*dx \
  - v*div(n*driftvel)*dx \
  + driftvel_n('-')*(n('-') - n('+'))*v('-')*dS \
  + driftvel_n('+')*(n('+') - n('-'))*v('+')*dS \

# matching params I've used before in advection work ...
luparams = {"mat_type": "aij",
            "ksp_type": "preonly",
            "pc_type": "lu"}

stepper = TimeStepper(F, butcher_tableau, t, dt, n, solver_parameters=luparams)

outfile = File("simple_advect_periodic_DG.pvd")

cnt=0
start = time.time()

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))
    if(cnt % 1 == 0):
       print("outputting data ...\n")
       outfile.write(n)
    cnt=cnt+1
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))

end = time.time()
wall_time = end-start

ns = n
L2_error_n_DG = norm(ns-n_exact, "L2")

print("done.")
print("\n")
print("L2 error norm n DG: %.2e" % L2_error_n_DG)
print("wall time:"+str(wall_time)+"\n")


