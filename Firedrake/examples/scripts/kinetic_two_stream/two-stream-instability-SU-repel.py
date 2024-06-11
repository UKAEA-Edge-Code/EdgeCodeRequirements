
# two-stream_instability_SU_repel.py
# Attempt at time-dependent solver of 1D1V two-stream instability
# non-concurrent Poisson solve (linear)
# SU correction as at
# https://github.com/volpatto/firedrake_scripts/blob/master/scripts/2D/adv-diff_supg_2D.py

from firedrake import *
import math
from irksome import Dt, MeshConstant, TimeStepper, RadauIIA
import time
import numpy

meshres = 64
mesh1d = PeriodicIntervalMesh(meshres, 2)
mesh2d = ExtrudedMesh(mesh1d, layers=meshres*4, layer_height=2/meshres, extrusion_type='uniform')
V1 = FunctionSpace(mesh2d, "CG", 2) # for f
V2 = FunctionSpace(mesh2d, "CG", 2, vfamily='R') # for phi
V3 = VectorFunctionSpace(mesh2d, "CG", 2)  # for phase-space advection velocity
V = V1*V2

# time parameters e.g. 40/4000
T = 40.0
timeres = 4000

t = Constant(0.0)
dt = Constant(T/timeres)

# parameters for irksome
butcher_tableau = RadauIIA(1)

# model parameters
omegaPsq = 10.0
sigma = 0.2

x, y = SpatialCoordinate(mesh2d)

f = Function(V1)
phi = TrialFunction(V2)
v1 = TestFunction(V1)
v2 = TestFunction(V2)

# Gaussian init data
f.interpolate(1.0*(1/sqrt(8*pi*sigma**2))*(exp(-((y-4-1)**2)/(2*sigma**2))+exp(-((y-4+1)**2)/(2*sigma**2)))*(1.0+0.01*sin(pi*(x-1.0))))

# try different initial data
#f.interpolate(conditional(le(((omegaPsq/(4*pi*pi))-cos(2*pi*(x-1.0))-0.5*(y-4.0)**2),0),0,sqrt(((omegaPsq/(4*pi*pi))-cos(2*pi*(x-1.0))-0.5*(y-4.0)**2))))
#f.interpolate(conditional(le(y-4.0,0),0,y-4.0))

driftvel = Function(V3)

#project charge density
U = FunctionSpace(mesh2d, 'CG', 1, vfamily='R', vdegree=0)
g = Function(U, name='g')
g.project(f)

g_norm = Function(U)
g_norm.interpolate(g-0.125)*8.0

# TRIALCODE check normalization
norm = assemble(g_norm*dx)
print("norm="+str(float(norm))+"\n")

# TRIALCODE check init data
File("two-stream_instability_SU_repel_init.pvd").write(f, g_norm)
#quit()

Lphi = inner(grad(phi),grad(v2))*dx
Rphi = omegaPsq*g*v2*dx
phi_s = Function(V2)

h = 2.0/meshres

# is just advection + SU term
F = Dt(f)*v1*dx \
  - v1*div(f*driftvel)*dx \
  + 0.5*h*(dot(driftvel, grad(f)))*dot(driftvel, grad(v1))*(1/sqrt((driftvel[0])**2+(driftvel[1])**2))*dx \

# params taken from Cahn-Hilliard example
# https://www.firedrakeproject.org/Irksome/demos/demo_cahnhilliard.py.html
params = {'snes_monitor': None, 'snes_max_it': 100,
          'snes_linesearch_type': 'l2',
          'ksp_type': 'preonly',
          'pc_type': 'lu', 'mat_type': 'aij',
          'pc_factor_mat_solver_type': 'mumps'}

stepper = TimeStepper(F, butcher_tableau, t, dt, f, solver_parameters=params)

# this is intended to be direct solver
linparams = {"mat_type": "aij",
          "snes_type": "ksponly",
          "ksp_type": "preonly",
          "pc_type": "lu"}

outfile = File("two_stream_instability_SU_repel.pvd")
energydatfile = open("two_stream_instability_SU_repel_energy.txt", "w")

cnt=0
start = time.time()

while float(t) < float(T):
    if (float(t) + float(dt)) >= T:
        dt.assign(T - float(t))
    g.project(f)
    Rphi = omegaPsq*(g-0.125)*8.0*v2*dx

    solve(Lphi==Rphi, phi_s, solver_parameters=params)
    driftvel.interpolate(as_vector([-(y-4.0), grad(phi_s)[0]]))
    energy = assemble(0.0625*(inner(grad(phi_s),grad(phi_s)))*dx)  # 1/2 cos energy growth rate is twice mode growth rate; 1/8 to cancel extent of velocity space
    energydatfile.write(str(float(t))+", "+str(energy)+"\n")

    if(cnt % 40 == 0):
       print("outputting data ...\n")
       outfile.write(f, phi_s, g)
    cnt=cnt+1
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))

end = time.time()
wall_time = end-start

energydatfile.close()

print("done.")
print("\n")
print("wall time:"+str(wall_time)+"\n")
