# Simple tests

# 1. Simple advection (simple_advect_CG_fenicsx.py)

Advection of scalar field at unit velocity on a periodic domain $40 \times 10$ length units.  Advect left for $40$ time units.\
Discretization: $64 \times 16$ quadrilaterals, order-3 Lagrange polynomial basis functions, continuous Galerkin.\
Boundary conditions: all-Neumann (I think).\
Time-stepper: First-order backward Euler.\
Initial data: either Gaussian $n=e^{-\frac{r^2}{s^2}}$ centered in the domain.\
Analytic solution: obviously just the uniformly-advected initial data.\
Other notes: TODOs include adding better time-stepper, L2 error evaluation, DG version.
