# Simple tests

# 1. Simple advection (simple_advect_CG_fenicsx.py, simple_advect_periodic_CG_fenicsx.py)

Advection of scalar field at unit velocity on a periodic domain $40 \times 10$ length units.  Advect left for $40$ time units.\
Discretization: $64 \times 16$ quadrilaterals, order-3 Lagrange polynomial basis functions, continuous Galerkin.\
Boundary conditions: all-Neumann (I think) or Dirichlet (top/bot) + periodic (left/right).\
Time-stepper: First-order backward Euler.\
Initial data: either Gaussian $n=e^{-\frac{r^2}{s^2}}$ centered in the domain.\
Analytic solution: obviously just the uniformly-advected initial data.\
Other notes: TODOs include adding better time-stepper, L2 error evaluation, DG version.

# Slab anistropic diffusion - Deluzet-Narski singular perturbation problem (anisotropic_diffusion_DeluzetNarski_fenicsx.py)

Anisotropic diffusion problem in square domain, see notes for the Firedrake version.  The script is a straight port of the Firedrake one. \

Discretization: triangle meshes generated by gmsh, edit square.geo to set refinement $h=0.1, 0.05, ... 0.00078125$ as used in the paper (finest mesh is 3.8M elements, fine on laptop), then "gmsh -2 square.geo".  Can also edit script to generate its own quad mesh instead.  First-order CG Lagrange elements.\
Boundary conditions: homogeneous Dirichlet top / bottom, homogeneous Neumann left / right. \
Time-stepper / initial data: none (elliptic problem).\
Analytic solution: $n=\sin(\pi y + \alpha (y^2-y) \cos m \pi x) + \epsilon \cos 2 \pi x \sin \pi y$.\
Other notes: TODO add L2 error evaluation, improve plotting of output (limited to first order so looks unnecessarily poor on a low res mesh at higher-order).
