# aniso_diffusion_DeluzetNarski.py

# What this does: tests straight- and curved-fieldlines cases of anisotropic diffusion
# after "A two field iterated AP method for highly anisotropic elliptic equations"
# by Deluzet and Narski.
# This shows the issues explored in the paper:
# 1. Null space problem for eps -> 0.
# 2. "Locking" problem for order-1 (comes if basis functions have no cpt in the direction of the B-field): solution amplitude comes out much smaller than should be.
# NOTE modified to use actual MMS solution from the paper
# NOTE is set up to use unstructured triangles mesh from file, can also choose Cartesian quad mesh but need to change the BC labels for this
# NOTE triangle mesh version tested to highest res used in the paper, h = 0.00078125 / 3.8M elements, fine on laptop

from firedrake import *

# MODEL PARAMETERS
alpha = Constant(0.0)  # parameter for varying B-field away from x-direction; 10 gives reasonable bent fieldlines, 0 is x-direction only
eps = Constant(1.0e-1)  # anisotropy parameter, 1.0 is isotropic, intended to be decreased to small positive values
m = 0
# Note if changing mesh need to make sure indices in bcB and bcT are correct - see comment.
# END OF MODEL PARAMETERS

#square mesh
#meshres=50
#import numpy as np
#xcoords = np.linspace(-0.5, 0.5, meshres + 1, dtype=np.double)
#ycoords = np.linspace(-0.5, 0.5, meshres + 1, dtype=np.double)
#mesh = TensorRectangleMesh(xcoords, ycoords, quadrilateral="True")

mesh = Mesh("square.msh")  # triangle mesh of unit square [0,1]x[0,1]

V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

x, y = SpatialCoordinate(mesh)

bcB = DirichletBC(V, 0.0, 11)   # y = 0 homogeneous Dirichlet, 11 on unit_square.msh, 3 on TensorRectangleMesh
bcT = DirichletBC(V, 0.0, 13)   # y = 1 homogeneous Dirichlet, 13 on unit_square.msh, 4 on TensorRectangleMesh

solution = Function(V)
solution.interpolate(sin(pi*y+alpha*(y**2-y)*cos(m*pi*x))+eps*cos(2*pi*x)*sin(pi*y))
bcLR = DirichletBC(V, solution, (1,2))  # optional Dirichlet BC for lam=0.0 case: removes null space problem

# magnetic field lines direction unit vector
norm_denom = sqrt((alpha*(2*y-1)*cos(m*pi*x)+pi)**2+(pi*alpha*m*(y**2-y)*sin(m*pi*x))**2)
bhat = as_vector([(alpha*(2*y-1)*cos(m*pi*x)+pi)/norm_denom, (pi*alpha*m*(y**2-y)*sin(m*pi*x))/norm_denom])

k_par = Constant(1.0)
k_per = eps  # THIS IS EPSILON IN ANISO DIFFUSION TENSOR

flux = k_par * bhat * dot(bhat, grad(u)) + k_per * (grad(u) - bhat * dot(bhat, grad(u)))

a = inner(flux, grad(v))*dx
f = Function(V)

# got to work out f from MMS given the solution in "solution" above
# algebra is a total PITA so let's use a high-order function space to do it numerically ...
# can easily be tested for straight fieldlines case (MMS is easy algebra in that case)
V_S = FunctionSpace(mesh, "CG", 10)
solution_S = Function(V_S)
solution_S.interpolate(sin(pi*y+alpha*(y**2-y)*cos(m*pi*x))+eps*cos(2*pi*x)*sin(pi*y))
flux_S = k_par * bhat * dot(bhat, grad(solution_S)) + k_per * (grad(solution_S) - bhat * dot(bhat, grad(solution_S)))
f.interpolate(-div(flux_S))

L = inner(f,v)*dx

T = Function(V)

solve( a==L, T, bcs=[bcB, bcT])  # optionally add bcLR to remove null space problem and locking problem

# output magnetic field lines
VB=VectorFunctionSpace(mesh,"CG", 1)
BField = Function(VB)
BField.interpolate(bhat)
File("aniso_diffusion_DeluzetNarski_BField.pvd").write(BField)

T_benchmark = sin(pi*y+alpha*(y**2-y)*cos(m*pi*x))+eps*cos(2*pi*x)*sin(pi*y)
L2_error = norm(T-T_benchmark, "L2")
print("L2 error norm: %.2e" % L2_error)

T.rename("numerical solution")
solution.rename("analytic solution")

File("aniso_diffusion_DeluzetNarski.pvd").write(T, solution)
