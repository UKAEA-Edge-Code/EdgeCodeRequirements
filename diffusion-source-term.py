#!/usr/bin/env python
#
# Obtain \nabla \cdot (\kappa \nabla u^{\epsilon}) symbolically and convert to callable function.
# Assumed \kappa_{\para} = 1, \kappa_{\perp} = \epsilon (here "eps").
from sympy import *
x, y, alpha, m, epsilon = symbols('x y alpha m epsilon') # a => α, m => m, eps => ϵ in Deluzet and Narski.
init_printing(use_unicode=True)

## unit vector components and u^eps (symbolic)
bxh = (alpha*(2*y-1)*cos(m*pi*x) + pi)  / sqrt( (alpha*(2*y-1)*cos(m*pi*x)+ pi)**2 + (pi*alpha*m*(y**2 - y)*sin(m*pi*x))**2 )
byh = pi*alpha*m*(y**2 - y)*sin(m*pi*x) / sqrt( (alpha*(2*y-1)*cos(m*pi*x)+ pi)**2 + (pi*alpha*m*(y**2 - y)*sin(m*pi*x))**2 )
ueps  = sin(pi*y + alpha*(y**2 - y)*cos(m*pi*x)) + epsilon*cos(2*pi*x)*sin(pi*y)

## \vec{F} = \kappa \nabla u^{\epsilon} (symbolic)
Fx = ((1-epsilon)*bxh**2 + epsilon)*diff(ueps,x) +        (1-epsilon)*bxh*byh*diff(ueps,y) # x-component of \kappa grad(u^eps)
Fy =        (1-epsilon)*bxh*byh*diff(ueps,x) + ((1-epsilon)*byh**2 + epsilon)*diff(ueps,y) # y-component " " "
SourceTerm = diff(Fx,x) + diff(Fy,y) # divergence of F (still symbolic)

## Callable function
source_term_DN   = lambdify((x,y,alpha,m,epsilon),SourceTerm, 'numpy') # convert symbols to (numpy flavoured) lambda function


## Slight rearrangement alternative
bbdotgradu = bxh*diff(ueps,x) + byh*diff(ueps,y)
Fx_2 = bxh*bbdotgradu + epsilon * (diff(ueps,x) - bxh*bbdotgradu)
Fy_2 = byh*bbdotgradu + epsilon * (diff(ueps,y) - byh*bbdotgradu)
SourceTerm_2 = diff(Fx_2,x) + diff(Fy_2,y)
source_term_DN_2 = lambdify((x,y,alpha,m,epsilon),SourceTerm_2, 'numpy')

print(cxxcode(SourceTerm)) # not so pretty

## Application in plotting the source term on unit square (40,000 points)
import numpy as np
import matplotlib.pyplot as plt
x=np.linspace(0,1,200)
X,Y=np.meshgrid(x,x)
plt.contourf(source_term_DN(X,Y,2.0,10.0,1e-3),35)
# plt.contourf(source_term_DN(X,Y,2.0,10.0,1e-5)-source_term_DN_2(X,Y,2.0,10.0,1e-5),35)
plt.colorbar()
plt.show()
