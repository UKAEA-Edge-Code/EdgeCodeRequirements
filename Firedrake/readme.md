[Home](../readme.md)
# Firedrake

## Overview
Description from the [Firedrake project page](https://www.firedrakeproject.org/):

> "Firedrake is an automated system for the solution of partial differential equations using the finite element method (FEM). Firedrake uses sophisticated code generation to provide mathematicians, scientists, and engineers with a very high productivity way to create sophisticated high performance simulations."


## Notes

- Tight integration with PETSc - effectively write callback functions for it in C, auto-generated from user Python
- No postproc tools at first glance, just generates vtk and relies on ParaView.
- Surprisingly difficult to infer from the docs that it's doing code gen...
- Parallelism all handled under-the-hood; 

### Installation

System install of several dependencies...

- PETSc

## Links

- [firedrakeproject/ufl](https://github.com/firedrakeproject/ufl): A fork of [FEniCS/UFL](https://github.com/fenics/ufl).
- [GitHub page](https://github.com/firedrakeproject/firedrake)
- [Installation](https://www.firedrakeproject.org/download.html)
- [Irksome](https://github.com/firedrakeproject/Irksome): Time-stepping package used to generate Runge Kutta methods for Firedrake

## Example scripts - benchmark problems

Please see repo (https://github.com/ethrelfall/FEM-framework-comparison/tree/main).
