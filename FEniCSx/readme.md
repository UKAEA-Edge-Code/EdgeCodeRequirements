[Home](../readme.md)
# FEniCSx

## Overview
Description from the [FEniCSx project page](https://fenicsproject.org/):

> "FEniCS is a popular open-source computing platform for solving partial differential equations (PDEs) with the finite element method (FEM). FEniCS enables users to quickly translate scientific models into efficient finite element code. With the high-level Python and C++ interfaces to FEniCS, it is easy to get started, but FEniCS offers also powerful capabilities for more experienced programmers. FEniCS runs on a multitude of platforms ranging from laptops to high-performance computers."

 
## Notes

### Installation

- Conda

```bash
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista
```

## Links

- [Instructions](https://fenicsproject.org/download/)


FEniCSx components:
- [basix](https://github.com/FEniCS/basix): "a finite element definition and tabulation runtime library"
- [dolfinx](https://github.com/fenics/dolfinx): "the computational environment of FEniCSx... implements the FEniCS Problem Solving Environment in C++ and Python"
- [ffcx](https://github.com/fenics/ffcx):  "a compiler for finite element variational forms" (UFL => C)
- [UFL](https://github.com/fenics/ufl) : "a flexible interface for choosing finite element spaces and defining expressions for weak forms in a notation close to mathematical notation"