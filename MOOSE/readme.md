[Home](../readme.md)
# MOOSE

## Overview

Description from the [MOOSE project webpage](https://mooseframework.inl.gov/):

> "The Multiphysics Object-Oriented Simulation Environment (MOOSE) is a finite-element, multiphysics framework primarily developed by Idaho National Laboratory. It provides a high-level interface to some of the most sophisticated nonlinear solver technology on the planet. MOOSE presents a straightforward API that aligns well with the real-world problems scientists and engineers need to tackle."


## Notes

Terrifying statement from the docs: "MOOSE does not use traditional versioning, is under heavy development, and is being updated continuously. Therefore, it is important that you continue to update MOOSE as you use it to develop your application(s); weekly updates are recommended."

Code is C++17, Python only used for postprocessing, data manipulation, docs generation etc.
Intel compilers NOT supported; gcc (>7.5) and clang (>10.0.1) only.


### Dependencies

- PETSc (numerous deps - BLAS, LAPACK, PARMETIS, PTSCOTCH, STRUMPACK , HYPRE, SLEPC, etc.)
- libMesh
- WASP

### Installation

- Has conda install (root not required), but probably not suitable for dev
  - Install miniForge (900 MB)
  - Add INL conda channel
  - Install MOOSE
- From src (create spack package later?)
  - Need to provide MPI, Boost
  - Seems to build petsc,libmesh, WASP ok, but test build failed with missing pyyaml

## Links

- [Examples and Tutorials](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/index.html)
- [GitHub](https://github.com/idaholab/moose)
- [Installation](https://mooseframework.inl.gov/getting_started/installation/index.html)