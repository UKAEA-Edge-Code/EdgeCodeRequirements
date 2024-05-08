[Home](../readme.md)
# MOOSE

## Overview

Description from the [MOOSE project webpage](https://mooseframework.inl.gov/):

> "The Multiphysics Object-Oriented Simulation Environment (MOOSE) is a finite-element, multiphysics framework primarily developed by Idaho National Laboratory. It provides a high-level interface to some of the most sophisticated nonlinear solver technology on the planet. MOOSE presents a straightforward API that aligns well with the real-world problems scientists and engineers need to tackle."

## Test problems

### Advection

Can be done via config file (see ```advection_dg.i```) - the default Moose app has all required features.

Moose doesn't support periodic bounday conditions for Elemental variables, hence
extended domain doubled size of domain in x direction from [-20,20] to [-60,20]
but element sizes are the same.

![Results at time 0 and time 40](Start_and_End.png)


![Error at time 40](error.png "Error")

#### Caveats
- Second order basis is highest available with quad elements

## Notes

Terrifying statement from the docs: "MOOSE does not use traditional versioning, is under heavy development, and is being updated continuously. Therefore, it is important that you continue to update MOOSE as you use it to develop your application(s); weekly updates are recommended."

Code is C++17, Python only used for postprocessing, data manipulation, docs generation etc.
Intel compilers NOT supported; gcc (>7.5) and clang (>10.0.1) only.

### Parallelism

- MPI + Threading supported
- GPU via PETSc only.
- Supports (in theory) CUDA, KOKKOS, HIP, OpenCL. SYCL ("Not yet supported")
- Have looked at this for another project will share when written up but..
    - Works on NVIDIA (CUDA or Kokkos), can't make PETSc work on Intel GPUs
      (uses SYCL via KOKKOS), never tried with AMD 
    - Seems to be that not enough work is offloaded to GPU to get speed-up with
      A100 (vs sapphire rapids).


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
