[Home](../readme.md)
# DUNE

## Overview

Description from the [DUNE project page](https://www.dune-project.org/):

> "DUNE, the Distributed and Unified Numerics Environment is a modular toolbox for solving partial differential equations (PDEs) with grid-based methods. It supports the easy implementation of methods like Finite Elements (FE), Finite Volumes (FV), and also Finite Differences (FD)."


## Review

The comments in this section apply to the DUNE-fem module and its Python bindings (`dune-fempy`) unless otherwise specified.
Version 2.10.dev20240316 was reviewed, which was the latest version that installed successfully with pip at the time of testing.

Overall:
- Lacks some key features, including support for GPUs, particles and higher dimensions.
- Documentation isn't great.
- Small development community is a concern.
- Encountered lots of quirks and hard-to-track-down bugs when coding up the examples.
- Overall impression is that doing anything more complex would uncover similar difficulties, tempering the advantage of the UFL interface.

### Community
Summary: The DUNE-fem community is small (a handful of developers) and makes relatively few releases. The DUNE core community is larger (a few tens of people?), with ~1-2 releases per year.

Pros:
- Module-specific mailing lists exist.
- There's an annual user meeting (covers all DUNE modules).

Cons:
- Couldn't find much evidence of user support online (forum responses etc.).
- dune-fem development team appears to be relatively small: 3 regular contributors over the last year or so (all are also contributors to DUNE core modules).
- Not aware of any experienced DUNE users in-house or among UKAEA collaborators.

### Ease of development
Summary: UFL-based interface, so users get the usual headstart with assembling matrices. Many simple things didn't work perfectly in practice.  e.g. Periodic BCs were almost undocumented, had weird quirks; couldn't get anything working with DG (via `dune-fem-dg`) at all.

Pros:
- The "development-with-Python, production-with-C++" philosophy is potentially a good one (Python bindings follow C++ code stucture).
- UFL-based (rapid development).
- Relatively small list of dependencies.

Cons:
- UFL-based (hard to hack?)
- Writing solvers for the example problems was far from straightforward in practice.
- spack package exists, but hasn't been updated for several years.
- Handling multiple DUNE modules as dependencies could add some build complexity.

### Documentation
Summary: The documentation and examples that exist for dune-fem are useful and well written, but coverage is patchy. Similarly, tutorials/examples are detailed, but there aren't many of them. Documentation for the DUNE core modules seems to be limited to the [Doxygen](https://dune-project.org/doxygen/2.9.1/index.html), which isn't that easy to navigate.

Pros:
-  dune-fem examples are detailed and focussed.

Cons:
- No tutorials for DUNE core modules, only for dune-fem and dune-pdelab.
- Very hard to track down guidance for some functionality (e.g. periodic BCs).
- Not always obvious which module contains which features, and therefore where to look for help.

### Performance
Summary: No GPU support, matrix-free methods look to be restricted to preconditioners.

Pros:
- Matrix-free preconditioners.

Cons:
- No GPU support.
- No matrix-free operators?


### Numerics and Specific features

Summary: Claims to handle CG and DG (but the latter was difficult to get working in practice). Lacks some other important features, particularly higher dimension support, particles. Range of timesteppers might be quite limited.

Pros:
- Just-in-time compiling of C++ code in Python apps allows easy testing of production code in productivity framework.
- On-the-fly plotting of fields helps with debugging.
- Supports CG and DG.
- Supports elliptic solves.
- The dune-grid interface allegedly supports non-conforming meshes (e.g. section 2. of [this pdf](https://wrap.warwick.ac.uk/146797/2/WRAP-The-DUNE-framework-basic-concepts-recent-development-Dedner-2020.pdf)).

Cons:
- Proprietary mesh format [DuneGridFormat(.dgf)](https://dune-project.org/doxygen/master/group__DuneGridFormatParser.html). Conversion from gmsh format is possible though.
- No support for higher dimensional problems.
- List of timesteppers looks fairly limited.
- No particles.

## Examples

For problem descriptions and equations see the [Firedrake examples readme](../Firedrake/examples/readme.md).

### CG advection ([simple_advect_periodic_CG.py](examples/simple_advection/simple_advect_periodic_CG.py))

<img src="examples/simple_advection/postproc/DUNE_advection_CG.gif" width="600">

More-or-less working (at least with Gaussian ICs), but:
- Result is out by half a unit (!?)
- Needs a much smaller timestep than used in implementations based on other frameworks
- Bizarre bug with the periodic boundary if run **without** MPI - appears to be a gap or ghost region between low-x and high-x boundaries. Can't find a solution in the docs (Or much about PBCs at all...).
### DG advection ([simple_advect_periodic_DG.py](examples/simple_advection/simple_advect_periodic_DG.py))

Unfinished.

### Anisotropic Diffusion ([aniso_diffusion_DeluzetNarski.py](examples/aniso_diffusion/aniso_diffusion_DeluzetNarski.py))

Decent results

<img src="examples/aniso_diffusion/output/aniso_diffa0_m0_4bcs_160x160.png" width="400" style="margin-right: 1.5rem">
<img src="examples/aniso_diffusion/output/aniso_diffa2_m1_4bcs_160x160.png" width="400" style="margin-right: 1.5rem">
<img src="examples/aniso_diffusion/output/aniso_diffa2_m10_4bcs_160x160.png" width="400" style="margin-right: 1.5rem">


### 1D outflow isothermal compressible Euler ([SOL1D_DG.py](examples/1doutflow/SOL1D_DG.py))

Unfinished.

### Non-conforming mesh ([non-conforming.py](examples/non-conformal_mesh/non-conforming.py))

Attempt to read a non-conforming mesh and run a simple problem. Couldn't get this to work (but see review comments re. non-conforming meshes).



## Installation (Ubuntu)

**pip**
- Using dune-fem module 
- core modules easy to install, no sudo needed
- additional (dune-fem, dune-fem-dg etc.) modules easy to install
- *Some* tutorials work out-of-the-box once you have the correct version

<!-- Dependencies
- SciPy
- Eigen (optional)
- PAPI (optional)
- PETSc (optional)
- SIONlib (optional)
- SuiteSparse (optional) -->

### Gripes
- Suggested (apt-based) installation on Ubuntu didn't work, at least not easily - couldn't find the dune-pdelab package, which seemed to be essential
- dune.fem tutorials require the latest development version, didn't work with `pip install dune.fem`

## Links

- [Installation](https://www.dune-project.org/doc/installation/)
- [Core modules on GitLab](https://gitlab.dune-project.org/core/)
- [DUNE-fem on GitLab](https://gitlab.dune-project.org/dune-fem/dune-fem/)