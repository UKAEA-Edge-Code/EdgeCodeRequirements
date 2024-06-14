[Home](../readme.md)
# MFEM

## Overview

Description from the [MFEM project page](https://mfem.org/):

> MFEM is a free, lightweight, scalable C++ library for finite element methods.

It provides fairly low-level tools to construct solvers once the week
form has been determined. A decent amount of code will be required to
assemble the solver, but it is not overly complicated. Python wrappers
are available.

MFEM is released under the BSD-3 license.

## Advection Example
An example of an advection solver using PyMFEM is provided in
`[advection.py](./advection.py)`. The simulation can be varied using
the following command-line options:

```
usage: advection.py [-h] [-vel VELOCITY] [-a ANGLE] [-rad RADIUS] [-mag MAGNITUDE]
                    [-fr FRANKENSTEIN_PARAM] [-vis] [-o ORDER] [-d DEVICE]
                    [-tf T_FINAL] [-dt TIME_STEP] [-nx X_RESOLUTION]
                    [-ny Y_RESOLUTION] [-sx X_SIZE] [-sy Y_SIZE]
                    [-vs VISUALIZATION_STEPS] [-paraview]

Simple advection app

options:
  -h, --help            show this help message and exit
  -vel VELOCITY, --velocity VELOCITY
                        The magnitude of the advection velocity.
  -a ANGLE, --angle ANGLE
                        The angle of the velocity field relative to the x-axis, in
                        degrees.
  -rad RADIUS, --radius RADIUS
                        Radius of the perturbation being advected.
  -mag MAGNITUDE, --magnitude MAGNITUDE
                        Magnitude of the perturbation being advected.
  -fr FRANKENSTEIN_PARAM, --frankenstein-param FRANKENSTEIN_PARAM
                        The 'a' parameter in the 'Frankenstein's monster' type
                        function for the perturbation. If 0 then a guassian will
                        be used instead.
  -vis, --visualization
                        Enable GLVis visualization
  -o ORDER, --order ORDER
                        Order (degree) of the finite elements.
  -d DEVICE, --device DEVICE
                        Device configuration string, see Device::Configure().
  -tf T_FINAL, --t-final T_FINAL
                        Final time; start time is 0.
  -dt TIME_STEP, --time-step TIME_STEP
                        Time step.
  -nx X_RESOLUTION, --x-resolution X_RESOLUTION
                        Number of elements in the x-direction.
  -ny Y_RESOLUTION, --y-resolution Y_RESOLUTION
                        Number of elements in the y-direction.
  -sx X_SIZE, --x-size X_SIZE
                        The length of the domain in the x-direction.
  -sy Y_SIZE, --y-size Y_SIZE
                        The length of the domain in the y-direction.
  -vs VISUALIZATION_STEPS, --visualization-steps VISUALIZATION_STEPS
                        Visualize every n-th timestep.
  -paraview, --paraview-datafiles
                        Save data files for ParaView (paraview.org) visualization.
```

## Diffusion Solver
An anisotropic diffusion solver was written in C++ with MFEM in
[diffusion.cpp](./diffusion.cpp). It requires MFEM to already have
been built and installed, but for the source still to be available as
well. It can be compiled with the command `make
MFEM_SOURCE_DIR=/path/mfem/source`. The simulation can be varied using
the following command-line options:

```
Usage: ./diffusion [options] ...
Options:
   -h, --help
	Print this help message and exit.
   -m <string>, --mesh <string>, current value: square.msh
	Mesh file to use.
   -o <int>, --order <int>, current value: 1
	Finite element order (polynomial degree) or -1 for isoparametric space.
   -sc, --static-condensation, -no-sc, --no-static-condensation, current option: --no-static-condensation
	Enable static condensation.
   -pa, --partial-assembly, -no-pa, --no-partial-assembly, current option: --no-partial-assembly
	Enable Partial Assembly.
   -fa, --full-assembly, -no-fa, --no-full-assembly, current option: --no-full-assembly
	Enable Full Assembly.
   -eps <double>, --epsilon <double>, current value: 0.1
	Ratio between perpendicular and parallel diffusivity.
   -a <double>, --alpha <double>, current value: 2
	Alpha parameter for magnetic field.
   -mp <double>, --m-param <double>, current value: 1
	m parameter for magnetic field.
   -d <string>, --device <string>, current value: cpu
	Device configuration string, see Device::Configure().
   -vis, --visualization, -no-vis, --no-visualization, current option: --visualization
	Enable or disable GLVis visualization.
   -paraview, --paraview-datafiles, -no-paraview, --no-paraview-datafiles, current option: --no-paraview-datafiles
	Save data files for ParaView visualization.
```

## General Impressions

The documentation, tutorials, and examples for MFEM are pretty
good. It is much easier to understand how to write a solver than for
Nektar++. API docs could be better. PyMFEM documentation is poorer,
consisting only of text files in the repository.

Implementing linear systems of equations is pretty straightforward, as
there are extensive bilinear operators already implemented. Nonlinear
operators look like they require more work; only a few nonlinear terms
are included in MFEM so the author is likely to need to write these
themselves. I have not tried this and don't know how difficult it is.

The basic design of the code seems sound, although there are some
areas that could be tidied up for greater consistency. The design is
fairly typical of C++, with lots of getter/setter methods and
mutability.

MFEM appears to be under active development, with pull requests being
opened and approved regularly. However, I've been told that it often
takes about a month for a PR to be approved.

There is currently support for CUDA and HIP. A PR exists to add SYCL
support, but it was not merged in time for the most recent release and
is now very far behind the `main` branch.

Hanging/nonconformal nodes are supported, but this is in the context
of mesh refinement. It appears the **hanging nodes need to be on an edge
which terminates in non-hanging nodes**. This means the field-aligned
non-conforming meshes generated by FAME will not be supported. It is
unclear how difficult it would be to add support for this.

MFEM only seems to support up to 3 dimensions. It does not look like
velocity space is supported.

There is a "miniapp" providing automatic differentiation. However, it
is not particularly easy to use and is not integrated into the core of
MFEM. When the miniapps are installed by `spack`, the automatic
differentiation headers end up containing broken include paths. There
are likely better 3rd party library for this.

There are unit tests but no information on coverage that I can
see. There do not appear to be integration tests. The examples could
in principle act as regression tests but they don't appear to be run
as part of CI.

The PyMFEM bindings are fairly low-level and not entirely
Pythonic. It's possible to end up with dangling pointers if the owner
of an object gets garbage-collected, leading to segfaults. This can
lead to difficult bugs.

Working with MFEM in C++ comes with the usual headaches of that
language. The containers provided by MFEM do not provide
bounds-checking, so memory errors can happen easily and lead to
strange/unreproducible behaviour. Raw pointers are used surprisingly
often, which could lead to problems.

## Installation

MFEM can be installed straightforwardly with `spack`. Installing the
Python wrappers can be done through pip but they will only work if you
have the same version of NumPy as the pre-compiled package was built
against. In principle it should be possible to have pip build this
itself, although this hasn't been tested. I have written a first draft
of a `spack` package to install it (`[package.py](./package.py)`),
which takes around 15 minutes. It seems to be very particular about
which versions of MFEM it builds against and many of the PyMFEM
releases build against particular commits of MFEM that are not a
tagged release.

## Links

- [GitHub Page](https://github.com/mfem/mfem)
- [PyMFEM](https://github.com/mfem/PyMFEM): Python wrapper
- [Installation](https://github.com/mfem/mfem/blob/master/INSTALL)
