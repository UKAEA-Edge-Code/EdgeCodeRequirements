[Home](../readme.md)
# Trixi

## Overview

Description from the [Trixi project page](https://trixi-framework.github.io/Trixi.jl/stable/overview/):

> "Trixi.jl is designed as a library of components for discretizations of
> hyperbolic conservation laws. Trixi.jl uses the method of lines, i.e., the
> full space-time discretization is separated into two steps; the spatial
> semidiscretization is performed at first and the resulting ODE system is
> solved numerically using a suitable time integration method. "

## Notes

- Visualisation with Plots.jl, Makie.jl ("experimental") and Paraview, via Trixi2Vtk.
- Support for various static/adaptive(AMR) meshes in 1D,2D and 3D
   - DG\[SEM\] possible with all.
   - gmsh import only for straight-sided elements or cells
   - MPI limited to some.
- Multi-threading and MPI should be possible (but should test for suitable problems).
- High order integration schemes with various callback types e.g. solution analysis, timestep control, AMR, visualisation.
- Easy to couple with other frameworks, python, fortran, C++ etc. (compiled/partially_compiled julia possible)
- Additional static variables, e.g. magnetic field need to be kept together solution with present semidiscretizations.

### Installation
```julia
julia> using Pkg

julia> Pkg.add(["Trixi", "Trixi2Vtk", "OrdinaryDiffEq", "Plots"])
````
## Links

- [GitHub page](https://github.com/trixi-framework/Trixi.jl)(std)

## Experience with simple test problems

Trixi is primarily designed to solve equations of the form
$$
\frac{\partial {\bf U}}{\partial t} + \nabla \cdot {\bf F_{\rm adv}} = {\bf S}
$$
where ${\bf F_{\rm adv}}$ are the advective fluxes and ${\bf S}$ contains sources. The framework is designed to be used at a high level, so for a new set of equations one does not need to do much more than defining the fluxes and sources if adapting one of the many examples. Parabolic terms can be included with a hybrid semidiscretization that treats hyperbolic and parabolic differently. A semidiscretization
provides a system of ODEs that are solved using the [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl) package that is an interface to a large collection of integration schemes. At present, elliptic equations are solved as time-dependent problems (perhaps a new type of semidiscretization could be added to couple with various solvers).

### Advection test

The advection equations are already implemented as ``LinearScalarAdvectionEquation2D`` and provided as convergence tests. You can provide an advection velocity and construct a structured mesh with ``StructuredMesh``. Unfortunately ``RadauIIA3()`` (the integration scheme provided by [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl)) crashes the code, so `BS3()` is used instead. With purely horizontal advection and a Gaussian of width $2$, at $t=40$ the solution has $L_2$ error of $2.28576555e-04$ and $L_\infty$ error of $2.13218567e-03$.
Changing the advection angle hardly affects the result. The test completes in under $0.5s$ on my desktop. A script is [here](./Trixi/simple_advect_periodic_DG/simple_advect_periodic_DG.jl).



I was able to use the equation abstraction to have 3 scalar variables, so as to include a 2D-vector field. I could then advect along a spatially-varying static field. This is helpful for the diffusion test.

### Diffusion test

**TODO**
- Can't use the gmsh mesh as curved edges and forcing them straight is messy (but might be able to use a similar mesh).
- Need to work out BCs.
  
