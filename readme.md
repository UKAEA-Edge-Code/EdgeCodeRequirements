# Edge code requirements

There are two broad goals of the project:
1. Create a SOLPS-EIRENE replacement, which involves a neutral particle solver coupled to a Braginskii code.
    - This will be a modern replacement of SOLPS that is more productive by being a combination of: faster; better engineered and easier to develop; algorithmically capable; and having more physics.
2. Build a framework for running non-linear drift- and gyro-kinetic simulations based upon the moment-kinetics numerical method.
    - This will enable simulations of large high power tokamaks in H-mode, for the calculation of fluxes that are subsequently fed into the SOLPS replacement.

## Framework evaluation directories

- [deal.II, life<sup>X</sup>](./dealII/readme.md) 
- [DUNE](./DUNE/readme.md)
- [FELTOR](FELTOR/readme.md)
- [FEniCSx](FEniCSx/readme.md)
- [Firedrake](./Firedrake/readme.md)
- [jax, jax-cfd, jax-fem](./jax/readme.md)
- [MFEM, pyMFEM](./MFEM/readme.md)
- [Moose](./MOOSE/readme.md)
- [Nektar++](./nektar++/readme.md)
- [NGSolve](./NGSolve/readme.md)
- [Trixi.jl](./Trixi/readme.md)

## Overview tables

[General info](docs/info.md)

[Overall ratings](docs/ratings.md)

## Framework requirements

Some requirements of the productivity framework are:

| General                                                         | Reason                                                                                 | Firedrake  | FENICs, MOOSE, AMRpX, jax   |
| --------------------------------------------------------------- | -------------------------------------------------------------------------------------- | --------- | ------------------------- |
| Is enjoyable to use                                             | People have to like what they do                                                       | Is not enjoyable to install. Is reasonable once actually installed. | Text |
| Well engineered                                                 | Makes it easier to learn the codebase                                                  | Based on capturing the mathematics in the abstractions/DSLs then the user composes the maths. | Text                      |
| Large user base                                                 | Accessible experience                                                                  | This is hard to quantify.      | Text                      |
| Mature and not fading                                           | Can't choose abandonware                                                               | Is funded by PostDoc/PhD money which comes with inherit instability over funding and quality of candidates.      | Text                      |
| Well documented                                                 | Sensibly hosted online documentation                                                   | copyleft but LGPL allows redistribution under terms  | Text                      |
| Permissive (non copyleft) open source                           | Make it easier to commercialise                                                        | Multiple tutorials available - good slack support.  | Text                      |
| Good tutorials available                                        | Make it easier to learn asynchronously                                                 | The examples are usually very good - usually do not need to dive into tests.      | Text                      |
| Tests are good learning resource                                | Sign of good code                                                                      | Tests exist - and in parallel.      | Text                      |
| Test coverage is good                                           | People have thought about correctness                                                  | Many tests, actual coverage is hard to quantify.  | Text                      |
| In a sensible language(s)                                       | Need to choose a language that can access different                                    | Python DSL frontend - generates C via several intermediate stages.      | Text                      |
| Has an escape hatch to low level interfaces (if applicable)     | Code generation techniques can be a big barrier to implementing unsupported capability | Intermediate levels accessible. Solvers are provided by PETSc. | Text                      |
| Is performance portable                                         | Got to run on different vendors' hardware                                              | Text      | Text                      |
| Low level code is understandable (if applicable)                | C++ or code generation must be legible and understandable                              | Similar level of complexity to other low-level representations.      | Text                      |
| Code generation is not insurmountably difficult (if applicable) | We need to understand how to work with it                                              | AST/Expression tree approach with successive lowering (standard approach).       | Text                      |
| Is differentiable                                               | Nice-to-have                                                                           | yes, see Pyadjoint.      | Text                      |
| Follows versioning best-practices                               | Need a stable API, reproducible builds                                                 | UFL interface stable, lower interfaces subject to change.  | Text                      |
| Has relatively few dependencies                                 | Avoid complicated installation on HPC                                                  | Typical dependencies for FEM frameworks, installing is non-trivial if the installation script fails. Currently does not have spack.      | Text                      |
| Dependencies build quickly                                      | Speed up development, installation on new platforms                                    | Most of the dependency install is PETSc build.  | Text                      |
| Pull requests are welcome and dealt with quickly                                   | Our changes will be made upstream                                    | Text  | Text                      |

| Field solvers                                           | Reason                                                                                   | Firedrake |
| ------------------------------------------------------- | ---------------------------------------------------------------------------------------- | -------- |
| The periodic table of finite elements                   | Need to support FEEC                                                                     | Large selection of elements, including HDiv/HCurl spaces.   |
| Discontinuous and Continuous Galerkin methods           | Need to support CG and DG                                                                | Yes     |
| UFL interface or similar                                | Must have a pleasant API for rapid prototyping                                           | UFL     |
| Distributed MPI-X                                       | Must be distributed by MPI and support threading / GPU etc                               | MPI+Serial with MPI+Loopy (Opencl) in developement. Uses PETSc GPU support for solvers.     |
| Performance portable field solve GPU support            | Run on many vendor GPUs                                                                  | Loopy provides portable loop abstraction, uses PETSc GPU (in development).     |
| Supports velocity space                                 | Have 3D1V and 3D2V i.e. integrals over some dimensions collapse down to the spatial ones | Possibly - maybe ET knows.     |
| Velocity space support possible                         | There is an accessible roadmap to velocity space support, if not supported already       | See above.     |
| Supports arbitrary explicit / implicit time integrators | e.g. irksome                                                                             | Yes Irksome uses Firedrake, supports arbitrary  Butcher tableau |

| Particle support                                    | Reason                                                                             | Firedrake |
| --------------------------------------------------- | ---------------------------------------------------------------------------------- | -------- |
| Already has GPU particle capability                 | Future FLOPS will increasingly come from GPUs                                      | No     |
| GPU particle capability does not require a re-write | GPU support is available in a short amount of time i.e. 1-2 years, not 5-10 years. | Probably     |
| Can find cell belonging to point                    | Basic requirement for particles                                                    | Depends, linear mesh probably. |
| Can evaluate basis functions in cell                | Basic requirement for particles                                                    | Almost certainly |
| Projection code exists                              | e.g. particle charge is converted to a charge density field                        | Yes, but not fast for particle use case, see VertexOnlyMesh  |


| Performance                                              | Reason                                                                                      | Option 1 |
| -------------------------------------------------------- | ------------------------------------------------------------------------------------------- | -------- |
| Scales well                                              | We don't want to spend money getting it into the middle of tha pack in terms of scalability | Scaling is mostly determined by PETSc solver choice which will be common to FEM using PETSc.     |
| Uses well known parallelisation techniques               | We want to use what others use                                                              | Standard domain decomposition.     |
| Has access to well known and easy to use preconditioners | Important for accuracy and speed                                                            | Uses PETSc, other Firedrake users develop preconditioners in Firedrake. |

| Gridding                                     | Reason                                                 | Firedrake        |
| -------------------------------------------- | ------------------------------------------------------ | --------------- |
| Supports multiple gridding inputs            | Flexibility and likelihood to continue growing support | Has some capability to read standard formats natively. Also has a netgen interface. |
| Supports unstructured meshes                 | To have realistic geometries                           | Yes            |
| Supports hanging nodes for DG                | To support FCI                                         | TBC            |
| Can use different elements in one simulation | Tets+hexes+prism often found together                  | TBC (probably) |
| h-refinement                                 | Nice-to-have                                           | TBC           |
| p-refinement                                 | Nice-to-have                                           | In development?     |
| r-refinement                                 | Nice-to-have                                           | Yes - the mesh is a function and r-refinement is actively used, see FireShape.  |


| IO                             | Reason                                      | Firedrake |
| ------------------------------ | ------------------------------------------- | -------- |
| Checkpoints function           | Req. for HPC                                | Yes     |
| In-situ visualisation          | Extract more physics without saving to disk | TBC     |
| Supports standard file formats | HDF5, NetCDF                                | HDF5     |


| Equation terms                                            | Reason | Option 1 |
| --------------------------------------------------------- | ------ | -------- |
| Source terms                                              | Text   | Text     |
| Divergence                                                | Text   | Text     |
| Curl                                                      | Text   | Text     |
| Advection                                                 | Text   | Text     |
| Diffusion                                                 | Text   | Text     |
| Anisotropic diffusion                                     | Text   | Text     |
| Parallel diffusion                                        | Text   | Text     |
| Perpendicular diffusion                                   | Text   | Text     |
| Hyper diffusion                                           | Text   | Text     |
| Poisson bracket                                           | Text   | Text     |
| Time evolving BCs                                         | Text   | Text     |
| BCs that depend on solver variables' values and gradients | Text   | Text     |
| BCs that depend on geometry  | Gradients in arbitrary directions wrt grid  | Text     |

| Use cases                                      | Reason | Option 1 |
| ---------------------------------------------- | ------ | -------- |
| Has a history of UQ applications               | Text   | Text     |
| Adjoint methods have been successfully applied | Text   | Text     |
| Has provenance capture built in for data, config etc | Text   | Text     |
