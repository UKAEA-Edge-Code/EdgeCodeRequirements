# Edge code requirements

There are two broad goals of the project:
1. Create a SOLPS-EIRENE replacement, which involves a neutral particle solver coupled to a Braginskii code.
    - This will be a modern replacement of SOLPS that is faster, better engineered, and ultimately more flexible offering more physical interpretation.
3. Build a framework for running non-linear drift- and gyro-kinetic simulations based upon the moment-kinetics numerical method.
    - This will enable simulations of large high power tokamaks in H-mode, for the calculation of fluxes that are subsequently fed into the SOLPS replacement.


## Solver framework

| General | Reason | Option 1 |
| -------- | -------- | -------- |
| Is enjoyable to use | People have to like what they do    | FireDrake, FENICs, MOOSE, AMRpX |
| Well engineered | Makes it easier to learn the codebase     | Text     |
| Large user base | Accessible experience     | Text     |
| Mature and not fading | Can't choose abandonware     | Text     |
| Well documented | Sensibly hosted online documentation     | Text     |
| Permissive (non copyleft) open source | Make it easier to commercialise  | Text     |
| Good tutorials available | Make it easier to learn asynchronously   | Text     |
| Tests are good learning resource | Sign of good code     | Text     |
| Test coverage is good | People have thought about correctness     | Text     |
| In a sensible language(s) | Need to choose a language that can access different     | Text     |
| Has an escape hatch to low level interfaces (if applicable) | Code generation techniques can be a big barrier to implementing unsupported capability     | Text     |
| Is performance portable | Got to run on different vendors' hardware     | Text     |
| Low level code is understandable (if applicable) | C++ or code generation must be legible and undertandable  | Text     |
| Code generation is not insurmountably difficult (if applicable) | We need to understand how to work with it    | Text     |
| Is differentiable | Nice-to-have  | Text     |

| Field solvers | Reason | Option 1 |
| -------- | -------- | -------- |
| The periodic table of finite elements | Need to support FEEC   | Text     |
| Discontinuous and Continuous Galerkin methods | Need to support CG and DG  | Text     |
| UFL interface or similar | Must have a pleasant API for rapid prototyping  | Text     |
| Distributed MPI-X | Must be distributed by MPI and support threading / GPU etc   | Text     |
| Performance portable field solve GPU support | Run on many vendor GPUs     | Text     |
| Supports velocity space | Have 3D1V and 3D2V i.e. integrals over some dimensions collapse down to the spatial ones    | Text     |
| Velocity space support possible | There is an accessible roadmap to velocity space support, if not supported already     | Text     |
| Supports arbitrary explicit / implicit time integrators   | e.g. irksome    | text |

| Particle support | Reason | Option 1 |
| -------- | -------- | -------- |
| Already has GPU particle capability | Future FLOPS will increasingly come from GPUs     | Text     |
| GPU particle capability does not require a re-write | GPU support is available in a short amount of time i.e. 1-2 years, not 5-10 years.      | Text     |
| Can find cell belonging to point | Basic requirement for particles     | Text     |
| Can evaluate basis functions in cell | Basic requirement for particles     | Text     |
| Projection code exists | e.g. particle charge is converted to a charge density field  | Text     |


| Performance | Reason | Option 1 |
| -------- | -------- | -------- |
| Scales well | We don't want to spend money getting it into the middle of tha pack in terms of scalability     | Text     |
| Uses well known parallelisation techniques | We want to use what others use  | Text     |

| Gridding | Reason | Option 1 |
| -------- | -------- | -------- |
| Supports multiple gridding inputs | Flexibility and likelihood to continue growing support   | For ease of use    | Text |
| Supports unstructured meshes   | To have realistic geometries     | Text |
| Supports hanging nodes for DG   | To support FCI     | Text |
| Can use different elements in one simulation | Tets+hexes+prism often found together   | Text     | Text |
| h-refinement | Nice-to-have     | Speed     |
| p-refinement | Nice-to-have     | Speed / acc    |
| r-refinement | Nice-to-have     | Text     |


| IO | Reason | Option 1 |
| -------- | -------- | -------- |
| Checkpoints function | Req. for HPC    | Text     |
| In-situ visualisation | Extract more physics without saving to disk     | Text     |
| Supports standard file formats | HDF5, NetCDF     | Text     |


| Equation terms | Reason | Option 1 |
| -------- | -------- | -------- |
| Source terms | Text    | Text     |
| Divergence | Text    | Text     |
| Curl | Text    | Text     |
| Advection | Text    | Text     |
| Diffusion | Text    | Text     |
| Anistropic diffusion | Text    | Text     |
| Parallel diffusion | Text    | Text     |
| Perpendicular diffusion | Text    | Text     |
| Hyper diffusion | Text    | Text     |
| Poisson bracket | Text    | Text     |
| Time evolving BCs | Text    | Text     |
| BCs that depend on solver variables' values and gradients | Text    | Text     |

| Use cases | Reason | Option 1 |
| -------- | -------- | -------- |
| Has a history of UQ applications | Text    | Text     |
| Adjoint methods have been successfully applied | Text    | Text     |



