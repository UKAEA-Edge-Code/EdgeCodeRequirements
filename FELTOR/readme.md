[Home](../readme.md)
# FELTOR

## Overview

Description from the [FELTOR GitHub repo page](https://www.nektar.info/):

> "FELTOR (Full-F ELectromagnetic code in TORoidal geometry) is both a numerical library and a scientific software package built on top of it. Its main physical target [sic] are plasma edge and scrape-off layer (gyro-)fluid simulations. The numerical methods centre around discontinuous Galerkin methods on structured grids. Our core level functions are parallelized for a variety of hardware from multi-core cpu to hybrid MPI+GPU"

## Notes
Writing a solver is fairly straightforward if it is based off an existing one - writing one from scratch would be quite an involved process.  They have implemented their own BLAS routines which seems somewhat unnecessary.  There is a large selection of functors provided which can be used to form initial conditions.  Most functions are templated so can be used with their own container types or standard library ones.
Several choices of preconditioners and time-stepping tableaux are available.

GPUs are supported, but only CUDA compatible ones.

Input is done through a json file which specifies grid details, initial conditions, time-stepper etc., and the output is a NetCDF file.

There seems to be only one or two people involved in the development and little to no active community.
Documentation exists, but does not go into much detail about the library itself.  The code itself does have doxygen comments.

### Pros
 - Easy and quick to install
 - Good performance
 - Few dependencies
 - Tokamak specific

### Cons
 - Only structured grids
 - Small developer team
 - GPU acceleration limited to NVIDIA devices
 - Not all of the tutorials run correctly with the supplied input files
 - Lack of desired features: no UFL, mesh refinement or element variety

### Installation
Follow instructions in repository.

```
cd ~
git clone https://www.github.com/feltor-dev/feltor
# Fetch dependencies
git clone https://www.github.com/nvidia/thrust
git clone https://www.github.com/cusplibrary/cusplibrary
git clone https://www.github.com/vectorclass/version1 vcl
git clone https://www.github.com/feltor-dev/draw

sudo apt-get install libnetcdf-dev
sudo apt-get install libhdf5-dev
sudo apt-get install libjsoncpp-dev
sudo apt-get install libglfw3-dev

# We need to checkout an older version of thrust compatible with cusp
cd thrust
git checkout 1.9.3

# Link to dependencies
cd ~
mkdir include
cd include
ln -s ~/thrust/thrust
ln -s ~/cusplibrary/cusp
ln -s ~/vcl
ln -s ~/draw
ln -s /usr/include/jsoncpp/json

```

## Advection Example
Advection of a Gaussian potential on a 64x16 quadrilateral mesh with a unit x velocity and periodic boundary conditions

![feltor](https://github.com/UKAEA-Edge-Code/EdgeCodeRequirements/assets/16880076/83dc8728-ddd9-433c-aab0-dec156e88c97)

Move the ```advection-diffusion``` directory into the ```feltor/src``` directory, then navigate to the ```advection-diffusion``` directory
To build, run make:

```
make advection-diffusion device=cpu
```

To run the simulation:

```
./advection-diffusion inputfile.json outputfile.nc
```

## Links

- [GitHub repo](https://github.com/feltor-dev/feltor)
- [Python bindings](https://github.com/feltor-dev/pyFeltor)
