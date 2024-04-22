[Home](../readme.md)
# FELTOR

## Overview

Description from the [FELTOR GitHub repo page](https://www.nektar.info/):

> "FELTOR (Full-F ELectromagnetic code in TORoidal geometry) is both a numerical library and a scientific software package built on top of it. Its main physical target [sic] are plasma edge and scrape-off layer (gyro-)fluid simulations. The numerical methods centre around discontinuous Galerkin methods on structured grids. Our core level functions are parallelized for a variety of hardware from multi-core cpu to hybrid MPI+GPU"

## Notes
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

## Links

- [GitHub repo](https://github.com/feltor-dev/feltor)
- [Python bindings](https://github.com/feltor-dev/pyFeltor)
