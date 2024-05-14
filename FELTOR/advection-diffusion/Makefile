# Makefile
device=omp
FELTOR_PATH=../../feltor

#configure machine
include ../../config/default.mk
include ../../config/version.mk
include ../../config/*.mk
include ../../config/devices/devices.mk

INCLUDE+=-I../../inc/

all: advection-diffusion advection-diffusion_hpc advection-diffusion_mpi

# only necessary if you use the draw library
#advection-diffusion: advection-diffusion.cpp equations.h init.h diag.h
#    $(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(JSONLIB) -g
advection-diffusion: advection-diffusion.cpp equations.h init.h diag.h
	$(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -DWITH_GLFW -g

advection-diffusion_hpc: advection-diffusion.cpp equations.h init.h diag.h $(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -g

advection-diffusion_mpi: advection-diffusion.cpp equations.h init.h diag.h $(MPICC) $(OPT) $(MPICFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -DWITH_MPI $(VERSION_FLAGS) -g

.PHONY: clean

clean: rm -rf advection-diffusion advection-diffusion_hpc advection-diffusion_mpi