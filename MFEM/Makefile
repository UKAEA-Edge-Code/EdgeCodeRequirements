MFEM_SOURCE_DIR := ~/mfem
AD_DIR := $(MFEM_SOURCE_DIR)/miniapps/autodiff

include $(MFEM_SOURCE_DIR)/config/config.mk

diffusion: diffusion.o
	$(MFEM_CXX) -o $@ $<  $(MFEM_LIBS) 

diffusion.o: diffusion.cpp
	$(MFEM_CXX) $(MFEM_FLAGS) -O0 -g -I$(AD_DIR) -o $@ -c $^
