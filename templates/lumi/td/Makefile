
# COMPILERS
CXX=cc
# NVCC=nvcc
HIPCC=hipcc
PYTHON=python3

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
CFLAGS=-D__HIP_PLATFORM_AMD__ -I$(HIP_PATH)/include
# NVCCFLAGS= 
HIPCCFLAGS=-std=c++14 --offload-arch=gfx90a 

# LIBRARIES
LIBS= -L$(HIP_LIB_PATH) -lamdhip64 -L$(EBROOTHIPFFT)/lib -lhipfft -lm


# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.td
