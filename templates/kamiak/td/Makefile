
# COMPILERS
CXX=mpicc
NVCC=nvcc
PYTHON=python

# DIRECTORY SETTINGS 
# W-SLDA installation folder (for read-only)
WSLDADIR=$(WSLDA)
# folder where temporary files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
# C compiler flags
CFLAGS= -std=gnu99 -O3
# CUDA compiler flags, make sure you select the correct compute capability for your gpus (-arch flag)
NVCCFLAGS = -arch sm_37 -O3 -std=c++14 -Wno-deprecated-gpu-targets
# prescription how to link with: CUDA, FFTW
LIBS=-lcudart -lcufft -lfftw3 -lstdc++ -lm


# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.td
