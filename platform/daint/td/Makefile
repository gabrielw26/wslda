
# COMPILERS
CXX=cc
NVCC=nvcc
PYTHON=python

# DIRECTORY SETTINGS (must end with /)
# W-SLDA instalation folder (for read only)
WSLDADIR=$(WSLDA)
# folder where temporarry files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS & LIBRARIES
# C compiler flags
CFLAGS= -std=gnu99 -O3 -DTARGET_MACHINE=TITAN
# CUDA compiler flags
NVCCFLAGS = -arch sm_60 -O3 
# prescription how to link with: FFTW, CUDA
LIBS=-lcudart -lcufft -lm -lstdc++


# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.td


