
# COMPILERS
CXX=mpicc
NVCC=nvcc
PYTHON=python

# DIRECTORY SETTINGS 
# W-SLDA installation folder (for read-only)
WSLDADIR=../
# folder where temporary files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS & LIBRARIES
# C compiler flags
CFLAGS= -std=gnu99 -O3 -I/usr/local/cuda/include/ -DTARGET_MACHINE=TITAN
# CUDA compiler flags
NVCCFLAGS = -arch sm_30 -O3 
# prescription how to link with: FFTW, CUDA
LIBS=-L/usr/local/cuda/lib64 -lcudart -lcufft -I/usr/local/cuda/include/ -lfftw3 -lstdc++ -lm


# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.td

