# LOAD MODULES BEFORE COMPILATION
# 
# module load openmpi/3.1.4_gcc620
# module unload cudatoolkit
# module load cuda/11.2.1_460.32.03

# COMPILERS
CXX=mpicc
NVCC=nvcc
PYTHON=python

# DIRECTORY SETTINGS (must end with /)
# W-SLDA instalation folder (for read only)
WSLDADIR=/mnt/storage_2/project_data/grant_518/wslda/
# folder where temporarry files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS & LIBRARIES
# C compiler flags
CFLAGS= -std=gnu99 -O3 
# CUDA compiler flags
NVCCFLAGS = -arch sm_60 -O3 
# prescription how to link with: FFTW, CUDA
LIBS=-lcudart -lcufft -lfftw3 -lstdc++ -lm


# ----- DO NOT MODIFY -----
include $(WSLDADIR)hpc-engine/mk.td

