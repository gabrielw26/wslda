# LOAD BEFORE COMPILATION
# 
# source env.sh

# COMPILERS
CXX=mpicc
NVCC=nvcc
PYTHON=python

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
CFLAGS= -O3
NVCCFLAGS= -O3

# LIBRARIES
LIBS=-lcudart -lcufft -lm -lfftw3 -lstdc++


# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.td
