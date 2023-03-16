# LOAD MODULES BEFORE COMPILATION
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
CFLAGS= -std=gnu99 -O3
NVCCFLAGS = -arch sm_80 -O3 

# LIBRARIES
LIBS=-lcudart -lcufft -lfftw3 -lstdc++ -lm


# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.td
