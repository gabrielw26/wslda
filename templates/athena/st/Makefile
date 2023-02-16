# LOAD MODULES BEFORE COMPILATION
# 
# source scl_source enable devtoolset-7 python27
# module load cuda/9.0
# module load openmpi-gcc721-Cuda90/3.1.1
# module load lapack/391
# module load scalapack-gcc721-cuda90-openmpi311/210
# module load elpa-gcc721/202005

# COMPILER
CXX=mpicc
PYTHON=python

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
CFLAGS= -std=gnu99 -O3 -I${EBROOTELPA}/include/elpa-2022.05.001
LIBS=-lfftw3 -lm -lopenblas -lscalapack -lelpa

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
