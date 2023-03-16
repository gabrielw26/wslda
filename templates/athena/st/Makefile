# LOAD MODULES BEFORE COMPILATION
# source env.sh

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
 
