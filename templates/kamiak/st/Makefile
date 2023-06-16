
# COMPILER
CXX=mpicc
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
CFLAGS= -std=gnu99 -O3 -I/opt/apps/elpa/2018.11.001/gcc/7.3.0/openmpi/4.0.3/include/elpa-2018.11.001
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS=-lfftw3 -lm -lblas -L/opt/apps/lapack/3.8.0/gcc/7.3.0/lib/ -llapack -L/opt/apps/scalapack/2.1.0/gcc/7.3.0/openmpi/4.0.3/lib/ -lscalapack -L/opt/apps/elpa/2018.11.001/gcc/7.3.0/openmpi/4.0.3/lib -lelpa

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
