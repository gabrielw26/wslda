
# COMPILER
CXX=mpicc
PYTHON=python

# DIRECTORY SETTINGS 
# W-SLDA installation folder (for read-only)
WSLDADIR=../
# folder where temporary files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
# C compiler flags
CFLAGS= -std=gnu99 -O3
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS=-lfftw3 -lm -llapack -L/usr/lib64/openmpi/lib -lscalapack 

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
