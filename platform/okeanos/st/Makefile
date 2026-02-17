
# COMPILER
CXX=cc
PYTHON=python

# DIRECTORY SETTINGS
# W-BSK installation folder (for read-only)
WSLDADIR=../
# folder where temporary files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
# C compiler flags
CFLAGS= -O3
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS=-lfftw3 -lm -L/opt/cray/diag/lib

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
