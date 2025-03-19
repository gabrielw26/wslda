
# COMPILER
CXX=cc
PYTHON=python3

# DIRECTORY SETTINGS 
# W-SLDA installation folder (for read-only)
WSLDADIR=$(WSLDA)
# folder where temporary files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
# C compiler flags
CFLAGS= -std=gnu99 -O3 -I$(ELPAINCLUDE)
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS=-lfftw3 -lm -L$(ELPALIB) -lelpa

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
