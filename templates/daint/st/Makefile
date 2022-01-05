
# COMPILER
CXX=cc
PYTHON=python

# DIRECTORY SETTINGS (must end with /)
# W-SLDA instalation folder (for read only)
WSLDADIR=$(WSLDA)
# folder where temporarry files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
# C compiler flags
CFLAGS= -std=gnu99 -O3  -I/project/pr125/share/elpa-2020.11.001/include/elpa-2020.11.001
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS=-lfftw3 -L/project/pr125/share/elpa-2020.11.001/lib -lelpa

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
