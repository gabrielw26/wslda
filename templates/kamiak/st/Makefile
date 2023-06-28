
# COMPILER
CXX=mpiicc
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
CFLAGS= -std=c11 -O3 -DMKL_LP64 -I$MKL_HOME/include/intel64/lp64 -I$ELPA_HOME/include/elpa-2023.05.001
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS=-lfftw3 -L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -I$MKL_HOME/include/intel64/lp64 -L$ELPA_HOME/lib -lelpa

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
