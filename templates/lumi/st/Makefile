
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
CFLAGS= -std=gnu99 -O3 -I/project/project_465000708/share/elpa-libsci/PrgEnv-gnu-8.3.3/LUMI-22.06/rocm-5.1.4/include/elpa-2022.05.001
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS=-lfftw3 -lm -L/project/project_465000708/share/elpa-libsci/PrgEnv-gnu-8.3.3/LUMI-22.06/rocm-5.1.4/lib -lelpa 

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
