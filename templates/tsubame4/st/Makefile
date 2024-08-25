# LOAD BEFORE COMPILATION
# 
# source env.sh

# COMPILER
CXX=mpiicc
PYTHON=python

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
# C compiler flags
CFLAGS= -std=c11 -O3 -qmkl=sequential -I/gs/bs/hp190063/share/libs/include/elpa-2024.05.001
# prescription how to link with: FFTW, Lapack, ScaLapack, ELPA (optionally)
LIBS= -lfftw3 -L${MKLROOT}/lib -qmkl=sequential  -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -L/gs/bs/hp190063/share/libs/lib -lelpa
 
# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st

