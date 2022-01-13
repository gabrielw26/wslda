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
CFLAGS= -std=c11 -O3  -DMKL_ILP64 -mkl=cluster -I/gs/hs1/hp190063/share/elpa/include/elpa-2020.05.001 -DM_PI=3.14159265358979323846 -DM_PI_2=1.570796326794896558
LIBS=-lfftw3  -lmkl_scalapack_ilp64 -lmkl_blacs_intelmpi_ilp64 -lpthread -ldl  -DMKL_ILP64 -mkl=cluster -L/gs/hs1/hp190063/share/elpa/lib -lelpa

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
