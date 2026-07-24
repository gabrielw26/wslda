# LOAD MODULES BEFORE COMPILATION
# 
# module load impi/2020.4.912
# module load mkl/2020.0.4
# module load cudatoolkit/8.0.61-21551265

# COMPILER
CXX=mpiicc
PYTHON=python

# DIRECTORY SETTINGS 
WSLDADIR=/mnt/storage_2/project_data/grant_518/wslda/
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
CFLAGS= -O3 -I/mnt/storage_2/project_data/grant_518/share/elpa/include/elpa-2021.05.002/
LIBS=-lfftw3 -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lstdc++ -ldl -L/mnt/storage_2/project_data/grant_518/share/elpa/lib/ -lelpa 

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
