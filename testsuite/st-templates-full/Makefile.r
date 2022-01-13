# LOAD MODULES BEFORE COMPILATION
# 
# source scl_source enable devtoolset-7 python27
# module load cuda/9.0
# module load openmpi-gcc721-Cuda90/3.1.1
# module load lapack/391
# module load scalapack-gcc721-cuda90-openmpi311/210
# module load elpa-gcc721/202005

# COMPILER
CXX=mpicc

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
CFLAGS= -std=gnu99 -O3 -I/usr/local/elpa202005-openmpi311-gcc721-cuda90-lapack391/include/elpa-2020.05.001/ \
	-DTESTSUITE -DMATRIX_IS_REAL -DREGULARIZATION_SCHEME=SPHERICAL_CUTOFF
LIBS=-lfftw3 -lm -llapack -lscalapack -lelpa

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.st
 
