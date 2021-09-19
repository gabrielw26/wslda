# LOAD BEFORE COMPILATION
# 
# module load cuda/9.0

# COMPILERS
CXX=mpicc
NVCC=nvcc 

# DIRECTORY SETTINGS (must end with /)
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
CFLAGS= -std=gnu99 -O3 -I/usr/local/cuda-9.0/include/ -DTARGET_MACHINE=DWARF -DTESTSUITE
NVCCFLAGS = -arch sm_35 -O3 -std=c++03 -DTESTSUITE

# LIBRARIES
LIBS=-L/usr/local/cuda-9.0/lib64 -lcudart -lcufft -lfftw3 -lm


# ----- DO NOT MODIFY -----
include $(WSLDADIR)hpc-engine/mk.td
