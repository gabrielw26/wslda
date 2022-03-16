# LOAD BEFORE COMPILATION
# 
# module load cuda/9.0

# COMPILERS
NVCC=nvcc

# ARCHITECTURE
SM=sm_35

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
NVCCFLAGS = -arch $(SM) -O3 -std=c++03

# LIBRARIES
LIBS=-lcudart -lcufft -lm

# DIMENSION
NX=1024
NY=128
NZ=128

# VERSION
# REAL - DEPEND ON STEP
# IMAGINARY - DEPEND ON TIME

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/gpe/make.gpe