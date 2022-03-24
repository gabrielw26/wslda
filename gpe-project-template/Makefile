# LOAD BEFORE COMPILATION
# 
# module load cuda/9.0

# DIMENSION
NX=256
NY=32
NZ=32

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

# PARAMETERS
PARAMS=-DNX=$(NX) -DNY=$(NY) -DNZ=$(NZ)

# VERSION
# REAL - DEPEND ON STEP
# IMAGINARY - DEPEND ON TIME

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/gpe/make.gpe
