# LOAD BEFORE COMPILATION
# 
# module load cuda/9.0

# DIMENSION
NX=256
NY=32
NZ=32

# COMPILERS
NVCC=nvcc
CXX=gcc

# ARCHITECTURE
SM=sm_35

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj
# folder where executable binary will be placed (will be created automatically)
BINDIR=.

# COMPILER FLAGS
NVCCFLAGS = -arch $(SM) -O3
CFLAGS= -std=gnu99 -O3 -I/usr/local/cuda-9.0/include/ -DTARGET_MACHINE=DWARF

# LIBRARIES
LIBS=-lcudart -lcufft -lm

# PARAMETERS
PARAMS= -DNX=$(NX) -DNY=$(NY) -DNZ=$(NZ)

# VERSION
# REAL - DEPEND ON STEP
# IMAG - DEPEND ON TIME

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/gpe/make.gpe
