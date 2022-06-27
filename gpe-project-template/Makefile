# LOAD BEFORE COMPILATION
# module load cuda/9.0

# COMPILERS
NVCC=nvcc
CXX=mpicc
PYTHON=python

# ARCHITECTURE
SM=sm_35

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj
BINDIR=.

# COMPILER FLAGS
NVCCFLAGS = -arch $(SM) -O3
CFLAGS=-std=gnu99 -O3 -I/usr/local/cuda-9.0/include/ -DTARGET_MACHINE=DWARF -Wl,--no-as-needed

# LIBRARIES
LIBS=-L/usr/local/cuda-9.0/lib64 -lcudart -lcufft -lfftw3 -lstdc++ -lm

# VERSION
# REAL - MODE 1 IN INPUT.TXT FILE
# IMAG - MODE 0 IN INPUT.TXT FILE

# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/gpe/make.gpe