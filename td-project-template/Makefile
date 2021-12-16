
# COMPILERS
CXX=mpicc
NVCC=nvcc
PYTHON=python

# DIRECTORY SETTINGS (must end with /)
# W-SLDA instalation folder (for read only)
WSLDADIR=../
# folder where temporarry files will be stored (will be created automatically)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
# C compiler flags
CFLAGS= -std=gnu99 -O3 -I/usr/local/cuda/include/
# CUDA compiler flags, make sure you select the correct compute capability for your gpus (-arch flag)
NVCCFLAGS = -arch sm_35 -O3 -std=c++03
# prescription how to link with: CUDA, FFTW
LIBS=-L/usr/local/cuda/lib64 -lcudart -lcufft -lfftw3 -lstdc++ -lm


# ----- DO NOT MODIFY -----
include $(WSLDADIR)hpc-engine/mk.td
