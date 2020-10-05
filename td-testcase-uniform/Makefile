
# COMPILERS
CXX=mpicc
NVCC = nvcc 

# DIRECTORY SETTINGS (must end with /)
WSLDADIR=../
OBJDIR=./obj/
BINDIR=./

# COMPILER FLAGS
CFLAGS= -std=gnu99 -O3 -I/usr/local/cuda/include/ -DTARGET_MACHINE=TITAN
NVCCFLAGS = -arch sm_30 -O3 

# LIBRARIES
LIBS=-L/usr/local/cuda/lib64 -lcudart -lcufft -I/usr/local/cuda/include/ -lm


include mk.td

