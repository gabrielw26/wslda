# COMPILERS
CXX=cc
HIPCC=hipcc
PYTHON=python3

# DIRECTORY SETTINGS 
WSLDADIR=$(WSLDA)
OBJDIR=./obj/
# folder where executable binary will be placed (will be created automatically)
BINDIR=./

# COMPILER FLAGS
CFLAGS= -D__HIP_PLATFORM_AMD__ -I$(ROCM_PATH)/include -O3
HIPCCFLAGS=-std=c++14 --offload-arch=gfx90a -xhip -O3

# LIBRARIES
LIBS= -L$(ROCM_PATH)/lib -lamdhip64 -lhipfft -L$(CRAYLIBS_X86_64)/lib -lcraymath


# ----- DO NOT MODIFY -----
include $(WSLDADIR)/hpc-engine/mk.td
