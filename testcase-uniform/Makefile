
# COMPILER
CXX=mpicc

# DIRECTORY SETTINGS
OBJDIR=./obj/
BINDIR=./

# COMPILER FLAGS
CFLAGS= -std=gnu99 -O3  
LIBS=-lfftw3 -lm -llapack -L/usr/lib64/openmpi/lib -lscalapack 



# =====================================================================
# =================== DO NOT MODIFY THIS PART =========================
# =====================================================================
SRCDIR=../pca/
INCLUDEDIR=../pca/
DEPS = pca_io.h pca_macro.h pca_settings.h pca_uniform.h pca_utils.h pca_edf.h pca_logger.h s2dpca_edf.h s3dpca_grid.h s2dpca_densities.h  s2dpca_me.h s2dpca_fft.h s2dpca_uext.h sxdpca_broyden.h s3dpca_me.h s3dpca_fft.h s3dpca_testfun.h s3dpca_densities.h s3dpca_uext.h s3dpca_inittype4.h
SRCS = sxdpca_broyden.o
SRCS2D = s2dpca.o s2dpca_fft.o s2dpca_me.o s2dpca_densities.o 
SRCS3D = s3dpca.o s3dpca_fft.o s3dpca_densities.o s3dpca_me.o s3dpca_inittype4.o
CFLAGS := $(CFLAGS) -I$(INCLUDEDIR) -I./

OBJS := $(addprefix $(OBJDIR),$(SRCS:.c=.o))
SRCS := $(addprefix $(SRCDIR),$(SRCS))
OBJS2D := $(addprefix $(OBJDIR),$(SRCS2D:.c=.o))
SRCS2D := $(addprefix $(SRCDIR),$(SRCS2D))
OBJS3D := $(addprefix $(OBJDIR),$(SRCS3D:.c=.o))
SRCS3D := $(addprefix $(SRCDIR),$(SRCS3D))
DEPS := $(addprefix $(INCLUDEDIR),$(DEPS))
	
all: 2d 3d

2d: $(OBJS2D) $(OBJS) 
	$(CXX) $(CFLAGS) $(LIBS) -DWSLDA -o $(BINDIR)wslda-2d $^  
	
3d: $(OBJS3D) $(OBJS) 
	$(CXX) $(CFLAGS) $(LIBS) -DWSLDA -o $(BINDIR)wslda-3d $^

$(OBJDIR)%.o: $(SRCDIR)%.c $(DEPS) dirs
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $< 
	

dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

.PHONY: clean all ver

clean:
	rm -f $(OBJDIR)*
	rm -f $(BINDIR)wslda-2d
	rm -f $(BINDIR)wslda-3d
	
ver:
	$(CXX) --version

 
