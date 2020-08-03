# =====================================================================
# =================== DO NOT MODIFY THIS PART =========================
# =====================================================================
SRCDIR=$(WSLDADIR)
INCLUDEDIR=$(WSLDADIR)
WDATAINC=$(WDATADIR)
WDATASRC=$(WDATADIR)
DEPS = pca_io.h pca_macro.h pca_settings.h pca_uniform.h pca_utils.h pca_edf.h pca_logger.h s2dpca_edf.h s3dpca_grid.h s2dpca_densities.h  s2dpca_me.h s2dpca_fft.h s2dpca_uext.h sxdpca_broyden.h s3dpca_me.h s3dpca_fft.h s3dpca_testfun.h s3dpca_densities.h s3dpca_uext.h s3dpca_inittype4.h wslda_writevars.h netlib-lapack.h
WDATADEPS = wdata.h
DEPSLOCAL = predefines.h  problem-definition.h
SRCS = sxdpca_broyden.o pca_utils.o wslda_writevars.o
SRCS2D = s2dpca.o s2dpca_fft.o s2dpca_me.o s2dpca_densities.o 
SRCS3D = s3dpca.o s3dpca_fft.o s3dpca_densities.o s3dpca_me.o s3dpca_inittype4.o
WDATALIB = wdata.o
CFLAGS := $(CFLAGS) -I$(INCLUDEDIR) -I$(WDATAINC) -I./

OBJS := $(addprefix $(OBJDIR),$(SRCS:.c=.o))
SRCS := $(addprefix $(SRCDIR),$(SRCS))
OBJS2D := $(addprefix $(OBJDIR),$(SRCS2D:.c=.o))
SRCS2D := $(addprefix $(SRCDIR),$(SRCS2D))
OBJS3D := $(addprefix $(OBJDIR),$(SRCS3D:.c=.o))
SRCS3D := $(addprefix $(SRCDIR),$(SRCS3D))
DEPS := $(addprefix $(INCLUDEDIR),$(DEPS))
WDATADEPS := $(addprefix $(WDATAINC),$(WDATADEPS))
WDATAOBJS := $(addprefix $(OBJDIR),$(WDATALIB:.c=.o))
WDATALIB := $(addprefix $(WDATASRC),$(WDATALIB))


all: 2d 3d

2d: $(OBJS2D) $(OBJS) $(WDATAOBJS)
	$(CXX) $(CFLAGS) $(LIBS) -DWSLDA -o $(BINDIR)st-wslda-2d $^  
	
3d: $(OBJS3D) $(OBJS) $(WDATAOBJS)
	$(CXX) $(CFLAGS) $(LIBS) -DWSLDA -o $(BINDIR)st-wslda-3d $^

$(OBJDIR)%.o: $(SRCDIR)%.c $(DEPS) $(WDATADEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $< 
	
$(OBJDIR)%.o: $(WDATASRC)%.c $(DEPS) $(WDATADEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $< 
	

dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

.PHONY: clean ver

clean:
	rm -f $(OBJDIR)*
	rm -f $(BINDIR)st-wslda-2d
	rm -f $(BINDIR)st-wslda-3d
	
ver:
	$(CXX) --version
 
