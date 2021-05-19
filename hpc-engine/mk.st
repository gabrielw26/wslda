# =====================================================================
# =================== DO NOT MODIFY THIS PART =========================
# =====================================================================
SRCDIR=$(WSLDADIR)hpc-engine/
INCLUDEDIR=$(WSLDADIR)hpc-engine/
WDATAINC=$(WSLDADIR)lib-wdata/
WDATASRC=$(WSLDADIR)lib-wdata/
DEPS = pca_io.h pca_macro.h pca_settings.h pca_uniform.h pca_utils.h pca_edf.h pca_logger.h s2dpca_edf.h s3dpca_grid.h s2dpca_densities.h  s2dpca_me.h s2dpca_fft.h s2dpca_uext.h sxdpca_broyden.h s3dpca_me.h s3dpca_fft.h s3dpca_testfun.h s3dpca_densities.h s3dpca_uext.h wslda_writevars.h netlib-lapack.h wslda_potdens.h wslda_functionals.h aslda_edf.h wslda_reproducibility.h wslda_wavevectors.h wslda_errors.h s1dpca_me.h s1dpca_fft.h s1dpca_densities.h wslda_st_checkpoint.h wslda_interpolation.h wslda_resize.h jdb.h ../VERSION.h
WDATADEPS = wdata.h
DEPSLOCAL = predefines.h  problem-definition.h logger.h
SRCS = wslda_functionals.o wslda_potdens.o sxdpca_broyden.o pca_utils.o wslda_writevars.o wslda_reproducibility.o wslda_errors.o wslda_st_checkpoint.o
SRCS1D = s1dpca.o s1dpca_fft.o s1dpca_me.o s1dpca_densities.o
SRCS2D = s2dpca.o s2dpca_fft.o s2dpca_me.o s2dpca_densities.o 
SRCS3D = s3dpca.o s3dpca_fft.o s3dpca_densities.o s3dpca_me.o 
WDATALIB = wdata.o
CFLAGS := $(CFLAGS) -I$(INCLUDEDIR) -I$(WDATAINC) -I./

OBJS := $(addprefix $(OBJDIR),$(SRCS:.c=.o))
SRCS := $(addprefix $(SRCDIR),$(SRCS))
OBJS1D := $(addprefix $(OBJDIR),$(SRCS1D:.c=.o))
SRCS1D := $(addprefix $(SRCDIR),$(SRCS1D))
OBJS2D := $(addprefix $(OBJDIR),$(SRCS2D:.c=.o))
SRCS2D := $(addprefix $(SRCDIR),$(SRCS2D))
OBJS3D := $(addprefix $(OBJDIR),$(SRCS3D:.c=.o))
SRCS3D := $(addprefix $(SRCDIR),$(SRCS3D))
DEPS := $(addprefix $(INCLUDEDIR),$(DEPS))
WDATADEPS := $(addprefix $(WDATAINC),$(WDATADEPS))
WDATAOBJS := $(addprefix $(OBJDIR),$(WDATALIB:.c=.o))
WDATALIB := $(addprefix $(WDATASRC),$(WDATALIB))

all: 1d 2d 3d

1d: dirs d1d

2d: dirs d2d

3d: dirs d3d

d1d: $(OBJS1D) $(OBJS) $(WDATAOBJS)
	cat $(WSLDADIR)VERSION.h
	$(CXX) -o $(BINDIR)st-wslda-1d $(CFLAGS) $(LIBS) -DWSLDA $^ 

d2d: $(OBJS2D) $(OBJS) $(WDATAOBJS)
	cat $(WSLDADIR)VERSION.h
	$(CXX) -o $(BINDIR)st-wslda-2d $(CFLAGS) $(LIBS) -DWSLDA $^  

d3d: $(OBJS3D) $(OBJS) $(WDATAOBJS)
	cat $(WSLDADIR)VERSION.h
	$(CXX) -o $(BINDIR)st-wslda-3d $(CFLAGS) $(LIBS) -DWSLDA $^

ext: $(WDATAOBJS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) $(WSLDADIR)extensions/solve-uniform.c -o $(BINDIR)solve-uniform -lm -lfftw3
	$(CXX) $(CFLAGS) $(WDATAOBJS) $(WSLDADIR)extensions/st-test-external-potentials.c -o $(BINDIR)test-external-potentials -lm -lfftw3
	
$(OBJDIR)wslda_reproducibility.o: $(SRCDIR)wslda_reproducibility.c $(SRCDIR)wslda_reproducibility.h $(SRCDIR)pca_settings.h predefines.h problem-definition.h logger.h
	python $(SRCDIR)reproducibility.py $(OBJDIR)
	$(CXX) $(CFLAGS) -DWSLDA -I$(OBJDIR) -c -o $@ $<

$(OBJDIR)s1dpca.o: $(SRCDIR)s2dpca.c $(DEPS) $(WDATADEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -DWSLDA -DCODEDIM=1 -c -o $@ $<
	
$(OBJDIR)s2dpca.o: $(SRCDIR)s2dpca.c $(DEPS) $(WDATADEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -DWSLDA -DCODEDIM=2 -c -o $@ $<
	
$(OBJDIR)%.o: $(SRCDIR)%.c $(DEPS) $(WDATADEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $< 
	
$(OBJDIR)%.o: $(WDATASRC)%.c $(DEPS) $(WDATADEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $< 
	
.PHONY: clean all ver dirs 

dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

clean:
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)reproducibility.h
	rm -f $(BINDIR)st-wslda-1d
	rm -f $(BINDIR)st-wslda-2d
	rm -f $(BINDIR)st-wslda-3d
	rm -f $(BINDIR)test-external-potentials
	rm -f $(BINDIR)solve-uniform
	
ver:
	$(CXX) --version
 
