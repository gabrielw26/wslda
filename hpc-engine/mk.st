# =====================================================================
# =================== DO NOT MODIFY THIS PART =========================
# =====================================================================

# sets it to empty by default (the ?= means "assign only if not already set")
# and now you can use: 
#   make EXTRA_CFLAGS="-DTESTSUITE"
EXTRA_CFLAGS ?= 

WSLDADIR:=$(WSLDADIR)/
OBJDIR:=$(OBJDIR)/
BINDIR:=$(BINDIR)/

SRCDIR=$(WSLDADIR)hpc-engine/
INCLUDEDIR=$(WSLDADIR)hpc-engine/

WDATAINC=$(WSLDADIR)lib/wdata/c/
WDATASRC=$(WSLDADIR)lib/wdata/c/

WDERIVINC=$(WSLDADIR)lib/wderiv/c/
WDERIVSRC=$(WSLDADIR)lib/wderiv/c/

WINTERPINC=$(WSLDADIR)lib/winterp/c/
WINTERPSRC=$(WSLDADIR)lib/winterp/c/

DEPS = pca_io.h pca_macro.h pca_settings.h pca_uniform.h pca_utils.h pca_edf.h pca_logger.h s2dpca_edf.h s3dpca_grid.h s2dpca_densities.h  s2dpca_me.h s2dpca_fft.h s2dpca_uext.h sxdpca_broyden.h s3dpca_me.h s3dpca_fft.h s3dpca_testfun.h s3dpca_densities.h s3dpca_uext.h wslda_writevars.h netlib-lapack.h wslda_potdens.h wslda_functionals.h aslda_edf.h sldae_functional.h wslda_reproducibility.h wslda_wavevectors.h wslda_errors.h s1dpca_me.h s1dpca_fft.h s1dpca_densities.h wslda_st_checkpoint.h wslda_resize.h wslda_api_version.h jdb.h wslda_predefines_test.h ../VERSION.h

WDATADEPS = wdata.h
WDERIVDEPS = wderiv.h
WINTERPDEPS = winterp.h

DEPSLOCAL = predefines.h
SRCS = wslda_functionals.o wslda_potdens.o sldae_functional.o sxdpca_broyden.o pca_utils.o wslda_writevars.o wslda_reproducibility.o wslda_errors.o wslda_st_checkpoint.o
SRCS1D = s1dpca.o s1dpca_fft.o s1dpca_me.o s1dpca_densities.o
SRCS2D = s2dpca.o s2dpca_fft.o s2dpca_me.o s2dpca_densities.o
SRCS3D = s3dpca.o s3dpca_fft.o s3dpca_densities.o s3dpca_me.o

WDATALIB = wdata.o
WDERIVLIB = wderiv.o
WINTERPLIB = winterp.o

CFLAGS := $(CFLAGS) -I$(INCLUDEDIR) -I$(WDATAINC) -I$(WDERIVINC) -I$(WINTERPINC) -I./ $(EXTRA_CFLAGS)

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

WDERIVDEPS := $(addprefix $(WDERIVINC),$(WDERIVDEPS))
WDERIVOBJS := $(addprefix $(OBJDIR),$(WDERIVLIB:.c=.o))
WDERIVLIB := $(addprefix $(WDERIVSRC),$(WDERIVLIB))

WINTERPDEPS := $(addprefix $(WINTERPINC),$(WINTERPDEPS))
WINTERPOBJS := $(addprefix $(OBJDIR),$(WINTERPLIB:.c=.o))
WINTERPLIB := $(addprefix $(WINTERPSRC),$(WINTERPLIB))

# Default values for variables
ifeq ($(strip $(PYTHON)),)
PYTHON=python
endif

# Default values for variables
ifneq ($(strip $(WSLDA_MACHINE)),)
CFLAGS := $(CFLAGS) -I$(WSLDA_MACHINE)
endif

ifneq (,$(wildcard machine.h))
DEPSLOCAL := $(DEPSLOCAL) machine.h
endif

all: 1d 2d 3d

1d: dirs d1d

2d: dirs d2d

3d: dirs d3d

0d: dirs d0d

ext: dirs \
     (BINDIR)test-external-potentials-1d \
     (BINDIR)test-external-potentials-2d \
     (BINDIR)test-external-potentials-3d \
     (BINDIR)interpolate-dataset

d1d: $(OBJS1D) $(OBJS) $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS)
	$(CXX) -o $(BINDIR)st-wslda-1d $(CFLAGS) -DWSLDA $^ $(LIBS)
	$(PYTHON) $(SRCDIR)review-settings.py $(WSLDADIR) $(BINDIR)st-wslda-1d

d2d: $(OBJS2D) $(OBJS) $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS)
	$(CXX) -o $(BINDIR)st-wslda-2d $(CFLAGS) -DWSLDA $^ $(LIBS)
	$(PYTHON) $(SRCDIR)review-settings.py $(WSLDADIR) $(BINDIR)st-wslda-2d

d3d: $(OBJS3D) $(OBJS) $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS)
	$(CXX) -o $(BINDIR)st-wslda-3d $(CFLAGS) -DWSLDA $^ $(LIBS)
	$(PYTHON) $(SRCDIR)review-settings.py $(WSLDADIR) $(BINDIR)st-wslda-3d

d0d: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL)
	$(CXX) -o $(BINDIR)st-wslda-0d $(CFLAGS) $(WSLDADIR)hpc-engine/s0dpca.c $(CFLAGS) $(LIBS)
	cat $(WSLDADIR)VERSION.h

(BINDIR)test-external-potentials-1d: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL)
	$(CXX) -o $(BINDIR)test-external-potentials-1d $(CFLAGS) $(WDATAOBJS) $(WSLDADIR)extensions/st-test-external-potentials.c -DCODEDIM=1 -lm -lfftw3

(BINDIR)test-external-potentials-2d: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL)
	$(CXX) -o $(BINDIR)test-external-potentials-2d $(CFLAGS) $(WDATAOBJS) $(WSLDADIR)extensions/st-test-external-potentials.c -DCODEDIM=2 -lm -lfftw3

(BINDIR)test-external-potentials-3d: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL)
	$(CXX) -o $(BINDIR)test-external-potentials-3d $(CFLAGS) $(WDATAOBJS) $(WSLDADIR)extensions/st-test-external-potentials.c -DCODEDIM=3 -lm -lfftw3

(BINDIR)interpolate-dataset: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL)
	$(CXX) -o $(BINDIR)interpolate-dataset $(CFLAGS) $(WDATAOBJS) $(WINTERPOBJS) $(WSLDADIR)extensions/interpolate-dataset.c -lm -lfftw3

$(OBJDIR)wslda_reproducibility.o: $(SRCDIR)wslda_reproducibility.c $(SRCDIR)wslda_reproducibility.h $(SRCDIR)pca_settings.h $(DEPSLOCAL) problem-definition.h logger.h
	$(PYTHON) $(SRCDIR)reproducibility.py $(OBJDIR)
	$(CXX) $(CFLAGS) -DWSLDA -I$(OBJDIR) -c -o $@ $<

$(OBJDIR)pca_utils.o: $(SRCDIR)pca_utils.c $(SRCDIR)pca_utils.h $(DEPS) $(WDATADEPS) $(WDERIVDEPS) $(DEPSLOCAL) problem-definition.h logger.h
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

$(OBJDIR)wslda_functionals.o: $(SRCDIR)wslda_functionals.c $(DEPS) $(WDATADEPS) $(WDERIVDEPS) $(DEPSLOCAL) problem-definition.h
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

$(OBJDIR)s1dpca.o: $(SRCDIR)s2dpca.c $(DEPS) $(WDATADEPS) $(WDERIVDEPS) $(DEPSLOCAL) logger.h
	$(CXX) $(CFLAGS) -DWSLDA -DCODEDIM=1 -c -o $@ $<

$(OBJDIR)s2dpca.o: $(SRCDIR)s2dpca.c $(DEPS) $(WDATADEPS) $(WDERIVDEPS) $(DEPSLOCAL) logger.h
	$(CXX) $(CFLAGS) -DWSLDA -DCODEDIM=2 -c -o $@ $<

$(OBJDIR)s3dpca.o: $(SRCDIR)s3dpca.c $(DEPS) $(WDATADEPS) $(WDERIVDEPS) $(DEPSLOCAL) logger.h
	$(CXX) $(CFLAGS) -DWSLDA -DCODEDIM=3 -c -o $@ $<

$(OBJDIR)%.o: $(SRCDIR)%.c $(DEPS) $(WDATADEPS) $(WDERIVDEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

$(OBJDIR)%.o: $(WDATASRC)%.c $(DEPS) $(WDATADEPS)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

$(OBJDIR)%.o: $(WDERIVSRC)%.c $(DEPS) $(WDERIVDEPS)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

$(OBJDIR)%.o: $(WINTERPSRC)%.c $(DEPS) $(WINTERPDEPS)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

.PHONY: clean all ver dirs

dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

clean:
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)reproducibility.h
	rm -f $(BINDIR)st-wslda-0d
	rm -f $(BINDIR)st-wslda-1d
	rm -f $(BINDIR)st-wslda-2d
	rm -f $(BINDIR)st-wslda-3d
	rm -f $(BINDIR)test-external-potentials-1d
	rm -f $(BINDIR)test-external-potentials-2d
	rm -f $(BINDIR)test-external-potentials-3d
	rm -f $(BINDIR)interpolate-dataset

ver:
	$(CXX) --version
