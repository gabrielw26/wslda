# =====================================================================
# =================== DO NOT MODIFY THIS PART =========================
# =====================================================================
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

DEPS = pca_io.h pca_kernels.h pca_macro.h pca_settings.h pca_uniform.h pca_utils.h pca_derivative.h pca_edf.h pca_dens.h pca_logger.h \
       ccpca_kernels.h ccpca_dens.h cpca_checkpoint.h \
       cpca_kernels.h cpca_derivative.h cpca_dens.h wslda_potdens.h \
       wslda_reproducibility.h wslda_errors.h wslda_wavevectors.h tdwslda_static_vars.h jdb.h sldae_functional.h wslda_predefines_test.h \
       wslda_api_version.h ../VERSION.h
NVCCDEPS = pca_macro.h pca_settings.h pca_uext.h pca_derivative.h pca_edf.h cpca_derivative.h wslda_cuda_utils.h tdwslda_functionals.h \
        reduce_many.h tdwslda_energy_generic.h

WDATADEPS = wdata.h
WDERIVDEPS = wderiv.h
WINTERPDEPS = winterp.h

DEPSLOCAL = predefines.h

SRCS = wslda_potdens.o pca_utils.o wslda_writevars.o wslda_reproducibility.o wslda_errors.o sldae_functional.o
NVCCSRCS = reduce_many.o

SRCS1D = ccpca.o ccpca_checkpoint.o
NVCCSRCS1D = ccpca_derivative.o ccpca_kernels.o ccpca_dens.o

SRCS2D = cpca.o cpca_checkpoint.o
NVCCSRCS2D = cpca_derivative.o cpca_kernels.o cpca_dens.o

SRCS3D = pca.o pca_checkpoint.o
NVCCSRCS3D = pca_derivative.o pca_kernels.o pca_dens.o

WDATALIB = wdata.o
WDERIVLIB = wderiv.o
WINTERPLIB = winterp.o

CFLAGS := $(CFLAGS) -I$(INCLUDEDIR) -I$(WDATAINC) -I$(WDERIVINC) -I$(WINTERPINC) -I./ -DTDWSLDA
NVCCFLAGS := $(NVCCFLAGS) -I$(INCLUDEDIR) -I./ -DTDWSLDA

OBJS := $(addprefix $(OBJDIR),$(SRCS:.c=.o))
SRCS := $(addprefix $(SRCDIR),$(SRCS))
NVCCOBJS := $(addprefix $(OBJDIR),$(NVCCSRCS:.cu=.o))
NVCCSRCS := $(addprefix $(SRCDIR),$(NVCCSRCS))
OBJS1D := $(addprefix $(OBJDIR),$(SRCS1D:.c=.o))
SRCS1D := $(addprefix $(SRCDIR),$(SRCS1D))
NVCCOBJS1D := $(addprefix $(OBJDIR),$(NVCCSRCS1D:.cu=.o))
NVCCSRCS1D := $(addprefix $(SRCDIR),$(NVCCSRCS1D))
OBJS2D := $(addprefix $(OBJDIR),$(SRCS2D:.c=.o))
SRCS2D := $(addprefix $(SRCDIR),$(SRCS2D))
NVCCOBJS2D := $(addprefix $(OBJDIR),$(NVCCSRCS2D:.cu=.o))
NVCCSRCS2D := $(addprefix $(SRCDIR),$(NVCCSRCS2D))
OBJS3D := $(addprefix $(OBJDIR),$(SRCS3D:.c=.o))
SRCS3D := $(addprefix $(SRCDIR),$(SRCS3D))
NVCCOBJS3D := $(addprefix $(OBJDIR),$(NVCCSRCS3D:.cu=.o))
NVCCSRCS3D := $(addprefix $(SRCDIR),$(NVCCSRCS3D))
DEPS := $(addprefix $(INCLUDEDIR),$(DEPS))
NVCCDEPS := $(addprefix $(INCLUDEDIR),$(NVCCDEPS))

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
NVCCFLAGS := $(NVCCFLAGS) -I$(WSLDA_MACHINE)
endif

ifneq (,$(wildcard machine.h))
DEPSLOCAL := $(DEPSLOCAL) machine.h
endif

ifneq ($(strip $(HIPCC)),)
HIPCCFLAGS := $(HIPCCFLAGS) -I$(INCLUDEDIR) -I./ -DTDWSLDA -DHIPMODE
CFLAGS := $(CFLAGS) -DHIPMODE
ifneq ($(strip $(WSLDA_MACHINE)),)
HIPCCFLAGS := $(HIPCCFLAGS) -I$(WSLDA_MACHINE)
endif
endif

all: 1d 2d 3d

1d: dirs d1d

2d: dirs d2d

3d: dirs d3d

ext: dirs \
     (BINDIR)solve-uniform \
     (BINDIR)interpolate-dataset

d1d: $(OBJS1D) $(OBJS) $(NVCCOBJS1D) $(NVCCOBJS) $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS)
	cat $(WSLDADIR)VERSION.h
	$(CXX) -o $(BINDIR)td-wslda-1d $(CFLAGS) $^ $(LIBS)

d2d: $(OBJS2D) $(OBJS) $(NVCCOBJS2D) $(NVCCOBJS) $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS)
	cat $(WSLDADIR)VERSION.h
	$(CXX) -o $(BINDIR)td-wslda-2d $(CFLAGS) $^ $(LIBS)

d3d: $(OBJS3D) $(OBJS) $(NVCCOBJS3D) $(NVCCOBJS) $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS)
	cat $(WSLDADIR)VERSION.h
	$(CXX) -o $(BINDIR)td-wslda-3d $(CFLAGS) $^ $(LIBS)

(BINDIR)solve-uniform: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL)
	$(CXX) -o $(BINDIR)solve-uniform $(CFLAGS) $(WSLDADIR)extensions/solve-uniform.c -lm -lfftw3

(BINDIR)interpolate-dataset: $(WDATAOBJS) $(WDERIVOBJS) $(WINTERPOBJS) $(DEPSLOCAL)
	$(CXX) -o $(BINDIR)interpolate-dataset $(CFLAGS) $(WDATAOBJS) $(WINTERPOBJS) $(WSLDADIR)extensions/interpolate-dataset.c -lm -lfftw3

$(OBJDIR)ccpca.o: $(SRCDIR)ccpca.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS) $(WDERIVDEPS) logger.h
	$(CXX) $(CFLAGS) -DCODEDIM=1 -c -o $@ $<

$(OBJDIR)cpca.o: $(SRCDIR)cpca.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS) $(WDERIVDEPS) logger.h
	$(CXX) $(CFLAGS) -DCODEDIM=2 -c -o $@ $<
	
$(OBJDIR)pca.o: $(SRCDIR)pca.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS) $(WDERIVDEPS) logger.h
	$(CXX) $(CFLAGS) -DCODEDIM=3 -c -o $@ $<

$(OBJDIR)cpca_checkpoint.o: $(SRCDIR)cpca_checkpoint.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS) $(WDERIVDEPS)
	$(CXX) $(CFLAGS) -DCODEDIM=2 -c -o $@ $<

$(OBJDIR)ccpca_checkpoint.o: $(SRCDIR)cpca_checkpoint.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS) $(WDERIVDEPS)
	$(CXX) $(CFLAGS) -DCODEDIM=1 -c -o $@ $<

$(OBJDIR)wslda_reproducibility.o: $(SRCDIR)wslda_reproducibility.c $(SRCDIR)wslda_reproducibility.h $(SRCDIR)pca_settings.h $(DEPSLOCAL) problem-definition.h logger.h
	$(PYTHON) $(SRCDIR)reproducibility.py $(OBJDIR)
	$(CXX) $(CFLAGS) -I$(OBJDIR) -c -o $@ $<

$(OBJDIR)pca_utils.o: $(SRCDIR)pca_utils.c $(SRCDIR)pca_utils.h $(SRCDIR)pca_settings.h $(DEPSLOCAL) problem-definition.h logger.h
	$(CXX) $(CFLAGS) -c -o $@ $<

$(OBJDIR)%.o: $(SRCDIR)%.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS) $(WDERIVDEPS)
	$(CXX) $(CFLAGS) -c -o $@ $<

$(OBJDIR)%.o: $(WDATASRC)%.c $(DEPS) $(WDATADEPS)
	$(CXX) $(CFLAGS) -c -o $@ $<

$(OBJDIR)%.o: $(WDERIVSRC)%.c $(DEPS) $(WDERIVDEPS)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

$(OBJDIR)%.o: $(WINTERPSRC)%.c $(DEPS) $(WINTERPDEPS)
	$(CXX) $(CFLAGS) -DWSLDA -c -o $@ $<

ifneq ($(strip $(HIPCC)),)

$(OBJDIR)ccpca_kernels.o: $(SRCDIR)hip/ccpca_kernels.cpp $(NVCCDEPS) $(DEPSLOCAL) problem-definition.h
	$(HIPCC) $(HIPCCFLAGS) -c -o $@ $<

$(OBJDIR)cpca_kernels.o: $(SRCDIR)hip/cpca_kernels.cpp $(NVCCDEPS) $(DEPSLOCAL) problem-definition.h
	$(HIPCC) $(HIPCCFLAGS) -c -o $@ $<

$(OBJDIR)pca_kernels.o: $(SRCDIR)hip/pca_kernels.cpp $(NVCCDEPS) $(DEPSLOCAL) problem-definition.h
	$(HIPCC) $(HIPCCFLAGS) -c -o $@ $<

$(OBJDIR)%.o: $(SRCDIR)hip/%.cpp $(NVCCDEPS) $(DEPSLOCAL)
	$(HIPCC) $(HIPCCFLAGS) -c -o $@ $<

else

$(OBJDIR)ccpca_kernels.o: $(SRCDIR)ccpca_kernels.cu $(NVCCDEPS) $(DEPSLOCAL) problem-definition.h
	$(NVCC) $(NVCCFLAGS) -Xcudafe "--diag_suppress=declared_but_not_referenced" -Xcudafe "--diag_suppress=set_but_not_used" -c -o $@ $<

$(OBJDIR)cpca_kernels.o: $(SRCDIR)cpca_kernels.cu $(NVCCDEPS) $(DEPSLOCAL) problem-definition.h
	$(NVCC) $(NVCCFLAGS) -Xcudafe "--diag_suppress=declared_but_not_referenced" -Xcudafe "--diag_suppress=set_but_not_used" -c -o $@ $<

$(OBJDIR)pca_kernels.o: $(SRCDIR)pca_kernels.cu $(NVCCDEPS) $(DEPSLOCAL) problem-definition.h
	$(NVCC) $(NVCCFLAGS) -Xcudafe "--diag_suppress=declared_but_not_referenced" -Xcudafe "--diag_suppress=set_but_not_used" -c -o $@ $<

$(OBJDIR)%.o: $(SRCDIR)%.cu $(NVCCDEPS) $(DEPSLOCAL)
	$(NVCC) $(NVCCFLAGS) -Xcudafe "--diag_suppress=declared_but_not_referenced" -Xcudafe "--diag_suppress=set_but_not_used" -c -o $@ $<

endif

.PHONY: clean all ver dirs

dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

clean:
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)reproducibility.h
	rm -f $(BINDIR)td-wslda-1d
	rm -f $(BINDIR)td-wslda-2d
	rm -f $(BINDIR)td-wslda-3d
	rm -f $(BINDIR)solve-uniform

ver:
	$(CXX) --version
