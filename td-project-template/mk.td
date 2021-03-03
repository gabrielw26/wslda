
# =====================================================================
# =================== DO NOT MODIFY THIS PART =========================
# =====================================================================
SRCDIR=$(WSLDADIR)hpc-engine/
INCLUDEDIR=$(WSLDADIR)hpc-engine/
WDATAINC=$(WSLDADIR)lib-wdata/
WDATASRC=$(WSLDADIR)lib-wdata/
DEPS = pca_io.h pca_kernels.h pca_macro.h pca_settings.h pca_uniform.h pca_utils.h pca_derivative.h pca_edf.h pca_dens.h pca_logger.h \
       ccpca_kernels.h ccpca_dens.h cpca_checkpoint.h \
       cpca_kernels.h cpca_derivative.h cpca_dens.h wslda_potdens.h \
       wslda_reproducibility.h wslda_errors.h wslda_wavevectors.h ../VERSION.h
NVCCDEPS = pca_macro.h pca_settings.h pca_uext.h pca_derivative.h pca_edf.h cpca_derivative.h wslda_cuda_utils.h
WDATADEPS = wdata.h
DEPSLOCAL = predefines.h  problem-definition.h logger.h

SRCS = wslda_potdens.o pca_utils.o wslda_writevars.o wslda_reproducibility.o wslda_errors.o

SRCS1D = ccpca.o ccpca_checkpoint.o
NVCCSRCS1D = ccpca_derivative.o ccpca_kernels.o ccpca_dens.o

SRCS2D = cpca.o cpca_checkpoint.o
NVCCSRCS2D = cpca_derivative.o cpca_kernels.o cpca_dens.o

SRCS3D = pca.o pca_checkpoint.o
NVCCSRCS3D = pca_derivative.o pca_kernels.o pca_dens.o

WDATALIB = wdata.o

CFLAGS := $(CFLAGS) -I$(INCLUDEDIR) -I$(WDATAINC) -I./ -DTDWSLDA
NVCCFLAGS := $(NVCCFLAGS) -I$(INCLUDEDIR) -I./ -DTDWSLDA

OBJS := $(addprefix $(OBJDIR),$(SRCS:.c=.o))
SRCS := $(addprefix $(SRCDIR),$(SRCS))
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

all: 1d 2d 3d

1d: dirs d1d

2d: dirs d2d

3d: dirs d3d

d1d: $(OBJS1D) $(OBJS) $(NVCCOBJS1D) $(NVCCOBJS) $(WDATAOBJS)
	$(CXX) $(CFLAGS) $(LIBS) -o $(BINDIR)td-wslda-1d $^ 

d2d: $(OBJS2D) $(OBJS) $(NVCCOBJS2D) $(NVCCOBJS) $(WDATAOBJS)
	$(CXX) $(CFLAGS) $(LIBS) -o $(BINDIR)td-wslda-2d $^ 

d3d: $(OBJS3D) $(OBJS) $(NVCCOBJS3D) $(NVCCOBJS) $(WDATAOBJS)
	$(CXX) $(CFLAGS) $(LIBS) -o $(BINDIR)td-wslda-3d $^
	
ext: $(WDATAOBJS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) $(WSLDADIR)extensions/solve-uniform.c -o $(BINDIR)solve-uniform -lm -lfftw3

$(OBJDIR)ccpca.o: $(SRCDIR)ccpca.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS)
	$(CXX) $(CFLAGS) -DCODEDIM=1 -c -o $@ $<
	
$(OBJDIR)cpca.o: $(SRCDIR)cpca.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS)
	$(CXX) $(CFLAGS) -DCODEDIM=2 -c -o $@ $<
	
$(OBJDIR)cpca_checkpoint.o: $(SRCDIR)cpca_checkpoint.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS)
	$(CXX) $(CFLAGS) -DCODEDIM=2 -c -o $@ $<

$(OBJDIR)ccpca_checkpoint.o: $(SRCDIR)cpca_checkpoint.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS)
	$(CXX) $(CFLAGS) -DCODEDIM=1 -c -o $@ $<

$(OBJDIR)wslda_reproducibility.o: $(SRCDIR)wslda_reproducibility.c $(SRCDIR)wslda_reproducibility.h $(SRCDIR)pca_settings.h predefines.h problem-definition.h logger.h
	python $(SRCDIR)reproducibility.py $(OBJDIR)
	$(CXX) $(CFLAGS) -I$(OBJDIR) -c -o $@ $<
	
$(OBJDIR)%.o: $(SRCDIR)%.c $(DEPS) $(DEPSLOCAL) $(WDATADEPS)
	$(CXX) $(CFLAGS) -c -o $@ $< 

$(OBJDIR)%.o: $(WDATASRC)%.c $(DEPS) $(WDATADEPS) $(DEPSLOCAL)
	$(CXX) $(CFLAGS) -c -o $@ $< 

$(OBJDIR)%.o: $(SRCDIR)%.cu $(NVCCDEPS) $(DEPSLOCAL)
	$(NVCC) $(NVCCFLAGS) -c -o $@ $< 

.PHONY: clean all ver dirs

dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

clean:
	rm -f $(OBJDIR)*
	rm -f $(BINDIR)td-wslda-1d
	rm -f $(BINDIR)td-wslda-2d
	rm -f $(BINDIR)td-wslda-3d
	
ver:
	$(CXX) --version 
