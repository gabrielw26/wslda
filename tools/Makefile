# C compiler settings - system dependent
CXX=gcc
MPICPP=mpic++
CFLAGS= -std=gnu99 -O3

WINTERP=-I../lib/winterp/c -L../lib/winterp/ -lwinterp -lfftw3
WDERIV=-I../lib/wderiv/c -L../lib/wderiv/ -lwderiv -lfftw3
WDATA=-I../lib/wdata/c -L../lib/wdata/ -lwdata 

# -------------------------------------------------------------------------------------
# ------------------------------- DO NOT MODIFY ---------------------------------------
# -------------------------------------------------------------------------------------

all: bindir wdata2checkpoint dpca2wdata dpca2wdata-hard wslda-hdc-3d

bindir:
	mkdir -p ./bin

dpca2wdata: bindir 
	$(CXX) $(CFLAGS) dpca2wdata.c -o ./bin/dpca2wdata -I../hpc-engine/ $(WDATA) -lm 
	
dpca2wdata-hard: bindir 
	$(CXX) $(CFLAGS) dpca2wdata-hard.c -o ./bin/dpca2wdata-hard -I../hpc-engine/ $(WDATA) -lm 
	
wdata2checkpoint: bindir
	$(CXX) $(CFLAGS) wdata2checkpoint.c -o ./bin/wdata2checkpoint -L. $(WDATA) -lm 

wslda-hdc-3d: bindir
	$(MPICPP) hdc-3d.cpp -o ./bin/wslda-hdc-3d $(WDATA) -lfftw3 -lm
	
clean:
	rm -f *.o
	rm -f ./bin/dpca2wdata
	rm -f ./bin/dpca2wdata-hard
	rm -f ./bin/wdata2checkpoint
	rm -f ./bin/wslda-hdc-3d
