# C compiler settings - system dependent
CXX=gcc
CFLAGS= -std=gnu99 -O3

WINTERP=-I../lib/winterp/c -L../lib/winterp/ -lwinterp -lfftw3
WDERIV=-I../lib/wderiv/c -L../lib/wderiv/ -lwderiv -lfftw3
WDATA=-I../lib/wdata/c -L../lib/wdata/ -lwdata -lfftw3

# -------------------------------------------------------------------------------------
# ------------------------------- DO NOT MODIFY ---------------------------------------
# -------------------------------------------------------------------------------------

all: bindir wdata2checkpoint dpca2wdata dpca2wdata-hard

bindir:
	mkdir -p ./bin

dpca2wdata: bindir 
	$(CXX) $(CFLAGS) dpca2wdata.c -o ./bin/dpca2wdata -I../hpc-engine/ $(WDATA) -lm 
	
dpca2wdata-hard: bindir 
	$(CXX) $(CFLAGS) dpca2wdata-hard.c -o ./bin/dpca2wdata-hard -I../hpc-engine/ $(WDATA) -lm 
	
wdata2checkpoint: bindir
	$(CXX) $(CFLAGS) wdata2checkpoint.c -o ./bin/wdata2checkpoint -L. $(WDATA) -lm 

clean:
	rm *.o
	rm ./bin/dpca2wdata
	rm ./bin/dpca2wdata-hard
	rm ./bin/wdata2checkpoint
