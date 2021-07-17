# C compiler settings - system dependent
CXX=gcc
CFLAGS= -std=gnu99 -O3

# -------------------------------------------------------------------------------------
# ------------------------------- DO NOT MODIFY ---------------------------------------
# -------------------------------------------------------------------------------------

all: bindir libwdata libwinterp wdata-cut wdata-sample dpca2wdata wdata2checkpoint dpca2wdata-hard wdata-interpolate

bindir:
	mkdir -p ./bin

libwdata: 
	$(CXX) -O3 -c -I../lib/wdata/c ../lib/wdata/c/wdata.c -o ./wdata.o -fPIC
	ar crf libwdata.a wdata.o
	
libwinterp: 
	$(CXX) -O3 -c -I../lib/winterp/c ../lib/winterp/c/winterp.c -o ./winterp.o -fPIC
	ar crf libwinterp.a winterp.o

wdata-cut: bindir libwdata
	$(CXX) $(CFLAGS) wdata-cut.c -o ./bin/wdata-cut -I../lib/wdata/c -L. -lwdata -lm 

wdata-interpolate: bindir libwdata libwinterp
	$(CXX) $(CFLAGS) wdata-interpolate.c -o ./bin/wdata-interpolate -I../lib/wdata/c -I../lib/winterp/c -I../hpc-engine -L. -lwdata -lwinterp -lfftw3 -lm 
	
wdata-sample: bindir libwdata
	$(CXX) $(CFLAGS) wdata-sample.c -o ./bin/wdata-sample -I../lib/wdata/c -L. -lwdata -lm 

dpca2wdata: bindir libwdata
	$(CXX) $(CFLAGS) dpca2wdata.c -o ./bin/dpca2wdata -I../lib/wdata/c -I../hpc-engine/ -L. -lwdata -lm 
	
dpca2wdata-hard: bindir libwdata
	$(CXX) $(CFLAGS) dpca2wdata-hard.c -o ./bin/dpca2wdata-hard -I../lib/wdata/c -I../hpc-engine/ -L. -lwdata -lm 
	
wdata2checkpoint: bindir libwdata
	$(CXX) $(CFLAGS) wdata2checkpoint.c -o ./bin/wdata2checkpoint -I../lib/wdata/c -L. -lwdata -lm 

clean:
	rm *.o
	rm libwdata.a
	rm ./bin/wdata-cut
	rm ./bin/wdata-sample
	rm ./bin/dpca2wdata
	rm ./bin/dpca2wdata-hard
