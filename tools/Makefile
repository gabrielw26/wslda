# C compiler settings - system dependent
CXX=gcc
CFLAGS= -std=gnu99 -O3

# -------------------------------------------------------------------------------------
# ------------------------------- DO NOT MODIFY ---------------------------------------
# -------------------------------------------------------------------------------------

all: bindir lib wdata-cut wdata-sample dpca2wdata wdata2checkpoint dpca2wdata-hard

bindir:
	mkdir -p ./bin

lib: 
	$(CXX) -O3 -c -I../lib-wdata/ ../lib-wdata/wdata.c -o ./wdata.o -fPIC
	ar crf libwdata.a wdata.o

wdata-cut: bindir lib
	$(CXX) $(CFLAGS) wdata-cut.c -o ./bin/wdata-cut -I../lib-wdata/ -L. -lwdata -lm 

wdata-sample: bindir lib
	$(CXX) $(CFLAGS) wdata-sample.c -o ./bin/wdata-sample -I../lib-wdata/ -L. -lwdata -lm 

dpca2wdata: bindir lib
	$(CXX) $(CFLAGS) dpca2wdata.c -o ./bin/dpca2wdata -I../lib-wdata/ -I../hpc-engine/ -L. -lwdata -lm 
	
dpca2wdata-hard: bindir lib
	$(CXX) $(CFLAGS) dpca2wdata-hard.c -o ./bin/dpca2wdata-hard -I../lib-wdata/ -I../hpc-engine/ -L. -lwdata -lm 
	
wdata2checkpoint: bindir lib
	$(CXX) $(CFLAGS) wdata2checkpoint.c -o ./bin/wdata2checkpoint -I../lib-wdata/ -L. -lwdata -lm 

clean:
	rm *.o
	rm libwdata.a
	rm ./bin/wdata-cut
	rm ./bin/wdata-sample
	rm ./bin/dpca2wdata
	rm ./bin/dpca2wdata-hard
