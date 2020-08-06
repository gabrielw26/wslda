# C compiler
CXX=gcc

CFLAGS= -std=gnu99 -O3

all: lib wdata-cut dpca2wdata

lib: 
	$(CXX) -O3 -c -I../lib-wdata/ ../lib-wdata/wdata.c -o ./wdata.o -fPIC
	ar crf libwdata.a wdata.o

wdata-cut: lib
	$(CXX) $(CFLAGS) wdata-cut.c -o wdata-cut -I../lib-wdata/ -L. -lwdata -lm 
	
dpca2wdata: lib
	$(CXX) $(CFLAGS) dpca2wdata.c -o dpca2wdata -I../lib-wdata/ -I../hpc-engine/ -L. -lwdata -lm 

clean:
	rm *.o
	rm libwdata.a
	rm wdata-cut
	rm dpca2wdata
