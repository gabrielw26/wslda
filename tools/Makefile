# C compiler
CXX=gcc

CFLAGS= -std=gnu99 -O3

all: lib wdata-cut

lib: 
	$(CXX) -O3 -c -I../lib-wdata/ ../lib-wdata/wdata.c -o ./wdata.o -fPIC
	ar crf libwdata.a wdata.o

wdata-cut: lib
	$(CXX) $(CFLAGS) wdata-cut.c -o wdata-cut -I../lib-wdata/ -L. -lwdata -lm 
