# C compiler
ifeq ($(strip $(CC)),)
    CC=gcc -std=gnu99
endif

ifeq ($(strip $(MPICPP)),)
    MPICPP=mpic++ -std=gnu99 
endif

# fftw linking statment
FFTW=-lfftw3

# WDATA lib need only for compiling examples
WDATA=-L../wdata/ -I../wdata/c/ -lwdata

# lib: generates lib static and dynamic (using C compiler)
# examples: generates example codes

all: lib examples
lib: libwinterp.a libwinterp.so	

libwinterp.a: ./c/winterp.c ./c/winterp.h
	$(CC) -O3 -c ./c/winterp.c -fPIC
	ar crf libwinterp.a winterp.o 
	
libwinterp.so: ./c/winterp.c ./c/winterp.h
	$(CC) -O3 -c ./c/winterp.c -fPIC -shared -o libwinterp.so

examples: lib
	$(CC) -O3 ./c-examples/lattice2d-simple.c -o ./c-examples/lattice2d-simple -I./c/ -L. -lwinterp $(FFTW) -lm
	$(CC) -O3 ./c-examples/interp1d_simple.c -o ./c-examples/interp1d_simple -I./c/ -L. -lwinterp $(FFTW) -lm
	$(CC) -O3 ./c-examples/interp2d_simple.c -o ./c-examples/interp2d_simple -I./c/ -L. -lwinterp $(FFTW) -lm
	$(CC) -O3 ./c-examples/interp3d_simple.c -o ./c-examples/interp3d_simple -I./c/ -L. -lwinterp $(FFTW) -lm
	$(CC) -O3 ./c-examples/section2d_simple.c -o ./c-examples/section2d_simple -I./c/ -L. -lwinterp $(FFTW) -lm
	
clean:
	rm -f *.o
	rm -f *.a
	rm -f *.so
	rm -f ./c-examples/lattice2d-simple
	rm -f ./c-examples/interp1d_simple
	rm -f ./c-examples/interp2d_simple
	rm -f ./c-examples/interp3d_simple
	rm -f ./c-examples/section2d_simple

