
# C compiler
CC=gcc -std=c99

# fftw linking statment
FFTW=-lfftw3


# lib: generates lib static and dynamic (using C compiler)
# examples: generates example codes

all: lib examples tests
tests: test test1d test2d
lib: libwderiv.a libwderiv.so	

libwderiv.a: ./c/wderiv.c ./c/wderiv.h
	$(CC) -O3 -c ./c/wderiv.c -fPIC
	ar crf libwderiv.a wderiv.o 
	
libwderiv.so: ./c/wderiv.c ./c/wderiv.h
	$(CC) -O3 -c ./c/wderiv.c -fPIC -shared -o libwderiv.so

examples: lib
	$(CC) -O3 ./c-examples/simple-1d.c -o ./c-examples/simple-1d -I./c/ -L. -lwderiv $(FFTW) -lm 

test: lib
	$(CC) -std=c99 ./c/test-wderiv.c -o ./c/test-wderiv -I./c/ -L. -lwderiv $(FFTW) -lm 

test2d: lib
	$(CC) -std=c99 ./c/test-wderiv2D.c -o ./c/test-wderiv2D -I./c/ -L. -lwderiv $(FFTW) -lm

test1d: lib
	$(CC) -std=c99 ./c/test-wderiv1D.c -o ./c/test-wderiv1D -I./c/ -L. -lwderiv $(FFTW) -lm 
	
clean:
	rm -f *.o
	rm -f *.a
	rm -f *.so
	rm -f ./c-examples/simple-1d
	rm -f ./c/test-wderiv
	rm -f ./c/test-wderiv2D
	rm -f ./c/test-wderiv1D
