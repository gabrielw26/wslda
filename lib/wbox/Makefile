# C compiler (with C99 dialect)
ifeq ($(strip $(CC)),)
    CC=gcc -std=gnu99
endif

# WDATA lib need only for compiling examples
WDATA=-L../wdata/ -I../wdata/c/ -lwdata

# lib: generates lib static and dynamic (using C compiler)
# examples: generates example codes

all: lib examples
lib: libwbox.a libwbox.so	

libwbox.a: ./c/wbox.c ./c/wbox.h
	$(CC) -O3 -c ./c/wbox.c -fPIC
	ar crf libwbox.a wbox.o 
	
libwbox.so: ./c/wbox.c ./c/wbox.h
	$(CC) -O3 -c ./c/wbox.c -fPIC -shared -o libwbox.so

examples: lib
	$(CC) -O3 ./c-examples/lattice1d-insert.c -o ./c-examples/lattice1d-insert -I./c/ -L. -lwbox -lm
	$(CC) -O3 ./c-examples/lattice2d-insert.c -o ./c-examples/lattice2d-insert -I./c/ -L. -lwbox -lm
	$(CC) -O3 ./c-examples/lattice3d-insert.c -o ./c-examples/lattice3d-insert -I./c/ -L. -lwbox -lm
	# $(CC) -O3 ./c-examples/lattice2d-merge.c -o ./c-examples/lattice2d-merge -I./c/ -L. -lwbox -lm
	# $(CC) -O3 ./c-examples/lattice2d-rotate.c -o ./c-examples/lattice2d-rotatee -I./c/ -L. -lwbox -lm

	
clean:
	rm -f *.o
	rm -f *.a
	rm -f *.so
	rm -f ./c-examples/lattice1d-insert
	rm -f ./c-examples/lattice2d-insert
	rm -f ./c-examples/lattice3d-insert

