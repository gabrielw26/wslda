
# C compiler (with C99 dialect)
CC=gcc -std=gnu99

# To compile tools you need to install first winterp library
# https://gitlab.fizyka.pw.edu.pl/wtools/winterp
WINTERP=-I../winterp/c -L../winterp/ -lwinterp -lfftw3

# lib: generates lib static and dynamic (using C compiler)
# examples: generates example codes
# tools: generates set of auxiliary tools for manipulating wdata format

all: lib examples tools
lib: libwdata.a libwdata.so


libwdata.a: ./c/wdata.c ./c/wdata.h
	$(CC) -O3 -c ./c/wdata.c -fPIC
	ar crf libwdata.a wdata.o
	
libwdata.so: ./c/wdata.c ./c/wdata.h
	$(CC) -O3 -c ./c/wdata.c -fPIC -shared -o libwdata.so
	
examples: lib
	$(CC) ./c-examples/example-write.c -o ./c-examples/example-write -I./c/ -L. -lwdata -lm 
	$(CC) ./c-examples/example-write-many.c -o ./c-examples/example-write-many -I./c/ -L. -lwdata -lm 
	$(CC) ./c-examples/example-write-many-t_varying.c -o ./c-examples/example-write-many-t_varying -I./c/ -L. -lwdata -lm 
	$(CC) ./c-examples/example-read.c -o ./c-examples/example-read -I./c/ -L. -lwdata -lm 
	$(CC) ./c-examples/example-addvar.c -o ./c-examples/example-addvar -I./c/ -L. -lwdata -lm 
	
tools: lib
	mkdir -p ./bin/
	$(CC) ./tools/wdata-cut.c -o ./bin/wdata-cut -I./c/ -L. -lwdata -lm
	$(CC) ./tools/wdata-cut.c -o ./bin/wdata-stride -I./c/ -L. -lwdata -lm -DWDATA_STRIDE
	$(CC) ./tools/wdata-datadim-up.c -o ./bin/wdata-datadim-up -I./c/ -L. -lwdata -lm
	$(CC) ./tools/wdata-merge.c -o ./bin/wdata-merge -I./c/ -L. -lwdata -lm
	$(CC) ./tools/wdata-interpolate.c -o ./bin/wdata-interpolate -I./c/ -L. -lwdata $(WINTERP) -lm
	$(CC) ./tools/wdata-section.c -o ./bin/wdata-section -I./c/ -L. -lwdata $(WINTERP) -lm

clean:
	rm *.o
	rm *.a
	rm *.so
	rm ./c-examples/example-write
	rm ./c-examples/example-read
	rm ./c-examples/example-addvar
	rm ./c-examples/example-write-many
	rm ./c-examples/example-write-many-t_varying
	rm ./bin/wdata-cut
	rm ./bin/wdata-stride
	rm ./bin/wdata-interpolate
	rm ./bin/wdata-datadim-up
	rm ./bin/wdata-merge
	rm ./bin/wdata-section
