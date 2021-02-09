
# C++ compiler
CXX=g++
# C compiler
CC=gcc

# lib: generates static lib
# examples: generates example codes

all: libc lib examples
libs: lib libc

lib: 
	$(CXX) -O3 -c wdata.c -fPIC
	ar crf libwdata.a wdata.o
	
libc: 
	$(CC) -O3 -c wdata.c -fPIC
	ar crf libwdatac.a wdata.o

examples: lib example-write.c
	$(CXX) example-write.c -o example-write -L. -lwdata -lm 
	$(CXX) example-read.c -o example-read -L. -lwdata -lm 
	$(CXX) example-addvar.c -o example-addvar -L. -lwdata -lm 

clean:
	rm wdata.o
	rm libwdata.a
	rm example-write
	rm example-read
	rm example-addvar
