
# C++ compiler
CXX=g++

# C compiler
CC=gcc

# lib: generates lib static and dynamic (using C++ compiler)
# lic: generates lib static and dynamic (using C compiler)
# examples: generates example codes

all: libc lib examples
libs: lib libc

lib: 
	$(CXX) -O3 -c ./c/wdata.c -fPIC
	ar crf libwdata.a wdata.o
	$(CXX) -O3 -c ./c/wdata.c -fPIC -shared -o libwdata.so
	
libc: 
	$(CC) -O3 -c ./c/wdata.c -fPIC
	ar crf libwdatac.a wdata.o 
	$(CC) -O3 -c ./c/wdata.c -fPIC -shared -o libwdatac.so

examples: lib
	$(CXX) ./c-examples/example-write.c -o ./c-examples/example-write -I./c/ -L. -lwdata -lm 
	$(CXX) ./c-examples/example-read.c -o ./c-examples/example-read -I./c/ -L. -lwdata -lm 
	$(CXX) ./c-examples/example-addvar.c -o ./c-examples/example-addvar -I./c/ -L. -lwdata -lm 

clean:
	rm *.o
	rm *.a
	rm *.so
	rm ./c-examples/example-write
	rm ./c-examples/example-read
	rm ./c-examples/example-addvar
