
# C++ compiler
CXX=g++

# lib: generates static lib
# examples: generates example codes

all: lib examples


lib: 
	$(CXX) -O3 -c wdata.c -fPIC
	ar crf libwdata.a wdata.o

examples: lib example-write.c
	$(CXX) example-write.c -o example-write -L. -lwdata -lm 
	$(CXX) example-read.c -o example-read -L. -lwdata -lm 

