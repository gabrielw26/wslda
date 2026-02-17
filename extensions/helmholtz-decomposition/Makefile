# settings

MPICPP=mpic++
WDATA=-I$(WSLDA)/lib/wdata/c -L$(WSLDA)/lib/wdata/ -lwdata


all: wslda-hdc-3d

wslda-hdc-3d:
	$(MPICPP) hdc-3d.cpp -o wslda-hdc-3d $(WDATA) -lfftw3 -lm

clean:
	rm -f ./wslda-hdc-3d
