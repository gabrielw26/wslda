# Info - wderiv

C library for computation derivatives of functions defined on a lattice. The library uses spectral methods to achieve high accuracy.

# Compiling lib
Edit header of `Makefile` and:
```bash
make
```
*Note*: you need to have installed [FFTW](http://www.fftw.org/).  

The compilation process will produce:
* `libwderiv.a` - static library
* `libwderiv.so` - dyanamic library

Header file is located in `./c/` folder.

In addition, testcases and example codes will be compiled. In order to compile only the lib:
```bash
make lib
```

# Basic usage
The simples code:

```c
/**
 * Simple example showing usage of wderiv
 * 
 * Compile command
 * gcc -std=c99 -O3 simple-1d.c -o simple-1d -I./c/ -L. -lwderiv -lfftw3 -lm
 *   or 
 * g++ -O3 simple-1d.c -o simple-1d -I./c/ -L. -lwderiv -lfftw3 -lm
 **/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

// Include header
#include "wderiv.h"

// test function
double tfx(double x, double a, double x0)
{
    return exp(a*pow(x-x0,2));
}
// ... and its first derivative, (analytical result)
double tdfdx(double x, double a, double x0)
{
    return tfx(x,a,x0)*a*2.*(x-x0);
}
// ... and its second derivative (analytical result)
double td2fdx2(double x, double a, double x0)
{
    return tdfdx(x,a,x0)*a*2.*(x-x0) + tfx(x,a,x0)*a*2.0;
}

int main()
{
    // lattice
    int nx=256;    // number of points
    double dx=0.5; // lattice spacing
    
    double *fx = (double *)malloc(nx * sizeof(double));
    double *dfdx = (double *)malloc(nx * sizeof(double));
    double *d2fdx2 = (double *)malloc(nx * sizeof(double));
    
    int ix;
    double ax=-0.01;
    for(ix=0; ix<nx; ix++) 
        fx[ix] = tfx(dx*ix, ax, dx*nx/2); // create function
    
    // Initilize wderiv - 1D variant
    wderiv_init_1d(nx, dx);
    
    // compute derivatives
    wderiv_dfdx_1d_r(fx, dfdx);     // first derivative
    wderiv_d2fdx2_1d_r(fx, d2fdx2); // second derivative

    // print results
    printf("%12s %12s %12s %12s %12s %12s\n", " ", " ", "numerical", "analytical", "numerical", "analytical");
    printf("%12s %12s %12s %12s %12s %12s\n", "x", "f(x)", "df/dx", "df/dx", "d2f/dx2", "d2f/dx2");
    for(ix=0; ix<nx; ix++) 
        printf("%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",
               dx*ix, fx[ix], dfdx[ix], tdfdx(dx*ix, ax, dx*nx/2), d2fdx2[ix], td2fdx2(dx*ix, ax, dx*nx/2));
    
    return 0; 
}
```

# Documentation
See [wiki pages](https://gitlab.fizyka.pw.edu.pl/wtools/wderiv/-/wikis/home) for documentation.

# Contributors
* Ruszczak Bartosz, Warsaw University of Technology, Faculty of Physics
* Gabriel Wlazłowski, Warsaw University of Technology, Faculty of Physic



