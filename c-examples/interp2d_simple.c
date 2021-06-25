/**
 * W-interp Library 
 * @see https://gitlab.fizyka.pw.edu.pl/wtools/winterp
 * 
 * This example shows how to interpolate 2D data
 * 
 * */ 

// winterp lib
#include "winterp.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 

// test function
double tfx(double x, double a, double x0)
{
    return exp(a*pow(x-x0,2));
} 


int main()
{
    // initial lattice
    int inNX=32, inNY=28; 
    double inDX=0.9, inDY=1.1;
    
    // fill input array with some function
    double complex *inputc = (double complex*)malloc(inNX*inNY*sizeof(double complex));
    double         *inputr = (double        *)malloc(inNX*inNY*sizeof(double        ));
    
    int ix, iy;
    int ixyz=0;
    double ax=-0.05, ay=-0.08;
    for(ix=0; ix<inNX; ix++) for(iy=0; iy<inNY; iy++)
    {
        double _r = sqrt(pow(inDX*ix-inDX*inNX/2, 2) + pow(inDY*iy-inDY*inNY/2, 2));
        inputr[ixyz] = tfx(_r, ax, 0);
        inputc[ixyz] = tfx(_r, ax, 0) + I*tfx(_r, ay, 0); //gaussian
        ixyz++;
    }
    
    // create interpolator
    winterp_interpolator iinputr;
    winterp_create_interpolator_2d_r(inNX, inNY, inDX, inDY, inputr, &iinputr);
    
    for(ix=0; ix<inNX; ix++)
    {
        double vr;
        winterp_getvalue_2d_r(&iinputr, inDX*ix, inDY*inNY/2, &vr);
        ixyz = inNY/2 + ix*inNY;
        printf("%12.8f %16.12g %16.12g\n", inDX*ix, inputr[ixyz], vr);
    }
    
    winterp_destroy_interpolator(&iinputr);
    
    // create interpolator
    winterp_interpolator iinputc;
    winterp_create_interpolator_2d_c(inNX, inNY, inDX, inDY, inputc, &iinputc);
    
    for(ix=0; ix<inNX; ix++)
    {
        double complex vc;
        winterp_getvalue_2d_c(&iinputc, inDX*ix, inDY*inNY/2, &vc);
        ixyz = inNY/2 + ix*inNY;
        printf("%12.8f %16.12g %16.12g %16.12g %16.12g\n", inDX*ix, creal(inputc[ixyz]), cimag(inputc[ixyz]), creal(vc), cimag(vc));
    }
        
    winterp_destroy_interpolator(&iinputc);
    
    return 0; 
}
