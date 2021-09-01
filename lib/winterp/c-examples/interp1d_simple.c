/**
 * W-interp Library 
 * @see https://gitlab.fizyka.pw.edu.pl/wtools/winterp
 * 
 * This example shows how to interpolate 1D data
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
    int inNX=32; 
    double inDX=1.0;
    
    // fill input array with some function
    double complex *inputc = (double complex*)malloc(inNX*sizeof(double complex));
    double         *inputr = (double        *)malloc(inNX*sizeof(double        ));
    
    int ix;
    int ixyz=0;
    double ax=-0.02, ay=-0.05;
    for(ix=0; ix<inNX; ix++)
    {
        inputr[ixyz] = tfx(inDX*ix, ax, inDX*inNX/2);
        inputc[ixyz] = tfx(inDX*ix, ax, inDX*inNX/2) + I*tfx(inDX*ix, ay, inDX*inNX/2); //gaussian
        ixyz++;
    }
    
    // create interpolator
    winterp_interpolator iinputr;
    winterp_create_interpolator_1d_r(inNX,inDX, inputr, &iinputr);
    
    for(ix=0; ix<inNX; ix++)
    {
        double vr;
        winterp_getvalue_1d_r(&iinputr, inDX*ix, &vr);
        printf("%12.8f %16.12g %16.12g\n", inDX*ix, inputr[ix], vr);
    }
    
    winterp_destroy_interpolator(&iinputr);
    
    // create interpolator
    winterp_interpolator iinputc;
    winterp_create_interpolator_1d_c(inNX,inDX, inputc, &iinputc);
    
    for(ix=0; ix<inNX; ix++)
    {
        double complex vc;
        winterp_getvalue_1d_c(&iinputc, inDX*ix, &vc);
        printf("%12.8f %16.12g %16.12g %16.12g %16.12g\n", inDX*ix, creal(inputc[ix]), cimag(inputc[ix]), creal(vc), cimag(vc));
    }
        
    winterp_destroy_interpolator(&iinputc);
    
    return 0; 
}
