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
    int inNX=32, inNY=32; 
    double inDX=1.0, inDY=1.0;
    
    // fill input array with some function
    double *input = (double *)malloc(inNX*inNY * sizeof(double));
    
    int ix, iy, ixy=0;
    double ax=-0.02, ay=-0.05;
    for(ix=0; ix<inNX; ix++) for(iy=0; iy<inNY; iy++) 
    {
        input[ixy] = tfx(inDX*ix, ax, inDX*inNX/2)*tfx(inDY*iy, ay, inDY*inNY/2); // 2D gaussian
        ixy++;
    }
    
    // print function along x, for fixed y=inDY*inNY/2
    for(ix=0; ix<inNX; ix++) printf("%12.8f %16.12g\n", inDX*ix, input[inNY/2 + ix*inNY]);
    printf("\n\n");
    
    // create interpolator
    winterp_interpolator iinput;
    winterp_create_interpolator_2d_r(inNX, inNY, inDX, inDY, input, &iinput);

    // print interpolated function along x, for fixed y=inDY*inNY/2
    double x=0.0, vr; 
    while(x<inDX*inNX)
    {
        winterp_getvalue_2d_r(&iinput, x, inDY*inNY/2, &vr);
        printf("%12.8f %16.12g\n", x, vr);
        x+=0.1; // increment
    }
    
    winterp_destroy_interpolator(&iinput);

    return 0; 
}
