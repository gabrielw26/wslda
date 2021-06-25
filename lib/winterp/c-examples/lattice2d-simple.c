/**
 * W-interp Library 
 * @see https://gitlab.fizyka.pw.edu.pl/wtools/winterp
 * 
 * This example shows how to interpolate data defined in the 2D lattice to a new resolution. 
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
    
    // parameters of new lattice
    int outNX=48, outNY=64; 
    
    double *output = (double *)malloc(outNX*outNY * sizeof(double));
    
    // interpolate to new lattice with new resolution, while keeping fixed volume
    winterp_interpolation_2dr(inNX, inNY, input, outNX, outNY, output);
    
    
    // save arrays
    FILE *f;
    // save input
    f=fopen("input_f.wdat", "wb");
    fwrite (input , sizeof(double)*inNX*inNY, 1, f);
    fclose(f);
    // save output
    f=fopen("output_f.wdat", "wb");
    fwrite (output , sizeof(double)*outNX*outNY, 1, f);
    fclose(f);    
    
    return 0; 
}
