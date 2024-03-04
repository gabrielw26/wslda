/**
 * W-SLDA Toolkit
 *
 * This tests functions in wslda_utils.h
 *
 * gcc wslda_utils_test.c -o wslda_utils_test -lm
 * */

// Standard libraries
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "wslda_utils.h"

int main( int argc , char ** argv )
{
    double xmin=-20.0;
    double xmax=20.0;
    int    N=500;

    double dx=(xmax-xmin)/N;

    double x0=-5.0;
    double width = 10.0;

    double lvalue=-5.0;
    double rvalue= 7.0;

    int i;
    for(i=0; i<=N; i++)
    {
        double x = xmin + dx*i;
        double f;

        if(x<-10)
            f = smooth_fromatob(x,-15, -12, 0, -10);
        else if(x<5)
            f = smooth_fromatob(x,-2, 2, -10, 7);
        else
            f = smooth_fromatob(x,12, 15, 7, 0);

        // double f = smooth_step(x, 0.5, 1.0, 3.0, 4.0, 2.0);
        // double f = harmonic_oscillator_smooth_edges(x, 1.0, 10., 15., 2);

        printf("%16.12f %16.12f\n", x, f);
    }

    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
