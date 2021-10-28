// Author: Gabriel Wlazlowski
// This is simple code for tracking the vortex position in 2d
// 
// compile:
//    gcc -std=gnu99 -c vdetect-simple-2d.c -o vdetect1.o
//    g++ -std=c++98 -c vdetect-simple-2d.cpp -o vdetect2.o -I$WSLDA/lib/wdata/c
//    g++ vdetect1.o vdetect2.o -o vdetect-simple-2d -L$WSLDA/lib/wdata -lwdata -lfftw3 -lgsl -lgslcblas -lm

// standard libraries
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "gsl/gsl_multimin.h"


inline double fbeta(double E, double beta)
{
    double bE=beta*E;
    if(bE>50.) return 0.0; // to avoid numerical problems
    else if(bE<-50.) return 1.0; // to avoid numerical problems
    else return 1.0/(exp(bE)+1.0);
}

// #include "pca_settings.h"
// #include "pca_macro.h"
// #include "pca_utils.h"
// #include "pca_io.h"

#define file_operationl( cmd )                                                  \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "FILE ERROR:: cannot execute: %s\n" , #cmd);          \
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        return( EXIT_FAILURE ) ;                                                \
    } }
    
// allocation of memory, not involving MPI
#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "ERROR: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        something_to_cheer_you_up_pid0(stdout);                                 \
        return -1 ;                                                             \
    } 
    

/**
 * Funtion returns inerpolated value of function for coordinate (x,y,z).
 * Fourier transform coefficients are reqired
 * */
double complex interpolate(double complex *fcoeffs, double x, double y, int nx, int ny)
{
    int ix, iy, ixyz;
    double kx, ky;
    double complex r = 0.0 + I*0.0;
    
    // multiply by momentum
    ixyz=0;
    for(ix=0; ix<nx; ix++)
    {
        for(iy=0; iy<ny; iy++)
        {
            // extract momentum
            if(ix<nx/2) kx=2.*M_PI/(( double )nx) * ( double )(ix   );
            else        kx=2.*M_PI/(( double )nx) * ( double )(ix-nx);

            
            if(iy<ny/2) ky=2.*M_PI/(( double )ny) * ( double )(iy   );
            else        ky=2.*M_PI/(( double )ny) * ( double )(iy-ny); 
            
            r += fcoeffs[ixyz]*cexp(I*( kx*x + ky*y ));
            
            ixyz++;
        }
    }
    
    return r;
}

/** 
 * Function computes phase difference between two argumenents
 * To be consistent with computation of velocity field by formula v=j/rho and v=grad Phi
 * return value is from range (-pi,+pi)
 * */
double pdiff(double arg1, double arg2)
{
    
    double p;
    if(arg1>=arg2) p=arg1-arg2;
    else           p=2.0*M_PI - (arg2-arg1);
    
    if(p>M_PI) p-=2.*M_PI;
    return p;
    
}

/**
 * @return 1 - pahase decreas monotnicaly 
 * */
int check_if_monotonic_down(const int size, double * args)
{
    int i;
    double argdiffs[size];
    for(i=0; i<size; i++) argdiffs[i]=pdiff(args[i],args[(i+1)%size]);
//     for(i=0; i<size; i++) printf("check_if_monotonic_down: %d %f\n", i, argdiffs[i]);
    for(i=0; i<size; i++) if(argdiffs[i]<0.0) return 0;
        
    return 1; 
}


// ----------------------- GSL MINIMIZATION ------------------------
extern int _gslnx;
extern int _gslny;

double my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double complex *p = (double complex *)params;

  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  double complex r = interpolate(p, x, y, _gslnx, _gslny);
  #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
  return cnorm(r);
}
