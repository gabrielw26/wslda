/**
 * This code was used to test correctness of computation of matrix elemnts for derivative operators
 * Author: Gabriel Wlazlowski, WUT
 *
 * Compilation command:
 * gcc test_compute_matrix_elements_of_momentum_operator.c -o test_compute_matrix_elements_of_momentum_operator -lfftw3 -lm
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// See for reference:
// https://en.cppreference.com/w/c/numeric/complex
#include <complex.h>

// See for reference:
// http://www.fftw.org/fftw3_doc/
#include <fftw3.h>
/**
 * Test function
 * */
double function_x(double x)
{
    #define A -0.03
    return exp(A*x*x);
}

/**
 * Its first derivative...
 * */
double dx_function_x(double x)
{
    return 2.*A*x*function_x(x);
}

/**
 * ...and its second derivative
 * */
double dx2_function_x(double x)
{
    return 2.*A*function_x(x) + 2.*A*x*dx_function_x(x);
    #undef A
}

/**
 * You can use this function to check diff between two arrays
 * */
void test_array_diff(int N, double *a, double *b)
{
    int ixyz=0;

    double d,d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ixyz=0; ixyz<N; ixyz++)
    {
        d = a[ixyz]-b[ixyz];

        d2=d*d;
        sumd2+=d2;
        if(d2>maxd2) maxd2=d2;
    }

    printf("#    COMPARISON RESULTS:\n");
    printf("#           |max[a-b]| : %16.8g\n", sqrt(maxd2));
    printf("#         SUM[(a-b)^2] : %16.8g\n", sumd2);
    printf("# SQRT(SUM[(a-b)^2])/N : %16.8g\n", sqrt(sumd2)/N);
}

#define wprintf printf
// allocation of memory, not involving MPI
#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "ERROR: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }

/**
 * TAKEN FROM s2dpca_me.c
 * Function that computes matrix elements of -i*d/dx operator
 * @param nx lattice size (INPUT)
 * @param dx lattice constant
 * @param h matrix of size [nx*nx] with computed matrix elements (OUTPUT)
 *          <x1| -id/dx |x2> is given in matrix element h[x1 + nx*x2]
 * */
int compute_matrix_elements_of_momentum_operator(int nx, double dx, double complex *me)
{
    int ix, ci, ri; // iterator
    double kx;

    // allocate memory
    double complex *fft1;
    cppmallocl(fft1, nx, double complex); // for working area of plan

    // fftw plan
    fftw_plan plan_f_1d;
    fftw_plan plan_b_1d;

    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    plan_f_1d = fftw_plan_dft_1d(nx, fft1, fft1, FFTW_FORWARD, USE_FFTW_PLANNER);
    plan_b_1d = fftw_plan_dft_1d(nx, fft1, fft1, FFTW_BACKWARD, USE_FFTW_PLANNER);

    // reset matrix elements
    for(ix=0; ix<nx*nx; ix++) me[ix] = 0.0 + I*0.0;

    for(ci=0; ci<nx; ci++) // column-major iteration fashion, for each column do:
    {
        for(ix=0; ix<nx; ix++) fft1[ix]=0.0+I*0.0;
        fft1[ci] = 1.0/(nx); // normalization factor already included
        fftw_execute(plan_f_1d); // to momentum space

        // multiply by momentum
        for(ix=0; ix<nx; ix++)
        {
            // extract momentum
            if(ix<nx/2) kx=2.*M_PI/( ( double )nx * dx ) * ( double )(ix   );
            else        kx=2.*M_PI/( ( double )nx * dx ) * ( double )(ix-nx);

            if(ix==nx/2) kx=0.0;

            fft1[ix]*=kx; // -i*(i*kx) = kx
        }

        fftw_execute(plan_b_1d); // to coordinate space

        // copy as matrix element
        for(ri=0; ri<nx; ri++) // for each row in column
            me[ri + nx*ci] = fft1[ri];
    }

    // just in case - check if hermitian
    double complex x1x2, x2x1, diff;
    int hermitian_violated=0;
    for(ri=0; ri<nx; ri++) for(ci=0; ci<nx; ci++)
    {
        x1x2 =      me[ri + nx*ci] ;
        x2x1 = conj(me[ci + nx*ri]);
        diff = x1x2 - x2x1;

        if(fabs(creal(diff))>1.0e-9 || fabs(cimag(diff))>1.0e-9)
        {
            wprintf("# ERROR[compute_matrix_elements_of_momentum_operator]: hermitian problem: %6d %6d (%f,%f) <=> (%f,%f)\n",
                ri, ci, creal(x1x2), cimag(x1x2), creal(x2x1), cimag(x2x1)
            );
            hermitian_violated=1;
        }
    }

    // TEST check with the formula
    // Both the formula and numerical approch should give EXACTLY the same result
    for(ri=0; ri<nx; ri++) for(ci=0; ci<nx; ci++)
    {
        double complex f=0.0;
        if(ri!=ci) f= M_PI/(nx*dx) * pow(-1.,ri-ci)*cos(M_PI*(ri-ci)/nx)/sin(M_PI*(ri-ci)/nx);
        f*=-1.0*I;

        if(cabs(me[ri + nx*ci]-f)>1.0e-6)
        {
            wprintf("ERROR: %d %d (%f,%f) == (%f,%f)\n", ri, ci, creal(me[ri + nx*ci]),cimag(me[ri + nx*ci]), creal(f), cimag(f));
            hermitian_violated=2;
        }
    }

    // clear
    free(fft1);
    fftw_destroy_plan(plan_f_1d);
    fftw_destroy_plan(plan_b_1d);

    return hermitian_violated;
}


int main()
{

    // Settings
    int nx = 200; // number of points in x-direction
    double Lx = 100.0; // width in x-direction

    double x0 = -Lx/2;
    double dx = Lx/nx;

    double *fx = (double *) malloc(nx*sizeof(double)); // function
    double *formula_dx_fx = (double *) malloc(nx*sizeof(double)); // derivative of function
    double *formula_dx2_fx = (double *) malloc(nx*sizeof(double)); // second derivative of function

    // Fill arrays with data
    int ix;
    double x;
    for(ix=0; ix<nx; ix++)
    {
        x = x0 + dx*ix;
        fx[ix] = function_x(x);
        formula_dx_fx[ix] = dx_function_x(x);
        formula_dx2_fx[ix] = dx2_function_x(x);
    }

    // derivative operator
    double complex *op_dx = (double complex*) malloc(nx*nx*sizeof(double complex));
    compute_matrix_elements_of_momentum_operator(nx, dx, op_dx);
    for(ix=0; ix<nx*nx; ix++) op_dx[ix]*=I; // -id/dx * i = d/dx

    // matrix vector multiplication
    double complex *dx_fx = (double complex *) malloc(nx*sizeof(double complex)); // derivative computed numerically
    int ri, ci;
    for(ri=0; ri<nx; ri++)
    {
        dx_fx[ri] = 0.0 + I*0.0;
        for(ci=0; ci<nx; ci++) dx_fx[ri]+=op_dx[ri + nx*ci]*fx[ci];
    }

    // print results
    for(ix=0; ix<nx; ix++)
    {
        x = x0 + dx*ix;
        printf("%12.6f %12.6g %12.6g %12.6g %12.6g %12.6g\n",
               x, fx[ix], formula_dx_fx[ix], formula_dx2_fx[ix], creal(dx_fx[ix]), cimag(dx_fx[ix]));
    }


    return 1;
}
