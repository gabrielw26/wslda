/**
 * W-SLDA Toolkit
 * Warsaw University of Technology, Faculty of Physics (2021)
 * https://wslda.fizyka.pw.edu.pl/
 * 
 * Tester of wderiv lib for 1D cases
 * */ 

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include "wderiv.h"

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "ERROR: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return WDERIV_ERR_CANNOT_ALLOCATE_MEMORY;                               \
    } 
    
    
int test_array_diff_r(int N, double *a, double *b, double epsilon)
{
    int ix=0;
    
    double d,d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ix=0; ix<N; ix++)
    {
        d = a[ix]-b[ix];
        
        d2=d*d;
        sumd2+=d2;
        if(d2>maxd2) maxd2=d2;
    }
    
    printf("#    COMPARISON RESULTS:\n");
    printf("#           |max[a-b]| : %16.8g\n", sqrt(maxd2));
    printf("#         SUM[(a-b)^2] : %16.8g\n", sumd2);
    printf("# SQRT(SUM[(a-b)^2])/N : %16.8g\n", sqrt(sumd2)/N);
    if(sqrt(sumd2)/N<epsilon)
    {
        printf("#          TEST RESULT : PASSED\n");
        return 1;
    }
    else
    {
        printf("#          TEST RESULT : FAILED\n");
        return 0;
    }
    
}
int test_array_diff_c(int N, double complex *a, double complex *b, double epsilon)
{
    int ix=0;
    
    double complex d;
    double d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ix=0; ix<N; ix++)
    {
        d = a[ix]-b[ix];
        d2=cabs(d)*cabs(d);
        sumd2+=d2;
        if(fabs(d2)>fabs(maxd2)) maxd2=d2;
    }
    
    printf("#    COMPARISON RESULTS:\n");
    printf("#           |max[a-b]| : %16.8g\n", sqrt(fabs(maxd2)));
    printf("#         SUM[(a-b)^2] : %16.8g\n", fabs(sumd2));
    printf("# SQRT(SUM[(a-b)^2])/N : %16.8g\n", sqrt(fabs(sumd2))/N);
    if(sqrt(fabs(sumd2))/N<epsilon)
    {
        printf("#          TEST RESULT : PASSED\n");
        return 1;
    }
    else
    {
        printf("#          TEST RESULT : FAILED\n");
        return 0;
    }
    
}

// test function
double ffx(double x, double a, double x0)
{
    return exp(a*pow(x-x0,2));
}

double dfdx(double x, double a, double x0)
{
    return ffx(x,a,x0)*a*2.*(x-x0);
}

double d2fdx2(double x, double a, double x0)
{
    return dfdx(x,a,x0)*a*2.*(x-x0) + ffx(x,a,x0)*a*2.0;
}

double d7fdx7(double x, double a, double x0)
{
    return ffx(x,a,x0)*pow(a,4)*(x-x0)*(pow(a,3)*pow(x-x0,6)*128 + pow(a,2)*pow(x-x0,4)*1344 + pow(x-x0,2)*a*3360 + 1680);
}
// complex test function
double complex ccx(double x, double a, double x0)//b, x1 for complex
{
    return exp(a*pow(x-x0,2)) + I*exp(a*pow(x-x0,2));
}

double complex dcdx(double x, double a, double x0)
{
    return creal(ccx(x,a,x0))*a*2.*(x-x0) + I*cimag(ccx(x,a,x0))*a*2.*(x-x0);
}

double complex d2cdx2(double x, double a, double x0)
{
    return creal(dcdx(x,a,x0))*a*2.*(x-x0) + creal(ccx(x,a,x0))*a*2.0 + I*(cimag(dcdx(x,a,x0))*a*2.*(x-x0) + cimag(ccx(x,a,x0))*a*2.0);
}

int main()
{
    // lattice
    int nx=128;
    
    double dx=1;
    
    // Initilize 1D lib
    wderiv_init_1d(nx, dx);
    
    double *fx, *fxDX, *fxDX2, *fxDXn, *fxGRADx, *fxGRAD2, *fxLAP, *fxDIV, *fxDXgen, *fxDX2gen, *fxDXngen, *fxGRADxgen, *fxGRAD2gen, *fxLAPgen, *fxDIVgen;
    double complex *cx, *cxDX, *cxDX2, *cxGRADx, *cxGRAD2, *cxLAP, *cxDIV, *cxDXgen, *cxDX2gen, *cxGRADxgen, *cxGRAD2gen, *cxLAPgen, *cxDIVgen;
    double *test_fxDX, *test_fxDX2, *test_fxDXn, *test_fxGRADx, *test_fxGRAD2, *test_fxLAP, *test_fxDIV;
    double complex *test_cxDX, *test_cxDX2, *test_cxGRADx, *test_cxGRAD2, *test_cxLAP, *test_cxDIV;
   

    
    //memory allocation for an output of fourier transform
    cppmallocl(fx,      nx, double);
    cppmallocl(fxDX,    nx, double);
    cppmallocl(fxDX2,   nx, double);
    cppmallocl(fxDXn,   nx, double);
    cppmallocl(fxGRADx, nx, double);
    cppmallocl(fxGRAD2, nx, double);
    cppmallocl(fxLAP,   nx, double);
    cppmallocl(fxDIV,   nx, double);

    cppmallocl(fxDXgen,    nx, double);
    cppmallocl(fxDX2gen,   nx, double);
    cppmallocl(fxDXngen,   nx, double);
    cppmallocl(fxGRADxgen, nx, double);
    cppmallocl(fxGRAD2gen, nx, double);
    cppmallocl(fxLAPgen,   nx, double);
    cppmallocl(fxDIVgen,   nx, double);

    cppmallocl(cx,      nx, double complex);
    cppmallocl(cxDX,    nx, double complex);
    cppmallocl(cxDX2,   nx, double complex);
    cppmallocl(cxGRADx, nx, double complex);
    cppmallocl(cxGRAD2, nx, double complex);
    cppmallocl(cxLAP,   nx, double complex);
    cppmallocl(cxDIV,   nx, double complex);
    
    cppmallocl(cxDXgen,    nx, double complex);
    cppmallocl(cxDX2gen,   nx, double complex);
    cppmallocl(cxGRADxgen, nx, double complex);
    cppmallocl(cxGRAD2gen, nx, double complex);
    cppmallocl(cxLAPgen,   nx, double complex);
    cppmallocl(cxDIVgen,   nx, double complex);

    //memory allocation for an output of analytical derivatives
    cppmallocl(test_fxDX,    nx, double);
    cppmallocl(test_fxDX2,   nx, double);
    cppmallocl(test_fxDXn,   nx, double);
    cppmallocl(test_fxGRADx, nx, double);
    cppmallocl(test_fxGRAD2, nx, double);
    cppmallocl(test_fxLAP,   nx, double);
    cppmallocl(test_fxDIV,   nx, double);

    cppmallocl(test_cxDX,    nx, double complex);
    cppmallocl(test_cxDX2,   nx, double complex);
    cppmallocl(test_cxGRADx, nx, double complex);
    cppmallocl(test_cxGRAD2, nx, double complex);
    cppmallocl(test_cxLAP,   nx, double complex);
    cppmallocl(test_cxDIV,   nx, double complex);
    
    int ix=0;
    const int n_deriv=7;
    double ax=-0.1;
    for(int ix=0; ix<nx; ix++)
    {
        fx[ix] = ffx(dx*ix, ax, dx*nx/2); // domain [0,nx*dx]
        cx[ix] = ccx(dx*ix, ax, dx*nx/2);

        // dx, dy, dz derivatives analytical
        test_fxDX[ix] = dfdx(dx*ix, ax, dx*nx/2); 
        test_fxDX2[ix] = d2fdx2(dx*ix, ax, dx*nx/2); 
        test_fxDXn[ix] = d7fdx7(dx*ix, ax, dx*nx/2); 
        test_fxGRADx[ix] = test_fxDX[ix];
        test_fxGRAD2[ix] = test_fxDX[ix] * test_fxDX[ix];
        test_fxLAP[ix] = test_fxDX2[ix];
        test_fxDIV[ix] = test_fxDX[ix];

        
        test_cxDX[ix] = dcdx(dx*ix, ax, dx*nx/2); 
        test_cxDX2[ix] = d2cdx2(dx*ix, ax, dx*nx/2); 
        test_cxGRADx[ix] = test_cxDX[ix];
        test_cxGRAD2[ix] = test_cxDX[ix] * conj(test_cxDX[ix]);
        test_cxLAP[ix] = test_cxDX2[ix];
        test_cxDIV[ix] = test_cxDX[ix];

    }
    
    printf("\n1D REAL\n");
    printf("df/dx\n");
    wderiv_dfdx_1d_r(fx, fxDX);
    test_array_diff_r(nx, fxDX, test_fxDX, 1.0e-12);
    wderiv_dfdx(1, 'r', fx,fxDXgen);
    test_array_diff_r(nx, fxDXgen, test_fxDX, 1.0e-12);
    
    printf("\nd2f/dx2\n");
    wderiv_d2fdx2_1d_r(fx, fxDX2);
    test_array_diff_r(nx, fxDX2, test_fxDX2, 1.0e-12);
    wderiv_d2fdx2(1, 'r', fx,fxDX2gen);
    test_array_diff_r(nx, fxDX2gen, test_fxDX2, 1.0e-12);

    printf("\ndnf/dxn\n");
    wderiv_dnfdxn_1d_r(n_deriv, fx, fxDXn);
    test_array_diff_r(nx, fxDXn, test_fxDXn, 1.0e-12);
    wderiv_dnfdxn(1, 'r', 7, fx,fxDXngen);
    test_array_diff_r(nx, fxDXngen, test_fxDXn, 1.0e-12);

    printf("\nGradient\n");
    wderiv_gradient_1d_r(fx, fxGRADx);
    test_array_diff_r(nx, fxGRADx, test_fxGRADx, 1.0e-12);
    wderiv_gradient(1, 'r', fx, fxGRADxgen, NULL, NULL);
    test_array_diff_r(nx, fxGRADxgen, test_fxGRADx, 1.0e-12);

    printf("\nGradient squared\n");
    wderiv_gradient2_1d_r(fx, fxGRAD2);
    test_array_diff_r(nx, fxGRAD2, test_fxGRAD2, 1.0e-12);
    wderiv_gradient2(1, 'r', fx,fxGRAD2gen);
    test_array_diff_r(nx, fxGRAD2gen, test_fxGRAD2, 1.0e-12);

    printf("\nLaplacian\n");
    wderiv_laplace_1d_r(fx, fxLAP);
    test_array_diff_r(nx, fxLAP, test_fxLAP, 1.0e-12);
    wderiv_laplace(1, 'r', fx, fxLAPgen);
    test_array_diff_r(nx, fxLAPgen, test_fxLAP, 1.0e-12);

    printf("\nDivergence\n");
    wderiv_divergence_1d_r(fx, fxDIV);
    test_array_diff_r(nx, fxDIV, test_fxDIV, 1.0e-12);
    wderiv_divergence(1, 'r', fx, fxDIVgen);
    test_array_diff_r(nx, fxDIVgen, test_fxDIV, 1.0e-12);

    printf("\n/----------1D COMPLEX---------/\n");
    printf("dc/dx\n");
    wderiv_dfdx_1d_c(cx, cxDX);
    test_array_diff_c(nx, cxDX, test_cxDX, 1.0e-12);
    wderiv_dfdx(1, 'c', cx,cxDXgen);
    test_array_diff_c(nx, cxDXgen, test_cxDX, 1.0e-12);
    
    printf("\nd2c/dx2\n");
    wderiv_d2fdx2_1d_c(cx, cxDX2);
    test_array_diff_c(nx, cxDX2, test_cxDX2, 1.0e-12);
    wderiv_d2fdx2(1, 'c', cx,cxDX2gen);
    test_array_diff_c(nx, cxDX2gen, test_cxDX2, 1.0e-12);

    printf("\nGradient\n");
    wderiv_gradient_1d_c(cx, cxGRADx);
    test_array_diff_c(nx, cxGRADx, test_cxGRADx, 1.0e-12);
    wderiv_gradient(1, 'c', cx, cxGRADxgen, NULL, NULL);
    test_array_diff_c(nx, cxGRADxgen, test_cxGRADx, 1.0e-12);

    printf("\nGradient squared\n");
    wderiv_gradient2_1d_c(cx, cxGRAD2);
    test_array_diff_c(nx, cxGRAD2, test_cxGRAD2, 1.0e-12);
    wderiv_gradient2(1, 'c', cx,cxGRAD2gen);
    test_array_diff_c(nx, cxGRAD2gen, test_cxGRAD2, 1.0e-12);

    printf("\nLaplacian\n");
    wderiv_laplace_1d_c(cx, cxLAP);
    test_array_diff_c(nx, cxLAP, test_cxLAP, 1.0e-12);
    wderiv_laplace(1, 'c', cx, cxLAPgen);
    test_array_diff_c(nx, cxLAPgen, test_cxLAP, 1.0e-12);

    printf("\nDivergence\n");
    wderiv_divergence_1d_c(cx, cxDIV);
    test_array_diff_c(nx, cxDIV, test_cxDIV, 1.0e-12);
    wderiv_divergence(1, 'c', cx, cxDIVgen);
    test_array_diff_c(nx, cxDIVgen, test_cxDIV, 1.0e-12);
   
    wderiv_clean();
    free(fx); free(fxDX); free(fxDX2); free(fxDXn); free(fxGRAD2); free(fxLAP);
    free(cx); free(cxDX); free(cxDX2); free(cxGRAD2); free(cxLAP);
    free(test_fxDX); free(test_fxDX2); free(test_fxDXn); free(test_fxGRAD2); free(test_fxLAP);
    free(test_cxDX); free(test_cxDX2); free(test_cxGRAD2); free(test_cxLAP);
    return 0; 
}

