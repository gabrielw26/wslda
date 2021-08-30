/**
 * W-SLDA Toolkit
 * Warsaw University of Technology, Faculty of Physics (2021)
 * https://wslda.fizyka.pw.edu.pl/
 * 
 * Tester of wderiv lib for 2D cases
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
    int ixy=0;
    
    double d,d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ixy=0; ixy<N; ixy++)
    {
        d = a[ixy]-b[ixy];
        
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
    int ixy=0;
    
    double complex d;
    double d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ixy=0; ixy<N; ixy++)
    {
        d = a[ixy]-b[ixy];
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
double fx(double x, double a, double x0)
{
    return exp(a*pow(x-x0,2));
}

double dfdx(double x, double a, double x0)
{
    return fx(x,a,x0)*a*2.*(x-x0);
}

double d2fdx2(double x, double a, double x0)
{
    return dfdx(x,a,x0)*a*2.*(x-x0) + fx(x,a,x0)*a*2.0;
}

double d7fdx7(double x, double a, double x0)
{
    return fx(x,a,x0)*pow(a,4)*(x-x0)*(pow(a,3)*pow(x-x0,6)*128 + pow(a,2)*pow(x-x0,4)*1344 + pow(x-x0,2)*a*3360 + 1680);
}
// complex test function
double complex cx(double x, double a, double x0)//b, x1 for complex
{
    return exp(a*pow(x-x0,2)) + I*exp(a*pow(x-x0,2));
}

double complex dcdx(double x, double a, double x0)
{
    return creal(cx(x,a,x0))*a*2.*(x-x0) + I*cimag(cx(x,a,x0))*a*2.*(x-x0);
}

double complex d2cdx2(double x, double a, double x0)
{
    return creal(dcdx(x,a,x0))*a*2.*(x-x0) + creal(cx(x,a,x0))*a*2.0 + I*(cimag(dcdx(x,a,x0))*a*2.*(x-x0) + cimag(cx(x,a,x0))*a*2.0);
}

int main()
{
    // lattice
    int nx=128;
    int ny=170;
    
    double dx=0.9;
    double dy=1.1;
    
    // Initilize 2D lib
    wderiv_init_2d(nx, ny, dx, dy);
    
    double *fxy, *fxyDX, *fxyDY, *fxyDX2, *fxyDY2, *fxyDXn, *fxyDYn, *fxyGRADx, *fxyGRADy, *fxyGRAD2, *fxyLAP, *fxyDIV, *fxyDXgen, *fxyDYgen, *fxyDX2gen, *fxyDY2gen, *fxyDXngen, *fxyDYngen, *fxyGRADxgen, *fxyGRADygen, *fxyGRAD2gen, *fxyLAPgen, *fxyDIVgen;
    double complex *cxy, *cxyDX, *cxyDY, *cxyDX2, *cxyDY2, *cxyGRADx, *cxyGRADy, *cxyGRAD2, *cxyLAP, *cxyDIV, *cxyDXgen, *cxyDYgen, *cxyDX2gen, *cxyDY2gen, *cxyGRADxgen, *cxyGRADygen, *cxyGRAD2gen, *cxyLAPgen, *cxyDIVgen;
    double *test_fxyDX, *test_fxyDY, *test_fxyDX2, *test_fxyDY2, *test_fxyDXn, *test_fxyDYn, *test_fxyGRAD2, *test_fxyLAP, *test_fxyDIV;
    double complex *test_cxyDX, *test_cxyDY, *test_cxyDX2, *test_cxyDY2, *test_cxyGRAD2, *test_cxyLAP, *test_cxyDIV;
   

    
    //memory allocation for an output of fourier transform
    cppmallocl(fxy,      nx*ny, double);
    cppmallocl(fxyDX,    nx*ny, double);
    cppmallocl(fxyDY,    nx*ny, double);
    cppmallocl(fxyDX2,   nx*ny, double);
    cppmallocl(fxyDY2,   nx*ny, double);
    cppmallocl(fxyDXn,   nx*ny, double);
    cppmallocl(fxyDYn,   nx*ny, double);
    cppmallocl(fxyGRADx, nx*ny, double);
    cppmallocl(fxyGRADy, nx*ny, double);
    cppmallocl(fxyGRAD2, nx*ny, double);
    cppmallocl(fxyLAP,   nx*ny, double);
    cppmallocl(fxyDIV,   nx*ny, double);
    
    cppmallocl(fxyDXgen,    nx*ny, double);
    cppmallocl(fxyDYgen,    nx*ny, double);
    cppmallocl(fxyDX2gen,   nx*ny, double);
    cppmallocl(fxyDY2gen,   nx*ny, double);
    cppmallocl(fxyDXngen,   nx*ny, double);
    cppmallocl(fxyDYngen,   nx*ny, double);
    cppmallocl(fxyGRADxgen, nx*ny, double);
    cppmallocl(fxyGRADygen, nx*ny, double);
    cppmallocl(fxyGRAD2gen, nx*ny, double);
    cppmallocl(fxyLAPgen,   nx*ny, double);
    cppmallocl(fxyDIVgen,   nx*ny, double);

    cppmallocl(cxy,      nx*ny, double complex);
    cppmallocl(cxyDX,    nx*ny, double complex);
    cppmallocl(cxyDY,    nx*ny, double complex);
    cppmallocl(cxyDX2,   nx*ny, double complex);
    cppmallocl(cxyDY2,   nx*ny, double complex);
    cppmallocl(cxyGRADx, nx*ny, double complex);
    cppmallocl(cxyGRADy, nx*ny, double complex);
    cppmallocl(cxyGRAD2, nx*ny, double complex);
    cppmallocl(cxyLAP,   nx*ny, double complex);
    cppmallocl(cxyDIV,   nx*ny, double complex);
    
    cppmallocl(cxyDXgen,    nx*ny, double complex);
    cppmallocl(cxyDYgen,    nx*ny, double complex);
    cppmallocl(cxyDX2gen,   nx*ny, double complex);
    cppmallocl(cxyDY2gen,   nx*ny, double complex);
    cppmallocl(cxyGRADxgen, nx*ny, double complex);
    cppmallocl(cxyGRADygen, nx*ny, double complex);
    cppmallocl(cxyGRAD2gen, nx*ny, double complex);
    cppmallocl(cxyLAPgen,   nx*ny, double complex);
    cppmallocl(cxyDIVgen,   nx*ny, double complex);

    //memory allocation for an output of analytical derivatives
    cppmallocl(test_fxyDX,    nx*ny, double);
    cppmallocl(test_fxyDY,    nx*ny, double);
    cppmallocl(test_fxyDX2,   nx*ny, double);
    cppmallocl(test_fxyDY2,   nx*ny, double);
    cppmallocl(test_fxyDXn,   nx*ny, double);
    cppmallocl(test_fxyDYn,   nx*ny, double);
    cppmallocl(test_fxyGRAD2, nx*ny, double);
    cppmallocl(test_fxyLAP,   nx*ny, double);
    cppmallocl(test_fxyDIV,   nx*ny, double);

    cppmallocl(test_cxyDX,    nx*ny, double complex);
    cppmallocl(test_cxyDY,    nx*ny, double complex);
    cppmallocl(test_cxyDX2,   nx*ny, double complex);
    cppmallocl(test_cxyDY2,   nx*ny, double complex);
    cppmallocl(test_cxyGRAD2, nx*ny, double complex);
    cppmallocl(test_cxyLAP,   nx*ny, double complex);
    cppmallocl(test_cxyDIV,   nx*ny, double complex);
    
    int check=0;
    int ixy=0;
    const int n_deriv=7;
    double ax=-0.1, ay=-0.09;
    for(int ix=0; ix<nx; ix++) for(int iy=0; iy<ny; iy++)
    {
        fxy[ixy] = fx(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2); // domain [0,nx*dx] x [0,ny*dy] x [0,nz*dz] 
        cxy[ixy] = cx(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2);

        // dx, dy, dz derivatives analytical
        test_fxyDX[ixy] = dfdx(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2); 
        test_fxyDY[ixy] = fx(dx*ix, ax, dx*nx/2)*dfdx(dy*iy, ay, dy*ny/2); 
        test_fxyDX2[ixy] = d2fdx2(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2); 
        test_fxyDY2[ixy] = fx(dx*ix, ax, dx*nx/2)*d2fdx2(dy*iy, ay, dy*ny/2); 
        test_fxyDXn[ixy] = d7fdx7(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2); 
        test_fxyDYn[ixy] = fx(dx*ix, ax, dx*nx/2)*d7fdx7(dy*iy, ay, dy*ny/2); 
        test_fxyGRAD2[ixy] = test_fxyDX[ixy] * test_fxyDX[ixy] + test_fxyDY[ixy] * test_fxyDY[ixy];
        test_fxyLAP[ixy] = test_fxyDX2[ixy] + test_fxyDY2[ixy];
        test_fxyDIV[ixy] = test_fxyDX[ixy] + test_fxyDY[ixy];

        
        test_cxyDX[ixy] = dcdx(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2); 
        test_cxyDY[ixy] = cx(dx*ix, ax, dx*nx/2)*dcdx(dy*iy, ay, dy*ny/2); 
        test_cxyDX2[ixy] = d2cdx2(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2); 
        test_cxyDY2[ixy] = cx(dx*ix, ax, dx*nx/2)*d2cdx2(dy*iy, ay, dy*ny/2); 
        test_cxyGRAD2[ixy] = test_cxyDX[ixy] * conj(test_cxyDX[ixy]) + test_cxyDY[ixy] * conj(test_cxyDY[ixy]);
        test_cxyLAP[ixy] = test_cxyDX2[ixy] + test_cxyDY2[ixy];
        test_cxyDIV[ixy] = test_cxyDX[ixy] + test_cxyDY[ixy];

        ixy++; 
    }
    
    printf("\n2D REAL\n");
    printf("df/dx,df/dy\n");
    wderiv_dfdx_2d_r(fxy, fxyDX);
    test_array_diff_r(nx*ny, fxyDX, test_fxyDX, 1.0e-12);
    wderiv_dfdx(2, 'r', fxy, fxyDXgen);
    test_array_diff_r(nx*ny, fxyDXgen, test_fxyDX, 1.0e-12);
    wderiv_dfdy_2d_r(fxy, fxyDY);
    test_array_diff_r(nx*ny, fxyDY, test_fxyDY, 1.0e-12);
    wderiv_dfdy(2, 'r', fxy, fxyDYgen);
    test_array_diff_r(nx*ny, fxyDYgen, test_fxyDY, 1.0e-12);
    
    printf("\nd2f/dx2,d2f/dy2\n");
    wderiv_d2fdx2_2d_r(fxy, fxyDX2);
    test_array_diff_r(nx*ny, fxyDX2, test_fxyDX2, 1.0e-12);
    wderiv_d2fdx2(2, 'r', fxy, fxyDX2gen);
    test_array_diff_r(nx*ny, fxyDX2gen, test_fxyDX2, 1.0e-12);
    wderiv_d2fdy2_2d_r(fxy, fxyDY2);
    test_array_diff_r(nx*ny, fxyDY2, test_fxyDY2, 1.0e-12);
    wderiv_d2fdy2(2, 'r', fxy, fxyDY2gen);
    test_array_diff_r(nx*ny, fxyDY2gen, test_fxyDY2, 1.0e-12);

    printf("\ndnf/dxn,dnf/dyn\n");
    wderiv_dnfdxn_2d_r(n_deriv, fxy, fxyDXn);
    test_array_diff_r(nx*ny, fxyDXn, test_fxyDXn, 1.0e-12);
    wderiv_dnfdxn(2, 'r', 7, fxy, fxyDXngen);
    test_array_diff_r(nx*ny, fxyDXngen, test_fxyDXn, 1.0e-12);
    wderiv_dnfdyn_2d_r(n_deriv, fxy, fxyDYn);
    test_array_diff_r(nx*ny, fxyDYn, test_fxyDYn, 1.0e-12);
    wderiv_dnfdyn(2, 'r', 7, fxy, fxyDYngen);
    test_array_diff_r(nx*ny, fxyDYngen, test_fxyDYn, 1.0e-12);

    printf("\nGradient\n");
    wderiv_gradient_2d_r(fxy, fxyGRADx, fxyGRADy);
    test_array_diff_r(nx*ny, fxyGRADx, test_fxyDX, 1.0e-12);
    test_array_diff_r(nx*ny, fxyGRADy, test_fxyDY, 1.0e-12);
    wderiv_gradient(2, 'r', fxy, fxyGRADxgen, fxyGRADygen, NULL);
    test_array_diff_r(nx*ny, fxyGRADxgen, test_fxyDX, 1.0e-12);
    test_array_diff_r(nx*ny, fxyGRADygen, test_fxyDY, 1.0e-12);

    printf("\nGradient squared\n");
    wderiv_gradient2_2d_r(fxy, fxyGRAD2);
    test_array_diff_r(nx*ny, fxyGRAD2, test_fxyGRAD2, 1.0e-12);
    wderiv_gradient2(2, 'r', fxy, fxyGRAD2gen);
    test_array_diff_r(nx*ny, fxyGRAD2gen, test_fxyGRAD2, 1.0e-12);

    printf("\nLaplacian\n");
    wderiv_laplace_2d_r(fxy, fxyLAP);
    test_array_diff_r(nx*ny, fxyLAP, test_fxyLAP, 1.0e-12);
    wderiv_laplace(2, 'r', fxy, fxyLAPgen);
    test_array_diff_r(nx*ny, fxyLAPgen, test_fxyLAP, 1.0e-12);

    printf("\nDivergence\n");
    wderiv_divergence_2d_r(fxy, fxyDIV);
    test_array_diff_r(nx*ny, fxyDIV, test_fxyDIV, 1.0e-12);
    wderiv_divergence(2, 'r', fxy, fxyDIVgen);
    test_array_diff_r(nx*ny, fxyDIVgen, test_fxyDIV, 1.0e-12);

    printf("\n/----------2D COMPLEX---------/\n");
    printf("dc/dx,dc/dy\n");
    wderiv_dfdx_2d_c(cxy, cxyDX);
    test_array_diff_c(nx*ny, cxyDX, test_cxyDX, 1.0e-12);
    wderiv_dfdx(2, 'c', cxy, cxyDXgen);
    test_array_diff_c(nx*ny, cxyDXgen, test_cxyDX, 1.0e-12);
    wderiv_dfdy_2d_c(cxy, cxyDY);
    test_array_diff_c(nx*ny, cxyDY, test_cxyDY, 1.0e-12);
    wderiv_dfdy(2, 'c', cxy, cxyDYgen);
    test_array_diff_c(nx*ny, cxyDYgen, test_cxyDY, 1.0e-12);
    
    printf("\nd2c/dx2,d2c/dy2,d2c/dz2\n");
    wderiv_d2fdx2_2d_c(cxy, cxyDX2);
    test_array_diff_c(nx*ny, cxyDX2, test_cxyDX2, 1.0e-12);
    wderiv_d2fdx2(2, 'c', cxy, cxyDX2gen);
    test_array_diff_c(nx*ny, cxyDX2gen, test_cxyDX2, 1.0e-12);
    wderiv_d2fdy2_2d_c(cxy, cxyDY2);
    test_array_diff_c(nx*ny, cxyDY2, test_cxyDY2, 1.0e-12);
    wderiv_d2fdy2(2, 'c', cxy, cxyDY2gen);
    test_array_diff_c(nx*ny, cxyDY2gen, test_cxyDY2, 1.0e-12);

    printf("\nGradient\n");
    wderiv_gradient_2d_c(cxy, cxyGRADx, cxyGRADy);
    test_array_diff_c(nx*ny, cxyGRADx, test_cxyDX, 1.0e-12);
    test_array_diff_c(nx*ny, cxyGRADy, test_cxyDY, 1.0e-12);
    wderiv_gradient(2, 'c', cxy, cxyGRADxgen, cxyGRADygen, NULL);
    test_array_diff_c(nx*ny, cxyGRADxgen, test_cxyDX, 1.0e-12);
    test_array_diff_c(nx*ny, cxyGRADygen, test_cxyDY, 1.0e-12);

    printf("\nGradient squared\n");
    wderiv_gradient2_2d_c(cxy, cxyGRAD2);
    test_array_diff_c(nx*ny, cxyGRAD2, test_cxyGRAD2, 1.0e-12);
    wderiv_gradient2(2, 'c', cxy, cxyGRAD2gen);
    test_array_diff_c(nx*ny, cxyGRAD2gen, test_cxyGRAD2, 1.0e-12);

    printf("\nLaplacian\n");
    wderiv_laplace_2d_c(cxy, cxyLAP);
    test_array_diff_c(nx*ny, cxyLAP, test_cxyLAP, 1.0e-12);
    wderiv_laplace(2, 'c', cxy, cxyLAPgen);
    test_array_diff_c(nx*ny, cxyLAPgen, test_cxyLAP, 1.0e-12);

    printf("\nDivergence\n");
    wderiv_divergence_2d_c(cxy, cxyDIV);
    test_array_diff_c(nx*ny, cxyDIV, test_cxyDIV, 1.0e-12);
    wderiv_divergence(2, 'c', cxy, cxyDIVgen);
    test_array_diff_c(nx*ny, cxyDIVgen, test_cxyDIV, 1.0e-12);
   
    //    printf("x analytical_dcdx numerical_dcdx\n");
    // ixy=0;
    // for(int ix=0; ix<nx; ix++) for(int iy=0; iy<ny; iy++)
    // {
    //     if(iy==ny/2) 
    //         // printf("%d\t%f\t%f\t%f\t%f\n", ix, creal(test_cxyLAP[ixy]), cimag(test_cxyLAP[ixy]), creal(cxyLAP[ixy]), cimag(cxyLAP[ixy]));
    //         printf("%d\t%f\t%f\n", ix, test_fxyDX[ixy], fxyDX[ixy]);
    //     ixy++; 
    // }
    // destroy lib
    wderiv_clean();
    free(fxy); free(fxyDX); free(fxyDY); free(fxyDX2); free(fxyDY2); free(fxyDXn); free(fxyDYn); free(fxyGRAD2); free(fxyLAP);
    free(cxy); free(cxyDX); free(cxyDY); free(cxyDX2); free(cxyDY2); free(cxyGRAD2); free(cxyLAP);
    free(test_fxyDX); free(test_fxyDY); free(test_fxyDX2); free(test_fxyDY2); free(test_fxyDXn); free(test_fxyDYn); free(test_fxyGRAD2); free(test_fxyLAP);
    free(test_cxyDX); free(test_cxyDY); free(test_cxyDX2); free(test_cxyDY2); free(test_cxyGRAD2); free(test_cxyLAP);
    return 0; 
}

