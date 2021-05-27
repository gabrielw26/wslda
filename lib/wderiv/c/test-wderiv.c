/**
 * W-SLDA Toolkit
 * Warsaw University of Technology, Faculty of Physics (2021)
 * https://wslda.fizyka.pw.edu.pl/
 * 
 * Tester of wderiv lib
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
    int ixyz=0;
    
    double complex d;
    double d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ixyz=0; ixyz<N; ixyz++)
    {
        d = a[ixyz]-b[ixyz];
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
    int nx=80;
    int ny=70;
    int nz=90;
    
    double dx=1.2;
    double dy=1.1;
    double dz=1.3;
    
    // Initilize 3D lib
    wderiv_init_3d(nx, ny, nz, dx, dy, dz);
    
    double *fxyz, *fxyzDX, *fxyzDY, *fxyzDZ, *fxyzDX2, *fxyzDY2, *fxyzDZ2, *fxyzDXn, *fxyzDYn, *fxyzDZn, *fxyzGRADx, *fxyzGRADy, *fxyzGRADz, *fxyzGRAD2, *fxyzLAP, *fxyzDIV, *fxyzCURLx, *fxyzCURLy, *fxyzCURLz, *fxyzDXgen, *fxyzDYgen, *fxyzDZgen, *fxyzDX2gen, *fxyzDY2gen, *fxyzDZ2gen, *fxyzDXngen, *fxyzDYngen, *fxyzDZngen, *fxyzGRADxgen, *fxyzGRADygen, *fxyzGRADzgen, *fxyzGRAD2gen, *fxyzLAPgen, *fxyzDIVgen, *fxyzCURLxgen, *fxyzCURLygen, *fxyzCURLzgen;
    double complex *cxyz, *cxyzDX, *cxyzDY, *cxyzDZ, *cxyzDX2, *cxyzDY2, *cxyzDZ2, *cxyzGRADx, *cxyzGRADy, *cxyzGRADz, *cxyzGRAD2, *cxyzLAP, *cxyzDIV, *cxyzCURLx, *cxyzCURLy, *cxyzCURLz, *cxyzDXgen, *cxyzDYgen, *cxyzDZgen, *cxyzDX2gen, *cxyzDY2gen, *cxyzDZ2gen, *cxyzGRADxgen, *cxyzGRADygen, *cxyzGRADzgen, *cxyzGRAD2gen, *cxyzLAPgen, *cxyzDIVgen, *cxyzCURLxgen, *cxyzCURLygen, *cxyzCURLzgen;
    double *test_fxyzDX, *test_fxyzDY, *test_fxyzDZ, *test_fxyzDX2, *test_fxyzDY2, *test_fxyzDZ2, *test_fxyzDXn, *test_fxyzDYn, *test_fxyzDZn, *test_fxyzGRAD2, *test_fxyzLAP, *test_fxyzDIV, *test_fxyzCURLx, *test_fxyzCURLy, *test_fxyzCURLz;
    double complex *test_cxyzDX, *test_cxyzDY, *test_cxyzDZ, *test_cxyzDX2, *test_cxyzDY2, *test_cxyzDZ2, *test_cxyzGRAD2, *test_cxyzLAP, *test_cxyzDIV, *test_cxyzCURLx, *test_cxyzCURLy, *test_cxyzCURLz;
   

    
    //memory allocation for an output of fourier transform
    cppmallocl(fxyz,      nx*ny*nz, double);
    cppmallocl(fxyzDX,    nx*ny*nz, double);
    cppmallocl(fxyzDY,    nx*ny*nz, double);
    cppmallocl(fxyzDZ,    nx*ny*nz, double);
    cppmallocl(fxyzDX2,   nx*ny*nz, double);
    cppmallocl(fxyzDY2,   nx*ny*nz, double);
    cppmallocl(fxyzDZ2,   nx*ny*nz, double);
    cppmallocl(fxyzDXn,   nx*ny*nz, double);
    cppmallocl(fxyzDYn,   nx*ny*nz, double);
    cppmallocl(fxyzDZn,   nx*ny*nz, double);
    cppmallocl(fxyzGRADx, nx*ny*nz, double);
    cppmallocl(fxyzGRADy, nx*ny*nz, double);
    cppmallocl(fxyzGRADz, nx*ny*nz, double);
    cppmallocl(fxyzGRAD2, nx*ny*nz, double);
    cppmallocl(fxyzLAP,   nx*ny*nz, double);
    cppmallocl(fxyzDIV,   nx*ny*nz, double);
    cppmallocl(fxyzCURLx, nx*ny*nz, double);
    cppmallocl(fxyzCURLy, nx*ny*nz, double);
    cppmallocl(fxyzCURLz, nx*ny*nz, double);

    cppmallocl(fxyzDXgen,    nx*ny*nz, double);
    cppmallocl(fxyzDYgen,    nx*ny*nz, double);
    cppmallocl(fxyzDZgen,    nx*ny*nz, double);
    cppmallocl(fxyzDX2gen,   nx*ny*nz, double);
    cppmallocl(fxyzDY2gen,   nx*ny*nz, double);
    cppmallocl(fxyzDZ2gen,   nx*ny*nz, double);
    cppmallocl(fxyzDXngen,   nx*ny*nz, double);
    cppmallocl(fxyzDYngen,   nx*ny*nz, double);
    cppmallocl(fxyzDZngen,   nx*ny*nz, double);
    cppmallocl(fxyzGRADxgen, nx*ny*nz, double);
    cppmallocl(fxyzGRADygen, nx*ny*nz, double);
    cppmallocl(fxyzGRADzgen, nx*ny*nz, double);
    cppmallocl(fxyzGRAD2gen, nx*ny*nz, double);
    cppmallocl(fxyzLAPgen,   nx*ny*nz, double);
    cppmallocl(fxyzDIVgen,   nx*ny*nz, double);
    cppmallocl(fxyzCURLxgen, nx*ny*nz, double);
    cppmallocl(fxyzCURLygen, nx*ny*nz, double);
    cppmallocl(fxyzCURLzgen, nx*ny*nz, double);

    cppmallocl(cxyz,      nx*ny*nz, double complex);
    cppmallocl(cxyzDX,    nx*ny*nz, double complex);
    cppmallocl(cxyzDY,    nx*ny*nz, double complex);
    cppmallocl(cxyzDZ,    nx*ny*nz, double complex);
    cppmallocl(cxyzDX2,   nx*ny*nz, double complex);
    cppmallocl(cxyzDY2,   nx*ny*nz, double complex);
    cppmallocl(cxyzDZ2,   nx*ny*nz, double complex);
    cppmallocl(cxyzGRADx, nx*ny*nz, double complex);
    cppmallocl(cxyzGRADy, nx*ny*nz, double complex);
    cppmallocl(cxyzGRADz, nx*ny*nz, double complex);
    cppmallocl(cxyzGRAD2, nx*ny*nz, double complex);
    cppmallocl(cxyzLAP,   nx*ny*nz, double complex);
    cppmallocl(cxyzDIV,   nx*ny*nz, double complex);
    cppmallocl(cxyzCURLx, nx*ny*nz, double complex);
    cppmallocl(cxyzCURLy, nx*ny*nz, double complex);
    cppmallocl(cxyzCURLz, nx*ny*nz, double complex);

    cppmallocl(cxyzDXgen,    nx*ny*nz, double complex);
    cppmallocl(cxyzDYgen,    nx*ny*nz, double complex);
    cppmallocl(cxyzDZgen,    nx*ny*nz, double complex);
    cppmallocl(cxyzDX2gen,   nx*ny*nz, double complex);
    cppmallocl(cxyzDY2gen,   nx*ny*nz, double complex);
    cppmallocl(cxyzDZ2gen,   nx*ny*nz, double complex);
    cppmallocl(cxyzGRADxgen, nx*ny*nz, double complex);
    cppmallocl(cxyzGRADygen, nx*ny*nz, double complex);
    cppmallocl(cxyzGRADzgen, nx*ny*nz, double complex);
    cppmallocl(cxyzGRAD2gen, nx*ny*nz, double complex);
    cppmallocl(cxyzLAPgen,   nx*ny*nz, double complex);
    cppmallocl(cxyzDIVgen,   nx*ny*nz, double complex);
    cppmallocl(cxyzCURLxgen, nx*ny*nz, double complex);
    cppmallocl(cxyzCURLygen, nx*ny*nz, double complex);
    cppmallocl(cxyzCURLzgen, nx*ny*nz, double complex);

    //memory allocation for an output of analytical derivatives
    cppmallocl(test_fxyzDX,    nx*ny*nz, double);
    cppmallocl(test_fxyzDY,    nx*ny*nz, double);
    cppmallocl(test_fxyzDZ,    nx*ny*nz, double);
    cppmallocl(test_fxyzDX2,   nx*ny*nz, double);
    cppmallocl(test_fxyzDY2,   nx*ny*nz, double);
    cppmallocl(test_fxyzDZ2,   nx*ny*nz, double);
    cppmallocl(test_fxyzDXn,   nx*ny*nz, double);
    cppmallocl(test_fxyzDYn,   nx*ny*nz, double);
    cppmallocl(test_fxyzDZn,   nx*ny*nz, double);
    cppmallocl(test_fxyzGRAD2, nx*ny*nz, double);
    cppmallocl(test_fxyzLAP,   nx*ny*nz, double);
    cppmallocl(test_fxyzDIV,   nx*ny*nz, double);
    cppmallocl(test_fxyzCURLx,  nx*ny*nz, double);
    cppmallocl(test_fxyzCURLy,  nx*ny*nz, double);
    cppmallocl(test_fxyzCURLz,  nx*ny*nz, double);

    cppmallocl(test_cxyzDX,    nx*ny*nz, double complex);
    cppmallocl(test_cxyzDY,    nx*ny*nz, double complex);
    cppmallocl(test_cxyzDZ,    nx*ny*nz, double complex);
    cppmallocl(test_cxyzDX2,   nx*ny*nz, double complex);
    cppmallocl(test_cxyzDY2,   nx*ny*nz, double complex);
    cppmallocl(test_cxyzDZ2,   nx*ny*nz, double complex);
    cppmallocl(test_cxyzGRAD2, nx*ny*nz, double complex);
    cppmallocl(test_cxyzLAP,   nx*ny*nz, double complex);
    cppmallocl(test_cxyzDIV,   nx*ny*nz, double complex);
    cppmallocl(test_cxyzCURLx,  nx*ny*nz, double complex);
    cppmallocl(test_cxyzCURLy,  nx*ny*nz, double complex);
    cppmallocl(test_cxyzCURLz,  nx*ny*nz, double complex);
    
    int check=0;
    int ixyz=0;
    const int n_deriv=7;
    double ax=-0.1, ay=-0.09, az=-0.07;
    for(int ix=0; ix<nx; ix++) for(int iy=0; iy<ny; iy++) for(int iz=0; iz<nz; iz++)
    {
        fxyz[ixyz] = fx(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2)*fx(dz*iz, az, dz*nz/2); // domain [0,nx*dx] x [0,ny*dy] x [0,nz*dz] 
        cxyz[ixyz] = cx(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2)*cx(dz*iz, az, dz*nz/2);

        // dx, dy, dz derivatives analytical
        test_fxyzDX[ixyz] = dfdx(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2)*fx(dz*iz, az, dz*nz/2); 
        test_fxyzDY[ixyz] = fx(dx*ix, ax, dx*nx/2)*dfdx(dy*iy, ay, dy*ny/2)*fx(dz*iz, az, dz*nz/2); 
        test_fxyzDZ[ixyz] = fx(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2)*dfdx(dz*iz, az, dz*nz/2); 
        test_fxyzDX2[ixyz] = d2fdx2(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2)*fx(dz*iz, az, dz*nz/2); 
        test_fxyzDY2[ixyz] = fx(dx*ix, ax, dx*nx/2)*d2fdx2(dy*iy, ay, dy*ny/2)*fx(dz*iz, az, dz*nz/2); 
        test_fxyzDZ2[ixyz] = fx(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2)*d2fdx2(dz*iz, az, dz*nz/2); 
        test_fxyzDXn[ixyz] = d7fdx7(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2)*fx(dz*iz, az, dz*nz/2); 
        test_fxyzDYn[ixyz] = fx(dx*ix, ax, dx*nx/2)*d7fdx7(dy*iy, ay, dy*ny/2)*fx(dz*iz, az, dz*nz/2); 
        test_fxyzDZn[ixyz] = fx(dx*ix, ax, dx*nx/2)*fx(dy*iy, ay, dy*ny/2)*d7fdx7(dz*iz, az, dz*nz/2); 
        test_fxyzGRAD2[ixyz] = test_fxyzDX[ixyz] * test_fxyzDX[ixyz] + test_fxyzDY[ixyz] * test_fxyzDY[ixyz] + test_fxyzDZ[ixyz] * test_fxyzDZ[ixyz];
        test_fxyzLAP[ixyz] = test_fxyzDX2[ixyz] + test_fxyzDY2[ixyz] + test_fxyzDZ2[ixyz];
        test_fxyzDIV[ixyz] = test_fxyzDX[ixyz] + test_fxyzDY[ixyz] + test_fxyzDZ[ixyz];
        test_fxyzCURLx[ixyz] = test_fxyzDY[ixyz] - test_fxyzDZ[ixyz]; 
        test_fxyzCURLy[ixyz] = test_fxyzDZ[ixyz] - test_fxyzDX[ixyz]; 
        test_fxyzCURLz[ixyz] = test_fxyzDX[ixyz] - test_fxyzDY[ixyz];
        
        test_cxyzDX[ixyz] = dcdx(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2)*cx(dz*iz, az, dz*nz/2); 
        test_cxyzDY[ixyz] = cx(dx*ix, ax, dx*nx/2)*dcdx(dy*iy, ay, dy*ny/2)*cx(dz*iz, az, dz*nz/2); 
        test_cxyzDZ[ixyz] = cx(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2)*dcdx(dz*iz, az, dz*nz/2);  
        test_cxyzDX2[ixyz] = d2cdx2(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2)*cx(dz*iz, az, dz*nz/2); 
        test_cxyzDY2[ixyz] = cx(dx*ix, ax, dx*nx/2)*d2cdx2(dy*iy, ay, dy*ny/2)*cx(dz*iz, az, dz*nz/2); 
        test_cxyzDZ2[ixyz] = cx(dx*ix, ax, dx*nx/2)*cx(dy*iy, ay, dy*ny/2)*d2cdx2(dz*iz, az, dz*nz/2); 
        test_cxyzGRAD2[ixyz] = test_cxyzDX[ixyz] * conj(test_cxyzDX[ixyz]) + test_cxyzDY[ixyz] * conj(test_cxyzDY[ixyz]) + test_cxyzDZ[ixyz] * conj(test_cxyzDZ[ixyz]);
        test_cxyzLAP[ixyz] = test_cxyzDX2[ixyz] + test_cxyzDY2[ixyz] + test_cxyzDZ2[ixyz];
        test_cxyzDIV[ixyz] = test_cxyzDX[ixyz] + test_cxyzDY[ixyz] + test_cxyzDZ[ixyz];
        test_cxyzCURLx[ixyz] = test_cxyzDY[ixyz] - test_cxyzDZ[ixyz]; 
        test_cxyzCURLy[ixyz] = test_cxyzDZ[ixyz] - test_cxyzDX[ixyz]; 
        test_cxyzCURLz[ixyz] = test_cxyzDX[ixyz] - test_cxyzDY[ixyz]; 

        ixyz++; 
    }
    
    printf("\n3D REAL\n");
    printf("df/dx,df/dy,df/dz\n");
    wderiv_dfdx_3d_r(fxyz, fxyzDX);
    test_array_diff_r(nx*ny*nz, fxyzDX, test_fxyzDX, 1.0e-12);
    wderiv_dfdx(3, 'r', fxyz, fxyzDXgen);
    test_array_diff_r(nx*ny*nz, fxyzDXgen, test_fxyzDX, 1.0e-12);
    wderiv_dfdy_3d_r(fxyz, fxyzDY);
    test_array_diff_r(nx*ny*nz, fxyzDY, test_fxyzDY, 1.0e-12);
    wderiv_dfdy(3, 'r', fxyz, fxyzDYgen);
    test_array_diff_r(nx*ny*nz, fxyzDYgen, test_fxyzDY, 1.0e-12);
    wderiv_dfdz_3d_r(fxyz, fxyzDZ);
    test_array_diff_r(nx*ny*nz, fxyzDZ, test_fxyzDZ, 1.0e-12);
    wderiv_dfdz(3, 'r', fxyz, fxyzDZgen);
    test_array_diff_r(nx*ny*nz, fxyzDZgen, test_fxyzDZ, 1.0e-12);
    
    printf("\nd2f/dx2,d2f/dy2,d2f/dz2\n");
    wderiv_d2fdx2_3d_r(fxyz, fxyzDX2);
    test_array_diff_r(nx*ny*nz, fxyzDX2, test_fxyzDX2, 1.0e-12);
    wderiv_d2fdx2(3, 'r', fxyz, fxyzDX2gen);
    test_array_diff_r(nx*ny*nz, fxyzDX2gen, test_fxyzDX2, 1.0e-12);
    wderiv_d2fdy2_3d_r(fxyz, fxyzDY2);
    wderiv_d2fdy2(3, 'r', fxyz, fxyzDY2gen);
    test_array_diff_r(nx*ny*nz, fxyzDY2gen, test_fxyzDY2, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzDY2, test_fxyzDY2, 1.0e-12);
    wderiv_d2fdz2_3d_r(fxyz, fxyzDZ2);
    test_array_diff_r(nx*ny*nz, fxyzDZ2, test_fxyzDZ2, 1.0e-12);
    wderiv_d2fdz2(3, 'r', fxyz, fxyzDZ2gen);
    test_array_diff_r(nx*ny*nz, fxyzDZ2gen, test_fxyzDZ2, 1.0e-12);
    
    printf("\ndnf/dxn,dnf/dyn,dnf/dzn\n");
    wderiv_dnfdxn_3d_r(n_deriv, fxyz, fxyzDXn);
    test_array_diff_r(nx*ny*nz, fxyzDXn, test_fxyzDXn, 1.0e-12);
    wderiv_dnfdxn(3, 'r', 7, fxyz, fxyzDXngen);
    test_array_diff_r(nx*ny*nz, fxyzDXngen, test_fxyzDXn, 1.0e-12);
    wderiv_dnfdyn_3d_r(n_deriv, fxyz, fxyzDYn);
    test_array_diff_r(nx*ny*nz, fxyzDYn, test_fxyzDYn, 1.0e-12);
    wderiv_dnfdyn(3, 'r', 7, fxyz, fxyzDYngen);
    test_array_diff_r(nx*ny*nz, fxyzDYngen, test_fxyzDYn, 1.0e-12);
    wderiv_dnfdzn_3d_r(n_deriv, fxyz, fxyzDZn);
    test_array_diff_r(nx*ny*nz, fxyzDZn, test_fxyzDZn, 1.0e-12);
    wderiv_dnfdzn(3, 'r', 7, fxyz, fxyzDZngen);
    test_array_diff_r(nx*ny*nz, fxyzDZngen, test_fxyzDZn, 1.0e-12);
    
    printf("\nGradient x, y, z\n");
    wderiv_gradient_3d_r(fxyz, fxyzGRADx, fxyzGRADy, fxyzGRADz);
    wderiv_gradient(3, 'r', fxyz, fxyzGRADxgen, fxyzGRADygen, fxyzGRADzgen);
    test_array_diff_r(nx*ny*nz, fxyzGRADx, test_fxyzDX, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzGRADxgen, test_fxyzDX, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzGRADy, test_fxyzDY, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzGRADygen, test_fxyzDY, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzGRADz, test_fxyzDZ, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzGRADzgen, test_fxyzDZ, 1.0e-12);

    printf("\nGradient squared\n");
    wderiv_gradient2_3d_r(fxyz, fxyzGRAD2);
    test_array_diff_r(nx*ny*nz, fxyzGRAD2, test_fxyzGRAD2, 1.0e-12);
    wderiv_gradient2(3, 'r', fxyz, fxyzGRAD2gen);
    test_array_diff_r(nx*ny*nz, fxyzGRAD2gen, test_fxyzGRAD2, 1.0e-12);

    printf("\nLaplacian\n");
    wderiv_laplace_3d_r(fxyz, fxyzLAP);
    test_array_diff_r(nx*ny*nz, fxyzLAP, test_fxyzLAP, 1.0e-12);
    wderiv_laplace(3, 'r', fxyz, fxyzLAPgen);
    test_array_diff_r(nx*ny*nz, fxyzLAPgen, test_fxyzLAP, 1.0e-12);
    
    printf("\nDivergence\n");
    wderiv_divergence_3d_r(fxyz, fxyzDIV);
    test_array_diff_r(nx*ny*nz, fxyzDIV, test_fxyzDIV, 1.0e-12);
    wderiv_divergence(3, 'r', fxyz, fxyzDIVgen);
    test_array_diff_r(nx*ny*nz, fxyzDIVgen, test_fxyzDIV, 1.0e-12);

    printf("\nCURLx\n");
    wderiv_curl_3d_r(fxyz,fxyz,fxyz, fxyzCURLx, fxyzCURLy, fxyzCURLz);
    test_array_diff_r(nx*ny*nz, fxyzCURLx, test_fxyzCURLx, 1.0e-12);
    wderiv_curl(3, 'r', fxyz,fxyz,fxyz, fxyzCURLxgen, fxyzCURLygen, fxyzCURLzgen);
    test_array_diff_r(nx*ny*nz, fxyzCURLxgen, test_fxyzCURLx, 1.0e-12);
    printf("\nCURLy\n");
    test_array_diff_r(nx*ny*nz, fxyzCURLy, test_fxyzCURLy, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzCURLygen, test_fxyzCURLy, 1.0e-12);
    printf("\nCURLz\n");
    test_array_diff_r(nx*ny*nz, fxyzCURLz, test_fxyzCURLz, 1.0e-12);
    test_array_diff_r(nx*ny*nz, fxyzCURLzgen, test_fxyzCURLz, 1.0e-12);

    printf("\n/----------3D COMPLEX---------/\n");
    printf("dc/dx,dc/dy,dc/dz\n");
    wderiv_dfdx_3d_c(cxyz, cxyzDX);
    test_array_diff_c(nx*ny*nz, cxyzDX, test_cxyzDX, 1.0e-12);
    wderiv_dfdx(3, 'c', cxyz, cxyzDXgen);
    test_array_diff_c(nx*ny*nz, cxyzDXgen, test_cxyzDX, 1.0e-12);
    wderiv_dfdy_3d_c(cxyz, cxyzDY);
    test_array_diff_c(nx*ny*nz, cxyzDY, test_cxyzDY, 1.0e-12);
    wderiv_dfdy(3, 'c', cxyz, cxyzDYgen);
    test_array_diff_c(nx*ny*nz, cxyzDYgen, test_cxyzDY, 1.0e-12);
    wderiv_dfdz_3d_c(cxyz, cxyzDZ);
    test_array_diff_c(nx*ny*nz, cxyzDZ, test_cxyzDZ, 1.0e-12);
    wderiv_dfdz(3, 'c', cxyz, cxyzDZgen);
    test_array_diff_c(nx*ny*nz, cxyzDZgen, test_cxyzDZ, 1.0e-12);

    printf("\nd2c/dx2,d2c/dy2,d2c/dz2\n");
    wderiv_d2fdx2_3d_c(cxyz, cxyzDX2);
    test_array_diff_c(nx*ny*nz, cxyzDX2, test_cxyzDX2, 1.0e-12);
    wderiv_d2fdx2(3, 'c', cxyz, cxyzDX2gen);
    test_array_diff_c(nx*ny*nz, cxyzDX2gen, test_cxyzDX2, 1.0e-12);
    wderiv_d2fdy2_3d_c(cxyz, cxyzDY2);
    test_array_diff_c(nx*ny*nz, cxyzDY2, test_cxyzDY2, 1.0e-12);
    wderiv_d2fdy2(3, 'c', cxyz, cxyzDY2gen);
    test_array_diff_c(nx*ny*nz, cxyzDY2gen, test_cxyzDY2, 1.0e-12);
    wderiv_d2fdz2_3d_c(cxyz, cxyzDZ2);
    test_array_diff_c(nx*ny*nz, cxyzDZ2, test_cxyzDZ2, 1.0e-12);
    wderiv_d2fdz2(3, 'c', cxyz, cxyzDZ2gen);
    test_array_diff_c(nx*ny*nz, cxyzDZ2gen, test_cxyzDZ2, 1.0e-12);
    
    printf("\nGradient x, y, z\n");
    wderiv_gradient_3d_c(cxyz, cxyzGRADx, cxyzGRADy, cxyzGRADz);
    wderiv_gradient(3, 'c', cxyz, cxyzGRADxgen, cxyzGRADygen, cxyzGRADzgen);
    test_array_diff_c(nx*ny*nz, cxyzGRADx, test_cxyzDX, 1.0e-12);
    test_array_diff_c(nx*ny*nz, cxyzGRADxgen, test_cxyzDX, 1.0e-12);
    test_array_diff_c(nx*ny*nz, cxyzGRADy, test_cxyzDY, 1.0e-12);
    test_array_diff_c(nx*ny*nz, cxyzGRADygen, test_cxyzDY, 1.0e-12);
    test_array_diff_c(nx*ny*nz, cxyzGRADz, test_cxyzDZ, 1.0e-12);
    test_array_diff_c(nx*ny*nz, cxyzGRADzgen, test_cxyzDZ, 1.0e-12);

    printf("\nGradient squared\n");
    wderiv_gradient2_3d_c(cxyz, cxyzGRAD2);
    test_array_diff_c(nx*ny*nz, cxyzGRAD2, test_cxyzGRAD2, 1.0e-12);
    wderiv_gradient2(3, 'c', cxyz, cxyzGRAD2gen);
    test_array_diff_c(nx*ny*nz, cxyzGRAD2gen, test_cxyzGRAD2, 1.0e-12);

    printf("\nLaplacian\n");
    wderiv_laplace_3d_c(cxyz, cxyzLAP);
    test_array_diff_c(nx*ny*nz, cxyzLAP, test_cxyzLAP, 1.0e-12);
    wderiv_laplace(3, 'c', cxyz, cxyzLAPgen);
    test_array_diff_c(nx*ny*nz, cxyzLAPgen, test_cxyzLAP, 1.0e-12);

    printf("\nDivergence\n");
    wderiv_divergence_3d_c(cxyz, cxyzDIV);
    test_array_diff_c(nx*ny*nz, cxyzDIV, test_cxyzDIV, 1.0e-12);
    wderiv_divergence(3, 'c', cxyz, cxyzDIVgen);
    test_array_diff_c(nx*ny*nz, cxyzDIVgen, test_cxyzDIV, 1.0e-12);

    printf("\nCURLx\n");
    wderiv_curl_3d_c(cxyz,cxyz,cxyz, cxyzCURLx, cxyzCURLy, cxyzCURLz);
    test_array_diff_c(nx*ny*nz, cxyzCURLx, test_cxyzCURLx, 1.0e-12);
    wderiv_curl(3, 'c', cxyz,cxyz,cxyz, cxyzCURLxgen, cxyzCURLygen, cxyzCURLzgen);
    test_array_diff_c(nx*ny*nz, cxyzCURLxgen, test_cxyzCURLx, 1.0e-12);
    printf("\nCURLy\n");
    test_array_diff_c(nx*ny*nz, cxyzCURLy, test_cxyzCURLy, 1.0e-12);
    test_array_diff_c(nx*ny*nz, cxyzCURLygen, test_cxyzCURLy, 1.0e-12);
    printf("\nCURLz\n");
    test_array_diff_c(nx*ny*nz, cxyzCURLz, test_cxyzCURLz, 1.0e-12);
    test_array_diff_c(nx*ny*nz, cxyzCURLygen, test_cxyzCURLy, 1.0e-12);
   
    // destroy lib
    wderiv_clean();
    free(fxyz); free(fxyzDX); free(fxyzDY); free(fxyzDZ); free(fxyzDX2); free(fxyzDY2); free(fxyzDZ2); free(fxyzDXn); free(fxyzDYn); free(fxyzDZn); free(fxyzGRAD2); free(fxyzLAP); free(fxyzCURLx); free(fxyzCURLy); free(fxyzCURLz);
    free(cxyz); free(cxyzDX); free(cxyzDY); free(cxyzDZ); free(cxyzDX2); free(cxyzDY2); free(cxyzDZ2); free(cxyzGRAD2); free(cxyzLAP); free(cxyzCURLx); free(cxyzCURLy); free(cxyzCURLz);
    free(test_fxyzDX); free(test_fxyzDY); free(test_fxyzDZ); free(test_fxyzDX2); free(test_fxyzDY2); free(test_fxyzDZ2); free(test_fxyzDXn); free(test_fxyzDYn); free(test_fxyzDZn); free(test_fxyzGRAD2); free(test_fxyzLAP); free(test_fxyzCURLx); free(test_fxyzCURLy); free(test_fxyzCURLz);
    free(test_cxyzDX); free(test_cxyzDY); free(test_cxyzDZ); free(test_cxyzDX2); free(test_cxyzDY2); free(test_cxyzDZ2); free(test_cxyzGRAD2); free(test_cxyzLAP); free(test_cxyzCURLx); free(test_cxyzCURLy); free(test_cxyzCURLz);
    return 0; 
}

