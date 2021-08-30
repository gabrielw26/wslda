/**
 * W-SLDA Toolkit
 * Warsaw University of Technology, Faculty of Physics (2021)
 * https://wslda.fizyka.pw.edu.pl/
 * 
 * This library provides set of functions for computation derivatives using sepectral methods. 
 * The lib depends on FFTW library.
 * The library is compatible with W-SLDA Toolkit however, it can be used as standalone lib. 
 * It is written in C99 standard.
 * 
 * For more information see: https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/wderiv-library
 * */ 


#ifndef __W_DERIV_LIB__
#define __W_DERIV_LIB__

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifdef __cplusplus
#include <complex>
#define Complex std::complex<double>
#else
#include <complex.h>
#define Complex double complex
#endif

#include <fftw3.h>

/**
 * List of errors
 * */
#define WDERIV_OK 20001
#define WDERIV_ERR_CANNOT_ALLOCATE_MEMORY 20002
#define WDERIV_ERR_NOINIT 20003
#define WDERIV_ERR_INCORRECT_DIRECTION 20004
#define WDERIV_ERR_INCORRECT_DATADIM 20005
#define WDERIV_ERR_INCORRECT_DATATYPE 20006
#define WDERIV_ERR_UNDEFINED 20007

#define WDERIV_DX 0
#define WDERIV_DY 1
#define WDERIV_DZ 2

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Function initilizes wderiv for 3d calcs
 * @param nx lattice size, x-coordinate
 * @param ny lattice size, y-coordinate
 * @param nz lattice size, z-coordinate
 * @param dx lattice spacing, x-coordinate
 * @param dy lattice spacing, y-coordinate
 * @param dz lattice spacing, z-coordinate
 * @return error code
 * */
int wderiv_init_3d(int nx, int ny, int nz, double dx, double dy, double dz);
int wderiv_init_2d(int nx, int ny, double dx, double dy);
int wderiv_init_1d(int nx, double dx);

/**
 * Function cleans wdderiv lib structures
 * @return error code
 * */
int wderiv_clean_3d();
int wderiv_clean_2d();
int wderiv_clean_1d();
int wderiv_clean();

/**
 * Function computes n-th derivative with respect to direction
 * of 3D real function f(x,y,z)   
 * @param direction WDERIV_DX or WDERIV_DY or WDERIV_DZ
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param deriv pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *              NOTE: deriv can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_derivative_3d_r(int direction, int n, double *f, double *deriv);


/**
 * Function computes first derivative with respect to x (df/dx)
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dfdx pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dfdx can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdx_3d_r(double *f, double *dfdx);

/**
 * Function computes first derivative with respect to y (df/dy)
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dfdy pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dfdy can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdy_3d_r(double *f, double *dfdy);

/**
 * Function computes first derivative with respect to z (df/dz)
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dfdz pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dfdz can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdz_3d_r(double *f, double *dfdz);

/**
 * Function computes second derivative with respect to x (d^2f/dx^2)
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param d2fdx2 pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: d2fdx2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdx2_3d_r(double *f, double *d2fdx2);

/**
 * Function computes second derivative with respect to y (d^2f/dy^2)
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dnfdyn pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dnfdyn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdy2_3d_r(double *f, double *d2fdy2);

/**
 * Function computes second derivative with respect to z (d^2f/dz^2)
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dnfdzn pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dnfdzn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdz2_3d_r(double *f, double *d2fdz2);

/**
 * Function computes n-th derivative with respect to x (d^nf/dx^n)
 * of 3D real function f(x,y,z)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dnfdxn pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dnfdxn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdxn_3d_r(int n, double *f, double *dnfdxn);

/**
 * Function computes n-th derivative with respect to y (d^nf/dy^n)
 * of 3D real function f(x,y,z)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dnfdyn pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dnfdyn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdyn_3d_r(int n, double *f, double *dnfdyn);

/**
 * Function computes n-th derivative with respect to z (d^nf/dz^n)
 * of 3D real function f(x,y,z)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dnfdzn pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dnfdzn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdzn_3d_r(int n, double *f, double *dnfdzn);

/**
 * Function computes gradient
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dfdx pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @param dfdy pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @param dfdz pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient_3d_r(double *f, double *dfdx, double *dfdy, double *dfdz);

/**
 * Function computes gradient_square = |grad f|^2 = |f/dx|^2 + |df/dy|^2 + |df/dz|^2
 * of 3D real function f(x,y,z)
 * @param f pointer to vector function, array of size [nx*ny*nz] (INPUT)
 * @param gradient_square pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient2_3d_r(double *f, double *gradient_square);

/**
 * Function computes laplace = d^2f/dx^2 + d^2f/dy^2 + d^2f/dz^2
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param laplace pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_laplace_3d_r(double *f, double *laplace);

/**
 * Function computes divergence = df/dx + df/dy + df/dz
 * of 3D real function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param divergence pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_divergence_3d_r(double *f, double *divergence);

/**
 * Function computes curl = [dfy/dz-dfz/dy, dfz/dx-dfx/dz, dfx/dy-dfy/dx]
 * of 3D real vector function [fx(x,y,z), fy(x,y,z), fz(x,y,z)]
 * @param fx pointer to x component of vector function, array of size [nx*ny*nz] (INPUT)
 * @param fy pointer to y component of vector function, array of size [nx*ny*nz] (INPUT)
 * @param fz pointer to z component of vector function, array of size [nx*ny*nz] (INPUT)
 * @param curl_x pointer to x component of computed curl, array of size [nx*ny*nz] (OUTPUT)
 * @param curl_y pointer to y component of computed curl, array of size [nx*ny*nz] (OUTPUT)
 * @param curl_z pointer to z component of computed curl, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_curl_3d_r(double *fx, double *fy, double *fz, double *curl_x, double *curl_y, double *curl_z);

/**
 * Function computes n-th derivative with respect to direction
 * of 2D real function f(x,y)   
 * @param direction WDERIV_DX or WDERIV_DY
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param deriv pointer to function, array of size [nx*ny] (OUTPUT)
 *              NOTE: deriv can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_derivative_2d_r(int direction, int n, double *f, double *deriv);

/**
 * Function computes first derivative with respect to x (df/dx)
 * of 2D real function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dfdx pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dfdx can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdx_2d_r(double *f, double *dfdx);

/**
 * Function computes first derivative with respect to y (df/dy)
 * of 2D real function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dfdy pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dfdy can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdy_2d_r(double *f, double *dfdy);

/**
 * Function computes second derivative with respect to x (d^2f/dx^2)
 * of 2D real function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param d2fdx2 pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: d2fdx2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdx2_2d_r(double *f, double *d2fdx2);

/**
 * Function computes second derivative with respect to y (d^2f/dy^2)
 * of 2D real function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param d2fdy2 pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: d2fdy2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdy2_2d_r(double *f, double *d2fdy2);

/**
 * Function computes n-th derivative with respect to x (d^nf/dx^n)
 * of 2D real function f(x,y)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dnfdxn pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dnfdxn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdxn_2d_r(int n, double *f, double *dnfdxn);

/**
 * Function computes n-th derivative with respect to y (d^nf/dy^n)
 * of 2D real function f(x,y)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dnfdyn pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dnfdyn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdyn_2d_r(int n, double *f, double *dnfdyn);

/**
 * Function computes gradient
 * of 2D real function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dfdx pointer to function, array of size [nx*ny] (OUTPUT)
 * @param dfdy pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient_2d_r(double *f, double *dfdx, double *dfdy);

/**
 * Function computes gradient_square = |grad f|^2 = |f/dx|^2 + |df/dy|^2
 * of 2D real function f(x,y)
 * @param f pointer to vector function, array of size [nx*ny] (INPUT)
 * @param gradient_square pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient2_2d_r(double *f, double * gradient_square);

/**
 * Function computes laplace = d^2f/dx^2 + d^2f/dy^2
 * of 2D real function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param laplace pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */
int wderiv_laplace_2d_r(double *f, double *laplace);

/**
 * Function computes divergence = df/dx + df/dy + df/dz
 * of 2D real function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param divergence pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */
int wderiv_divergence_2d_r(double *f, double *divergence);

/**
 * Function computes n-th derivative
 * @param n order of derivative
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param deriv pointer to function, array of size [nx] (OUTPUT)
 *              NOTE: deriv can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_derivative_1d_r(int n, double *f, double *deriv);

/**
 * Function computes first derivative with respect to x (df/dx)
 * of 1D real function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param dfdx pointer to function, array of size [nx] (OUTPUT)
 *                NOTE: dfdx can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdx_1d_r(double *f, double *dfdx);

/**
 * Function computes second derivative with respect to x (d^2f/dx^2)
 * of 1D real function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param d2fdx2 pointer to function, array of size [nx] (OUTPUT)
 *                NOTE: d2fdx2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdx2_1d_r(double *f, double *d2fdx2);

/**
 * Function computes n-th derivative with respect to x (d^nf/dx^n)
 * of 1D real function f(x)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param dnfdxn pointer to function, array of size [nx] (OUTPUT)
 *                NOTE: dnfdxn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdxn_1d_r(int n, double *f, double *dnfdxn);

/**
 * Function computes gradient = dfdx
 * of 1D real function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param dfdx pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient_1d_r(double *f, double *dfdx);

/**
 * Function computes gradient_square = |grad f|^2 = |f/dx|^2
 * of 1D real function f(x)
 * @param f pointer to vector function, array of size [nx] (INPUT)
 * @param gradient_square pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient2_1d_r(double *f, double * gradient_square);

/**
 * Function computes laplace = d^2f/dx^2
 * of 1D real function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param laplace pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_laplace_1d_r(double *f, double *laplace);

/**
 * Function computes divergence = df/dx
 * of 1D complex function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param divergence pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_divergence_1d_r(double *f, double *divergence);

//COMPLEX

/**
 * Function computes n-th derivative with respect to direction
 * of 3D complex function f(x,y,z)   
 * @param direction WDERIV_DX or WDERIV_DY or WDERIV_DZ
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param deriv pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *              NOTE: deriv can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_derivative_3d_c(int direction, int n, Complex *f, Complex *deriv);

/**
 * Function computes first derivative with respect to x (df/dx)
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nz*ny*nz] (INPUT)
 * @param dfdx pointer to function, array of size [nz*ny*nz] (OUTPUT)
 *                NOTE: dfdx can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdx_3d_c(Complex *f, Complex *dfdx);

/**
 * Function computes first derivative with respect to y (df/dy)
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dfdy pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dfdy can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdy_3d_c(Complex *f, Complex *dfdy);

/**
 * Function computes first derivative with respect to z (df/dz)
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dfdz pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dfdz can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdz_3d_c(Complex *f, Complex *dfdz);

/**
 * Function computes second derivative with respect to x (d^2f/dx^2)
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nz*ny*nz] (INPUT)
 * @param d2fdx2 pointer to function, array of size [nz*ny*nz] (OUTPUT)
 *                NOTE: d2fdx2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdx2_3d_c(Complex *f, Complex *d2fdx2);

/**
 * Function computes second derivative with respect to y (d^2f/dy^2)
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param d2fdy2 pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: d2fdy2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdy2_3d_c(Complex *f, Complex *d2fdy2);

/**
 * Function computes second derivative with respect to z (d^2f/dz^2)
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param d2fdz2 pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: d2fdz2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdz2_3d_c(Complex *f, Complex *d2fdz2);


/**
 * Function computes n-th derivative with respect to x (d^nf/dx^n)
 * of 3D complex function f(x,y,z)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nz*ny*nz] (INPUT)
 * @param dnfdxn pointer to function, array of size [nz*ny*nz] (OUTPUT)
 *                NOTE: dnfdxn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdxn_3d_c(int n, Complex *f, Complex *dnfdxn);

/**
 * Function computes n-th derivative with respect to y (d^nf/dx^y)
 * of 3D complex function f(x,y,z)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dnfdyn pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dnfdyn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdyn_3d_c(int n, Complex *f, Complex *dnfdyn);

/**
 * Function computes n-th derivative with respect to z (d^nf/dx^z)
 * of 3D complex function f(x,y,z)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dnfdzn pointer to function, array of size [nx*ny*nz] (OUTPUT)
 *                NOTE: dnfdzn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdzn_3d_c(int n, Complex *f, Complex *dnfdzn);

/**
 * Function computes gradient
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param dfdx pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @param dfdy pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @param dfdz pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient_3d_c(Complex *f, Complex *dfdx, Complex *dfdy, Complex *dfdz);

/**
 * Function computes gradient_square = |grad f|^2 = |df/dx|^2 + |df/dy|^2 + |df/dz|^2 (where |c|^2 = c x c*)
 * of 3D complex function f(x,y,z)
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param gradient_square pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient2_3d_c(Complex *f, Complex *gradient_square);

/**
 * Function computes laplace = d^2f/dx^2 + d^2f/dy^2 + d^2f/dz^2
 * of 3D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param laplace pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */ 
int wderiv_laplace_3d_c(Complex *f, Complex *laplace);

/**
 * Function computes divergence = df/dx + df/dy + df/dz
 * of 3D complex function f(x,y,z)  
 * @param f pointer to function, array of size [nx*ny*nz] (INPUT)
 * @param divergence pointer to function, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_divergence_3d_c(Complex *f, Complex *divergence);


/**
 * Function computes curl = [dfy/dz-dfz/dy, dfz/dx-dfx/dz, dfx/dy-dfy/dx]
 * of 3D complex vector function [fx(x,y,z), fy(x,y,z), fz(x,y,z)]
 * @param fx pointer to x component of vector function, array of size [nx*ny*nz] (INPUT)
 * @param fy pointer to y component of vector function, array of size [nx*ny*nz] (INPUT)
 * @param fz pointer to z component of vector function, array of size [nx*ny*nz] (INPUT)
 * @param curl_x pointer to x component of computed curl, array of size [nx*ny*nz] (OUTPUT)
 * @param curl_y pointer to y component of computed curl, array of size [nx*ny*nz] (OUTPUT)
 * @param curl_z pointer to z component of computed curl, array of size [nx*ny*nz] (OUTPUT)
 * @return error code
 * */
int wderiv_curl_3d_c(Complex *fx, Complex *fy, Complex *fz, Complex *curl_x, Complex *curl_y, Complex *curl_z);

/**
 * Function computes n-th derivative with respect to direction
 * of 2D complex function f(x,y)   
 * @param direction WDERIV_DX or WDERIV_DY
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param deriv pointer to function, array of size [nx*ny] (OUTPUT)
 *              NOTE: deriv can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_derivative_2d_c(int direction, int n, Complex *f, Complex *deriv);

/**
 * Function computes first derivative with respect to x (df/dx)
 * of 2D complex function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dfdx pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dfdx can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdx_2d_c(Complex *f, Complex *dfdx);

/**
 * Function computes first derivative with respect to y (df/dy)
 * of 2D complex function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dfdy pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dfdy can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdy_2d_c(Complex *f, Complex *dfdy);

/**
 * Function computes second derivative with respect to x (d^2f/dx^2)
 * of 2D complex function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param d2fdx2 pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: d2fdx2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdx2_2d_c(Complex *f, Complex *d2fdx2);

/**
 * Function computes second derivative with respect to y (d^2f/dy^2)
 * of 2D complex function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param d2fdy2 pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: d2fdy2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdy2_2d_c(Complex *f, Complex *d2fdy2);

/**
 * Function computes n-th derivative with respect to x (d^nf/dx^n)
 * of 2D complex function f(x,y)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dnfdxn pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dnfdxn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdxn_2d_c(int n, Complex *f, Complex *dnfdxn);

/**
 * Function computes n-th derivative with respect to y (d^nf/dy^n)
 * of 2D complex function f(x,y)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dnfdyn pointer to function, array of size [nx*ny] (OUTPUT)
 *                NOTE: dnfdyn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdyn_2d_c(int n, Complex *f, Complex *dnfdyn);

/**
 * Function computes gradient
 * of 2D complex function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param dfdx pointer to function, array of size [nx*ny] (OUTPUT)
 * @param dfdy pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient_2d_c(Complex *f, Complex *dfdx, Complex *dfdy);

/**
 * Function computes gradient_square = |grad f|^2 = |df/dx|^2 + |df/dy|^2  (where |c|^2 = c x c*)
 * of 2D complex vector function f(x,y)
 * @param f pointer to vector function, array of size [nx*ny] (INPUT)
 * @param gradient_square pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient2_2d_c(Complex *f, Complex *gradient_square);

/**
 * Function computes laplace = d^2f/dx^2 + d^2f/dy^2
 * of 2D complex function f(x,y)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param laplace pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */ 
int wderiv_laplace_2d_c(Complex *f, Complex *laplace);

/**
 * Function computes divergence = df/dx + df/dy + df/dz
 * of 2D complex function f(x,y,z)   
 * @param f pointer to function, array of size [nx*ny] (INPUT)
 * @param divergence pointer to function, array of size [nx*ny] (OUTPUT)
 * @return error code
 * */
int wderiv_divergence_2d_c(Complex *f, Complex *divergence);

/**
 * Function computes n-th derivative of 1D complex function f(x)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param deriv pointer to function, array of size [nx] (OUTPUT)
 *              NOTE: deriv can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_derivative_1d_c(int n, Complex *f, Complex *deriv);

/**
 * Function computes first derivative with respect to x (df/dx)
 * of 1D complex function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param dfdx pointer to function, array of size [nx] (OUTPUT)
 *                NOTE: dfdx can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdx_1d_c(Complex *f, Complex *dfdx);

/**
 * Function computes second derivative with respect to x (d^2f/dx^2)
 * of 1D complex function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param d2fdx2 pointer to function, array of size [nx] (OUTPUT)
 *                NOTE: d2fdx2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdx2_1d_c(Complex *f, Complex *d2fdx2);

/**
 * Function computes n-th derivative with respect to x (d^nf/dx^n)
 * of 1D complex function f(x)   
 * @param n order of derivative
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param dnfdxn pointer to function, array of size [nx] (OUTPUT)
 *                NOTE: dnfdxn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdxn_1d_c(int n, Complex *f, Complex *dnfdxn);
/**
 * Function computes gradient
 * of 1D real function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param dfdx pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient_1d_c(Complex *f, Complex *dfdx);

/**
 * Function computes gradient_square = |grad f|^2 = |df/dx|^2  (where |c|^2 = c x c*)
 * of 1D complex function f(x)
 * @param f pointer to vector function, array of size [nx] (INPUT)
 * @param gradient_square pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_gradient2_1d_c(Complex *f, Complex *gradient_square);

/**
 * Function computes laplace = d^2f/dx^2
 * of 1D complex function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param laplace pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_laplace_1d_c(Complex *f, Complex *laplace);

/**
 * Function computes divergence = df/dx
 * of 1D complex function f(x)   
 * @param f pointer to function, array of size [nx] (INPUT)
 * @param divergence pointer to function, array of size [nx] (OUTPUT)
 * @return error code
 * */
int wderiv_divergence_1d_c(Complex *f, Complex *divergence);

//GENERIC FUNCTIONS

/**
 * GENERIC function that computes first derivative with respect to x (df/dx)
 * @param datadim dimensonality of function: 1->f(x), 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param dfdx pointer to function (OUTPUT)
 *                NOTE: dfdx can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdx(int datadim, char type, void *f, void *dfdx);

/**
 * GENERIC function that computes first derivative with respect to y (df/dy)
 * @param datadim dimensonality of function: 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param dfdy pointer to function (OUTPUT)
 *                NOTE: dfdy can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdy(int datadim, char type, void *f, void *dfdy);

/**
 * GENERIC function that computes first derivative with respect to z (df/dz)
 * @param datadim dimensonality of function: 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param dfdz pointer to function (OUTPUT)
 *                NOTE: dfdz can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dfdz(int datadim, char type, void *f, void *dfdz);

/**
 * GENERIC function that computes second derivative with respect to x (df/dx)
 * @param datadim dimensonality of function: 1->f(x), 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param d2fdx2 pointer to function (OUTPUT)
 *                NOTE: d2fdx2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdx2(int datadim, char type, void *f, void *d2fdx2);

/**
 * GENERIC function that computes second derivative with respect to y (df/dy)
 * @param datadim dimensonality of function: 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param df2dy2 pointer to function (OUTPUT)
 *                NOTE: d2fdy2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdy2(int datadim, char type, void *f, void *d2fdy2);

/**
 * GENERIC function that computes second derivative with respect to z (df/dz)
 * @param datadim dimensonality of function: 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param d2fdz2 pointer to function (OUTPUT)
 *                NOTE: d2fdz2 can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_d2fdz2(int datadim, char type, void *f, void *d2fdz2);

/**
 * GENERIC function that computes n-th derivative with respect to x (dnf/dxn)
 * @param n order of derivative
 * @param datadim dimensonality of function: 1->f(x), 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param dnfdxn pointer to function (OUTPUT)
 *                NOTE: dnfdxn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdxn(int datadim, char type, int n,  void *f, void *dnfdxn);

/**
 * GENERIC function that computes n-th derivative with respect to y (dnf/dyn)
 * @param n order of derivative
 * @param datadim dimensonality of function: 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param dnfdyn pointer to function (OUTPUT)
 *                NOTE: dnfdyn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdyn(int datadim, char type, int n,  void *f, void *dnfdyn);
/**
 * GENERIC function that computes n-th derivative with respect to z (dnf/dzn)
 * @param n order of derivative
 * @param datadim dimensonality of function: 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param dnfdzn pointer to function (OUTPUT)
 *                NOTE: dnfdzn can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_dnfdzn(int datadim, char type, int n,  void *f, void *dnfdzn);

/**
 * GENERIC function that computes gradient vector = [df/dx, df/dy, df/dz]
 * if datadim==3 function returns [df/dx, df/dy, df/dz]
 * if datadim==2 function returns [df/dx, df/dy], pointer dfdz is not used (can be NULL)
 * if datadim==1 function returns [df/dx], pointers dfdy and dfdz are not used (can be NULL)
 * @param datadim dimensonality of function: 1->f(x), 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function                                          (INPUT)
 * @param dfdx pointer to function                                      (OUTPUT)
 * @param dfdy pointer to function, ignored if datadim==1               (OUTPUT)
 * @param dfdz pointer to function, ignored if datadim==1 or datadim==2 (OUTPUT)
 * @return error code
 * */
int wderiv_gradient(int datadim, char type, void *f, void *dfdx, void *dfdy, void *dfdz);

/**
 * GENERIC function that computes laplacian = d2f/dx2 + d2f/dy2 + d2f/dz2
 * @param datadim dimensonality of function: 1->f(x), 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param laplace pointer to function (OUTPUT)
 *                NOTE: laplace can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_laplace(int datadim, char type, void *f, void *laplace);

/**
 * GENERIC function that computes gradient squared = |df/dx|^2 + |df/dy|^2 + |df/dz|^2
 * @param datadim dimensonality of function: 1->f(x), 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param gradient_squared pointer to function (OUTPUT)
 *                NOTE: gradient_squared can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_gradient2(int datadim, char type, void *f, void *gradient_squared);

/**
 * GENERIC function that computes divergence = df/dx + df/dy + df/dz
 * @param datadim dimensonality of function: 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param f pointer to function (INPUT)
 * @param divergence pointer to function (OUTPUT)
 *                NOTE: divergence can be the same pointer as f, then result overwrites input 
 * @return error code
 * */
int wderiv_divergence(int datadim, char type, void *f, void *divergence);

/**
 * GENERIC function that computes curl = [dfy/dz-dfz/dy, dfz/dx-dfx/dz, dfx/dy-dfy/dx]
 * @param datadim dimensonality of function: 2->f(x,y), 3->f(x,y,z)
 * @param type r-real function, c-complex function
 * @param fx pointer to function (INPUT)
 * @param fy pointer to function (INPUT)
 * @param fz pointer to function (INPUT)
 * @param curl_x pointer to function (OUTPUT)
 * @param curl_y pointer to function (OUTPUT)
 * @param curl_z pointer to function (OUTPUT)
 * @return error code
 * */
int wderiv_curl(int datadim, char type, void *fx, void *fy, void *fz, void *curl_x, void *curl_y, void *curl_z);

#ifdef __cplusplus
} // closing brace for extern "C"
#endif
#endif
