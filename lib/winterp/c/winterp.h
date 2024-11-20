/**
 * W-interp Library 
 * @see https://gitlab.fizyka.pw.edu.pl/wtools/winterp
 * 
 * Set of functions for the spectral interpolation of data defined on a lattice
 * 
 * Warsaw University of Technology,
 * Nuclear Theory Group
 * 2021
 * 
 * @author Gabriel Wlazlowski
 * */  


#ifndef W_INTERPOLATION_H_ 
#define W_INTERPOLATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#ifdef __cplusplus
#include <complex>
#define Complex std::complex<double>
#else
#include <complex.h>
#define Complex double complex
#endif

#include <fftw3.h>
#include <math.h>
#include <stdbool.h>

// ERROR CODES
#define WINTERP_OK 0
#define WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED 10015
#define WINTERP_ERR_INTRISTIC_ERROR 10022
#define WINTERP_ERR_CANNOT_ALLOCATE_MEMORY 20002

#ifdef __cplusplus
extern "C" {
#endif

/**********************************************************************/
/***************************** DEFINITIONS ****************************/
/**********************************************************************/


/**
 * Function changing the lattice resolution, while keeping fixed volume
 * @see README.md#lattice-interpolation-to-a-new-resolution
 * */ 

int winterp_interpolation_1dr(int nxi, double *funIn, int nxo, double *funOut);
int winterp_interpolation_1dv(int nxi, double *funIn, int nxo, double *funOut);
int winterp_interpolation_1dc(int nxi, Complex *funIn, int nxo, Complex *funOut);
int winterp_interpolation_1d(char type, int nxi, void *funIn, int nxo, void *funOut);

int winterp_interpolation_2dr(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut);
int winterp_interpolation_2dv(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut);
int winterp_interpolation_2dc(int nxi, int nyi, Complex *funIn, int nxo, int nyo, Complex *funOut);
int winterp_interpolation_2d(char type, int nxi, int nyi, void *funIn, int nxo, int nyo, void *funOut);

int winterp_interpolation_3dr(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut);
int winterp_interpolation_3dv(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut);
int winterp_interpolation_3dc(int nxi, int nyi, int nzi, Complex *funIn, int nxo, int nyo, int nzo, Complex *funOut);  
int winterp_interpolation_3d(char type, int nxi, int nyi, int nzi, void *funIn, int nxo, int nyo, int nzo, void *funOut);
 


/**
 * Functions for extracting values for given point
 * @see README.md#interpolation-for-arbitrary-point-in-space
 * */

/**
 * Structure representing interpolator for arbitrary data type
 * */
typedef struct
{
    char datatype; // r-real, c-complex, v-vector
    int datadim; /// dimensonality of data
    int nx; /// lattice size
    int ny;
    int nz; 
    double dx; /// lattice spacing
    double dy;
    double dz; 
    Complex *datak; /// Fourier transform of data
} winterp_interpolator;

/** 
 * Function creates interpolator representing data
 * @param datatype r-real, c-complex, v-vector (INPUT)
 * @param datadim 1,2 or 3 (INPUT)
 * @param dims array of size datadim with dimensions (INPUT)
 * @param spacings array of size datadim with lattice spacings (INPUT)
 * @param data pointer to data to be interpolated (INPUT)
 * @param interp interpolator (OUTPUT)
 * @return error code
 * */
int winterp_create_interpolator(char datatype, int datadim, int *dims, double *spacings, void *data, winterp_interpolator *interp);
int winterp_create_interpolator_1d_c(int nx, double dx, Complex *data, winterp_interpolator *interp);
int winterp_create_interpolator_1d_r(int nx, double dx, double *data, winterp_interpolator *interp);
int winterp_create_interpolator_1d_v(int nx, double dx, double *data, winterp_interpolator *interp);
int winterp_create_interpolator_2d_c(int nx, int ny, double dx, double dy, Complex *data, winterp_interpolator *interp);
int winterp_create_interpolator_2d_r(int nx, int ny, double dx, double dy, double *data, winterp_interpolator *interp);
int winterp_create_interpolator_2d_v(int nx, int ny, double dx, double dy, double *data, winterp_interpolator *interp);
int winterp_create_interpolator_3d_c(int nx, int ny, int nz, double dx, double dy, double dz, Complex *data, winterp_interpolator *interp);
int winterp_create_interpolator_3d_r(int nx, int ny, int nz, double dx, double dy, double dz, double *data, winterp_interpolator *interp);
int winterp_create_interpolator_3d_v(int nx, int ny, int nz, double dx, double dy, double dz, double *data, winterp_interpolator *interp);

/**
 * Returns value of function for arbitrary point in 1D
 * @param interp interpolator of function (INPUT)
 * @param x point (INPUT)
 * @param val value of function (OUTPUT)
 * @return error code
 * */
int winterp_getvalue_1d_c(winterp_interpolator *interp, double x, Complex *val);
int winterp_getvalue_1d_r(winterp_interpolator *interp, double x, double *val);
int winterp_getvalue_1d_v(winterp_interpolator *interp, double x, double *valx, double *valy, double *valz);

/**
 * Returns value of function for arbitrary point in 2D
 * @param interp interpolator of function (INPUT)
 * @param x point (INPUT)
 * @param y point (INPUT)
 * @param val value of function (OUTPUT)
 * @return error code
 * */
int winterp_getvalue_2d_c(winterp_interpolator *interp, double x, double y, Complex *val);
int winterp_getvalue_2d_r(winterp_interpolator *interp, double x, double y, double *val);
int winterp_getvalue_2d_v(winterp_interpolator *interp, double x, double y, double *valx, double *valy, double *valz);

/**
 * Returns value of function for arbitrary point in 3D
 * @param interp interpolator of function (INPUT)
 * @param x point (INPUT)
 * @param y point (INPUT)
 * @param z point (INPUT)
 * @param val value of function (OUTPUT)
 * @return error code
 * */
int winterp_getvalue_3d_c(winterp_interpolator *interp, double x, double y, double z, Complex *val);
int winterp_getvalue_3d_r(winterp_interpolator *interp, double x, double y, double z, double *val);
int winterp_getvalue_3d_v(winterp_interpolator *interp, double x, double y, double z, double *valx, double *valy, double *valz);
/**
 * Function destroys interolator and frees memory
 * */
int winterp_destroy_interpolator(winterp_interpolator *interp);

#ifdef __cplusplus
} // closing brace for extern "C"
#endif
#endif
