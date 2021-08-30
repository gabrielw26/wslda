/**
 * W-SLDA Toolkit
 * 
 * Functions for expanding array from xD dim to yD dim.
 * 
 * @see example: ...
 * 
 * @author Gabriel Wlazlowski
 * */ 

#ifndef __WSLDA_RESIZE_H__ 
#define __WSLDA_RESIZE_H__ 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>

/**********************************************************************/
/***************************** DEFINITIONS ****************************/
/**********************************************************************/

int wslda_resize_array_2d_to_3d_r(int nx, int ny, double *in_2d, int nz, double *out_3d); 
int wslda_resize_array_2d_to_3d_v(int nx, int ny, double *in_2d, int nz, double *out_3d);
int wslda_resize_array_2d_to_3d_c(int nx, int ny, double complex *in_2d, int nz, double complex *out_3d);
int wslda_resize_array_2d_to_3d(char type, int nx, int ny, void *in_2d, int nz, void *out_3d);

int wslda_resize_array_1d_to_3d_r(int nx, double *in_1d, int ny, int nz, double *out_3d);
int wslda_resize_array_1d_to_3d_v(int nx, double *in_1d, int ny, int nz, double *out_3d);
int wslda_resize_array_1d_to_3d_c(int nx, double complex *in_1d, int ny, int nz, double complex *out_3d);
int wslda_resize_array_1d_to_3d(char type, int nx, void *in_1d, int ny, int nz, void *out_3d);

int wslda_resize_array_1d_to_2d_r(int nx, double *in_1d, int ny, double *out_2d);
int wslda_resize_array_1d_to_2d_v(int nx, double *in_1d, int ny, double *out_2d);
int wslda_resize_array_1d_to_2d_c(int nx, double complex *in_1d, int ny, double complex *out_2d); 
int wslda_resize_array_1d_to_2d(char type, int nx, void *in_1d, int ny, void *out_2d);

/**********************************************************************/
/**************************** IMPLEMENTATION **************************/
/**********************************************************************/

/**
 * Function expands 2D array of size (nx,ny) to 3D version of size (nx,ny,nz)
 * assuming uniformity along z direction.
 * Real version. 
 **/ 
int wslda_resize_array_2d_to_3d_r(int nx, int ny, double *in_2d, int nz, double *out_3d) 
{

    int ix, iy, iz, ixyz;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        out_3d[ixyz]=in_2d[iy+ny*ix];
        ixyz++;
    }    
    
    return 0;
}

/**
 * Function expands 2D array of size (nx,ny) to 3D version of size (nx,ny,nz)
 * assuming uniformity along z direction.
 * Vector version. 
 **/ 
int wslda_resize_array_2d_to_3d_v(int nx, int ny, double *in_2d, int nz, double *out_3d)
{
    wslda_resize_array_2d_to_3d_r(nx, ny, in_2d+0*nx*ny, nz, out_3d+0*nx*ny*nz);
    wslda_resize_array_2d_to_3d_r(nx, ny, in_2d+1*nx*ny, nz, out_3d+1*nx*ny*nz);
    wslda_resize_array_2d_to_3d_r(nx, ny, in_2d+2*nx*ny, nz, out_3d+2*nx*ny*nz);
    return 0;
}

/**
 * Function expands 2D array of size (nx,ny) to 3D version of size (nx,ny,nz)
 * assuming uniformity along z direction.
 * Complex version. 
 **/ 
int wslda_resize_array_2d_to_3d_c(int nx, int ny, double complex *in_2d, int nz, double complex *out_3d) 
{

    int ix, iy, iz, ixyz;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        out_3d[ixyz]=in_2d[iy+ny*ix];
        ixyz++;
    }    
    
    return 0;
}

/**
 * Function expands 2D array of size (nx,ny) to 3D version of size (nx,ny,nz)
 * assuming uniformity along z direction.
 * Generic version.
 * @param type c-complex, v-vector, otherwise real 
 **/
int wslda_resize_array_2d_to_3d(char type, int nx, int ny, void *in_2d, int nz, void *out_3d)
{
    if(type=='c')      return wslda_resize_array_2d_to_3d_c(nx, ny, (double complex *)in_2d, nz, (double complex *)out_3d);
    else if(type=='v') return wslda_resize_array_2d_to_3d_v(nx, ny, (double *)in_2d, nz, (double *)out_3d);
    else               return wslda_resize_array_2d_to_3d_r(nx, ny, (double *)in_2d, nz, (double *)out_3d);
}


/**
 * Function expands 1D array of size (nx) to 3D version of size (nx,ny,nz)
 * assuming uniformity along y and z direction.
 * Real version. 
 **/ 
int wslda_resize_array_1d_to_3d_r(int nx, double *in_1d, int ny, int nz, double *out_3d) 
{

    int ix, iy, iz, ixyz;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        out_3d[ixyz]=in_1d[ix];
        ixyz++;
    }    
    
    return 0;
}

/**
 * Function expands 1D array of size (nx) to 3D version of size (nx,ny,nz)
 * assuming uniformity along y and z direction.
 * Vector version. 
 **/ 
int wslda_resize_array_1d_to_3d_v(int nx, double *in_1d, int ny, int nz, double *out_3d) 
{
    wslda_resize_array_1d_to_3d_r(nx, in_1d+0*nx, ny, nz, out_3d+0*nx*ny*nz); 
    wslda_resize_array_1d_to_3d_r(nx, in_1d+1*nx, ny, nz, out_3d+1*nx*ny*nz); 
    wslda_resize_array_1d_to_3d_r(nx, in_1d+2*nx, ny, nz, out_3d+2*nx*ny*nz); 
    return 0;
}

/**
 * Function expands 1D array of size (nx) to 3D version of size (nx,ny,nz)
 * assuming uniformity along y and z direction.
 * Complex version. 
 **/ 
int wslda_resize_array_1d_to_3d_c(int nx, double complex *in_1d, int ny, int nz, double complex *out_3d) 
{

    int ix, iy, iz, ixyz;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        out_3d[ixyz]=in_1d[ix];
        ixyz++;
    }    
    
    return 0;
}

/**
 * Function expands 1D array of size (nx) to 3D version of size (nx,ny,nz)
 * assuming uniformity along y and z direction.
 * Generic version.
 * @param type c-complex, v-vector, otherwise real  
 **/ 
int wslda_resize_array_1d_to_3d(char type, int nx, void *in_1d, int ny, int nz, void *out_3d)
{
    if(type=='c')      return wslda_resize_array_1d_to_3d_c(nx, (double complex *)in_1d, ny, nz, (double complex *)out_3d);
    else if(type=='v') return wslda_resize_array_1d_to_3d_v(nx, (double *)in_1d, ny, nz, (double *)out_3d);
    else               return wslda_resize_array_1d_to_3d_r(nx, (double *)in_1d, ny, nz, (double *)out_3d);
}

/**
 * Function expands 1D array of size (nx) to 2D version of size (nx,ny)
 * assuming uniformity along y direction.
 * Real version. 
 **/ 
int wslda_resize_array_1d_to_2d_r(int nx, double *in_1d, int ny, double *out_2d) 
{

    int ix, iy, ixyz;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++)
    {
        out_2d[ixyz]=in_1d[ix];
        ixyz++;
    }    
    
    return 0;
}

/**
 * Function expands 1D array of size (nx) to 2D version of size (nx,ny)
 * assuming uniformity along y direction.
 * Vector version. 
 **/ 
int wslda_resize_array_1d_to_2d_v(int nx, double *in_1d, int ny, double *out_2d) 
{
    wslda_resize_array_1d_to_2d_r(nx, in_1d+0*nx, ny, out_2d+0*nx*ny); 
    wslda_resize_array_1d_to_2d_r(nx, in_1d+1*nx, ny, out_2d+1*nx*ny); 
    wslda_resize_array_1d_to_2d_r(nx, in_1d+2*nx, ny, out_2d+2*nx*ny); 
    return 0;
}

/**
 * Function expands 1D array of size (nx) to 2D version of size (nx,ny)
 * assuming uniformity along y direction.
 * Complex version. 
 **/ 
int wslda_resize_array_1d_to_2d_c(int nx, double complex *in_1d, int ny, double complex *out_2d) 
{

    int ix, iy, ixyz;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++)
    {
        out_2d[ixyz]=in_1d[ix];
        ixyz++;
    }    
    
    return 0;
}

/**
 * Function expands 1D array of size (nx) to 2D version of size (nx,ny)
 * assuming uniformity along y direction.
 * Generic version.
 * @param type c-complex, v-vector, otherwise real  
 **/ 
int wslda_resize_array_1d_to_2d(char type, int nx, void *in_1d, int ny, void *out_2d)
{
    if(type=='c')      return wslda_resize_array_1d_to_2d_c(nx, (double complex *)in_1d, ny, (double complex *)out_2d);
    else if(type=='v') return wslda_resize_array_1d_to_2d_v(nx, (double *)in_1d, ny, (double *)out_2d);
    else               return wslda_resize_array_1d_to_2d_r(nx, (double *)in_1d, ny, (double *)out_2d);
}

#endif
