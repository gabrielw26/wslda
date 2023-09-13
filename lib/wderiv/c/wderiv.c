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

#include "wderiv.h"

#define WDERIV_USE_FFTW_PLANNER FFTW_ESTIMATE

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "ERROR: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return WDERIV_ERR_CANNOT_ALLOCATE_MEMORY;                               \
    } 
    
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
    /**
 * Global structure containing basic info about lattice
 * */
typedef struct
{
    int nx;
    int ny;
    int nz;
    double dx;
    double dy;
    double dz;
    int dim;
} wderiv_info;

// global variable, visible only within this file
static wderiv_info wd = {0,0,0,0.0,0.0,0.0,0};


// --- Initialization routines ---
int wderiv_init_3d(int nx, int ny, int nz, double dx, double dy, double dz)
{
    wd.dim=3;
    wd.nx=nx;
    wd.ny=ny;
    wd.nz=nz;
    wd.dx=dx;
    wd.dy=dy;
    wd.dz=dz;    
    
    return WDERIV_OK;
}
int wderiv_init_2d(int nx, int ny, double dx, double dy)
{
    wd.dim=2;
    wd.nx=nx;
    wd.ny=ny;
    wd.dx=dx;
    wd.dy=dy;
    
    return WDERIV_OK;
}
int wderiv_init_1d(int nx, double dx)
{
    wd.dim=1;
    wd.nx=nx;
    wd.dx=dx;
    
    return WDERIV_OK;
}

// --- Cleaning routines ---
int wderiv_clean()
{
    wd.dim=0;
    wd.nx= 0;
    wd.ny= 0;
    wd.nz= 0;
    wd.dx= 0;
    wd.dy= 0;
    wd.dz= 0;
    return WDERIV_OK;
}

int wderiv_clean_3d() {return wderiv_clean();}
int wderiv_clean_2d() {return wderiv_clean();}
int wderiv_clean_1d() {return wderiv_clean();}

// --- 3D main low level functions REAL ---
int wderiv_derivative_3d_r(int direction, int n, double *f, double *deriv)
{
    if(wd.dim!=3) return WDERIV_ERR_NOINIT;
    
    register int nx=wd.nx, ny=wd.ny, nz=wd.nz; // send to cache
    register int ix, iy, iz, ixyz=0;
    register int nxyz = nx*ny*nz;
    
    // for multiplication in momentum space
    register double deltakx = 2.*M_PI/(wd.dx*wd.nx); // send to cache
    register double deltaky = 2.*M_PI/(wd.dy*wd.ny); // send to cache
    register double deltakz = 2.*M_PI/(wd.dz*wd.nz); // send to cache
    
    register double kx, ky, kz; 
    
    // allocate memory
    double complex * wrk; // working array
    cppmallocl(wrk,nx*ny*nz,double complex); 
    
    // Create fftw plans
    fftw_plan plan_f = fftw_plan_dft_r2c_3d(nx, ny, nz, f, wrk, WDERIV_USE_FFTW_PLANNER); 
    fftw_plan plan_b = fftw_plan_dft_c2r_3d(nx, ny, nz, wrk, deriv, WDERIV_USE_FFTW_PLANNER);
    
    // f --FT--> wrk
    fftw_execute(plan_f);
    
    // multiply by momentum
    
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz/2+1; iz++) 
    {
        if(ix<nx/2) kx=deltakx*(ix   );
        else        kx=deltakx*(ix-nx);
        
        if(iy<ny/2) ky=deltaky*(iy   );
        else        ky=deltaky*(iy-ny);
        
        if(iz<nz/2) kz=deltakz*(iz   );
        else        kz=deltakz*(iz-nz);
        
        if(n%2==1) // odd derivative
        {
           if(ix==nx/2) kx=0.0;
           if(iy==ny/2) ky=0.0;
           if(iz==nz/2) kz=0.0;
        }
        
        // multiply
        if     (direction==WDERIV_DX) wrk[ixyz]*=cpow(I*kx, n)/nxyz;
        else if(direction==WDERIV_DY) wrk[ixyz]*=cpow(I*ky, n)/nxyz;
        else if(direction==WDERIV_DZ) wrk[ixyz]*=cpow(I*kz, n)/nxyz;
        else return WDERIV_ERR_INCORRECT_DIRECTION;
        
        ixyz++;
    }
    
    // wrk --IFT--> deriv
    fftw_execute(plan_b);
    
    
    // clear memory
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
    free(wrk);
    
    return WDERIV_OK;
}

int wderiv_dfdx_3d_r(double *f, double *dfdx) {return wderiv_derivative_3d_r(WDERIV_DX, 1, f, dfdx);}
int wderiv_dfdy_3d_r(double *f, double *dfdy) {return wderiv_derivative_3d_r(WDERIV_DY, 1, f, dfdy);}
int wderiv_dfdz_3d_r(double *f, double *dfdz) {return wderiv_derivative_3d_r(WDERIV_DZ, 1, f, dfdz);}
int wderiv_d2fdx2_3d_r(double *f, double *d2fdx2) {return wderiv_derivative_3d_r(WDERIV_DX, 2, f, d2fdx2);}
int wderiv_d2fdy2_3d_r(double *f, double *d2fdy2) {return wderiv_derivative_3d_r(WDERIV_DY, 2, f, d2fdy2);}
int wderiv_d2fdz2_3d_r(double *f, double *d2fdz2) {return wderiv_derivative_3d_r(WDERIV_DZ, 2, f, d2fdz2);}
int wderiv_dnfdxn_3d_r(int n, double *f, double *dnfdxn) {return wderiv_derivative_3d_r(WDERIV_DX, n, f, dnfdxn);}
int wderiv_dnfdyn_3d_r(int n, double *f, double *dnfdyn) {return wderiv_derivative_3d_r(WDERIV_DY, n, f, dnfdyn);}
int wderiv_dnfdzn_3d_r(int n, double *f, double *dnfdzn) {return wderiv_derivative_3d_r(WDERIV_DZ, n, f, dnfdzn);}
int wderiv_gradient_3d_r(double *f, double *dfdx, double *dfdy, double *dfdz)
{
    int ierr;
    ierr=wderiv_derivative_3d_r(WDERIV_DX, 1, f, dfdx); if(ierr!=WDERIV_OK) return ierr;
    ierr=wderiv_derivative_3d_r(WDERIV_DY, 1, f, dfdy); if(ierr!=WDERIV_OK) return ierr;
    ierr=wderiv_derivative_3d_r(WDERIV_DZ, 1, f, dfdz); if(ierr!=WDERIV_OK) return ierr;
    return WDERIV_OK;
}
int wderiv_gradient2_3d_r(double *f, double *gradient_square)
{
    int ierr;
    int i;
    int nxyz = wd.nx*wd.ny*wd.nz;
    double *tmp_dfdx2, *tmp_dfdy2, *tmp_dfdz2;
    cppmallocl(tmp_dfdx2, nxyz, double);
    cppmallocl(tmp_dfdy2, nxyz, double);
    cppmallocl(tmp_dfdz2, nxyz, double);
    
    ierr=wderiv_gradient_3d_r(f ,tmp_dfdx2 ,tmp_dfdy2 ,tmp_dfdz2); if(ierr!=WDERIV_OK) return ierr;

    for(i=0; i<nxyz; i++)
    {
        tmp_dfdx2[i]*=tmp_dfdx2[i];
        tmp_dfdy2[i]*=tmp_dfdy2[i];
        tmp_dfdz2[i]*=tmp_dfdz2[i];
        gradient_square[i] = tmp_dfdx2[i] + tmp_dfdy2[i] + tmp_dfdz2[i];
    }
    
    free(tmp_dfdx2);
    free(tmp_dfdy2);
    free(tmp_dfdz2);
    return WDERIV_OK;
}
int wderiv_laplace_3d_r(double *f, double *laplace)
{
    int ierr;
    int i;
    int nxyz = wd.nx*wd.ny*wd.nz;
    double *tmp_d2fdx2, *tmp_d2fdy2, *tmp_d2fdz2;
    cppmallocl(tmp_d2fdx2, nxyz, double);
    cppmallocl(tmp_d2fdy2, nxyz, double);
    cppmallocl(tmp_d2fdz2, nxyz, double);
    
    ierr = wderiv_d2fdx2_3d_r(f,tmp_d2fdx2); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_d2fdy2_3d_r(f,tmp_d2fdy2); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_d2fdz2_3d_r(f,tmp_d2fdz2); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxyz; i++){        laplace[i]=tmp_d2fdx2[i]+tmp_d2fdy2[i]+tmp_d2fdz2[i];    }

    free(tmp_d2fdx2);
    free(tmp_d2fdy2);
    free(tmp_d2fdz2);
    return WDERIV_OK;
}
int wderiv_divergence_3d_r(double *f, double *divergence) {
    int ierr;
    int i;
    int nxyz = wd.nx*wd.ny*wd.nz;
    double *tmp_dfdx, *tmp_dfdy, *tmp_dfdz;
    cppmallocl(tmp_dfdx, nxyz, double);
    cppmallocl(tmp_dfdy, nxyz, double);
    cppmallocl(tmp_dfdz, nxyz, double);

    ierr = wderiv_gradient_3d_r(f,tmp_dfdx,tmp_dfdy,tmp_dfdz); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxyz; i++){        divergence[i]=tmp_dfdx[i]+tmp_dfdy[i]+tmp_dfdz[i];    }

    free(tmp_dfdx);
    free(tmp_dfdy);
    free(tmp_dfdz);
    return WDERIV_OK;
}
int wderiv_curl_3d_r(double *fx, double *fy, double *fz, double *curl_x, double *curl_y, double *curl_z)
{    
    int ierr;
    int i;
    int nxyz = wd.nx*wd.ny*wd.nz;
    double *tmp_dfxdy, *tmp_dfxdz, *tmp_dfydx, *tmp_dfydz, *tmp_dfzdx, *tmp_dfzdy;
    cppmallocl(tmp_dfxdy, nxyz, double);
    cppmallocl(tmp_dfxdz, nxyz, double);
    cppmallocl(tmp_dfydx, nxyz, double);
    cppmallocl(tmp_dfydz, nxyz, double);
    cppmallocl(tmp_dfzdx, nxyz, double);
    cppmallocl(tmp_dfzdy, nxyz, double);

    ierr = wderiv_derivative_3d_r(WDERIV_DY,1,fx,tmp_dfxdy); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_r(WDERIV_DZ,1,fx,tmp_dfxdz); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_r(WDERIV_DX,1,fy,tmp_dfydx); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_r(WDERIV_DZ,1,fy,tmp_dfydz); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_r(WDERIV_DX,1,fz,tmp_dfzdx); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_r(WDERIV_DY,1,fz,tmp_dfzdy); if(ierr!=WDERIV_OK) return ierr;
    for (i=0; i < nxyz; i++){
        curl_x[i] = tmp_dfzdy[i] - tmp_dfydz[i];
        curl_y[i] = tmp_dfxdz[i] - tmp_dfzdx[i];
        curl_z[i] = tmp_dfydx[i] - tmp_dfxdy[i];
    }
    
    free(tmp_dfxdy);
    free(tmp_dfxdz);
    free(tmp_dfydx);
    free(tmp_dfydz);
    free(tmp_dfzdx);
    free(tmp_dfzdy);
    return WDERIV_OK;    

}

// --- 2D main low level functions REAL ---
int wderiv_derivative_2d_r(int direction, int n, double *f, double *deriv)
{
    if(wd.dim!=2) return WDERIV_ERR_NOINIT;
    
    register int nx=wd.nx, ny=wd.ny; // send to cache
    register int ix, iy, ixy=0;
    register int nxy = nx*ny;
    
    // for multiplication in momentum space
    register double deltakx = 2.*M_PI/(wd.dx*wd.nx); // send to cache
    register double deltaky = 2.*M_PI/(wd.dy*wd.ny); // send to cache
    
    register double kx, ky; 
    
    // allocate memory
    double complex * wrk; // working array
    cppmallocl(wrk,nx*ny,double complex); 
    
    // Create fftw plans
    fftw_plan plan_f = fftw_plan_dft_r2c_2d(nx, ny, f, wrk, WDERIV_USE_FFTW_PLANNER); 
    fftw_plan plan_b = fftw_plan_dft_c2r_2d(nx, ny, wrk, deriv, WDERIV_USE_FFTW_PLANNER);
    
    // f --FT--> wrk
    fftw_execute(plan_f);
    
    // multiply by momentum
    
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny/2+1; iy++)
    {
        if(ix<nx/2) kx=deltakx*(ix   );
        else        kx=deltakx*(ix-nx);
        
        if(iy<ny/2) ky=deltaky*(iy   );
        else        ky=deltaky*(iy-ny);
        
        if(n%2==1) // odd derivative
        {
           if(ix==nx/2) kx=0.0;
           if(iy==ny/2) ky=0.0;
        }
        
        // multiply
        if     (direction==WDERIV_DX) wrk[ixy]*=cpow(I*kx, n)/nxy;
        else if(direction==WDERIV_DY) wrk[ixy]*=cpow(I*ky, n)/nxy;
        else return WDERIV_ERR_INCORRECT_DIRECTION;
        
        ixy++;
    }
    
    // wrk --IFT--> deriv
    fftw_execute(plan_b);
    
    
    // clear memory
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
    free(wrk);
    
    return WDERIV_OK;
}

int wderiv_dfdx_2d_r(double *f, double *dfdx) {return wderiv_derivative_2d_r(WDERIV_DX, 1, f, dfdx);}
int wderiv_dfdy_2d_r(double *f, double *dfdy) {return wderiv_derivative_2d_r(WDERIV_DY, 1, f, dfdy);}
int wderiv_d2fdx2_2d_r(double *f, double *d2fdx2) {return wderiv_derivative_2d_r(WDERIV_DX, 2, f, d2fdx2);}
int wderiv_d2fdy2_2d_r(double *f, double *d2fdy2) {return wderiv_derivative_2d_r(WDERIV_DY, 2, f, d2fdy2);}
int wderiv_dnfdxn_2d_r(int n, double *f, double *dnfdxn) {return wderiv_derivative_2d_r(WDERIV_DX, n, f, dnfdxn);}
int wderiv_dnfdyn_2d_r(int n, double *f, double *dnfdyn) {return wderiv_derivative_2d_r(WDERIV_DY, n, f, dnfdyn);}
int wderiv_gradient_2d_r(double *f, double *dfdx, double *dfdy) 
{ 
    int ierr;
    ierr=wderiv_derivative_2d_r(WDERIV_DX, 1, f, dfdx); if(ierr!=WDERIV_OK) return ierr;
    ierr=wderiv_derivative_2d_r(WDERIV_DY, 1, f, dfdy); if(ierr!=WDERIV_OK) return ierr;
    return WDERIV_OK;
}
int wderiv_gradient2_2d_r(double *f, double *gradient_square)
{
    int ierr;
    int i;
    int nxy = wd.nx*wd.ny;
    double *tmp_dfdx2, *tmp_dfdy2;
    cppmallocl(tmp_dfdx2, nxy, double);
    cppmallocl(tmp_dfdy2, nxy, double);
    
    ierr=wderiv_gradient_2d_r(f ,tmp_dfdx2 ,tmp_dfdy2); if(ierr!=WDERIV_OK) return ierr;

    for(i=0; i<nxy; i++)
    {
        tmp_dfdx2[i]*=tmp_dfdx2[i];
        tmp_dfdy2[i]*=tmp_dfdy2[i];
        gradient_square[i] = tmp_dfdx2[i] + tmp_dfdy2[i];
    }
    
    free(tmp_dfdx2);
    free(tmp_dfdy2);
    return WDERIV_OK;
}
int wderiv_laplace_2d_r(double *f, double *laplace)
{
    int ierr;
    int i;
    double *tmp_d2fdx2, *tmp_d2fdy2;
    int nxyz = wd.nx*wd.ny;
    cppmallocl(tmp_d2fdx2, nxyz, double);
    cppmallocl(tmp_d2fdy2, nxyz, double);
    
    ierr = wderiv_d2fdx2_2d_r(f,tmp_d2fdx2); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_d2fdy2_2d_r(f,tmp_d2fdy2); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxyz; i++){ laplace[i]=tmp_d2fdx2[i]+tmp_d2fdy2[i]; }

    free(tmp_d2fdx2);
    free(tmp_d2fdy2);
    return WDERIV_OK;
}
int wderiv_divergence_2d_r(double *f, double *divergence) 
{
    int ierr;
    int i;
    double *tmp_dfdx, *tmp_dfdy;
    int nxy = wd.nx*wd.ny;
    cppmallocl(tmp_dfdx, nxy, double);
    cppmallocl(tmp_dfdy, nxy, double);
    
    ierr = wderiv_gradient_2d_r(f,tmp_dfdx,tmp_dfdy); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxy; i++){ divergence[i]=tmp_dfdx[i]+tmp_dfdy[i]; }
    
    free(tmp_dfdx);
    free(tmp_dfdy);
    return WDERIV_OK;
}



// --- 1D main low level functions REAL ---
int wderiv_derivative_1d_r(int n, double *f, double *deriv)
{
    if(wd.dim!=1) return WDERIV_ERR_NOINIT;
    
    register int nx=wd.nx; // send to cache
    register int ix;
    
    // for multiplication in momentum space
    register double deltakx = 2.*M_PI/(wd.dx*wd.nx); // send to cache
    
    register double kx; 
    
    // allocate memory
    double complex * wrk; // working array
    cppmallocl(wrk,nx,double complex); 
    
    // Create fftw plans
    fftw_plan plan_f = fftw_plan_dft_r2c_1d(nx, f, wrk, WDERIV_USE_FFTW_PLANNER); 
    fftw_plan plan_b = fftw_plan_dft_c2r_1d(nx, wrk, deriv, WDERIV_USE_FFTW_PLANNER);
    
    // f --FT--> wrk
    fftw_execute(plan_f);
    
    // multiply by momentum
    
    for(ix=0; ix<nx/2+1; ix++)
    {
        if(ix<nx/2) kx=deltakx*(ix   );
        else        kx=deltakx*(ix-nx);
        
        if(n%2==1) // odd derivative
        {
           if(ix==nx/2) kx=0.0;
        }
        
        // multiply
        wrk[ix]*=cpow(I*kx, n)/nx;

    }
    
    // wrk --IFT--> deriv
    fftw_execute(plan_b);
    
    
    // clear memory
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
    free(wrk);
    
    return WDERIV_OK;
}

int wderiv_dfdx_1d_r(double *f, double *dfdx) {return wderiv_derivative_1d_r(1, f, dfdx);}
int wderiv_d2fdx2_1d_r(double *f, double *d2fdx2) {return wderiv_derivative_1d_r(2, f, d2fdx2);}
int wderiv_dnfdxn_1d_r(int n, double *f, double *dnfdxn) {return wderiv_derivative_1d_r(n, f, dnfdxn);}
int wderiv_gradient_1d_r(double *f, double *dfdx) {return wderiv_dfdx_1d_r(f, dfdx);}
int wderiv_gradient2_1d_r(double *f, double *gradient_square) {
    int ierr;
    int i;
    int nx = wd.nx;
    double *tmp_dfdx2;
    cppmallocl(tmp_dfdx2, nx, double);
    
    ierr=wderiv_gradient_1d_r(f ,tmp_dfdx2); if(ierr!=WDERIV_OK) return ierr;

    for(i=0; i<nx; i++)
    {
        tmp_dfdx2[i]*=tmp_dfdx2[i];
        gradient_square[i] = tmp_dfdx2[i];
    }
    
    free(tmp_dfdx2);
    return WDERIV_OK;
}
int wderiv_laplace_1d_r(double *f, double *laplace) {return wderiv_d2fdx2_1d_r(f, laplace);}
int wderiv_divergence_1d_r(double *f, double *divergence) 
{
    int ierr;
    int i;
    double *tmp_dfdx;
    int nx = wd.nx;
    cppmallocl(tmp_dfdx, nx, double);
    
    ierr = wderiv_gradient_1d_r(f,tmp_dfdx); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nx; i++){ divergence[i]=tmp_dfdx[i]; }
    
    free(tmp_dfdx);
    return WDERIV_OK;
}

// --- 3D main low level functions COMPLEX ---
int wderiv_derivative_3d_c(int direction, int n, double complex *f, double complex *deriv)
{
    if(wd.dim!=3) return WDERIV_ERR_NOINIT;
    
    register int nx=wd.nx, ny=wd.ny, nz=wd.nz; // send to cache
    register int ix, iy, iz, ixyz=0;
    register int nxyz = nx*ny*nz;
    
    // for multiplication in momentum space
    register double deltakx = 2.*M_PI/(wd.dx*wd.nx); // send to cache
    register double deltaky = 2.*M_PI/(wd.dy*wd.ny); // send to cache
    register double deltakz = 2.*M_PI/(wd.dz*wd.nz); // send to cache
    
    register double kx, ky, kz; 
    
    // allocate memory
    double complex * wrk; // working array
    cppmallocl(wrk,nx*ny*nz,double complex); 
    
    // Create fftw plans
    fftw_plan plan_f = fftw_plan_dft_3d(nx, ny, nz, f, wrk, FFTW_FORWARD, WDERIV_USE_FFTW_PLANNER); 
    fftw_plan plan_b = fftw_plan_dft_3d(nx, ny, nz, wrk, deriv, FFTW_BACKWARD, WDERIV_USE_FFTW_PLANNER);
    
    // f --FT--> wrk
    fftw_execute(plan_f);
    
    // multiply by momentum
    
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++) 
    {
        if(ix<nx/2) kx=deltakx*(ix   );
        else        kx=deltakx*(ix-nx);
        
        if(iy<ny/2) ky=deltaky*(iy   );
        else        ky=deltaky*(iy-ny);
        
        if(iz<nz/2) kz=deltakz*(iz   );
        else        kz=deltakz*(iz-nz);
        
        if(n%2==1) // odd derivative
        {
           if(ix==nx/2) kx=0.0;
           if(iy==ny/2) ky=0.0;
           if(iz==nz/2) kz=0.0;
        }
        
        // multiply
        if     (direction==WDERIV_DX) wrk[ixyz]*=cpow(I*kx, n)/nxyz;
        else if(direction==WDERIV_DY) wrk[ixyz]*=cpow(I*ky, n)/nxyz;
        else if(direction==WDERIV_DZ) wrk[ixyz]*=cpow(I*kz, n)/nxyz;
        else return WDERIV_ERR_INCORRECT_DIRECTION;
        
        ixyz++;
    }
    
    // wrk --IFT--> deriv
    fftw_execute(plan_b);
    
    // clear memory
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
    free(wrk);
    
    return WDERIV_OK;
}

int wderiv_dfdx_3d_c(double complex *f, double complex *dfdx) {return wderiv_derivative_3d_c(WDERIV_DX, 1, f, dfdx);}
int wderiv_dfdy_3d_c(double complex *f, double complex *dfdy) {return wderiv_derivative_3d_c(WDERIV_DY, 1, f, dfdy);}
int wderiv_dfdz_3d_c(double complex *f, double complex *dfdz) {return wderiv_derivative_3d_c(WDERIV_DZ, 1, f, dfdz);}
int wderiv_d2fdx2_3d_c(double complex *f, double complex *d2fdx2) {return wderiv_derivative_3d_c(WDERIV_DX, 2, f, d2fdx2);}
int wderiv_d2fdy2_3d_c(double complex *f, double complex *d2fdy2) {return wderiv_derivative_3d_c(WDERIV_DY, 2, f, d2fdy2);}
int wderiv_d2fdz2_3d_c(double complex *f, double complex *d2fdz2) {return wderiv_derivative_3d_c(WDERIV_DZ, 2, f, d2fdz2);}
int wderiv_dnfdxn_3d_c(int n, double complex *f, double complex *dnfdxn) {return wderiv_derivative_3d_c(WDERIV_DX, n, f, dnfdxn);}
int wderiv_dnfdyn_3d_c(int n, double complex *f, double complex *dnfdyn) {return wderiv_derivative_3d_c(WDERIV_DY, n, f, dnfdyn);}
int wderiv_dnfdzn_3d_c(int n, double complex *f, double complex *dnfdzn) {return wderiv_derivative_3d_c(WDERIV_DZ, n, f, dnfdzn);}
int wderiv_gradient_3d_c(double complex *f, double complex *dfdx, double complex *dfdy, double complex *dfdz)
{
    int ierr;
    ierr=wderiv_derivative_3d_c(WDERIV_DX, 1, f, dfdx); if(ierr!=WDERIV_OK) return ierr;
    ierr=wderiv_derivative_3d_c(WDERIV_DY, 1, f, dfdy); if(ierr!=WDERIV_OK) return ierr;
    ierr=wderiv_derivative_3d_c(WDERIV_DZ, 1, f, dfdz); if(ierr!=WDERIV_OK) return ierr;
    return WDERIV_OK;
}
int wderiv_gradient2_3d_c(double complex *f, double complex *gradient_square)
{
    int ierr;
    int i;
    int nxyz = wd.nx*wd.ny*wd.nz;
    double complex *tmp_dfdx2, *tmp_dfdy2, *tmp_dfdz2;
    cppmallocl(tmp_dfdx2, nxyz, double complex);
    cppmallocl(tmp_dfdy2, nxyz, double complex);
    cppmallocl(tmp_dfdz2, nxyz, double complex);
    
    ierr=wderiv_gradient_3d_c(f ,tmp_dfdx2 ,tmp_dfdy2 ,tmp_dfdz2); if(ierr!=WDERIV_OK) return ierr;

    for(i=0; i<nxyz; i++)
    {
        tmp_dfdx2[i]=tmp_dfdx2[i] * conj(tmp_dfdx2[i]);
        tmp_dfdy2[i]=tmp_dfdy2[i] * conj(tmp_dfdy2[i]);
        tmp_dfdz2[i]=tmp_dfdz2[i] * conj(tmp_dfdz2[i]);
        gradient_square[i] = tmp_dfdx2[i] + tmp_dfdy2[i] + tmp_dfdz2[i];
    }
    
    free(tmp_dfdx2);
    free(tmp_dfdy2);
    free(tmp_dfdz2);    
    return WDERIV_OK;
}
int wderiv_laplace_3d_c(double complex *f, double complex *laplace)
{
    int ierr;
    int i;
    int nxyz = wd.nx*wd.ny*wd.nz;
    double complex *tmp_d2fdx2, *tmp_d2fdy2, *tmp_d2fdz2;
    cppmallocl(tmp_d2fdx2, nxyz, double complex);
    cppmallocl(tmp_d2fdy2, nxyz, double complex);
    cppmallocl(tmp_d2fdz2, nxyz, double complex);
    
    ierr = wderiv_d2fdx2_3d_c(f,tmp_d2fdx2); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_d2fdy2_3d_c(f,tmp_d2fdy2); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_d2fdz2_3d_c(f,tmp_d2fdz2); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxyz; i++){    laplace[i]=tmp_d2fdx2[i]+tmp_d2fdy2[i]+tmp_d2fdz2[i];    }

    free(tmp_d2fdx2);
    free(tmp_d2fdy2);
    free(tmp_d2fdz2);
    return WDERIV_OK;
}
int wderiv_divergence_3d_c(double complex *f, double complex *divergence) 
{
    int ierr;
    int i;
    double complex *tmp_dfdx, *tmp_dfdy, *tmp_dfdz;
    int nxyz = wd.nx*wd.ny*wd.nz;
    cppmallocl(tmp_dfdx, nxyz, double complex);
    cppmallocl(tmp_dfdy, nxyz, double complex);
    cppmallocl(tmp_dfdz, nxyz, double complex);
    
    ierr = wderiv_gradient_3d_c(f,tmp_dfdx,tmp_dfdy,tmp_dfdz); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxyz; i++){        divergence[i]=tmp_dfdx[i]+tmp_dfdy[i]+tmp_dfdz[i];    }
    
    free(tmp_dfdx);
    free(tmp_dfdy);
    free(tmp_dfdz);
    return WDERIV_OK;
}
int wderiv_curl_3d_c(double complex *fx, double complex *fy, double complex *fz, double complex *curl_x, double complex *curl_y, double complex *curl_z){
    
    int ierr;
    int i;
    int nxyz = wd.nx*wd.ny*wd.nz;
    double complex *tmp_dfxdy, *tmp_dfxdz, *tmp_dfydx, *tmp_dfydz, *tmp_dfzdx, *tmp_dfzdy;
    cppmallocl(tmp_dfxdy, nxyz, double complex);
    cppmallocl(tmp_dfxdz, nxyz, double complex);
    cppmallocl(tmp_dfydx, nxyz, double complex);
    cppmallocl(tmp_dfydz, nxyz, double complex);
    cppmallocl(tmp_dfzdx, nxyz, double complex);
    cppmallocl(tmp_dfzdy, nxyz, double complex);
    
    ierr = wderiv_derivative_3d_c(WDERIV_DY,1,fx,tmp_dfxdy); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_c(WDERIV_DZ,1,fx,tmp_dfxdz); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_c(WDERIV_DX,1,fy,tmp_dfydx); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_c(WDERIV_DZ,1,fy,tmp_dfydz); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_c(WDERIV_DX,1,fz,tmp_dfzdx); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_derivative_3d_c(WDERIV_DY,1,fz,tmp_dfzdy); if(ierr!=WDERIV_OK) return ierr;
    for (i=0; i < nxyz; i++)
    {
        curl_x[i] = tmp_dfzdy[i] - tmp_dfydz[i];
        curl_y[i] = tmp_dfxdz[i] - tmp_dfzdx[i];
        curl_z[i] = tmp_dfydx[i] - tmp_dfxdy[i];
    }
    
    free(tmp_dfxdy);
    free(tmp_dfxdz);
    free(tmp_dfydx);
    free(tmp_dfydz);
    free(tmp_dfzdx);
    free(tmp_dfzdy);
    return WDERIV_OK;    

}

// --- 2D main low level functions COMPLEX ---
int wderiv_derivative_2d_c(int direction, int n, double complex *f, double complex *deriv)
{
    if(wd.dim!=2) return WDERIV_ERR_NOINIT;
    
    register int nx=wd.nx, ny=wd.ny; // send to cache
    register int ix, iy, ixy=0;
    register int nxy = nx*ny;
    
    // for multiplication in momentum space
    register double deltakx = 2.*M_PI/(wd.dx*wd.nx); // send to cache
    register double deltaky = 2.*M_PI/(wd.dy*wd.ny); // send to cache
    
    register double kx, ky; 
    
    // allocate memory
    double complex * wrk; // working array
    cppmallocl(wrk,nx*ny,double complex); 
    
    // Create fftw plans
    fftw_plan plan_f = fftw_plan_dft_2d(nx, ny, f, wrk, FFTW_FORWARD, WDERIV_USE_FFTW_PLANNER); 
    fftw_plan plan_b = fftw_plan_dft_2d(nx, ny, wrk, deriv, FFTW_BACKWARD, WDERIV_USE_FFTW_PLANNER);
    
    // f --FT--> wrk
    fftw_execute(plan_f);
    
    // multiply by momentum
    
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) 
    {
        if(ix<nx/2) kx=deltakx*(ix   );
        else        kx=deltakx*(ix-nx);
        
        if(iy<ny/2) ky=deltaky*(iy   );
        else        ky=deltaky*(iy-ny);
        
        if(n%2==1) // odd derivative
        {
           if(ix==nx/2) kx=0.0;
           if(iy==ny/2) ky=0.0;
        }
        
        // multiply
        if     (direction==WDERIV_DX) wrk[ixy]*=cpow(I*kx, n)/nxy;
        else if(direction==WDERIV_DY) wrk[ixy]*=cpow(I*ky, n)/nxy;
        else return WDERIV_ERR_INCORRECT_DIRECTION;
        
        ixy++;
    }
    
    // wrk --IFT--> deriv
    fftw_execute(plan_b);
    
    // clear memory
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
    free(wrk);
    
    return WDERIV_OK;
}

int wderiv_dfdx_2d_c(double complex *f, double complex *dfdx) {return wderiv_derivative_2d_c(WDERIV_DX, 1, f, dfdx);}
int wderiv_dfdy_2d_c(double complex *f, double complex *dfdy) {return wderiv_derivative_2d_c(WDERIV_DY, 1, f, dfdy);}
int wderiv_d2fdx2_2d_c(double complex *f, double complex *d2fdx2) {return wderiv_derivative_2d_c(WDERIV_DX, 2, f, d2fdx2);}
int wderiv_d2fdy2_2d_c(double complex *f, double complex *d2fdy2) {return wderiv_derivative_2d_c(WDERIV_DY, 2, f, d2fdy2);}
int wderiv_dnfdxn_2d_c(int n, double complex *f, double complex *dnfdxn) {return wderiv_derivative_2d_c(WDERIV_DX, n, f, dnfdxn);}
int wderiv_dnfdyn_2d_c(int n, double complex *f, double complex *dnfdyn) {return wderiv_derivative_2d_c(WDERIV_DY, n, f, dnfdyn);}
int wderiv_gradient_2d_c(double complex *f, double complex *dfdx, double complex *dfdy) 
{ 
    int ierr;
    ierr=wderiv_derivative_2d_c(WDERIV_DX, 1, f, dfdx); if(ierr!=WDERIV_OK) return ierr;
    ierr=wderiv_derivative_2d_c(WDERIV_DY, 1, f, dfdy); if(ierr!=WDERIV_OK) return ierr;
    return WDERIV_OK;
}
int wderiv_gradient2_2d_c(double complex *f, double complex *gradient_square)
{
    int ierr;
    int i;
    int nxy = wd.nx*wd.ny;
    double complex *tmp_dfdx2, *tmp_dfdy2;
    cppmallocl(tmp_dfdx2, nxy, double complex);
    cppmallocl(tmp_dfdy2, nxy, double complex);
    
    ierr=wderiv_gradient_2d_c(f ,tmp_dfdx2 ,tmp_dfdy2); if(ierr!=WDERIV_OK) return ierr;

    for(i=0; i<nxy; i++)
    {
        tmp_dfdx2[i]=tmp_dfdx2[i] * conj(tmp_dfdx2[i]);
        tmp_dfdy2[i]=tmp_dfdy2[i] * conj(tmp_dfdy2[i]);
        gradient_square[i] = tmp_dfdx2[i] + tmp_dfdy2[i];
    }
    
    free(tmp_dfdx2);
    free(tmp_dfdy2);
    return WDERIV_OK;
}
int wderiv_laplace_2d_c(double complex *f, double complex *laplace)
{
    int ierr;
    int i;
    double complex *tmp_d2fdx2, *tmp_d2fdy2;
    int nxyz = wd.nx*wd.ny;
    cppmallocl(tmp_d2fdx2, nxyz, double complex);
    cppmallocl(tmp_d2fdy2, nxyz, double complex);
    
    ierr = wderiv_d2fdx2_2d_c(f,tmp_d2fdx2); if(ierr!=WDERIV_OK) return ierr;
    ierr = wderiv_d2fdy2_2d_c(f,tmp_d2fdy2); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxyz; i++){        laplace[i]=tmp_d2fdx2[i]+tmp_d2fdy2[i];    }

    free(tmp_d2fdx2);
    free(tmp_d2fdy2);
    return WDERIV_OK;

}
int wderiv_divergence_2d_c(double complex *f, double complex *divergence) 
{
    int ierr;
    int i;
    double complex *tmp_dfdx, *tmp_dfdy;
    int nxy = wd.nx*wd.ny;
    cppmallocl(tmp_dfdx, nxy, double complex);
    cppmallocl(tmp_dfdy, nxy, double complex);
    
    ierr = wderiv_gradient_2d_c(f,tmp_dfdx,tmp_dfdy); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nxy; i++){        divergence[i]=tmp_dfdx[i]+tmp_dfdy[i];    }
    
    free(tmp_dfdx);
    free(tmp_dfdy);
    return WDERIV_OK;
}


// --- 1D main low level functions COMPLEX ---
int wderiv_derivative_1d_c(int n, double complex *f, double complex *deriv)
{
    if(wd.dim!=1) return WDERIV_ERR_NOINIT;
    
    register int nx=wd.nx; // send to cache
    register int ix=0;
    
    // for multiplication in momentum space
    register double deltakx = 2.*M_PI/(wd.dx*wd.nx); // send to cache
    
    register double kx; 
    
    // allocate memory
    double complex * wrk; // working array
    cppmallocl(wrk,nx,double complex); 
    
    // Create fftw plans
    fftw_plan plan_f = fftw_plan_dft_1d(nx, f, wrk, FFTW_FORWARD, WDERIV_USE_FFTW_PLANNER); 
    fftw_plan plan_b = fftw_plan_dft_1d(nx, wrk, deriv, FFTW_BACKWARD, WDERIV_USE_FFTW_PLANNER);
    
    // f --FT--> wrk
    fftw_execute(plan_f);
    
    // multiply by momentum
    
    for(ix=0; ix<nx; ix++)
    {
        if(ix<nx/2) kx=deltakx*(ix   );
        else        kx=deltakx*(ix-nx);
        
        if(n%2==1) // odd derivative
        {
           if(ix==nx/2) kx=0.0;
        }
        
        // multiply
        wrk[ix]*=cpow(I*kx, n)/nx;
        
    }
    
    // wrk --IFT--> deriv
    fftw_execute(plan_b);
    
    // clear memory
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
    free(wrk);
    
    return WDERIV_OK;
}

int wderiv_dfdx_1d_c(double complex *f, double complex *dfdx) {return wderiv_derivative_1d_c(1, f, dfdx);}
int wderiv_d2fdx2_1d_c(double complex *f, double complex *d2fdx2) {return wderiv_derivative_1d_c(2, f, d2fdx2);}
int wderiv_dnfdxn_1d_c(int n, double complex *f, double complex *dnfdxn) {return wderiv_derivative_1d_c(n, f, dnfdxn);}
int wderiv_gradient_1d_c(double complex *f, double complex *dfdx) {return wderiv_dfdx_1d_c(f, dfdx);}
int wderiv_gradient2_1d_c(double complex *f, double complex *gradient_square) 
{
    int ierr;
    int i;
    int nx = wd.nx;
    double complex *tmp_dfdx2;
    cppmallocl(tmp_dfdx2, nx, double complex);
    
    ierr=wderiv_gradient_1d_c(f ,tmp_dfdx2); if(ierr!=WDERIV_OK) return ierr;

    for(i=0; i<nx; i++)
    {
        tmp_dfdx2[i]=tmp_dfdx2[i] * conj(tmp_dfdx2[i]);
        gradient_square[i] = tmp_dfdx2[i];
    }
    
    free(tmp_dfdx2);
    return WDERIV_OK;
}
int wderiv_laplace_1d_c(double complex *f, double complex *laplace) {return wderiv_d2fdx2_1d_c(f, laplace);}
int wderiv_divergence_1d_c(double complex *f, double complex *divergence) 
{
    int ierr;
    int i;
    double complex *tmp_dfdx;
    int nx = wd.nx;
    cppmallocl(tmp_dfdx, nx, double complex);
    
    ierr = wderiv_gradient_1d_c(f,tmp_dfdx); if(ierr!=WDERIV_OK) return ierr;
    for(i=0; i<nx; i++){        divergence[i]=tmp_dfdx[i];    }
    
    free(tmp_dfdx);
    return WDERIV_OK;
}


// --- GENERIC functions ---
int wderiv_dfdx(int datadim, char type, void *f, void *dfdx)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=1 && datadim!=2 && datadim!=3) return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==1 && type=='r') return wderiv_dfdx_1d_r((double         *)f, (double         *)dfdx);
    else if(datadim==1 && type=='c') return wderiv_dfdx_1d_c((double complex *)f, (double complex *)dfdx);
    else if(datadim==2 && type=='r') return wderiv_dfdx_2d_r((double         *)f, (double         *)dfdx);
    else if(datadim==2 && type=='c') return wderiv_dfdx_2d_c((double complex *)f, (double complex *)dfdx);
    else if(datadim==3 && type=='r') return wderiv_dfdx_3d_r((double         *)f, (double         *)dfdx);
    else if(datadim==3 && type=='c') return wderiv_dfdx_3d_c((double complex *)f, (double complex *)dfdx);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_dfdy(int datadim, char type, void *f, void *dfdy)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=2 && datadim!=3)               return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==2 && type=='r') return wderiv_dfdy_2d_r((double         *)f, (double         *)dfdy);
    else if(datadim==2 && type=='c') return wderiv_dfdy_2d_c((double complex *)f, (double complex *)dfdy);
    else if(datadim==3 && type=='r') return wderiv_dfdy_3d_r((double         *)f, (double         *)dfdy);
    else if(datadim==3 && type=='c') return wderiv_dfdy_3d_c((double complex *)f, (double complex *)dfdy);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_dfdz(int datadim, char type, void *f, void *dfdz)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=3)                             return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==3 && type=='r') return wderiv_dfdz_3d_r((double         *)f, (double         *)dfdz);
    else if(datadim==3 && type=='c') return wderiv_dfdz_3d_c((double complex *)f, (double complex *)dfdz);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_d2fdx2(int datadim, char type, void *f, void *d2fdx2)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=1 && datadim!=2 && datadim!=3) return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==1 && type=='r') return wderiv_d2fdx2_1d_r((double         *)f, (double         *)d2fdx2);
    else if(datadim==1 && type=='c') return wderiv_d2fdx2_1d_c((double complex *)f, (double complex *)d2fdx2);
    else if(datadim==2 && type=='r') return wderiv_d2fdx2_2d_r((double         *)f, (double         *)d2fdx2);
    else if(datadim==2 && type=='c') return wderiv_d2fdx2_2d_c((double complex *)f, (double complex *)d2fdx2);
    else if(datadim==3 && type=='r') return wderiv_d2fdx2_3d_r((double         *)f, (double         *)d2fdx2);
    else if(datadim==3 && type=='c') return wderiv_d2fdx2_3d_c((double complex *)f, (double complex *)d2fdx2);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_d2fdy2(int datadim, char type, void *f, void *d2fdy2)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=2 && datadim!=3)               return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==2 && type=='r') return wderiv_d2fdy2_2d_r((double         *)f, (double         *)d2fdy2);
    else if(datadim==2 && type=='c') return wderiv_d2fdy2_2d_c((double complex *)f, (double complex *)d2fdy2);
    else if(datadim==3 && type=='r') return wderiv_d2fdy2_3d_r((double         *)f, (double         *)d2fdy2);
    else if(datadim==3 && type=='c') return wderiv_d2fdy2_3d_c((double complex *)f, (double complex *)d2fdy2);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_d2fdz2(int datadim, char type, void *f, void *d2fdz2)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=3)                             return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==3 && type=='r') return wderiv_d2fdz2_3d_r((double         *)f, (double         *)d2fdz2);
    else if(datadim==3 && type=='c') return wderiv_d2fdz2_3d_c((double complex *)f, (double complex *)d2fdz2);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_dnfdxn(int datadim, char type, int n, void *f, void *dnfdxn)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=1 && datadim!=2 && datadim!=3) return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==1 && type=='r') return wderiv_dnfdxn_1d_r(n, (double         *)f, (double         *)dnfdxn);
    else if(datadim==1 && type=='c') return wderiv_dnfdxn_1d_c(n, (double complex *)f, (double complex *)dnfdxn);
    else if(datadim==2 && type=='r') return wderiv_dnfdxn_2d_r(n, (double         *)f, (double         *)dnfdxn);
    else if(datadim==2 && type=='c') return wderiv_dnfdxn_2d_c(n, (double complex *)f, (double complex *)dnfdxn);
    else if(datadim==3 && type=='r') return wderiv_dnfdxn_3d_r(n, (double         *)f, (double         *)dnfdxn);
    else if(datadim==3 && type=='c') return wderiv_dnfdxn_3d_c(n, (double complex *)f, (double complex *)dnfdxn);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_dnfdyn(int datadim, char type, int n, void *f, void *dnfdyn)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=2 && datadim!=3)               return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==2 && type=='r') return wderiv_dnfdyn_2d_r(n, (double         *)f, (double         *)dnfdyn);
    else if(datadim==2 && type=='c') return wderiv_dnfdyn_2d_c(n, (double complex *)f, (double complex *)dnfdyn);
    else if(datadim==3 && type=='r') return wderiv_dnfdyn_3d_r(n, (double         *)f, (double         *)dnfdyn);
    else if(datadim==3 && type=='c') return wderiv_dnfdyn_3d_c(n, (double complex *)f, (double complex *)dnfdyn);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_dnfdzn(int datadim, char type, int n,  void *f, void *dnfdzn)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=3)                             return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==3 && type=='r') return wderiv_dnfdzn_3d_r(n, (double         *)f, (double         *)dnfdzn);
    else if(datadim==3 && type=='c') return wderiv_dnfdzn_3d_c(n, (double complex *)f, (double complex *)dnfdzn);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_gradient(int datadim, char type, void *f, void *dfdx, void *dfdy, void *dfdz)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=1 && datadim!=2 && datadim!=3) return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==1 && type=='r') return wderiv_gradient_1d_r((double         *)f, (double         *)dfdx);
    else if(datadim==1 && type=='c') return wderiv_gradient_1d_c((double complex *)f, (double complex *)dfdx);
    else if(datadim==2 && type=='r') return wderiv_gradient_2d_r((double         *)f, (double         *)dfdx, (double         *)dfdy);
    else if(datadim==2 && type=='c') return wderiv_gradient_2d_c((double complex *)f, (double complex *)dfdx, (double complex *)dfdy);
    else if(datadim==3 && type=='r') return wderiv_gradient_3d_r((double         *)f, (double         *)dfdx, (double         *)dfdy, (double         *)dfdz);
    else if(datadim==3 && type=='c') return wderiv_gradient_3d_c((double complex *)f, (double complex *)dfdx, (double complex *)dfdy, (double complex *)dfdz);

    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_gradient2(int datadim, char type, void *f, void *gradient_squared)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=1 && datadim!=2 && datadim!=3) return WDERIV_ERR_INCORRECT_DATADIM;

    if     (datadim==1 && type=='r') return wderiv_gradient2_1d_r((double         *)f, (double         *)gradient_squared);
    else if(datadim==1 && type=='c') return wderiv_gradient2_1d_c((double complex *)f, (double complex *)gradient_squared);
    else if(datadim==2 && type=='r') return wderiv_gradient2_2d_r((double         *)f, (double         *)gradient_squared);
    else if(datadim==2 && type=='c') return wderiv_gradient2_2d_c((double complex *)f, (double complex *)gradient_squared);
    else if(datadim==3 && type=='r') return wderiv_gradient2_3d_r((double         *)f, (double         *)gradient_squared);
    else if(datadim==3 && type=='c') return wderiv_gradient2_3d_c((double complex *)f, (double complex *)gradient_squared);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_laplace(int datadim, char type, void *f, void *laplace)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=1 && datadim!=2 && datadim!=3) return WDERIV_ERR_INCORRECT_DATADIM;
      
    if     (datadim==1 && type=='r') return wderiv_laplace_1d_r((double         *)f, (double         *)laplace);
    else if(datadim==1 && type=='c') return wderiv_laplace_1d_c((double complex *)f, (double complex *)laplace);
    else if(datadim==2 && type=='r') return wderiv_laplace_2d_r((double         *)f, (double         *)laplace);
    else if(datadim==2 && type=='c') return wderiv_laplace_2d_c((double complex *)f, (double complex *)laplace);
    else if(datadim==3 && type=='r') return wderiv_laplace_3d_r((double         *)f, (double         *)laplace);
    else if(datadim==3 && type=='c') return wderiv_laplace_3d_c((double complex *)f, (double complex *)laplace);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_divergence(int datadim, char type, void *f, void *divergence)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=1 && datadim!=2 && datadim!=3) return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==1 && type=='r') return wderiv_divergence_1d_r((double         *)f, (double         *)divergence);
    else if(datadim==1 && type=='c') return wderiv_divergence_1d_c((double complex *)f, (double complex *)divergence);
    else if(datadim==2 && type=='r') return wderiv_divergence_2d_r((double         *)f, (double         *)divergence);
    else if(datadim==2 && type=='c') return wderiv_divergence_2d_c((double complex *)f, (double complex *)divergence);
    else if(datadim==3 && type=='r') return wderiv_divergence_3d_r((double         *)f, (double         *)divergence);
    else if(datadim==3 && type=='c') return wderiv_divergence_3d_c((double complex *)f, (double complex *)divergence);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}
int wderiv_curl(int datadim, char type, void *fx, void *fy, void *fz, void *curl_x, void *curl_y, void *curl_z)
{
    if(type!='r' && type!='c')                 return WDERIV_ERR_INCORRECT_DATATYPE;
    if(datadim!=3)                             return WDERIV_ERR_INCORRECT_DATADIM;
    
    if     (datadim==3 && type=='r') return wderiv_curl_3d_r((double         *)fx, (double         *)fy, (double         *)fz, (double         *)curl_x, (double         *)curl_y, (double         *)curl_z);
    else if(datadim==3 && type=='c') return wderiv_curl_3d_c((double complex *)fx, (double complex *)fy, (double complex *)fz, (double complex *)curl_x, (double complex *)curl_y, (double complex *)curl_z);
    else                             return WDERIV_ERR_UNDEFINED;
    return WDERIV_OK;
}