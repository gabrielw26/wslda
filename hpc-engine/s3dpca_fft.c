// Author: Gabriel Wlazlowski
// File that provides support for fftw operations

#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include "pca_settings.h"
#include "pca_macro.h"
#include "s3dpca_fft.h"

/**
 * Function creates plans and allocates memory
 * @param mdfft pointer to structure holding metadata for fft handling
 * @param batch number of vectors transformed by plan many  
 * */
int create_fft_plans(metadata_s3dpca_fft *mdfft, int batch)
{
    mdfft->batch=batch;
    
    // allocate memory for buffers
    cppmallocl(mdfft->fft3,NXYZ,double complex); 
//     cppmallocl(mdfft->fft3many,batch*NXYZ,double complex);
    cppmallocl(mdfft->fft3grad,2*3*NXYZ,double complex); // (u,v) * (dx,dy,dz)
    cppmallocl(mdfft->fft3uv,2*NXYZ,double complex); // (u,v)
    cppmallocl(mdfft->fft3rc,NXYZ,double); 
    
    
    // create plans
    // FFTW_ESTIMATE or FFTW_MEASURE
#define USE_FFTW_PLANNER FFTW_ESTIMATE
    mdfft->plan_f = fftw_plan_dft_3d(NX, NY, NZ, mdfft->fft3, mdfft->fft3, FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b = fftw_plan_dft_3d(NX, NY, NZ, mdfft->fft3, mdfft->fft3, FFTW_BACKWARD, USE_FFTW_PLANNER);  
    
    // FFTW many plans
    int fftwn[3] = {NX,NY,NZ};
    int idist, odist, istride, ostride;
    idist = odist = NXYZ; /* the distance in memory between the first element of the first array and the first element of the second array */
    istride = ostride = 1; /* array is contiguous in memory */
    int *inembed = fftwn, *onembed = fftwn;
//     mdfft->plan_f_many = fftw_plan_many_dft(2, fftwn, mdfft->batch,
//                                 mdfft->fft3many, inembed, istride, idist,
//                                 mdfft->fft3many, onembed, ostride, odist,
//                                 FFTW_FORWARD, USE_FFTW_PLANNER);
//     mdfft->plan_b_many = fftw_plan_many_dft(2, fftwn, mdfft->batch,
//                                 mdfft->fft3many, inembed, istride, idist,
//                                 mdfft->fft3many, onembed, ostride, odist,
//                                 FFTW_BACKWARD, USE_FFTW_PLANNER);
//     
    // gradients
    mdfft->plan_f_grad = fftw_plan_many_dft(3, fftwn, 2*3,
                                mdfft->fft3grad, inembed, istride, idist,
                                mdfft->fft3grad, onembed, ostride, odist,
                                FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b_grad = fftw_plan_many_dft(3, fftwn, 2*3,
                                mdfft->fft3grad, inembed, istride, idist,
                                mdfft->fft3grad, onembed, ostride, odist,
                                FFTW_BACKWARD, USE_FFTW_PLANNER);
//         
    // uv
    mdfft->plan_f_uv = fftw_plan_many_dft(3, fftwn, 2,
                                mdfft->fft3uv, inembed, istride, idist,
                                mdfft->fft3uv, onembed, ostride, odist,
                                FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b_uv = fftw_plan_many_dft(3, fftwn, 2,
                                mdfft->fft3uv, inembed, istride, idist,
                                mdfft->fft3uv, onembed, ostride, odist,
                                FFTW_BACKWARD, USE_FFTW_PLANNER);
    // r2c and c2r
    mdfft->plan_f_rc=fftw_plan_dft_r2c_3d(NX, NY, NZ, 
                               mdfft->fft3rc, mdfft->fft3,
                               USE_FFTW_PLANNER);
    mdfft->plan_b_cr=fftw_plan_dft_c2r_3d(NX, NY, NZ, 
                               mdfft->fft3, mdfft->fft3rc, 
                               USE_FFTW_PLANNER);
    return 0;
}


/**
 * Function clears memory
 * @param mdfft pointer to structure holding metadata for fft handling
 * */
int destroy_fft_plans(metadata_s3dpca_fft *mdfft)
{
    free(mdfft->fft3);
//     free(mdfft->fft3many);
    free(mdfft->fft3grad);
    free(mdfft->fft3uv);
    free(mdfft->fft3rc);
    return 0;
}

/**
 * Function computes laplace of real function laplace_f = d^2f/dx^2 + d^2f/dy^2 + d^2f/dz^2
 * @param f pointer to real function (INPUT)
 * @param laplace_f laplace of function, can be the same as f (OUTPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
int compute_laplace_real_f(double *f, double *laplace_f, metadata_s3dpca_fft *mdfft)
{
    int ixyz;
    int ix, iy, iz;
    double kx, ky, kz, k2;
    
    // copy data to working array
    for(ixyz=0; ixyz<NXYZ; ixyz++) mdfft->fft3rc[ixyz] = f[ixyz];
    
    fftw_execute(mdfft->plan_f_rc);
    
    // multiply by momentum
    ixyz=0;
    for(ix=0; ix<NX; ix++)
    {
        for(iy=0; iy<NY; iy++)
        {
            for(iz=0; iz<(NZ/2+1); iz++)
            {
                // extract momentum
                if(ix<NX/2) kx=2.*M_PI/(( double )NX * DX) * ( double )(ix   );
                else        kx=2.*M_PI/(( double )NX * DX) * ( double )(ix-NX);
                
                if(iy<NY/2) ky=2.*M_PI/(( double )NY * DY) * ( double )(iy   );
                else        ky=2.*M_PI/(( double )NY * DY) * ( double )(iy-NY);
                
                if(iz<NZ/2) kz=2.*M_PI/(( double )NZ * DZ) * ( double )(iz   );
                else        kz=2.*M_PI/(( double )NZ * DZ) * ( double )(iz-NZ);
                        
                k2 = -1.0*(kx*kx + ky*ky + kz*kz)/(NXYZ); // note: normalization factor is included
                                        
                mdfft->fft3[ixyz]*=k2;
                
                ixyz++;
            }
        }
    }
    
    fftw_execute(mdfft->plan_b_cr);
    
    // copy data to result table
    for(ixyz=0; ixyz<NXYZ; ixyz++) laplace_f[ixyz]=mdfft->fft3rc[ixyz];
 
    return 0;
}

int high_frequency_filter_d(double *in, double *out, double fd_mu, double fd_T, metadata_s3dpca_fft *mdfft)
{
    int ixyz;
    int ix, iy, iz;
    double kx, ky, kz, ek;
    
    // copy data to working array
    for(ixyz=0; ixyz<NXYZ; ixyz++) mdfft->fft3rc[ixyz] = in[ixyz];
    
    fftw_execute(mdfft->plan_f_rc);
    
    // multiply by momentum
    ixyz=0;
    for(ix=0; ix<NX; ix++)
    {
        for(iy=0; iy<NY; iy++)
        {
            for(iz=0; iz<(NZ/2+1); iz++)
            {
                // extract momentum
                if(ix<NX/2) kx=2.*M_PI/(( double )NX * DX) * ( double )(ix   );
                else        kx=2.*M_PI/(( double )NX * DX) * ( double )(ix-NX);
                
                if(iy<NY/2) ky=2.*M_PI/(( double )NY * DY) * ( double )(iy   );
                else        ky=2.*M_PI/(( double )NY * DY) * ( double )(iy-NY);
                
                if(iz<NZ/2) kz=2.*M_PI/(( double )NZ * DZ) * ( double )(iz   );
                else        kz=2.*M_PI/(( double )NZ * DZ) * ( double )(iz-NZ);
                        
                ek = 0.5*(kx*kx + ky*ky + kz*kz); 
                ek = 1.0/(exp((ek-fd_mu)/fd_T)+1.0) /NXYZ; // note: normalization factor is included
                                        
                mdfft->fft3[ixyz]*=ek;
                
                ixyz++;
            }
        }
    }
    
    fftw_execute(mdfft->plan_b_cr);
    
    // copy data to result table
    for(ixyz=0; ixyz<NXYZ; ixyz++) out[ixyz]=mdfft->fft3rc[ixyz];
    
    return 0;
}
