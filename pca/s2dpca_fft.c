// Author: Gabriel Wlazlowski
// File that provides support for fftw operations

#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include "pca_settings.h"
#include "pca_macro.h"
#include "s2dpca_fft.h"

/**
 * Function creates plans and allocates memory
 * @param mdfft pointer to structure holding metadata for fft handling
 * @param batch number of vectors transformed by plan many  
 * */
int create_fft_plans(metadata_s2dpca_fft *mdfft, int batch)
{
    mdfft->batch=batch;
    
    // allocate memory for buffers
    cppmallocl(mdfft->fft2,NX*NY,double complex); 
    cppmallocl(mdfft->fft2many,batch*NX*NY,double complex);
    cppmallocl(mdfft->fft2grad,2*2*NX*NY,double complex); // (u,v) * (dx,dy)
    cppmallocl(mdfft->fft2uv,2*NX*NY,double complex); // (u,v)
    cppmallocl(mdfft->fft2rc,NX*NY,double); 
    
    // create plans
    // FFTW_ESTIMATE or FFTW_MEASURE
#define USE_FFTW_PLANNER FFTW_ESTIMATE
    mdfft->plan_f = fftw_plan_dft_2d(NX, NY, mdfft->fft2, mdfft->fft2, FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b = fftw_plan_dft_2d(NX, NY, mdfft->fft2, mdfft->fft2, FFTW_BACKWARD, USE_FFTW_PLANNER);  
    
    // FFTW many plans
    int fftwn[2] = {NX,NY};
    int idist, odist, istride, ostride;
    idist = odist = NX*NY; /* the distance in memory between the first element of the first array and the first element of the second array */
    istride = ostride = 1; /* array is contiguous in memory */
    int *inembed = fftwn, *onembed = fftwn;
    mdfft->plan_f_many = fftw_plan_many_dft(2, fftwn, mdfft->batch,
                                mdfft->fft2many, inembed, istride, idist,
                                mdfft->fft2many, onembed, ostride, odist,
                                FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b_many = fftw_plan_many_dft(2, fftwn, mdfft->batch,
                                mdfft->fft2many, inembed, istride, idist,
                                mdfft->fft2many, onembed, ostride, odist,
                                FFTW_BACKWARD, USE_FFTW_PLANNER);
    
    // gradients
    mdfft->plan_f_grad = fftw_plan_many_dft(2, fftwn, 2*2,
                                mdfft->fft2grad, inembed, istride, idist,
                                mdfft->fft2grad, onembed, ostride, odist,
                                FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b_grad = fftw_plan_many_dft(2, fftwn, 2*2,
                                mdfft->fft2grad, inembed, istride, idist,
                                mdfft->fft2grad, onembed, ostride, odist,
                                FFTW_BACKWARD, USE_FFTW_PLANNER);
        
    // uv
    mdfft->plan_f_uv = fftw_plan_many_dft(2, fftwn, 2,
                                mdfft->fft2uv, inembed, istride, idist,
                                mdfft->fft2uv, onembed, ostride, odist,
                                FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b_uv = fftw_plan_many_dft(2, fftwn, 2,
                                mdfft->fft2uv, inembed, istride, idist,
                                mdfft->fft2uv, onembed, ostride, odist,
                                FFTW_BACKWARD, USE_FFTW_PLANNER);
    // r2c and c2r
    mdfft->plan_f_rc=fftw_plan_dft_r2c_2d(NX, NY,
                               mdfft->fft2rc, mdfft->fft2,
                               USE_FFTW_PLANNER);
    mdfft->plan_b_cr=fftw_plan_dft_c2r_2d(NX, NY,
                               mdfft->fft2, mdfft->fft2rc, 
                               USE_FFTW_PLANNER);
    return 0;
}


/**
 * Function clears memory
 * @param mdfft pointer to structure holding metadata for fft handling
 * */
int destroy_fft_plans(metadata_s2dpca_fft *mdfft)
{
    free(mdfft->fft2);
    free(mdfft->fft2many);
    free(mdfft->fft2grad);
    free(mdfft->fft2uv);
    free(mdfft->fft2rc);
    return 0;
}

/**
 * Function computes laplace of real function laplace_f = d^2f/dx^2 + d^2f/dy^2 + d^2f/dz^2
 * @param f pointer to real function (INPUT)
 * @param laplace_f laplace of function, can be the same as f (OUTPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
int compute_laplace_real_f(double *f, double *laplace_f, metadata_s2dpca_fft *mdfft)
{
    int ixyz;
    int ix, iy;
    double kx, ky, k2;
    
    // copy data to working array
    for(ixyz=0; ixyz<NX*NY; ixyz++) mdfft->fft2rc[ixyz] = f[ixyz];
    
    fftw_execute(mdfft->plan_f_rc);
    
    // multiply by momentum
    ixyz=0;
    for(ix=0; ix<NX; ix++)
    {
        for(iy=0; iy<(NY/2+1); iy++)
        {
            // extract momentum
            if(ix<NX/2) kx=2.*M_PI/(( double )NX*DX) * ( double )(ix   );
            else        kx=2.*M_PI/(( double )NX*DX) * ( double )(ix-NX);
            
            if(iy<NY/2) ky=2.*M_PI/(( double )NY*DY) * ( double )(iy   );
            else        ky=2.*M_PI/(( double )NY*DY) * ( double )(iy-NY);
                    
            k2 = -1.0*(kx*kx + ky*ky)/(NX*NY); // note: normalization factor is included
                                    
            mdfft->fft2[ixyz]*=k2;
            
            ixyz++;
        }
    }
    
    fftw_execute(mdfft->plan_b_cr);
    
    // copy data to result table
    for(ixyz=0; ixyz<NX*NY; ixyz++) laplace_f[ixyz]=mdfft->fft2rc[ixyz];
 
    return 0;
}

