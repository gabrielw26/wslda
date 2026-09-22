// Author: Gabriel Wlazlowski
// File that provides support for fftw operations

#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include "pca_settings.h"
#include "pca_macro.h"
#include "s1dpca_fft.h"

/**
 * Function creates plans and allocates memory
 * @param mdfft pointer to structure holding metadata for fft handling
 * @param batch number of vectors transformed by plan many  
 * */
int create_fft_plans(metadata_s1dpca_fft *mdfft, int batch)
{
    mdfft->batch=batch;
    
    // allocate memory for buffers
    cppmallocl(mdfft->fft1,NX,double complex); 
    cppmallocl(mdfft->fft1grad,2*1*NX,double complex); // (u,v) * (dx)
    cppmallocl(mdfft->fft1uv,2*NX,double complex); // (u,v)
    cppmallocl(mdfft->fft1rc,NX,double); 
    
    // create plans
    // FFTW_ESTIMATE or FFTW_MEASURE
#define USE_FFTW_PLANNER FFTW_ESTIMATE
    mdfft->plan_f = fftw_plan_dft_1d(NX, mdfft->fft1, mdfft->fft1, FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b = fftw_plan_dft_1d(NX, mdfft->fft1, mdfft->fft1, FFTW_BACKWARD, USE_FFTW_PLANNER);  
    
    // FFTW many plans
    int fftwn[1] = {NX};
    int idist, odist, istride, ostride;
    idist = odist = NX; /* the distance in memory between the first element of the first array and the first element of the second array */
    istride = ostride = 1; /* array is contiguous in memory */
    int *inembed = fftwn, *onembed = fftwn;
    
    // gradients
    mdfft->plan_f_grad = fftw_plan_many_dft(1, fftwn, 2*1,
                                mdfft->fft1grad, inembed, istride, idist,
                                mdfft->fft1grad, onembed, ostride, odist,
                                FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b_grad = fftw_plan_many_dft(1, fftwn, 2*1,
                                mdfft->fft1grad, inembed, istride, idist,
                                mdfft->fft1grad, onembed, ostride, odist,
                                FFTW_BACKWARD, USE_FFTW_PLANNER);
        
    // uv
    mdfft->plan_f_uv = fftw_plan_many_dft(1, fftwn, 2,
                                mdfft->fft1uv, inembed, istride, idist,
                                mdfft->fft1uv, onembed, ostride, odist,
                                FFTW_FORWARD, USE_FFTW_PLANNER);
    mdfft->plan_b_uv = fftw_plan_many_dft(1, fftwn, 2,
                                mdfft->fft1uv, inembed, istride, idist,
                                mdfft->fft1uv, onembed, ostride, odist,
                                FFTW_BACKWARD, USE_FFTW_PLANNER);
    // r2c and c2r
    mdfft->plan_f_rc=fftw_plan_dft_r2c_1d(NX,
                               mdfft->fft1rc, mdfft->fft1,
                               USE_FFTW_PLANNER);
    mdfft->plan_b_cr=fftw_plan_dft_c2r_1d(NX,
                               mdfft->fft1, mdfft->fft1rc, 
                               USE_FFTW_PLANNER);
    return 0;
}


/**
 * Function clears memory
 * @param mdfft pointer to structure holding metadata for fft handling
 * */
int destroy_fft_plans(metadata_s1dpca_fft *mdfft)
{
    free(mdfft->fft1);
    free(mdfft->fft1grad);
    free(mdfft->fft1uv);
    free(mdfft->fft1rc);
    return 0;
}

/**
 * Function computes laplace of real function laplace_f = d^2f/dx^2 + d^2f/dy^2 + d^2f/dz^2
 * @param f pointer to real function (INPUT)
 * @param laplace_f laplace of function, can be the same as f (OUTPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
int compute_laplace_real_f(double *f, double *laplace_f, metadata_s1dpca_fft *mdfft)
{
    int ixyz;
    int ix;
    double kx, k2;
    
    // copy data to working array
    for(ixyz=0; ixyz<NX; ixyz++) mdfft->fft1rc[ixyz] = f[ixyz];
    
    fftw_execute(mdfft->plan_f_rc);
    
    // multiply by momentum
    ixyz=0;
    for(ix=0; ix<NX; ix++)
    {
            if(ix<NX/2) kx=2.*M_PI/(( double )NX*DX) * ( double )(ix   );
            else        kx=2.*M_PI/(( double )NX*DX) * ( double )(ix-NX);
            
                    
            k2 = -1.0*(kx*kx)/(NX); // note: normalization factor is included
                                    
            mdfft->fft1[ixyz]*=k2;
            
            ixyz++;
    }
    
    fftw_execute(mdfft->plan_b_cr);
    
    // copy data to result table
    for(ixyz=0; ixyz<NX; ixyz++) laplace_f[ixyz]=mdfft->fft1rc[ixyz];
 
    return 0;
}

/**
 * The function removes high frequencies from a real signal.
 * The spectral filter is FD(mu,T)=1/(exp((k^2/2-mu)/T)+1).
 * @param in signal to be filtered (INPUT)
 * @param out filtered signal; may be the same as in (OUTPUT)
 * @param fd_mu chemical-potential parameter of the filter
 * @param fd_T temperature parameter of the filter
 * @param mdfft FFT metadata and workspace
 * @return 0
 * */
int high_frequency_filter_d(double *in, double *out, double fd_mu, double fd_T, metadata_s1dpca_fft *mdfft)
{
    int ix;
    double kx, ek;

    // copy data to working array
    for(ix=0; ix<NX; ix++) mdfft->fft1rc[ix] = in[ix];

    fftw_execute(mdfft->plan_f_rc);

    // multiply by the filter function in momentum space
    for(ix=0; ix<(NX/2+1); ix++)
    {
        if(ix<NX/2) kx=2.*M_PI/((double)NX*DX) * (double)(ix   );
        else        kx=2.*M_PI/((double)NX*DX) * (double)(ix-NX);

        ek = 0.5*kx*kx;
        ek = 1.0/(exp((ek-fd_mu)/fd_T)+1.0) /NX; // inverse FFT normalization included

        mdfft->fft1[ix]*=ek;
    }

    fftw_execute(mdfft->plan_b_cr);

    // copy data to result table
    for(ix=0; ix<NX; ix++) out[ix]=mdfft->fft1rc[ix];

    return 0;
}

/**
 * The function removes high frequencies from a complex signal.
 * The spectral filter is FD(mu,T)=1/(exp((k^2/2-mu)/T)+1).
 * @param in signal to be filtered (INPUT)
 * @param out filtered signal; may be the same as in (OUTPUT)
 * @param fd_mu chemical-potential parameter of the filter
 * @param fd_T temperature parameter of the filter
 * @param mdfft FFT metadata and workspace
 * @return 0
 * */
int high_frequency_filter_c(double complex *in, double complex *out, double fd_mu, double fd_T, metadata_s1dpca_fft *mdfft)
{
    int ix;
    double kx, ek;

    // copy data to working array
    for(ix=0; ix<NX; ix++) mdfft->fft1[ix] = in[ix];

    fftw_execute(mdfft->plan_f);

    // multiply by the filter function in momentum space
    for(ix=0; ix<NX; ix++)
    {
        if(ix<NX/2) kx=2.*M_PI/((double)NX*DX) * (double)(ix   );
        else        kx=2.*M_PI/((double)NX*DX) * (double)(ix-NX);

        ek = 0.5*kx*kx;
        ek = 1.0/(exp((ek-fd_mu)/fd_T)+1.0) /NX; // inverse FFT normalization included

        mdfft->fft1[ix]*=ek;
    }

    fftw_execute(mdfft->plan_b);

    // copy data to result table
    for(ix=0; ix<NX; ix++) out[ix]=mdfft->fft1[ix];

    return 0;
}
