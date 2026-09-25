// Author: Gabriel Wlazlowski
// File that provides support for fftw operations

// Structure holding metadata

#ifndef __WSLDA_FFT1D__
#define __WSLDA_FFT1D__

#include <math.h>
#include <complex.h>
#include <fftw3.h>

typedef struct
{
    double complex *fft1;
//     double complex *fft1many;
    double complex *fft1grad; // gradients for computing densities
    double complex *fft1uv; // for transforming (u,v) simultaniesly 
    fftw_plan plan_f;
    fftw_plan plan_b;
//     fftw_plan plan_f_many;
//     fftw_plan plan_b_many;
    fftw_plan plan_f_grad;
    fftw_plan plan_b_grad;
    fftw_plan plan_f_uv;
    fftw_plan plan_b_uv;
    
    fftw_plan plan_f_rc;
    fftw_plan plan_b_cr;
    double *fft1rc; // for rc and cr plans
    
    int batch;
    
} metadata_s1dpca_fft;

int create_fft_plans(metadata_s1dpca_fft *mdfft, int batch);
int destroy_fft_plans(metadata_s1dpca_fft *mdfft);
int compute_laplace_real_f(double *f, double *laplace_f, metadata_s1dpca_fft *mdfft);
int high_frequency_filter_d(double *in, double *out, double fd_mu, double fd_T, metadata_s1dpca_fft *mdfft);
int high_frequency_filter_c(double complex *in, double complex *out, double fd_mu, double fd_T, metadata_s1dpca_fft *mdfft);

#endif
