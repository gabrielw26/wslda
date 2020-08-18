// Author: Gabriel Wlazlowski
// File that provides support for fftw operations

// Structure holding metadata

#ifndef __KZPCA_FFT__
#define __KZPCA_FFT__

#include <math.h>
#include <complex.h>
#include <fftw3.h>

typedef struct
{
    double complex *fft2;
//     double complex *fft2many;
    double complex *fft2grad; // gradients for computing densities
    double complex *fft2uv; // for transforming (u,v) simultaniesly 
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
    double *fft2rc; // for rc and cr plans
    
    int batch;
    
} metadata_s2dpca_fft;

int create_fft_plans(metadata_s2dpca_fft *mdfft, int batch);
int destroy_fft_plans(metadata_s2dpca_fft *mdfft);
int compute_laplace_real_f(double *f, double *laplace_f, metadata_s2dpca_fft *mdfft);

#endif
