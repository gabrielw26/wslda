/**
 * W-SLDA Toolkit
 * 
 * Two kind of functions for the spectral interpolation that take into account real or complex numbers.
 * NOTE: All functions from this lib assume that spatial domain is fixed, while only lattice spacing changes!  
 * 
 * @see example: $WSLDA/extensions/interpolate-dataset.c
 * 
 * @author Wioleta Rzesa
 * @modified Gabriel Wlazlowski
 * */ 

#ifndef INTERPOLATION_H_ 
#define INTERPOLATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdbool.h>
#include <tgmath.h>
#include "wslda_errors.h"

/**********************************************************************/
/***************************** DEFINITIONS ****************************/
/**********************************************************************/

int wslda_interpolation_1dr(int nxi, double *funIn, int nxo, double *funOut);
int wslda_interpolation_1dc(int nxi, double complex *funIn, int nxo, double complex *funOut);

int wslda_interpolation_2dr(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut);
int wslda_interpolation_2dc(int nxi, int nyi, double complex *funIn, int nxo, int nyo, double complex *funOut);

int wslda_interpolation_3dr(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut);
int wslda_interpolation_3dc(int nxi, int nyi, int nzi, double complex *funIn, int nxo, int nyo, int nzo, double complex *funOut);  
 

/**********************************************************************/
/**********************************1D**********************************/
/**********************************************************************/
int wslda_interpolation_1dc(int nxi, double complex *funIn, int nxo, double complex *funOut){    

    double complex *out_forward, *in_backward;
    out_forward = (double complex *) malloc(sizeof(double complex)*nxi);
    in_backward = (double complex *) malloc(sizeof(double complex)*nxo);

    //set up fourier plans (creating plans for the direction of transform of a given array)
    fftw_plan plan_forward, plan_backward;
    plan_forward = fftw_plan_dft_1d(nxi, funIn, out_forward, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(nxo, in_backward, funOut, FFTW_BACKWARD, FFTW_ESTIMATE);

    //computing the forward fourier transform 
    fftw_execute(plan_forward);

    int ixyzi, ixi;
    int ixyzo, ixo;
    
    // reset
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) 
    {
        in_backward[ixyzo]=0.0+I*0.0;
        ixyzo++;
    }
    
    if (nxi <= nxo){//interpolation in 1 direction, zero padding 

        ixyzi=0;
        for (ixi = 0; ixi < nxi; ixi++)
        {
            if(ixi<nxi/2) ixo=ixi;
            else          ixo=ixi+(nxo-nxi);
            
            ixyzo = ixo;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzi++;
        }
    }
    else if(nxi >= nxo){//inverse interpolation in 1 direction
        
        ixyzo=0;
        for (ixo = 0; ixo < nxo; ixo++)
        {
            if(ixo<nxo/2) ixi=ixo;
            else          ixi=ixo+(nxi-nxo);
            
            ixyzi = ixi;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzo++;
        }
    }
    else{
        return WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED;
    }


    //computing the inverse fourier transform 
    fftw_execute(plan_backward);

    //normalisation
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) {

        funOut[ixyzo]/=(nxi);
        ixyzo++;
    }

    //memory deallocation 
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    free(in_backward);
    free(out_forward);
    
    return WSLDA_OK;
}

int wslda_interpolation_1dr(int nxi, double *funIn, int nxo, double *funOut){    

    double complex *out_forward, *in_backward;
    out_forward = (double complex *) malloc(sizeof(double complex)*(nxi/2+1));
    in_backward = (double complex *) malloc(sizeof(double complex)*(nxo/2+1));

    //set up fourier plans (creating plans for the direction of transform of a given array)
    fftw_plan plan_forward, plan_backward;
    plan_forward = fftw_plan_dft_r2c_1d(nxi, funIn, out_forward, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d(nxo, in_backward, funOut, FFTW_ESTIMATE);

    //computing the forward fourier transform 
    fftw_execute(plan_forward);

    int ixyzi, ixi;
    int ixyzo, ixo;
    
    // reset
    ixyzo=0;
    for (ixo = 0; ixo < (nxo/2+1); ixo++) 
    {
        in_backward[ixyzo]=0.0+I*0.0;
        ixyzo++;
    }
    
    if (nxi <= nxo){//interpolation in 1 direction, zero padding 

        ixyzi=0;
        for (ixi = 0; ixi < (nxi/2+1); ixi++)
        {
            if(ixi<nxi/2) ixo=ixi;
            else          ixo=ixi+(nxo-nxi);
            
            ixyzo = ixo;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzi++;
        }
    }
    else if(nxi >= nxo){//inverse interpolation in 1 direction
        
        ixyzo=0;
        for (ixo = 0; ixo < (nxo/2+1); ixo++)
        {
            if(ixo<nxo/2) ixi=ixo;
            else          ixi=ixo+(nxi-nxo);

            ixyzi = ixi;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzo++;
        }
    }
    else{
        return WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED;
    }

    //computing the inverse fourier transform 
    fftw_execute(plan_backward);

    //normalisation
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) {

        funOut[ixyzo]/=(nxi);
        ixyzo++;
    }

    //memory deallocation 
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    free(in_backward);
    free(out_forward);
    
    return WSLDA_OK;
}


/**********************************************************************/
/**********************************2D**********************************/
/**********************************************************************/
int wslda_interpolation_2dc(int nxi, int nyi, double complex *funIn, int nxo, int nyo, double complex *funOut){    

    double complex *out_forward, *in_backward;
    out_forward = (double complex *) malloc(sizeof(double complex)*nxi*nyi);
    in_backward = (double complex *) malloc(sizeof(double complex)*nxo*nyo);

    //set up fourier plans (creating plans for the direction of transform of a given array)
    fftw_plan plan_forward, plan_backward;
    plan_forward = fftw_plan_dft_2d(nxi, nyi, funIn, out_forward, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_2d(nxo, nyo, in_backward, funOut, FFTW_BACKWARD, FFTW_ESTIMATE);

    //computing the forward fourier transform 
    fftw_execute(plan_forward);

    int ixyzi, ixi, iyi;
    int ixyzo, ixo, iyo;
    
    // reset
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) 
    {
        in_backward[ixyzo]=0.0+I*0.0;
        ixyzo++;
    }
    
    if (nxi <= nxo && nyi <= nyo){//interpolation in 2 directions, zero padding 

        ixyzi=0;
        for (ixi = 0; ixi < nxi; ixi++) for (iyi = 0; iyi < nyi; iyi++)
        {
            if(ixi<nxi/2) ixo=ixi;
            else          ixo=ixi+(nxo-nxi);
            
            if(iyi<nyi/2) iyo=iyi;
            else          iyo=iyi+(nyo-nyi);
            
            ixyzo = iyo + nyo*ixo;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzi++;
        }
    }
    else if(nxi >= nxo && nyi >= nyo){//inverse interpolation in 2 directions
        
        ixyzo=0;
        for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++)
        {
            if(ixo<nxo/2) ixi=ixo;
            else          ixi=ixo+(nxi-nxo);
            
            if(iyo<nyo/2) iyi=iyo;
            else          iyi=iyo+(nyi-nyo);
            
            ixyzi = iyi + nyi*ixi;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzo++;
        }
    }
    else{
        return WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED;
    }


    //computing the inverse fourier transform 
    fftw_execute(plan_backward);

    //normalisation
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) {

        funOut[ixyzo]/=(nxi*nyi);
        ixyzo++;
    }

    //memory deallocation 
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    free(in_backward);
    free(out_forward);
    
    return WSLDA_OK;
}

int wslda_interpolation_2dr(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut){    

    double complex *out_forward, *in_backward;
    out_forward = (double complex *) malloc(sizeof(double complex)*nxi*(nyi/2+1));
    in_backward = (double complex *) malloc(sizeof(double complex)*nxo*(nyo/2+1));

    //set up fourier plans (creating plans for the direction of transform of a given array)
    fftw_plan plan_forward, plan_backward;
    plan_forward = fftw_plan_dft_r2c_2d(nxi, nyi, funIn, out_forward, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_2d(nxo, nyo, in_backward, funOut, FFTW_ESTIMATE);

    //computing the forward fourier transform 
    fftw_execute(plan_forward);

    int ixyzi, ixi, iyi;
    int ixyzo, ixo, iyo;
    
    // reset
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < (nyo/2+1); iyo++) 
    {
        in_backward[ixyzo]=0.0+I*0.0;
        ixyzo++;
    }
    
    if (nxi <= nxo && nyi <= nyo){//interpolation in 2 directions, zero padding 

        ixyzi=0;
        for (ixi = 0; ixi < nxi; ixi++) for (iyi = 0; iyi < (nyi/2+1); iyi++)
        {
            if(ixi<nxi/2) ixo=ixi;
            else          ixo=ixi+(nxo-nxi);
            
            if(iyi<nyi/2) iyo=iyi;
            else          iyo=iyi+(nyo-nyi);
            
            ixyzo = iyo + (nyo/2+1)*ixo;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzi++;
        }
    }
    else if(nxi >= nxo && nyi >= nyo){//inverse interpolation in 2 directions
        
        ixyzo=0;
        for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < (nyo/2+1); iyo++)
        {
            if(ixo<nxo/2) ixi=ixo;
            else          ixi=ixo+(nxi-nxo);
            
            if(iyo<nyo/2) iyi=iyo;
            else          iyi=iyo+(nyi-nyo);
            
            
            ixyzi = iyi + (nyi/2+1)*ixi;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzo++;
        }
    }
    else{
        return WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED;
    }

    //computing the inverse fourier transform 
    fftw_execute(plan_backward);

    //normalisation
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) {

        funOut[ixyzo]/=(nxi*nyi);
        ixyzo++;
    }

    //memory deallocation 
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    free(in_backward);
    free(out_forward);
    
    return WSLDA_OK;
}

/**********************************************************************/
/**********************************3D**********************************/
/**********************************************************************/
int wslda_interpolation_3dc(int nxi, int nyi, int nzi, double complex *funIn, int nxo, int nyo, int nzo, double complex *funOut){    

    double complex *out_forward, *in_backward;
    out_forward = (double complex *) malloc(sizeof(double complex)*nxi*nyi*nzi);
    in_backward = (double complex *) malloc(sizeof(double complex)*nxo*nyo*nzo);

    //set up fourier plans (creating plans for the direction of transform of a given array)
    fftw_plan plan_forward, plan_backward;
    plan_forward = fftw_plan_dft_3d(nxi, nyi, nzi, funIn, out_forward, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_3d(nxo, nyo, nzo, in_backward, funOut, FFTW_BACKWARD, FFTW_ESTIMATE);

    //computing the forward fourier transform 
    fftw_execute(plan_forward);


    int ixyzi, ixi, iyi, izi;
    int ixyzo, ixo, iyo, izo;
    
    // reset
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) for (izo = 0; izo < nzo; izo++) 
    {
        in_backward[ixyzo]=0.0+I*0.0;
        ixyzo++;
    }
    
    if (nxi <= nxo && nyi <= nyo && nzi <= nzo && nzi <= nzo){//interpolation in 3 directions, zero padding 

        ixyzi=0;
        for (ixi = 0; ixi < nxi; ixi++) for (iyi = 0; iyi < nyi; iyi++) for (izi = 0; izi < nzi; izi++)
        {
            if(ixi<nxi/2) ixo=ixi;
            else          ixo=ixi+(nxo-nxi);
            
            if(iyi<nyi/2) iyo=iyi;
            else          iyo=iyi+(nyo-nyi);
            
            if(izi<nzi/2) izo=izi;
            else          izo=izi+(nzo-nzi);
            
            ixyzo = izo + nzo*iyo + nzo*nyo*ixo;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzi++;
        }
    }
    else if(nxi >= nxo && nyi >= nyo && nzi >= nzo){//inverse interpolation in 3 directions
        
        ixyzo=0;
        for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) for (izo = 0; izo < nzo; izo++)
        {
            if(ixo<nxo/2) ixi=ixo;
            else          ixi=ixo+(nxi-nxo);
            
            if(iyo<nyo/2) iyi=iyo;
            else          iyi=iyo+(nyi-nyo);
            
            if(izo<nzo/2) izi=izo;
            else          izi=izo+(nzi-nzo);
            
            ixyzi = izi + nzi*iyi + nzi*nyi*ixi;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzo++;
        }
    }
    else{
        return WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED;
    }


    //computing the inverse fourier transform 
    fftw_execute(plan_backward);

    //normalisation
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) for (izo = 0; izo < nzo; izo++) {

        funOut[ixyzo]/=(nxi*nyi*nzi);
        ixyzo++;
    }

    //memory deallocation 
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);

    free(in_backward);
    free(out_forward);
    
    return WSLDA_OK;
}

int wslda_interpolation_3dr(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut){    

    double complex *out_forward, *in_backward;
    out_forward = (double complex *) malloc(sizeof(double complex)*nxi*nyi*(nzi/2+1));
    in_backward = (double complex *) malloc(sizeof(double complex)*nxo*nyo*(nzo/2+1));

    //set up fourier plans (creating plans for the direction of transform of a given array)
    fftw_plan plan_forward, plan_backward;
    plan_forward = fftw_plan_dft_r2c_3d(nxi, nyi, nzi, funIn, out_forward, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_3d(nxo, nyo, nzo, in_backward, funOut, FFTW_ESTIMATE);

    //computing the forward fourier transform 
    fftw_execute(plan_forward);


    int ixyzi, ixi, iyi, izi;
    int ixyzo, ixo, iyo, izo;
    
    // reset
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) for (izo = 0; izo < (nzo/2+1); izo++) 
    {
        in_backward[ixyzo]=0.0+I*0.0;
        ixyzo++;
    }
    
    if (nxi <= nxo && nyi <= nyo && nzi <= nzo && nzi <= nzo){//interpolation in 3 directions, zero padding 

        ixyzi=0;
        for (ixi = 0; ixi < nxi; ixi++) for (iyi = 0; iyi < nyi; iyi++) for (izi = 0; izi < (nzi/2+1); izi++)
        {
            if(ixi<nxi/2) ixo=ixi;
            else          ixo=ixi+(nxo-nxi);
            
            if(iyi<nyi/2) iyo=iyi;
            else          iyo=iyi+(nyo-nyi);
            
            if(izi<nzi/2) izo=izi;
            else          izo=izi+(nzo-nzi);
            
            ixyzo = izo + (nzo/2+1)*iyo + (nzo/2+1)*nyo*ixo;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzi++;
        }
    }
    else if(nxi >= nxo && nyi >= nyo && nzi >= nzo){//inverse interpolation in 3 directions
        
        ixyzo=0;
        for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) for (izo = 0; izo < (nzo/2+1); izo++)
        {
            if(ixo<nxo/2) ixi=ixo;
            else          ixi=ixo+(nxi-nxo);
            
            if(iyo<nyo/2) iyi=iyo;
            else          iyi=iyo+(nyi-nyo);
            
            if(izo<nzo/2) izi=izo;
            else          izi=izo+(nzi-nzo);
            
            ixyzi = izi + (nzi/2+1)*iyi + (nzi/2+1)*nyi*ixi;
            
            in_backward[ixyzo] =  out_forward[ixyzi];
            
            ixyzo++;
        }
    }
    else{
        return WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED;
    }


    //computing the inverse fourier transform 
    fftw_execute(plan_backward);

    //normalisation
    ixyzo=0;
    for (ixo = 0; ixo < nxo; ixo++) for (iyo = 0; iyo < nyo; iyo++) for (izo = 0; izo < nzo; izo++) {

        funOut[ixyzo]/=(nxi*nyi*nzi);
        ixyzo++;
    }

    //memory deallocation 
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    
    free(in_backward);
    free(out_forward);
    
    return WSLDA_OK;
}


#endif
