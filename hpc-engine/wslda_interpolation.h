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
int wslda_interpolation_1dv(int nxi, double *funIn, int nxo, double *funOut);
int wslda_interpolation_1dc(int nxi, double complex *funIn, int nxo, double complex *funOut);
int wslda_interpolation_1d(char type, int nxi, void *funIn, int nxo, void *funOut);

int wslda_interpolation_2dr(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut);
int wslda_interpolation_2dv(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut);
int wslda_interpolation_2dc(int nxi, int nyi, double complex *funIn, int nxo, int nyo, double complex *funOut);
int wslda_interpolation_2d(char type, int nxi, int nyi, void *funIn, int nxo, int nyo, void *funOut);

int wslda_interpolation_3dr(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut);
int wslda_interpolation_3dv(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut);
int wslda_interpolation_3dc(int nxi, int nyi, int nzi, double complex *funIn, int nxo, int nyo, int nzo, double complex *funOut);  
int wslda_interpolation_3d(char type, int nxi, int nyi, int nzi, void *funIn, int nxo, int nyo, int nzo, void *funOut);
 

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
            
            if(ixo>nxi/2) // special case - utilize symmetry
            {
                ixyzo = (ixo-(nxo-nxi));
                in_backward[ixyzo] =  conj(out_forward[ixyzi]);
            }
            else
            {
                ixyzo = ixo;
                in_backward[ixyzo] =  out_forward[ixyzi];
            }
            
            if(ixyzi>=(nxi/2+1)) return WSLDA_ERR_INTRISTIC_ERROR;
            if(ixyzo>=(nxo/2+1)) return WSLDA_ERR_INTRISTIC_ERROR;
            
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

int wslda_interpolation_1dv(int nxi, double *funIn, int nxo, double *funOut){ 
    wslda_interpolation_1dr(nxi, funIn+0*nxi, nxo, funOut+0*nxo);
    wslda_interpolation_1dr(nxi, funIn+1*nxi, nxo, funOut+1*nxo);
    wslda_interpolation_1dr(nxi, funIn+2*nxi, nxo, funOut+2*nxo);
    return WSLDA_OK;
} 

int wslda_interpolation_1d(char type, int nxi, void *funIn, int nxo, void *funOut){
    if(type=='c')      return wslda_interpolation_1dc(nxi, (double complex *)funIn, nxo, (double complex *)funOut);
    else if(type=='v') return wslda_interpolation_1dv(nxi, (double *)funIn, nxo, (double *)funOut);
    else               return wslda_interpolation_1dr(nxi, (double *)funIn, nxo, (double *)funOut);
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
            
            if(iyo>nyi/2) // special case - utilize symmetry
            {
                ixyzo = (iyo-(nyo-nyi)) + (nyo/2+1)*ixo;
                in_backward[ixyzo] =  conj(out_forward[ixyzi]);
            }
            else
            {
                ixyzo = iyo + (nyo/2+1)*ixo;
                in_backward[ixyzo] =  out_forward[ixyzi];
            }
            
            if(ixyzi>=nxi*(nyi/2+1)) return WSLDA_ERR_INTRISTIC_ERROR;
            if(ixyzo>=nxo*(nyo/2+1)) return WSLDA_ERR_INTRISTIC_ERROR;
            
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

int wslda_interpolation_2dv(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut){ 
    wslda_interpolation_2dr(nxi, nyi, funIn+0*nxi*nyi, nxo, nyo, funOut+0*nxo*nyo);
    wslda_interpolation_2dr(nxi, nyi, funIn+1*nxi*nyi, nxo, nyo, funOut+1*nxo*nyo);
    wslda_interpolation_2dr(nxi, nyi, funIn+2*nxi*nyi, nxo, nyo, funOut+2*nxo*nyo);
    return WSLDA_OK;
} 

int wslda_interpolation_2d(char type, int nxi, int nyi, void *funIn, int nxo, int nyo, void *funOut){
    if(type=='c')      return wslda_interpolation_2dc(nxi, nyi, (double complex *)funIn, nxo, nyo, (double complex *)funOut);
    else if(type=='v') return wslda_interpolation_2dv(nxi, nyi, (double *)funIn, nxo, nyo, (double *)funOut);
    else               return wslda_interpolation_2dr(nxi, nyi, (double *)funIn, nxo, nyo, (double *)funOut);
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
            
            if(izo>nzi/2) // special case - utilize symmetry
            {
                ixyzo = (izo-(nzo-nzi)) + (nzo/2+1)*iyo + (nzo/2+1)*nyo*ixo;
                in_backward[ixyzo] =  conj(out_forward[ixyzi]);
//                 printf("SC: %d --> %d\n", izo, (izo-(nzo-nzi)));
            }
            else
            {
                ixyzo = izo + (nzo/2+1)*iyo + (nzo/2+1)*nyo*ixo;
                in_backward[ixyzo] =  out_forward[ixyzi];
            }
            
            if(ixyzi>=nxi*nyi*(nzi/2+1)) return WSLDA_ERR_INTRISTIC_ERROR;
            if(ixyzo>=nxo*nyo*(nzo/2+1)) return WSLDA_ERR_INTRISTIC_ERROR;
            
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

int wslda_interpolation_3dv(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut){ 
    wslda_interpolation_3dr(nxi, nyi, nzi, funIn+0*nxi*nyi*nzi, nxo, nyo, nzo, funOut+0*nxo*nyo*nzo);
    wslda_interpolation_3dr(nxi, nyi, nzi, funIn+1*nxi*nyi*nzi, nxo, nyo, nzo, funOut+1*nxo*nyo*nzo);
    wslda_interpolation_3dr(nxi, nyi, nzi, funIn+2*nxi*nyi*nzi, nxo, nyo, nzo, funOut+2*nxo*nyo*nzo);
    return WSLDA_OK;
} 

int wslda_interpolation_3d(char type, int nxi, int nyi, int nzi, void *funIn, int nxo, int nyo, int nzo, void *funOut){
    if(type=='c')      return wslda_interpolation_3dc(nxi, nyi, nzi, (double complex *)funIn, nxo, nyo, nzo, (double complex *)funOut);
    else if(type=='v') return wslda_interpolation_3dv(nxi, nyi, nzi, (double *)funIn, nxo, nyo, nzo, (double *)funOut);
    else               return wslda_interpolation_3dr(nxi, nyi, nzi, (double *)funIn, nxo, nyo, nzo, (double *)funOut);
} 

#endif
