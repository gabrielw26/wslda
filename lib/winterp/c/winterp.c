/**
 * W-interp Library 
 * @see https://gitlab.fizyka.pw.edu.pl/wtools/winterp
 * */

#include "winterp.h"

#define WINTERP_USE_FFTW_PLANNER FFTW_ESTIMATE

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "ERROR: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return WINTERP_ERR_CANNOT_ALLOCATE_MEMORY;                              \
    }
    
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**********************************************************************/
/**********************************1D**********************************/
/**********************************************************************/
int winterp_interpolation_1dc(int nxi, double complex *funIn, int nxo, double complex *funOut){    

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
        return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED;
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
    
    return WINTERP_OK;
}

int winterp_interpolation_1dr(int nxi, double *funIn, int nxo, double *funOut){    

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
            
            if(ixyzi>=(nxi/2+1)) return WINTERP_ERR_INTRISTIC_ERROR;
            if(ixyzo>=(nxo/2+1)) return WINTERP_ERR_INTRISTIC_ERROR;
            
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
        return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED;
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
    
    return WINTERP_OK;
}

int winterp_interpolation_1dv(int nxi, double *funIn, int nxo, double *funOut){ 
    winterp_interpolation_1dr(nxi, funIn+0*nxi, nxo, funOut+0*nxo);
    winterp_interpolation_1dr(nxi, funIn+1*nxi, nxo, funOut+1*nxo);
    winterp_interpolation_1dr(nxi, funIn+2*nxi, nxo, funOut+2*nxo);
    return WINTERP_OK;
} 

int winterp_interpolation_1d(char type, int nxi, void *funIn, int nxo, void *funOut){
    if(type=='c')      return winterp_interpolation_1dc(nxi, (double complex *)funIn, nxo, (double complex *)funOut);
    else if(type=='v') return winterp_interpolation_1dv(nxi, (double *)funIn, nxo, (double *)funOut);
    else               return winterp_interpolation_1dr(nxi, (double *)funIn, nxo, (double *)funOut);
}    

/**********************************************************************/
/**********************************2D**********************************/
/**********************************************************************/
int winterp_interpolation_2dc(int nxi, int nyi, double complex *funIn, int nxo, int nyo, double complex *funOut){    

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
        return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED;
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
    
    return WINTERP_OK;
}

int winterp_interpolation_2dr(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut){    

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
            
            if(ixyzi>=nxi*(nyi/2+1)) return WINTERP_ERR_INTRISTIC_ERROR;
            if(ixyzo>=nxo*(nyo/2+1)) return WINTERP_ERR_INTRISTIC_ERROR;
            
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
        return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED;
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
    
    return WINTERP_OK;
}

int winterp_interpolation_2dv(int nxi, int nyi, double *funIn, int nxo, int nyo, double *funOut){ 
    winterp_interpolation_2dr(nxi, nyi, funIn+0*nxi*nyi, nxo, nyo, funOut+0*nxo*nyo);
    winterp_interpolation_2dr(nxi, nyi, funIn+1*nxi*nyi, nxo, nyo, funOut+1*nxo*nyo);
    winterp_interpolation_2dr(nxi, nyi, funIn+2*nxi*nyi, nxo, nyo, funOut+2*nxo*nyo);
    return WINTERP_OK;
} 

int winterp_interpolation_2d(char type, int nxi, int nyi, void *funIn, int nxo, int nyo, void *funOut){
    if(type=='c')      return winterp_interpolation_2dc(nxi, nyi, (double complex *)funIn, nxo, nyo, (double complex *)funOut);
    else if(type=='v') return winterp_interpolation_2dv(nxi, nyi, (double *)funIn, nxo, nyo, (double *)funOut);
    else               return winterp_interpolation_2dr(nxi, nyi, (double *)funIn, nxo, nyo, (double *)funOut);
}    


/**********************************************************************/
/**********************************3D**********************************/
/**********************************************************************/
int winterp_interpolation_3dc(int nxi, int nyi, int nzi, double complex *funIn, int nxo, int nyo, int nzo, double complex *funOut){    

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
        return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED;
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
    
    return WINTERP_OK;
}

int winterp_interpolation_3dr(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut){    
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
//                 wprintf("SC: %d --> %d\n", izo, (izo-(nzo-nzi)));
            }
            else
            {
                ixyzo = izo + (nzo/2+1)*iyo + (nzo/2+1)*nyo*ixo;
                in_backward[ixyzo] =  out_forward[ixyzi];
            }
            
            if(ixyzi>=nxi*nyi*(nzi/2+1)) return WINTERP_ERR_INTRISTIC_ERROR;
            if(ixyzo>=nxo*nyo*(nzo/2+1)) return WINTERP_ERR_INTRISTIC_ERROR;
            
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
        return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED;
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
    
    return WINTERP_OK;
}

int winterp_interpolation_3dv(int nxi, int nyi, int nzi, double *funIn, int nxo, int nyo, int nzo, double *funOut){ 
    winterp_interpolation_3dr(nxi, nyi, nzi, funIn+0*nxi*nyi*nzi, nxo, nyo, nzo, funOut+0*nxo*nyo*nzo);
    winterp_interpolation_3dr(nxi, nyi, nzi, funIn+1*nxi*nyi*nzi, nxo, nyo, nzo, funOut+1*nxo*nyo*nzo);
    winterp_interpolation_3dr(nxi, nyi, nzi, funIn+2*nxi*nyi*nzi, nxo, nyo, nzo, funOut+2*nxo*nyo*nzo);
    return WINTERP_OK;
} 

int winterp_interpolation_3d(char type, int nxi, int nyi, int nzi, void *funIn, int nxo, int nyo, int nzo, void *funOut){
    if(type=='c')      return winterp_interpolation_3dc(nxi, nyi, nzi, (double complex *)funIn, nxo, nyo, nzo, (double complex *)funOut);
    else if(type=='v') return winterp_interpolation_3dv(nxi, nyi, nzi, (double *)funIn, nxo, nyo, nzo, (double *)funOut);
    else               return winterp_interpolation_3dr(nxi, nyi, nzi, (double *)funIn, nxo, nyo, nzo, (double *)funOut);
} 


int winterp_create_interpolator(char datatype, int datadim, int *dims, double *spacings, void *data, winterp_interpolator *interp)
{
    interp->datadim=datadim;
    interp->datatype=datatype;
    if(datatype!='r' && datatype!='c' && datatype!='v') return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED;
    
    fftw_plan plan_f;
    int i,j,ik=1;
    size_t shift, shiftr;
    
    if(datadim==1)
    {
        interp->nx=dims[0]; interp->dx=spacings[0];
        interp->ny=1;       interp->dy=1.0;
        interp->nz=1;       interp->dz=1.0;
        
        shift=interp->nx*interp->ny*interp->nz;
        cppmallocl(interp->datak,shift,double complex);
        
        if(datatype=='c')
        {
            plan_f = fftw_plan_dft_1d(interp->nx, (double complex *)data, interp->datak, FFTW_FORWARD, WINTERP_USE_FFTW_PLANNER);
            fftw_execute(plan_f);
            fftw_destroy_plan(plan_f);
        }
        else // `r` or `v`
        {
            if(datatype=='v') ik=3;
            for(i=0; i<ik; i++)
            {
                for(j=0; j<shift; j++) *(interp->datak+i*shift+j) = *((double *)data+i*shift+j) + I*0.0;
                plan_f = fftw_plan_dft_1d(interp->nx, interp->datak+i*shift, interp->datak+i*shift, FFTW_FORWARD, WINTERP_USE_FFTW_PLANNER);
                fftw_execute(plan_f);
                fftw_destroy_plan(plan_f);
            }
        }
    }
    else if(datadim==2)
    {
        interp->nx=dims[0]; interp->dx=spacings[0];
        interp->ny=dims[1]; interp->dy=spacings[1];
        interp->nz=1;       interp->dz=1.0;

        shift=interp->nx*interp->ny*interp->nz;
        cppmallocl(interp->datak,shift,double complex);
        
        if(datatype=='c')
        {
            plan_f = fftw_plan_dft_2d(interp->nx, interp->ny, (double complex *)data, interp->datak, FFTW_FORWARD, WINTERP_USE_FFTW_PLANNER);
            fftw_execute(plan_f);
            fftw_destroy_plan(plan_f);
        }
        else // `r` or `v`
        {
            if(datatype=='v') ik=3;
            for(i=0; i<ik; i++)
            {
                for(j=0; j<shift; j++) *(interp->datak+i*shift+j) = *((double *)data+i*shift+j) + I*0.0;
                plan_f = fftw_plan_dft_2d(interp->nx, interp->ny, interp->datak+i*shift, interp->datak+i*shift, FFTW_FORWARD, WINTERP_USE_FFTW_PLANNER);
                fftw_execute(plan_f);
                fftw_destroy_plan(plan_f);
            }
        }
    }
    else if(datadim==3)
    {
        interp->nx=dims[0]; interp->dx=spacings[0];
        interp->ny=dims[1]; interp->dy=spacings[1];
        interp->nz=dims[2]; interp->dz=spacings[2];
        
        shift=interp->nx*interp->ny*interp->nz;
        cppmallocl(interp->datak,shift,double complex);
        
        if(datatype=='c')
        {
            plan_f = fftw_plan_dft_3d(interp->nx, interp->ny, interp->nz, (double complex *)data, interp->datak, FFTW_FORWARD, WINTERP_USE_FFTW_PLANNER);
            fftw_execute(plan_f);
            fftw_destroy_plan(plan_f);
        }
        else // `r` or `v`
        {
            if(datatype=='v') ik=3;
            for(i=0; i<ik; i++)
            {
                for(j=0; j<shift; j++) *(interp->datak+i*shift+j) = *((double *)data+i*shift+j) + I*0.0;
                plan_f = fftw_plan_dft_3d(interp->nx, interp->ny, interp->nz, interp->datak+i*shift, interp->datak+i*shift, FFTW_FORWARD, WINTERP_USE_FFTW_PLANNER);
                fftw_execute(plan_f);
                fftw_destroy_plan(plan_f);
            }
        }
    }
    else return WINTERP_ERR_INTERPOLATION_NOT_IMPLEMENTED; 
    
    // normalization
    for(i=0; i<shift*ik; i++) interp->datak[i]/=(interp->nx*interp->ny*interp->nz);
    
    return WINTERP_OK;
}

int winterp_create_interpolator_1d_c(int nx, double dx, Complex *data, winterp_interpolator *interp)
{
    int dims[1]={nx};
    double spacings[1]={dx};
    return winterp_create_interpolator('c', 1, dims, spacings, data, interp);
}

int winterp_create_interpolator_1d_r(int nx, double dx, double *data, winterp_interpolator *interp)
{
    int dims[1]={nx};
    double spacings[1]={dx};
    return winterp_create_interpolator('r', 1, dims, spacings, data, interp);
}

int winterp_create_interpolator_1d_v(int nx, double dx, double *data, winterp_interpolator *interp)
{
    int dims[1]={nx};
    double spacings[1]={dx};
    return winterp_create_interpolator('v', 1, dims, spacings, data, interp);
}

int winterp_create_interpolator_2d_c(int nx, int ny, double dx, double dy, Complex *data, winterp_interpolator *interp)
{
    int dims[2]={nx,ny};
    double spacings[2]={dx,dy};
    return winterp_create_interpolator('c', 2, dims, spacings, data, interp);
}

int winterp_create_interpolator_2d_r(int nx, int ny, double dx, double dy, double *data, winterp_interpolator *interp)
{
    int dims[2]={nx,ny};
    double spacings[2]={dx,dy};
    return winterp_create_interpolator('r', 2, dims, spacings, data, interp);
}

int winterp_create_interpolator_2d_v(int nx, int ny, double dx, double dy, double *data, winterp_interpolator *interp)
{
    int dims[2]={nx,ny};
    double spacings[2]={dx,dy};
    return winterp_create_interpolator('v', 2, dims, spacings, data, interp);
}

int winterp_create_interpolator_3d_c(int nx, int ny, int nz, double dx, double dy, double dz, Complex *data, winterp_interpolator *interp)
{
    int dims[3]={nx,ny,nz};
    double spacings[3]={dx,dy,dz};
    return winterp_create_interpolator('c', 3, dims, spacings, data, interp);
}

int winterp_create_interpolator_3d_r(int nx, int ny, int nz, double dx, double dy, double dz, double *data, winterp_interpolator *interp)
{
    int dims[3]={nx,ny,nz};
    double spacings[3]={dx,dy,dz};
    return winterp_create_interpolator('r', 3, dims, spacings, data, interp);
}

int winterp_create_interpolator_3d_v(int nx, int ny, int nz, double dx, double dy, double dz, double *data, winterp_interpolator *interp)
{
    int dims[3]={nx,ny,nz};
    double spacings[3]={dx,dy,dz};
    return winterp_create_interpolator('v', 3, dims, spacings, data, interp);
}

int winterp_destroy_interpolator(winterp_interpolator *interp)
{
    interp->datadim=0;
    free(interp->datak);
    return WINTERP_OK;
}

int __winterp_getvalue_1d_c(winterp_interpolator *interp, int shift, double x, Complex *val)
{
    register int ixyz=0;
    register int ix;
    register double kx;
    register double deltakx = 2.*M_PI/(interp->dx*interp->nx); // send to cache
    register double complex r = 0.0 + I*0.0;
    for(ix=0; ix<interp->nx; ix++)
    {
        if(ix<interp->nx/2) kx=deltakx*(ix           );
        else                kx=deltakx*(ix-interp->nx);
        
        r+=interp->datak[ixyz+shift]*cexp(I*( kx*x ));
        ixyz++;
    }
    
    *val = r; // save result
    return WINTERP_OK;
}

int winterp_getvalue_1d_c(winterp_interpolator *interp, double x, Complex *val)
{
    Complex r;
    int ierr=__winterp_getvalue_1d_c(interp, 0, x, &r); if(ierr!=WINTERP_OK) return ierr;
    *val = r;
    return ierr;
}

int winterp_getvalue_1d_r(winterp_interpolator *interp, double x, double *val)
{
    Complex r;
    int ierr=__winterp_getvalue_1d_c(interp, 0, x, &r); if(ierr!=WINTERP_OK) return ierr;
    *val = creal(r);
    return ierr;
}

int winterp_getvalue_1d_v(winterp_interpolator *interp, double x, double *valx, double *valy, double *valz)
{
    Complex r;
    int ierr;
    int shift=interp->nx;
    ierr=__winterp_getvalue_1d_c(interp, 0*shift, x, &r); if(ierr!=WINTERP_OK) return ierr;
    *valx = creal(r);
    ierr=__winterp_getvalue_1d_c(interp, 1*shift, x, &r); if(ierr!=WINTERP_OK) return ierr;
    *valy = creal(r);
    ierr=__winterp_getvalue_1d_c(interp, 2*shift, x, &r); if(ierr!=WINTERP_OK) return ierr;
    *valz = creal(r);
    return ierr;
}

int __winterp_getvalue_2d_c(winterp_interpolator *interp, int shift, double x, double y, Complex *val)
{
    register int ixyz=0;
    register int ix, iy;
    register double kx, ky;
    register double deltakx = 2.*M_PI/(interp->dx*interp->nx); // send to cache
    register double deltaky = 2.*M_PI/(interp->dy*interp->ny); // send to cache
    register double complex r = 0.0 + I*0.0;
    for(ix=0; ix<interp->nx; ix++)
    {
        if(ix<interp->nx/2) kx=deltakx*(ix           );
        else                kx=deltakx*(ix-interp->nx);
        
        for(iy=0; iy<interp->ny; iy++)
        {
            if(iy<interp->ny/2) ky=deltaky*(iy           );
            else                ky=deltaky*(iy-interp->ny);
        
            r+=interp->datak[ixyz+shift]*cexp(I*( kx*x + ky*y ));
            ixyz++;
        }
    }
    
    *val = r; // save result
    return WINTERP_OK;
}

int winterp_getvalue_2d_c(winterp_interpolator *interp, double x, double y, Complex *val)
{
    Complex r;
    int ierr=__winterp_getvalue_2d_c(interp, 0, x, y, &r); if(ierr!=WINTERP_OK) return ierr;
    *val = r;
    return ierr;
}

int winterp_getvalue_2d_r(winterp_interpolator *interp, double x, double y, double *val)
{
    Complex r;
    int ierr=__winterp_getvalue_2d_c(interp, 0, x, y, &r); if(ierr!=WINTERP_OK) return ierr;
    *val = creal(r);
    return ierr;
}

int winterp_getvalue_2d_v(winterp_interpolator *interp, double x, double y, double *valx, double *valy, double *valz)
{
    Complex r;
    int ierr;
    int shift=interp->nx*interp->ny;
    ierr=__winterp_getvalue_2d_c(interp, 0*shift, x, y, &r); if(ierr!=WINTERP_OK) return ierr;
    *valx = creal(r);
    ierr=__winterp_getvalue_2d_c(interp, 1*shift, x, y, &r); if(ierr!=WINTERP_OK) return ierr;
    *valy = creal(r);
    ierr=__winterp_getvalue_2d_c(interp, 2*shift, x, y, &r); if(ierr!=WINTERP_OK) return ierr;
    *valz = creal(r);
    return ierr;
}

int __winterp_getvalue_3d_c(winterp_interpolator *interp, int shift, double x, double y, double z, Complex *val)
{
    register int ixyz=0;
    register int ix, iy, iz;
    register double kx, ky, kz;
    register double deltakx = 2.*M_PI/(interp->dx*interp->nx); // send to cache
    register double deltaky = 2.*M_PI/(interp->dy*interp->ny); // send to cache
    register double deltakz = 2.*M_PI/(interp->dz*interp->nz); // send to cache
    register double complex r = 0.0 + I*0.0;
    for(ix=0; ix<interp->nx; ix++)
    {
        if(ix<interp->nx/2) kx=deltakx*(ix           );
        else                kx=deltakx*(ix-interp->nx);
        
        for(iy=0; iy<interp->ny; iy++)
        {
            if(iy<interp->ny/2) ky=deltaky*(iy           );
            else                ky=deltaky*(iy-interp->ny);
            
            for(iz=0; iz<interp->nz; iz++)
            {
                if(iz<interp->nz/2) kz=deltakz*(iz           );
                else                kz=deltakz*(iz-interp->nz);
            
                r+=interp->datak[ixyz+shift]*cexp(I*( kx*x + ky*y + kz*z));
                ixyz++;
            }
        }
    }
    
    *val = r; // save result
    return WINTERP_OK;
}

int winterp_getvalue_3d_c(winterp_interpolator *interp, double x, double y, double z, Complex *val)
{
    Complex r;
    int ierr=__winterp_getvalue_3d_c(interp, 0, x, y, z, &r); if(ierr!=WINTERP_OK) return ierr;
    *val = r;
    return ierr;
}

int winterp_getvalue_3d_r(winterp_interpolator *interp, double x, double y, double z, double *val)
{
    Complex r;
    int ierr=__winterp_getvalue_3d_c(interp, 0, x, y, z, &r); if(ierr!=WINTERP_OK) return ierr;
    *val = creal(r);
    return ierr;
}

int winterp_getvalue_3d_v(winterp_interpolator *interp, double x, double y, double z, double *valx, double *valy, double *valz)
{
    Complex r;
    int ierr;
    int shift=interp->nx*interp->ny*interp->nz;
    ierr=__winterp_getvalue_3d_c(interp, 0*shift, x, y, z, &r); if(ierr!=WINTERP_OK) return ierr;
    *valx = creal(r);
    ierr=__winterp_getvalue_3d_c(interp, 1*shift, x, y, z, &r); if(ierr!=WINTERP_OK) return ierr;
    *valy = creal(r);
    ierr=__winterp_getvalue_3d_c(interp, 2*shift, x, y, z, &r); if(ierr!=WINTERP_OK) return ierr;
    *valz = creal(r);
    return ierr;
}

