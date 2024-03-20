// Author: Gabriel Wlazlowski
// Date: 15-09-2016

// This file implements all functions needed for numerical computation of derivatives

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#include <stdio.h>

#include "pca_settings.h"
#include "pca_macro.h"

// structure that contains info needed to fast reconstruct the solution
typedef struct
{
    int batch_size;
    int nwfip;
    cufftHandle plans[CUFFT_NUMBER_OF_PLANS]; // handle to plans [BATCH Z2Z, ONE Z2Z, ONE D2Z, ONE Z2D]
    void *work_area; // work area of cufft plans 
} metadata_pca_cufftplans_t;

// Stores results;
metadata_pca_cufftplans_t __md_pca_cufftplans;
extern "C" void *pca_cufft_work_area; /* Declaration of the variable */
void *pca_cufft_work_area=NULL; /* Definition checked against declaration */

// *************** cuFFT ************************
/**
 * This function creates plans handle,
 * and provides size of workspace.
 * NOTE: you must call `set_workspace_for_cufftPlans` to set workspace for your plans
 * @param batch_size - for BATCH plans
 * @param nwfip - number of wave fuctions; 
 * @param workSize - requested worksize for this plan - OUTPUT
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int create_cufftPlans(int batch_size,  int nwfip, size_t *workSize)
{
    int fftwn[2] = {NX,NY};
    int idist, odist, istride, ostride;
    idist = odist = NXY; /* the distance in memory between the first element of the first array and the first element of the second array */
    istride = ostride = 1; /* array is contiguous in memory */
    int *inembed = fftwn, *onembed = fftwn;
    cufftResult cufft_result;
    int i;
    size_t max_work_Size =0;
    *workSize = 0;
    cufftHandle tplan;
    
    // store data
    if(batch_size>2*nwfip) batch_size=2*nwfip;
    __md_pca_cufftplans.batch_size=batch_size;
    __md_pca_cufftplans.nwfip=nwfip;
    
    for(i=0; i<CUFFT_NUMBER_OF_PLANS; i++)
    {
        cufft_result=cufftCreate(&tplan);
        __md_pca_cufftplans.plans[i]=tplan;
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;    
        
        cufft_result=cufftSetAutoAllocation(__md_pca_cufftplans.plans[i], 0);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    }
    
    // Create plans
    cufft_result=cufftMakePlanMany(__md_pca_cufftplans.plans[PLAN_Z2Z_BATCH], 2, fftwn,
                            inembed,
                            istride, idist,
                            onembed,
                            ostride, odist,
                            CUFFT_Z2Z, batch_size, &max_work_Size);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    if(max_work_Size>*workSize) *workSize=max_work_Size;
    
   
    cufft_result=cufftMakePlanMany(__md_pca_cufftplans.plans[PLAN_Z2Z_ONE], 2, fftwn,
                            inembed,
                            istride, idist,
                            onembed,
                            ostride, odist,
                            CUFFT_Z2Z, 1, &max_work_Size);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    if(max_work_Size>*workSize) *workSize=max_work_Size; 
    
    cufft_result=cufftMakePlan2d(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], NX, NY, CUFFT_D2Z, &max_work_Size);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    max_work_Size+=sizeof(cufftDoubleComplex)*NX*(NY/2+1)*3; // extra 3 arrays - see implementation of compute_gradient_real_f
#ifdef DERIVATIVE_COPY_DATA_MODE
    max_work_Size+=sizeof(double)*NXY;
#endif
    max_work_Size+=sizeof(cufftDoubleComplex)*NXY*PCA_WORKSPACE_SHIFT; // add extra PCA_WORKSPACE_SHIFT buffers as temporary data for hamiltonian execution (2*PCA_WORKSPACE_SHIFT buffers in double precision)
    if(max_work_Size>*workSize) *workSize=max_work_Size;
    
    cufft_result=cufftMakePlan2d(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], NX, NY, CUFFT_Z2D, &max_work_Size);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    max_work_Size+=sizeof(cufftDoubleComplex)*NX*(NY/2+1)*3; // extra 3 arrays - see implementation of compute_gradient_real_f
#ifdef DERIVATIVE_COPY_DATA_MODE
    max_work_Size+=sizeof(double)*NXY;
#endif
    max_work_Size+=sizeof(cufftDoubleComplex)*NXY*PCA_WORKSPACE_SHIFT; // add extra PCA_WORKSPACE_SHIFT buffers as temporary data for hamiltonian execution (2*PCA_WORKSPACE_SHIFT buffers in double precision)
    if(max_work_Size>*workSize) *workSize=max_work_Size;
        
    return 0;
}

/**
 * Function sets workspace for plans
 * NOTE: call create_cufftPlan before, to get plan handles and workSize of workArea
 * @param workArea - pointer to work area
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int set_workspace_for_cufftPlan(void *workArea)
{
    cufftResult cufft_result;
    int i;
    
    __md_pca_cufftplans.work_area = workArea; // save pointer to global structure
    pca_cufft_work_area = workArea; // save pointer to global structure
    cufftDoubleComplex *prt = (cufftDoubleComplex *)workArea;
    prt+=NX*(NY/2+1)*3;
#ifdef DERIVATIVE_COPY_DATA_MODE
    prt+=sizeof(double)*NXY;
#endif
    prt+=NXY*PCA_WORKSPACE_SHIFT;
    
    for(i=0; i<CUFFT_NUMBER_OF_PLANS; i++)
    {
        if(i==PLAN_D2Z_ONE || i==PLAN_Z2D_ONE)
        {
            cufft_result=cufftSetWorkArea(__md_pca_cufftplans.plans[i], (void *)prt);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;            
        }
        else
        {
            cufft_result=cufftSetWorkArea(__md_pca_cufftplans.plans[i], workArea);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        }
    }
    
    return 0;
}


// ================================================================================================
// =================================== compute_derivatives ========================================
// ================================================================================================
__global__ void kernel_compute_derivatives(int nwf, cufftDoubleComplex *in, cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz, cufftDoubleComplex *wf_laplace)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, i;
    double kx, ky, k2;
    size_t shift;
    cufftDoubleComplex z, zz;

    if(ixyz<NXY)
    {
        ixy2ixiy2d(ixyz,ix,iy); // decode cartesian coordinates
        
        // extract momentum
        if(ix<NX/2) kx=2.*M_PI/( double )LX * ( double )(ix   );
        else        kx=2.*M_PI/( double )LX * ( double )(ix-NX);
        
        if(iy<NY/2) ky=2.*M_PI/( double )LY * ( double )(iy   );
        else        ky=2.*M_PI/( double )LY * ( double )(iy-NY);
        
        
        k2 = -1.0*(kx*kx + ky*ky)/NXY; // note: normalization factor is included
        kx/=NXY; // note: normalization factor is included
        ky/=NXY; // note: normalization factor is included
        
        // corrections for gradient computation
        if(ix==NX/2) kx=0.0;
        if(iy==NY/2) ky=0.0;    
        
        // for each wave-function
        shift = ixyz;
        for(i=0; i<nwf; i++)
        {   
            // load data
            z = in[shift];
            
            // Multiply by momentum and save
            // zz=i*k*z = i*k *(z.x + i*z.y) =i*k*z.x - k*z.y             
            // dx
            zz.x=-1.0*kx*z.y;
            zz.y=     kx*z.x;
            wf_d_dx[shift]=zz;
            
            // dy 
            zz.x=-1.0*ky*z.y;
            zz.y=     ky*z.x;
            wf_d_dy[shift]=zz;
                        
            // laplace
            zz.x=k2*z.x;
            zz.y=k2*z.y;      
            wf_laplace[shift]=zz;

            // pointer to next wave-function
            shift+=(size_t)NXY;            
        }
    }
}

/**
 * @param n  number of wave-functions to process, where n=2*nwfip
 * @param wf array with wave-functions (INPUT)
 * @param wf_d_dx derivative with respect to dx (OUTPUT)
 * @param wf_d_dy derivative with respect to dy (OUTPUT)
 * @param wf_d_dz derivative with respect to dz (OUTPUT) --- NOT USED !!!!
 * @param wf_d_laplace laplace of wave-functions (OUTPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int compute_derivatives(int n, cufftDoubleComplex *wf, cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz, cufftDoubleComplex *wf_laplace, int nthreads)
{
    size_t shift;
    int iters1=n / __md_pca_cufftplans.batch_size;
    int iters2=n % __md_pca_cufftplans.batch_size;
    int i;
    cufftResult cufft_result;
    
    // number of blocks
    int nblocks = (int)ceil((float)NXY/nthreads);
    
//     printf("n=%d, iters1=%d, iters2=%d\n", n, iters1, iters2);
    
    // Step 1: Go to momentum space
    shift=0;
    for(i=0; i<iters1; i++) // transform with batches
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_BATCH], wf+shift, wf_laplace+shift, CUFFT_FORWARD);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)__md_pca_cufftplans.batch_size*NXY; // shift pointer by batch number of functions
    }
    for(i=0; i<iters2; i++) // remaining one by one
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_ONE], wf+shift, wf_laplace+shift, CUFFT_FORWARD);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)NXY; // shift pointer by one function
    }
        
    // Step 2: Multiply by momentum 
    kernel_compute_derivatives<<<nblocks, nthreads>>>(n, wf_laplace, wf_d_dx, wf_d_dy, NULL, wf_laplace);
        
    // Step 3: Go back to coordinate space
    shift=0;
    for(i=0; i<iters1; i++) // transform with batches
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_BATCH], wf_laplace+shift, wf_laplace+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_BATCH], wf_d_dx+shift, wf_d_dx+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_BATCH], wf_d_dy+shift, wf_d_dy+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)__md_pca_cufftplans.batch_size*NXY; // shift pointer by batch number of functions
    }
    for(i=0; i<iters2; i++) // remaining one by one
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_ONE], wf_laplace+shift, wf_laplace+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_ONE], wf_d_dx+shift, wf_d_dx+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_ONE], wf_d_dy+shift, wf_d_dy+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)NXY; // shift pointer by one function
    }

    
    return 0;
}

// ================================================================================================
// ==================================== compute_laplace ===========================================
// ================================================================================================
__global__ void kernel_compute_laplace(int nwf, cufftDoubleComplex *inout)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, i;
    double kx, ky, k2;
    size_t shift;
    cufftDoubleComplex z, zz;

    if(ixyz<NXY)
    {
        ixy2ixiy2d(ixyz,ix,iy); // decode cartesian coordinates
        
        // extract momentum
        if(ix<NX/2) kx=2.*M_PI/( double )LX * ( double )(ix   );
        else        kx=2.*M_PI/( double )LX * ( double )(ix-NX);
        
        if(iy<NY/2) ky=2.*M_PI/( double )LY * ( double )(iy   );
        else        ky=2.*M_PI/( double )LY * ( double )(iy-NY);
        
        k2 = -1.0*(kx*kx + ky*ky)/NXY; // note: normalization factor is included
        
        // for each wave-function
        shift = ixyz;
        for(i=0; i<nwf; i++)
        {
            // load data
            z = inout[shift];
            
            // laplace
            zz.x=k2*z.x;
            zz.y=k2*z.y;      
            inout[shift]=zz;

            // pointer to next wave-function
            shift+=(size_t)NXY;            
        }
    }
}

/**
 * @param n  number of wave-functions to process
 * @param wf array with wave-functions (INPUT)
 * @param wf_laplace laplace of wave-functions, can be the same is wf (OUTPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int compute_laplace(int n, cufftDoubleComplex *wf, cufftDoubleComplex *wf_laplace, int nthreads)
{
    size_t shift;
    int iters1=n / __md_pca_cufftplans.batch_size;
    int iters2=n % __md_pca_cufftplans.batch_size;
    int i;
    cufftResult cufft_result;
    
    // number of blocks
    int nblocks = (int)ceil((float)NXY/nthreads);
    
//     printf("n=%d, iters1=%d, iters2=%d\n", n, iters1, iters2);
    
    // Step 1: Go to momentum space
    shift=0;
    for(i=0; i<iters1; i++) // transform with batches
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_BATCH], wf+shift, wf_laplace+shift, CUFFT_FORWARD);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)__md_pca_cufftplans.batch_size*NXY; // shift pointer by batch number of functions
    }
    for(i=0; i<iters2; i++) // remaining one by one
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_ONE], wf+shift, wf_laplace+shift, CUFFT_FORWARD);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)NXY; // shift pointer by one function
    }
        
    // Step 2: Multiply by momentum 
    kernel_compute_laplace<<<nblocks, nthreads>>>(n, wf_laplace);
        
    // Step 3: Go back to coordinate space
    shift=0;
    for(i=0; i<iters1; i++) // transform with batches
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_BATCH], wf_laplace+shift, wf_laplace+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)__md_pca_cufftplans.batch_size*NXY; // shift pointer by batch number of functions
    }
    for(i=0; i<iters2; i++) // remaining one by one
    {
        cufft_result=cufftExecZ2Z(__md_pca_cufftplans.plans[PLAN_Z2Z_ONE], wf_laplace+shift, wf_laplace+shift, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        shift+=(size_t)NXY; // shift pointer by one function
    }

    
    return 0;
}

// ================================================================================================
// =================================== compute_gradient_real_f ====================================
// ================================================================================================
__global__ void kernel_compute_gradient_real_f(cufftDoubleComplex *in, cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy;
    double kx, ky/*, k2*/;
    cufftDoubleComplex z, zz;

    if(ixyz<NX*(NY/2+1))
    {
        ixy2ixiy2dD2Z(ixyz,ix,iy); // decode cartesian coordinates
        
        // extract momentum
        if(ix<NX/2) kx=2.*M_PI/( double )LX/( double )NXY * ( double )(ix   ); // note: normalization factor is included
        else        kx=2.*M_PI/( double )LX/( double )NXY * ( double )(ix-NX); // note: normalization factor is included
        
        if(iy<NY/2) ky=2.*M_PI/( double )LY/( double )NXY * ( double )(iy   ); // note: normalization factor is included
        else        ky=2.*M_PI/( double )LY/( double )NXY * ( double )(iy-NY); // note: normalization factor is included
        
        /*k2 = -1.0*(kx*kx + ky*ky)/NXY; // note: normalization factor is included */
        
        // load data
        z = in[ixyz];
            
        // Multiply by momentum and save
        // zz=i*k*z = i*k *(z.x + i*z.y) =i*k*z.x - k*z.y 
            
        // dx
        if(ix==NX/2) kx=0.0; 
        zz.x=-1.0*kx*z.y;
        zz.y=     kx*z.x;
        wf_d_dx[ixyz]=zz;
            
        // dy
        if(iy==NY/2) ky=0.0; 
        zz.x=-1.0*ky*z.y;
        zz.y=     ky*z.x;
        wf_d_dy[ixyz]=zz;
                        
        /*
        // laplace
        zz.x=k2*z.x;
        zz.y=k2*z.y;      
        wf_laplace[ixyz]=zz;
        */         
    }
}

/**
 * Function computes gradient of real function
 * @param f pointer to real function (INPUT)
 * @param df_dx derivative with respect to dx (OUTPUT)
 * @param df_dy derivative with respect to dy (OUTPUT)
 * @param df_dz derivative with respect to dz (OUTPUT) --- NOT USED !!!!
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int compute_gradient_real_f(double *f, double *df_dx, double *df_dy, double *df_dz, int nthreads)
{
    cufftResult cufft_result;
    
    // number of blocks
    int nblocks = (int)ceil((float)(NX*(NY/2+1))/nthreads);

    // get pointer to workspace
    cufftDoubleComplex * p_df_dx = (cufftDoubleComplex *)__md_pca_cufftplans.work_area;
    p_df_dx+=PCA_WORKSPACE_SHIFT*NXY;
    cufftDoubleComplex * p_df_dy = p_df_dx+NX*(NY/2+1);
    
    // Step 1: go to momentum space
#ifdef DERIVATIVE_COPY_DATA_MODE
    double * p_tmp = (double *) (p_df_dy+NX*(NY/2+1));
    if( cudaMemcpy( p_tmp , f , NXY*sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], p_tmp, p_df_dx);
#else
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], f, p_df_dx);
#endif
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;  
        
    // Step 2: Multiply by momentum
    kernel_compute_gradient_real_f<<<nblocks, nthreads>>>(p_df_dx, p_df_dx, p_df_dy, NULL);

    // Step 3: go back to coordinate space
    cufft_result=cufftExecZ2D(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], p_df_dx, df_dx);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result; 
    cufft_result=cufftExecZ2D(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], p_df_dy, df_dy);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;   
    
    return 0;
}

// ================================================================================================
// =============================== compute_derivative_real_vector_f ===============================
// ================================================================================================
__global__ void kernel_compute_derivative_real_vector_f(cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy;
    double kx, ky/*, k2*/;
    cufftDoubleComplex z, zz;

    if(ixyz<NX*(NY/2+1))
    {
        ixy2ixiy2dD2Z(ixyz,ix,iy); // decode cartesian coordinates
        
        // extract momentum
        if(ix<NX/2) kx=2.*M_PI/( double )LX/( double )NXY * ( double )(ix   ); // note: normalization factor is included
        else        kx=2.*M_PI/( double )LX/( double )NXY * ( double )(ix-NX); // note: normalization factor is included
        
        if(iy<NY/2) ky=2.*M_PI/( double )LY/( double )NXY * ( double )(iy   ); // note: normalization factor is included
        else        ky=2.*M_PI/( double )LY/( double )NXY * ( double )(iy-NY); // note: normalization factor is included
        
        /*k2 = -1.0*(kx*kx + ky*ky + kz*kz)/NXY; // note: normalization factor is included */
                                
        // dx
        if(ix==NX/2) kx=0.0; 
        z=wf_d_dx[ixyz];
        zz.x=-1.0*kx*z.y;
        zz.y=     kx*z.x;
        wf_d_dx[ixyz]=zz;
            
        // dy
        if(iy==NY/2) ky=0.0; 
        z=wf_d_dy[ixyz];
        zz.x=-1.0*ky*z.y;
        zz.y=     ky*z.x;
        wf_d_dy[ixyz]=zz;
                                 
    }
}

/**
 * Function computes derivatives of real vector function: dfx/dx,  dfy/dy , dfz/dz
 * @param fx pointer to real function, x coordinate (INPUT)
 * @param fy pointer to real function, y coordinate (INPUT)
 * @param fz pointer to real function, z coordinate (INPUT) --- NOT USED !!!!
 * @param dfx_dx derivative dfx/dx, can be the same as fx (OUTPUT)
 * @param dfy_dy derivative dfy/dy, can be the same as fy (OUTPUT)
 * @param dfz_dz derivative dfz/dz, can be the same as fz (OUTPUT) --- NOT USED !!!!
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int compute_derivative_real_vector_f(double *fx, double *fy, double *fz, 
                                                double *dfx_dx, double *dfy_dy, double *dfz_dz,
                                                int nthreads)
{
    cufftResult cufft_result;
    
    // number of blocks
    int nblocks = (int)ceil((float)(NX*(NY/2+1))/nthreads);

    // get pointer to workspace
    cufftDoubleComplex * p_fx = (cufftDoubleComplex *)__md_pca_cufftplans.work_area;
    p_fx+=PCA_WORKSPACE_SHIFT*NXY;
    cufftDoubleComplex * p_fy = p_fx+NX*(NY/2+1);
    
    // Step 1: go to momentum space
#ifdef DERIVATIVE_COPY_DATA_MODE
    double * p_tmp = (double *) (p_fy+NX*(NY/2+1));
    if( cudaMemcpy( p_tmp , fx , NXY*sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], p_tmp, p_fx);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result; 
    if( cudaMemcpy( p_tmp , fy , NXY*sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], p_tmp, p_fy);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result; 
#else
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], fx, p_fx);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;  
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], fy, p_fy);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;   
#endif
        
    // Step 2: Multiply by momentum
    kernel_compute_derivative_real_vector_f<<<nblocks, nthreads>>>(p_fx, p_fy, NULL);

    // Step 3: go back to coordinate space
    cufft_result=cufftExecZ2D(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], p_fx, dfx_dx);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result; 
    cufft_result=cufftExecZ2D(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], p_fy, dfy_dy);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;   
    
    return 0;
}

// ================================================================================================
// ==================================== compute_laplace_real_f ====================================
// ================================================================================================
__global__ void kernel_compute_laplace_real_f(cufftDoubleComplex *wf_d)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy;
    double kx, ky, k2;
    cufftDoubleComplex z, zz;

    if(ixyz<NX*(NY/2+1))
    {
        ixy2ixiy2dD2Z(ixyz,ix,iy); // decode cartesian coordinates
        
        // extract momentum
        if(ix<NX/2) kx=2.*M_PI/( double )LX * ( double )(ix   );
        else        kx=2.*M_PI/( double )LX * ( double )(ix-NX);
        
        if(iy<NY/2) ky=2.*M_PI/( double )LY * ( double )(iy   );
        else        ky=2.*M_PI/( double )LY * ( double )(iy-NY);
        
        
        k2 = -1.0*(kx*kx + ky*ky)/NXY; // note: normalization factor is included
                                
        // laplace
        z=wf_d[ixyz];
        zz.x=k2*z.x;
        zz.y=k2*z.y;      
        wf_d[ixyz]=zz;         
    }
}

/**
 * Function computes laplace of real function laplace_f = d^2f/dx^2 + d^2f/dy^2
 * @param f pointer to real function (INPUT)
 * @param laplace_f laplace of function, can be the same as f (OUTPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int compute_laplace_real_f(double *f, double *laplace_f, int nthreads)
{
    cufftResult cufft_result;
    
    // number of blocks
    int nblocks = (int)ceil((float)(NX*(NY/2+1))/nthreads);

    // get pointer to workspace
    cufftDoubleComplex * p_fx = (cufftDoubleComplex *)__md_pca_cufftplans.work_area;
    p_fx+=PCA_WORKSPACE_SHIFT*NXY;
    
    // Step 1: go to momentum space
#ifdef DERIVATIVE_COPY_DATA_MODE
    double * p_tmp = (double *) (p_fx+NX*(NY/2+1));
    if( cudaMemcpy( p_tmp , f , NXY*sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], p_tmp, p_fx);
#else
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], f, p_fx);
#endif
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;  
        
    // Step 2: Multiply by momentum
    kernel_compute_laplace_real_f<<<nblocks, nthreads>>>(p_fx);

    // Step 3: go back to coordinate space
    cufft_result=cufftExecZ2D(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], p_fx, laplace_f);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result; 
    
    return 0;
}

// ================================================================================================
// =============================== compute_divergence_real_vector_f ===============================
// ================================================================================================
__global__ void kernel_compute_divergence_real_vector_f(cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy;
    double kx, ky/*, k2*/;
    cufftDoubleComplex z, zz, res;

    if(ixyz<NX*(NY/2+1))
    {
        ixy2ixiy2dD2Z(ixyz,ix,iy); // decode cartesian coordinates
        
        // extract momentum
        if(ix<NX/2)  kx=2.*M_PI/( double )LX/( double )NXY * ( double )(ix   ); // note: normalization factor is included
        else         kx=2.*M_PI/( double )LX/( double )NXY * ( double )(ix-NX); // note: normalization factor is included
        
        if(iy<NY/2)  ky=2.*M_PI/( double )LY/( double )NXY * ( double )(iy   ); // note: normalization factor is included
        else         ky=2.*M_PI/( double )LY/( double )NXY * ( double )(iy-NY); // note: normalization factor is included
        
        /*k2 = -1.0*(kx*kx + ky*ky + kz*kz)/NXY; // note: normalization factor is included */
                                
        res.x=0.0; res.y=0.0;
        
        // dx
        if(ix==NX/2) kx=0.0; 
        z=wf_d_dx[ixyz];
        zz.x=-1.0*kx*z.y;
        zz.y=     kx*z.x;
        res.x+=zz.x; res.y+=zz.y;
            
        // dy
        if(iy==NY/2) ky=0.0; 
        z=wf_d_dy[ixyz];
        zz.x=-1.0*ky*z.y;
        zz.y=     ky*z.x;
        res.x+=zz.x; res.y+=zz.y;
            
        wf_d_dz[ixyz]=res; // save result
                     
    }
}

/**
 * Function computes derivatives of real vector function: dfx/dx,  dfy/dy , dfz/dz
 * @param fx pointer to real function, x coordinate (INPUT)
 * @param fy pointer to real function, y coordinate (INPUT)
 * @param fz pointer to real function, z coordinate (INPUT) --- NOT USED !!!!
 * @param divf divergence of vector (fx,fy,fz), ie divf = dfx/dx + dfy/dy (OUTPUT)
 *             can be the same as one of input pointers
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int compute_divergence_real_vector_f(double *fx, double *fy, double *fz, 
                                                double *divf,
                                                int nthreads)
{
    cufftResult cufft_result;
    
    // number of blocks
    int nblocks = (int)ceil((float)(NX*(NY/2+1))/nthreads);

    // get pointer to workspace
    cufftDoubleComplex * p_fx = (cufftDoubleComplex *)__md_pca_cufftplans.work_area;
    p_fx+=PCA_WORKSPACE_SHIFT*NXY;
    cufftDoubleComplex * p_fy = p_fx+NX*(NY/2+1);
    cufftDoubleComplex * p_fz = p_fy+NX*(NY/2+1);
    
    // Step 1: go to momentum space
#ifdef DERIVATIVE_COPY_DATA_MODE
    double * p_tmp = (double *) (p_fz+NX*(NY/2+1));
    if( cudaMemcpy( p_tmp , fx , NXY*sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], p_tmp, p_fx);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;  
    if( cudaMemcpy( p_tmp , fy , NXY*sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], p_tmp, p_fy);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;  
#else
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], fx, p_fx);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;  
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], fy, p_fy);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;   
#endif
        
    // Step 2: Multiply by momentum
    // (ikx*p_fx + iky*p_fy) -> p_fz
    kernel_compute_divergence_real_vector_f<<<nblocks, nthreads>>>(p_fx, p_fy, p_fz);

    // Step 3: go back to coordinate space
    cufft_result=cufftExecZ2D(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], p_fz, divf);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result; 
    
    return 0;
}

// ================================================================================================
// ====================================== high_frequency_filter ===================================
// ================================================================================================
__global__ void kernel_high_frequency_filter_d(cufftDoubleComplex *wf_d, double fd_mu, double fd_T)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, i;
    double kx, ky, ek;
    cufftDoubleComplex z, zz;

    if(ixyz<NX*(NY/2+1))
    {
        ixy2ixiy2dD2Z(ixyz,ix,iy); // decode cartesian coordinates

        // extract momentum
        if(ix<NX/2) kx=2.*M_PI/( double )LX * ( double )(ix   );
        else        kx=2.*M_PI/( double )LX * ( double )(ix-NX);

        if(iy<NY/2) ky=2.*M_PI/( double )LY * ( double )(iy   );
        else        ky=2.*M_PI/( double )LY * ( double )(iy-NY);

        ek = 0.5*(kx*kx + ky*ky);
        ek = 1.0/(exp((ek-fd_mu)/fd_T)+1.0) /NXY; // note: normalization factor is included

        // laplace
        z=wf_d[ixyz];
        zz.x=ek*z.x;
        zz.y=ek*z.y;
        wf_d[ixyz]=zz;
    }
}

/**
 * The function removes high frequencies from the signal represented by in array.
 * The method uses the spectral method:
 *     in ->fft->multiply by FD(mu,T)->ifft->out,
 * where FD(mu,T)=1.0/(exp((ek-mu)/T)+1.0) is the Fermi-Dirac function and ek=k^2/2.
 * @param in signal to be cleaned (INPUT)
 * @param out can be the same as input (OUTPUT)
 * @param mu parameter of FD function (INPUT)
 * @param T parameter of FD function (INPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int high_frequency_filter_d(double *in, double *out,
                                     double fd_mu, double fd_T,
                                     int nthreads)
{
    cufftResult cufft_result;

    // number of blocks
    int nblocks = (int)ceil((float)(NX*(NY/2+1))/nthreads);

    // get pointer to workspace
    cufftDoubleComplex * p_fx = (cufftDoubleComplex *)__md_pca_cufftplans.work_area;
    p_fx+=PCA_WORKSPACE_SHIFT*NXY;

    // Step 1: go to momentum space
#ifdef DERIVATIVE_COPY_DATA_MODE
    double * p_tmp = (double *)(p_fx+NX*(NY/2+1));
    if( cudaMemcpy( p_tmp , in , NXY*sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], p_tmp, p_fx);
#else
    cufft_result=cufftExecD2Z(__md_pca_cufftplans.plans[PLAN_D2Z_ONE], in, p_fx);
#endif
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;

    // Step 2: Multiply by filter function
    kernel_high_frequency_filter_d<<<nblocks, nthreads>>>(p_fx, fd_mu, fd_T);

    // Step 3: go back to coordinate space
    cufft_result=cufftExecZ2D(__md_pca_cufftplans.plans[PLAN_Z2D_ONE], p_fx, out);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;

    return 0;
}

/**
 * The function removes high frequencies from the signal represented by in array.
 * The method uses the spectral method:
 *     in ->fft->multiply by FD(mu,T)->ifft->out,
 * where FD(mu,T)=1.0/(exp((ek-mu)/T)+1.0) is the Fermi-Dirac function and ek=k^2/2.
 * @param n number of arrays (MASSIVE MODE)
 * @param in signal to be cleaned (INPUT)
 * @param out can be the same as input (OUTPUT)
 * @param mu parameter of FD function (INPUT)
 * @param T parameter of FD function (INPUT)
 * @return 0-OK, otherwise PROBLEM
 * */
extern "C" int high_frequency_filter_massive_d(int n, double *in, double *out,
                                     double fd_mu, double fd_T,
                                     int nthreads)
{
    // FIXME: it can be improved by executing FFT Many
    int i, ierr;
    for(i=0; i<n; i++)
    {
        ierr = high_frequency_filter_d(in+i*NXY, out+i*NXY, fd_mu,fd_T,nthreads);
        if(ierr!=0) return ierr;
    }

    return 0;
}
