// Author: Gabriel Wlazlowski
// Date: 23-04-2016

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>
#include <stdio.h>
typedef thrust::complex<double> Complex;

#include "pca_settings.h"
#include "pca_macro.h"
#include "pca_edf.h"
#include "wslda_cuda_utils.h"

// ===========================================================================
// ============================ CONSTANTS ====================================
// ===========================================================================
__constant__ double dc_mu_a; // chemical potentials 
__constant__ double dc_mu_b; // chemical potentials 
__constant__ double dc_ec; // energy cut-off 
__constant__ double dc_t0; // initial time, time=dc_t0 + dc_dt*it
__constant__ double dc_dt; // itegration time step
__constant__ double dc_kF; // reference kF
__constant__ double dc_eF; // reference eF (=kF^2/2)
__constant__ double dc_nF; // reference density nF (=kF^3 / (3*pi^2)) 

#ifdef BDG_MODE
__constant__ double dc_gBdG;
#endif

__constant__ void *dc_extra_data;
__constant__ size_t dc_extra_data_size;

extern "C" int memcopy_extra_data(size_t extra_data_size, void *extra_data)
{
    if( cudaMemcpyToSymbol(dc_extra_data,      &extra_data,      sizeof(void *))!= cudaSuccess ) return 1;
    if( cudaMemcpyToSymbol(dc_extra_data_size, &extra_data_size, sizeof(size_t))!= cudaSuccess ) return 2;
    return 0;
}


#ifdef TDWSLDA

__constant__ double dc_params[MAX_USER_PARAMS]; // array with params from input file
#include "problem-definition.h"

#ifdef ENABLE_V_EXT
#define u_ext(ix, iy, iz, it, spin) v_ext(ix, 0, 0, it, spin, dc_params, dc_extra_data_size, dc_extra_data)
#else 
#define u_ext(ix, iy, iz, it, spin) 0.0
#endif 

#ifdef ENABLE_DELTA_EXT
#define macro_delta_ext(ix, iy, iz, it, delta) delta_ext(ix, 0, 0, it, delta, dc_params, dc_extra_data_size, dc_extra_data)
#else 
#define macro_delta_ext(ix, iy, iz, it, delta) Complex(0.0,0.0)
#endif 

#else

#include "pca_uext.h"

#ifdef ENABLE_DELTA_EXT
#define macro_delta_ext(ix, iy, iz, it, delta) delta_ext(ix, 0, 0, it, delta)
#else 
#define macro_delta_ext(ix, iy, iz, it, delta) Complex(0.0,0.0)
#endif 

#endif

/**
 * This function copies data to constant memory buffers
 * */
extern "C" int memcopy_const(double mu_a, double mu_b, double ec, double t0, double dt, double kF)
{
    if( cudaMemcpyToSymbol(dc_mu_a, &mu_a, sizeof(double))!= cudaSuccess ) return 1;
    if( cudaMemcpyToSymbol(dc_mu_b, &mu_b, sizeof(double))!= cudaSuccess ) return 2;
    if( cudaMemcpyToSymbol(dc_ec, &ec, sizeof(double))!= cudaSuccess ) return 3;
    if( cudaMemcpyToSymbol(dc_t0, &t0, sizeof(double))!= cudaSuccess ) return 4;
    if( cudaMemcpyToSymbol(dc_dt, &dt, sizeof(double))!= cudaSuccess ) return 5;
    if( cudaMemcpyToSymbol(dc_kF, &kF, sizeof(double))!= cudaSuccess ) return 6;
    double eF = kF*kF/2.;
    if( cudaMemcpyToSymbol(dc_eF, &eF, sizeof(double))!= cudaSuccess ) return 7;
    double nF = kF*kF*kF / (3.*M_PI*M_PI);
    if( cudaMemcpyToSymbol(dc_nF, &nF, sizeof(double))!= cudaSuccess ) return 7;
    
    return 0;
}

/**
 * This function copies params to constant memory buffer
 * */
extern "C" int memcopy_const_params(double *params)
{
    if( cudaMemcpyToSymbol(dc_params, params, MAX_USER_PARAMS*sizeof(double))!= cudaSuccess ) return 1;
    
    return 0;
}

#ifdef BDG_MODE
/**
 * This function copies BdG functional data
 * */
extern "C" int memcopy_const_BdG(double aBdG)
{
    double gBdG = 4.0*M_PI*aBdG;
    if( cudaMemcpyToSymbol(dc_gBdG, &gBdG, sizeof(double))!= cudaSuccess ) return 1;
    
    return 0;
}
#endif

// ===========================================================================
// ============================ FUNCTIONS ====================================
// ===========================================================================

extern "C" int set_gpu(int device)
{
    cudaError err=cudaSetDevice( device );
    return (int)(err);    
}

extern "C" int gpu_malloc(size_t size, void ** pointer)
{
    // Allocate memory on GPU
    cudaError err=cudaMalloc( pointer , size );
    return (int)(err);
}

extern "C" int gpu_free(void * pointer)
{
    cudaError err=cudaFree(pointer );
    return (int)(err);
}

extern "C" int host_malloc_pl(size_t size, void ** pointer)
{
    // Allocate memory on GPU using page locked fashion
    cudaError err=cudaHostAlloc( pointer , size, cudaHostAllocDefault );
    return (int)(err);
}

extern "C" int host_free_pl(void * pointer)
{
    cudaError err=cudaFreeHost(pointer );
    return (int)(err);
}

extern "C" int memcopy_host2gpu(void * host, void * gpu, size_t size)
{
    cudaError err=cudaMemcpy( gpu , host , size, cudaMemcpyHostToDevice );    
    return (int)(err);
}

extern "C" int memcopy_gpu2host(void * gpu, void * host, size_t size)
{
    cudaError err=cudaMemcpy( host , gpu , size, cudaMemcpyDeviceToHost );   
    return (int)(err);
}

extern "C" int memcopy_gpu2gpu(void * gpusrc, void * gpudst, size_t size)
{
    cudaError err=cudaMemcpy( gpudst , gpusrc , size, cudaMemcpyDeviceToDevice );   
    return (int)(err);
}


////////////////////////////////////////////////////////////////////////////////
// LOCAL REDUCTIONS
////////////////////////////////////////////////////////////////////////////////
int opt_threads(int new_blocks,int threads, int current_size)
{
    int new_threads;
    if(new_blocks==1) 
    {
        new_threads=2; 
        while(new_threads<threads) 
        { 
            if(new_threads>=current_size) break;
            new_threads*=2;
        }
    }
    else new_threads=threads;
    return new_threads;
}

template <unsigned int blockSize>
__device__ void warpReduceR(volatile double *sdata, unsigned int tid) 
{
    if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
    if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
    if (blockSize >= 16) sdata[tid] += sdata[tid +  8];
    if (blockSize >=  8) sdata[tid] += sdata[tid +  4];
    if (blockSize >=  4) sdata[tid] += sdata[tid +  2];
    if (blockSize >=  2) sdata[tid] += sdata[tid +  1];
}

template <unsigned int blockSize>
__global__ void __reduce_kernelR__(double *g_idata, double *g_odata, int n, int mode)
{
    extern __shared__ double sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    unsigned int ishift=i+blockDim.x;
    
    // Loading data
//     if(mode==0) // sum of doubles
    {
        if(ishift<n) sdata[tid] = g_idata[i] + g_idata[ishift];
        else if(i<n) sdata[tid] = g_idata[i];
        else         sdata[tid] = 0.0;
    }
//     else // add here other modes
    

    __syncthreads();
    
    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    if (tid < 32) warpReduceR<blockSize>(sdata, tid);

    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

void call_reduction_kernelR(int dimGrid, int dimBlock, int size, double *d_idata, double *d_odata, int mode)
{
    int smemSize=dimBlock*sizeof(double);
    switch (dimBlock)
    {
        case 1024:
            __reduce_kernelR__<1024><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size, mode); break;
        case 512:
            __reduce_kernelR__< 512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size, mode); break;
        case 256:
            __reduce_kernelR__< 256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size, mode); break;
        case 128:
            __reduce_kernelR__< 128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size, mode); break;
        case 64:
            __reduce_kernelR__<  64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, size, mode); break;
    }   
}
/**
 * Function does fast reduction (sum of elements) of array. 
 * Result is located in partial_sums[0] element
 * If partial_sums==array then array will be destroyed
 * @param mode 0: add numbers (no transformation)
 * */ 
extern "C" int local_reductionR(double *array, int size, double *partial_sums, int threads, int mode)
{
    int blocks=(int)ceil((float)size/threads);
    unsigned int lthreads=threads/2; // Threads is always power of 2
    if(lthreads<64) lthreads=64; // at least 2*warp_size
    unsigned int new_blocks, current_size;

    // First reduction of the array
    call_reduction_kernelR(blocks, lthreads, size, array, partial_sums, mode);
    
    // Do iteratively reduction of partial_sums
    current_size=blocks;
    while(current_size>1)
    {
        new_blocks=(int)ceil((float)current_size/threads);
        lthreads=opt_threads(new_blocks,threads, current_size)/2;
        if(lthreads<64) lthreads=64; // at least 2*warp_size
        call_reduction_kernelR(new_blocks, lthreads, current_size, partial_sums, partial_sums, 0);
        current_size=new_blocks;
    }    
    
    return 0;
}

// =======================================================================================
// ================================ compute_potentials ===================================
// =======================================================================================
__global__ void kernel_compute_potentials(int it, 
                                          double *rho_a, double *rho_b, Complex *nu,
                                          double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                          double *tau_a, double *tau_b, 
                                          double *V_a, double *V_b, Complex *delta, double cccoeff
                                         )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, i;
    
    // registers
    double na, nb;
    double t1, t2, t3, t4, t5, t6, t7; // working buffers
    
    double alph_plus;
#ifdef CURRENT_CORRECTIONS
    double alph_minus;
#endif
    double dalphm_dna, dalphm_dnb, dalphp_dna, dalphp_dnb;
    double Va, Vb, Vanew, Vbnew, Va_const, Vb_const;
    Complex p0, kc, wz_0, Zone, lnu, ldelta;
    
    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        
        // load data to registers from global memory and form constant part of potentials
        Va_const=u_ext(ix,0,0,it,SPINA);
        Vb_const=u_ext(ix,0,0,it,SPINB);
        
        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        t1 = polarization(na, nb);
        alph_plus = alpha_plus(t1);
#ifdef CURRENT_CORRECTIONS
        alph_minus = alpha_minus(t1);
#endif
        
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
        dalphm_dna=0.0;
        dalphm_dnb=0.0;
        dalphp_dna=0.0;
        dalphp_dnb=0.0;
#else
        // term dalphm_dna*tau_m/2.0
        t1=tau_a[ixyz]; // tau_a
        t2=tau_b[ixyz]; // tau_b       
        t3=t1 - t2; // tau_m
        dalphm_dna=der_alpha_minus__der_na(na, nb);
        dalphm_dnb=der_alpha_minus__der_nb(na, nb);
        Va_const+=dalphm_dna*t3/2.0; // dalphm_dna*tau_m/2.0
        Vb_const+=dalphm_dnb*t3/2.0; // dalphm_dnb*tau_m/2.0
        
        // term dalphp_dna*tau_p/2.0
        t3=t1 + t2; // tau_p
        dalphp_dna=der_alpha_plus__der_na(na, nb);
        dalphp_dnb=der_alpha_plus__der_nb(na, nb);
        Va_const+=dalphp_dna*t3/2.0; // dalphp_dna*tau_p/2.0
        Vb_const+=dalphp_dnb*t3/2.0; // dalphp_dnb*tau_p/2.0
        // no other terms with tau, now I can resue t1 and t2
#endif

        // term dD_dna and dD_dnb
        Va_const += der_funD__der_na(na, nb);
        Vb_const += der_funD__der_nb(na, nb);
        
#ifdef CURRENT_CORRECTIONS
        // current terms
        t1=j_a_x[ixyz];
        t2=0.0;
        t3=0.0;
        t4=j_b_x[ixyz];
        t5=0.0;
        t6=0.0;
        // terms with ja^2
        t7 = p_regularization(na) * cccoeff;
        if(t7!=0.0)
        {
            t7 = t7*(t1*t1 + t2*t2 + t3*t3)/(2.0*na); // fr(na)*ja^2/2na
            Va_const+=( (alph_plus+alph_minus-1.0)/na - (dalphp_dna+dalphm_dna) ) *t7; // fr(na)*(alpha_a-1)*ja^2/2na^2 - fr(na)*dalpha_dna*ja^2/2na
            Va_const-=cccoeff*der_p_regularization(na)*(alph_plus+alph_minus-1.0)*(t1*t1 + t2*t2 + t3*t3)/(2.0*na); // derivative of regularization function
            Vb_const-=(dalphp_dnb+dalphm_dnb) *t7; // -dalpha_dnb*fr(na)*ja^2/2na            
        }

        // terms with jb^2
        t7 = p_regularization(nb) * cccoeff;
        if(t7!=0.0)
        {
            t7 = t7*(t4*t4 + t5*t5 + t6*t6)/(2.0*nb); // fr(nb)*jb^2/2nb
            Va_const-=(dalphp_dna-dalphm_dna) *t7; // -dalphb_dna*fr(nb)*jb^2/2nb
            Vb_const+=( (alph_plus-alph_minus-1.0)/nb - (dalphp_dnb-dalphm_dnb-1.0) ) *t7; // fr(nb)*(alpha_b-1)*jb^2/2nb^2 - fr(nb)*dalphb_dnb*jb^2/2nb     
            Vb_const-=cccoeff*der_p_regularization(nb)*(alph_plus-alph_minus-1.0)*(t4*t4 + t5*t5 + t6*t6)/(2.0*nb); // derivative of regularization function
        }
        
// //         // terms with j+^2
// //         t7 = p_regularization(na+nb)*((t1+t4)*(t1+t4) + (t2+t5)*(t2+t5) + (t3+t6)*(t3+t6))/(2.0*(na+nb)*(na+nb)); // j+^2/2n+^2
// //         Va_const-=t7;
// //         Vb_const-=t7;
// //         t7 = der_p_regularization(na+nb)*((t1+t4)*(t1+t4) + (t2+t5)*(t2+t5) + (t3+t6)*(t3+t6))/(2.0*(na+nb)); // derivative of regularization function
// //         Va_const+=t7;
// //         Vb_const+=t7;
        // NOTE: divergence terms include when applied hamiltonian - here not needed
#endif
        
        // prepare other variables for self-consistent process
        t1=dalphp_dna/alph_plus; 
        t2=dalphp_dnb/alph_plus;
        t3=der_tildeC__der_na(na, nb, 1.0) / alph_plus; // dtildeC_dna / alph_plus
        t4=der_tildeC__der_nb(na, nb, 1.0) / alph_plus; // dtildeC_dnb / alph_plus
        t5 = tildeC(na, nb, 1.0); // tC
        Va = V_a[ixyz]; // initial values
        Vb = V_b[ixyz]; // initial values
        lnu = nu[ixyz];
        Zone = Complex(1.0, 0.0);
  
        // computation of Va and Vb and delta
        for(i=0; i<UD_SCITERS; i++) // self-consistent loop
        {
            // pairing
#ifdef USE_CUBIC_CUTOFF
            wz_0=Complex(REGULARIZATION_SCHEME_K_CONST/(4.0*M_PI*DX), 0.0);
#else
            t7=(dc_mu_a-Va+dc_mu_b-Vb)/2.0;
            p0 = thrust::sqrt( Complex(2.0*t7/ alph_plus, 0.0) );
            if(p0.imag()<0.) p0 *= -1. ;
            kc = thrust::sqrt( Complex(2.0*(dc_ec+t7)/ alph_plus, 0.0) );
            if(kc.imag()<0.) kc *= -1. ;
            
            wz_0 = thrust::log( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( wz_0.imag() < 0. ) wz_0 += Complex(0.0, 2. * M_PI) ;    
            wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
#endif
            wz_0 = Zone*alph_plus / (Zone*t5 - wz_0);
            // g_eff = wz_0.real(); 
            ldelta = lnu*(-1.0*wz_0.real());
            
            // potential
            t6=(thrust::conj(ldelta)*lnu).real(); // delta^+ * nu 
            t7=thrust::norm(ldelta);
            Vanew = Va_const - t1*t6 - t3*t7;
            Vbnew = Vb_const - t2*t6 - t4*t7;
             
            // mixing of potentials
            Va = UD_MIX_COEFF*Vanew+(1.0-UD_MIX_COEFF)*Va;
            Vb = UD_MIX_COEFF*Vbnew+(1.0-UD_MIX_COEFF)*Vb;
        }
        
        // save results to global memory
        V_a[ixyz]=Va;
        V_b[ixyz]=Vb;
        delta[ixyz]=ldelta;
    }
}

#ifdef BDG_MODE
__global__ void kernel_compute_potentials_bdg(int it, 
                                          double *rho_a, double *rho_b, Complex *nu,
                                          double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                          double *tau_a, double *tau_b, 
                                          double *V_a, double *V_b, Complex *delta, double cccoeff
                                         )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    // registers
    double t5, t7; // working buffers
    
    double Va, Vb;
    Complex p0, kc, wz_0, Zone, lnu, ldelta;
    
    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        
        // load data to registers from global memory and form constant part of potentials
        Va=u_ext(ix,0,0,it,SPINA);
        Vb=u_ext(ix,0,0,it,SPINB);
        
        t5 = 1.0/ (dc_gBdG);
        lnu = nu[ixyz];
        Zone = Complex(1.0, 0.0);
  
        // pairing
#ifdef USE_CUBIC_CUTOFF
        wz_0=Complex(REGULARIZATION_SCHEME_K_CONST/(4.0*M_PI*DX), 0.0); // FIXME: account for effective mass
#else
        t7=(dc_mu_a-Va+dc_mu_b-Vb)/2.0;
        p0 = thrust::sqrt( Complex(2.0*t7, 0.0) );
        if(p0.imag()<0.) p0 *= -1. ;
        kc = thrust::sqrt( Complex(2.0*(dc_ec+t7), 0.0) );
        if(kc.imag()<0.) kc *= -1. ;
            
        wz_0 = thrust::log( ( kc + p0 ) / ( kc - p0 ) ) ;
        if ( wz_0.imag() < 0. ) wz_0 += Complex(0.0, 2. * M_PI) ;    
        wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
#endif
        wz_0 = Zone / (Zone*t5 - wz_0);
        // g_eff = wz_0.real(); 
        ldelta = lnu*(-1.0*wz_0.real());
                    
        // save results to global memory
        V_a[ixyz]=Va;
        V_b[ixyz]=Vb;
        delta[ixyz]=ldelta;
    }
}
#endif 

/**
 * Function computes potentials V_a, V_b and delta 
 * using formulas from section "9.3.2.2 Summary"
 * @param it index of time step, it is ised for proper evaluation of external potential
 * @param d_densites (INPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NX,
 *                          nu - double complex array of size NX,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NX
 *                   In total size of d_densites is 12*NX
 * @param d_potentials (INPUT/OUTPUT)
 *                     collective array with potentials [V_a, V_b, delta]
 *                     where: V_a, V_b - double arrays of size NX
 *                            delta - double complex array of size NX
 *                     In total size of d_potentials is 4*NX
 *                     NOTE: I assume that d_potentials contains potentials from previous iteration,
 *                           i.e. they are good starting point for self-consistent process.
 * @param cccoeff the current corrections coefficient
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/ 
extern "C" int compute_potentials(int it, double *d_densities, double *d_potentials, double cccoeff, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
    double *tau_a = (double *)(d_densities +  3*NX);
    double *j_a_x = (double *)(d_densities +  4*NX);
    double *j_a_y = (double *)(d_densities +  5*NX);
    double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
    double *tau_b = (double *)(d_densities +  8*NX);
    double *j_b_x = (double *)(d_densities +  9*NX);
    double *j_b_y = (double *)(d_densities + 10*NX);
    double *j_b_z = (double *)(d_densities + 11*NX);
    // pontentials
    double *V_a = (double *)(d_potentials +  0*NX);
    double *V_b = (double *)(d_potentials +  1*NX);
    Complex *delta = (Complex *)(d_potentials +  2*NX);    
    
#ifdef BDG_MODE
    kernel_compute_potentials_bdg<<<nblocks, nthreads>>>(it,
                                                     rho_a, rho_b, nu, 
                                                     j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                     tau_a, tau_b, 
                                                     V_a, V_b, delta, cccoeff
                                                    );
#else
    kernel_compute_potentials<<<nblocks, nthreads>>>(it,
                                                     rho_a, rho_b, nu, 
                                                     j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                     tau_a, tau_b, 
                                                     V_a, V_b, delta, cccoeff
                                                    );
#endif    
    return 0;
}

// =======================================================================================
// ====================================== zero_array =====================================
// =======================================================================================
__global__ void kernel_zero_array(int asize, double *array)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<asize) array[ixyz]=0.0;
}

/**
 * Fills array with zeros
 * @param asize size of array
 * @param array pointer to array
 * */
int zero_array(int asize, double *array)
{
    // number of blocks
    int nblocks = (int)ceil((float)asize/512);
    kernel_zero_array<<<nblocks, 512>>>(asize,array);
    return 0;
}

// =======================================================================================
// ================================== compute_energy =====================================
// =======================================================================================
__global__ void kernel_compute_energy_v_ext(int it, 
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b, 
                                      Complex *delta,
                                      double *E_ext
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    if(ixyz<NXY)
    {
        ix = ixyz; // decode cartesian coordinates
        
        // External potential energy
        E_ext[ixyz]=(rho_a[ixyz]*u_ext(ix,0,0,it,SPINA) + rho_b[ixyz]*u_ext(ix,0,0,it,SPINB))*DXYZ*NY*NZ;
    }
}

__global__ void kernel_compute_energy_delta_ext(int it, 
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b, 
                                      Complex *delta,
                                      double *E_ext
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    double na, nb;
    
    if(ixyz<NXY)
    {
        ix = ixyz; // decode cartesian coordinates
        
        // External potential energy
        E_ext[ixyz]=(thrust::conj(nu[ixyz])*macro_delta_ext(ix, 0, 0, it, delta[ixyz])).real()*(-1.0)*DXYZ*NY*NZ;
    }
}

__global__ void kernel_compute_energy_velocity_ext(int it, 
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b, 
                                      Complex *delta,
                                      double *E_ext
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    if(ixyz<NXY)
    {
        ix = ixyz; // decode cartesian coordinates
        
        // External potential energy
        E_ext[ixyz]=  (
                         j_a_x[ixyz]*velocity_ext(ix, 0, 0, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                        +j_b_x[ixyz]*velocity_ext(ix, 0, 0, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)
                       )*(-1.0)*DXYZ*NY*NZ;
        
    }
}

__global__ void kernel_compute_energy(int it, 
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b, 
                                      Complex *delta,
                                      double *E_kin, double *E_pot, double *E_pair, double *E_CM
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    double na, nb;
    double p;
    double taua, taub;
    double tx1, ty1, tz1, tx2, ty2, tz2;
    
    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        
        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        
        // kinetic energy
        p=polarization(na, nb);
        taua=tau_a[ixyz]; // tau_a
        taub=tau_b[ixyz]; // tau_b    
        
        // current corrections
        tx1=j_a_x[ixyz];
        ty1=0.0;
        tz1=0.0;
        tx2=j_b_x[ixyz];
        ty2=0.0;
        tz2=0.0;
        taua-=p_regularization(na)*(tx1*tx1+ty1*ty1+tz1*tz1)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(tx2*tx2+ty2*ty2+tz2*tz2)/nb; // -jb^2/nb: correction for tilde{tau}_b

        // galilean invariant contribution
        E_kin[ixyz]=0.5*(alpha_a(p)*taua + alpha_b(p)*taub)*(NY*NZ)*DXYZ;
        
        // potential energy
        E_pot[ixyz]=funD(na, nb)*(NY*NZ)*DXYZ;
        
        // pairing energy
        E_pair[ixyz]=(delta[ixyz]*thrust::conj(nu[ixyz])).real()*(-1.0)*(NY*NZ)*DXYZ;
        
        // flow energy
        E_CM[ixyz]=   p_regularization(na)*(tx1*tx1 + ty1*ty1 + tz1*tz1)/(2.*na)*(NY*NZ)*DXYZ  
                    + p_regularization(nb)*(tx2*tx2 + ty2*ty2 + tz2*tz2)/(2.*nb)*(NY*NZ)*DXYZ;  
    }
}

__global__ void kernel_compute_energy_bdg(int it, 
                                      double *rho_a, double *rho_b, Complex *nu,
                                      double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                      double *tau_a, double *tau_b, 
                                      Complex *delta,
                                      double *E_kin, double *E_pot, double *E_pair, double *E_CM
                                      )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    double na, nb;
    double taua, taub;
    double tx1, ty1, tz1, tx2, ty2, tz2;
    
    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        
        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        
        // kinetic energy
        taua=tau_a[ixyz]; // tau_a
        taub=tau_b[ixyz]; // tau_b  
        
        // current corrections
        tx1=j_a_x[ixyz];
        ty1=0.0;
        tz1=0.0;
        tx2=j_b_x[ixyz];
        ty2=0.0;
        tz2=0.0;
        taua-=p_regularization(na)*(tx1*tx1+ty1*ty1+tz1*tz1)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(tx2*tx2+ty2*ty2+tz2*tz2)/nb; // -jb^2/nb: correction for tilde{tau}_b

        // galilean invariant contribution
        E_kin[ixyz]=0.5*(taua + taub)*(NY*NZ)*DXYZ;
        
        // potential energy
        E_pot[ixyz]=0.0;
        
        // pairing energy
        E_pair[ixyz]=(delta[ixyz]*thrust::conj(nu[ixyz])).real()*(-1.0)*(NY*NZ)*DXYZ;
        
        // flow energy
        E_CM[ixyz]=   p_regularization(na)*(tx1*tx1 + ty1*ty1 + tz1*tz1)/(2.*na)*(NY*NZ)*DXYZ  
                    + p_regularization(nb)*(tx2*tx2 + ty2*ty2 + tz2*tz2)/(2.*nb)*(NY*NZ)*DXYZ;

    }
}

    
/**
 * Function computes energy and particle number
 * NOTE: this function executes cudaMemcpy, thus it is BLOCKING!
 * @param it index of time step, it is ised for proper evaluation of external potential
 * @param d_densites (INPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NX,
 *                          nu - double complex array of size NX,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NX
 *                   In total size of d_densites is 12*NX
 * @param d_potentials (INPUT)
 *                     collective array with potentials [V_a, V_b, delta]
 *                     where: V_a, V_b - double arrays of size NX
 *                            delta - double complex array of size NX
 *                     In total size of d_potentials is 4*NX
 * @param d_workarea (INPUT/OUTPUT)
 *                   working buffer of size 5*NX, 
 *                   On OUTPUT first 9 elements contain energies and particle number [E_kin, E_pot, E_pair, E_CM, E_ext, Na, Nb, Laz, Lbz]
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int compute_energy(int it, double *d_densities, double *d_potentials, double *d_workarea, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
    double *tau_a = (double *)(d_densities +  3*NX);
    double *j_a_x = (double *)(d_densities +  4*NX);
    double *j_a_y = (double *)(d_densities +  5*NX);
    double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
    double *tau_b = (double *)(d_densities +  8*NX);
    double *j_b_x = (double *)(d_densities +  9*NX);
    double *j_b_y = (double *)(d_densities + 10*NX);
    double *j_b_z = (double *)(d_densities + 11*NX);
    // pontentials
//     double *V_a = (double *)(d_potentials +  0*NX);
//     double *V_b = (double *)(d_potentials +  1*NX);
    Complex *delta = (Complex *)(d_potentials +  2*NX);  
    
    // buffers for energies
    double *E_kin = (double *)(d_workarea +  0*NX);
    double *E_pot = (double *)(d_workarea +  1*NX);
    double *E_pair= (double *)(d_workarea +  2*NX);
    double *E_CM  = (double *)(d_workarea +  3*NX);
    double *E_ext = (double *)(d_workarea +  4*NX);
    
//     // TODO - remove
//     zero_array(NX,j_a_x);
//     zero_array(NX,j_a_y);
//     zero_array(NX,j_a_z);
//     zero_array(NX,j_b_x);
//     zero_array(NX,j_b_y);
//     zero_array(NX,j_b_z);
//     zero_array(NX,j_b_z);
    
    // Step 1: prepare buffers for local reductions
#ifdef BDG_MODE
    kernel_compute_energy_bdg<<<nblocks, nthreads>>>(it,
                                                 rho_a, rho_b, nu, 
                                                 j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                 tau_a, tau_b, 
                                                 delta,
                                                 E_kin, E_pot, E_pair, E_CM
                                                 );

#else
    kernel_compute_energy<<<nblocks, nthreads>>>(it,
                                                 rho_a, rho_b, nu, 
                                                 j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                 tau_a, tau_b, 
                                                 delta,
                                                 E_kin, E_pot, E_pair, E_CM
                                                 );
#endif
    
    // Step 2: do local reductions
    int ierr, i;
    for(i=0; i<4; i++)
    {
        ierr = local_reductionR(d_workarea +  i*NX, NX, d_workarea +  i*NX, nthreads, 0);
        if(ierr!=0) return ierr;
    }
    
    // Step 3: copy data to correct elements of d_workarea
    if( cudaMemcpy( d_workarea+EKIN     , d_workarea +  0*NX , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -100;
    if( cudaMemcpy( d_workarea+EPOT     , d_workarea +  1*NX , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -101;
    if( cudaMemcpy( d_workarea+EPAIR    , d_workarea +  2*NX , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -102;
    if( cudaMemcpy( d_workarea+ECURRENT , d_workarea +  3*NX , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -103;
    
    // Step 4: external energies
#ifdef ENABLE_V_EXT
    kernel_compute_energy_v_ext<<<nblocks, nthreads>>>(it,
                                                    rho_a, rho_b, nu, 
                                                    j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                    tau_a, tau_b, 
                                                    delta,
                                                    E_ext
                                                    );
    ierr = local_reductionR(d_workarea +  4*NX, NX, d_workarea +  4*NX, nthreads, 0);
    if(ierr!=0) return ierr;
    if( cudaMemcpy( d_workarea+EPOTEXT , d_workarea +  4*NX , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -104;
#else
    cuda_set_array_elements(1, d_workarea+EPOTEXT,  0.0, 1); // set this contribution to 0.0
#endif
    
#ifdef ENABLE_DELTA_EXT
    kernel_compute_energy_delta_ext<<<nblocks, nthreads>>>(it,
                                                    rho_a, rho_b, nu, 
                                                    j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                    tau_a, tau_b, 
                                                    delta,
                                                    E_ext
                                                    );
    ierr = local_reductionR(d_workarea +  4*NX, NX, d_workarea +  4*NX, nthreads, 0);
    if(ierr!=0) return ierr;
    if( cudaMemcpy( d_workarea+EPAIREXT , d_workarea +  4*NX , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -105;
#else
    cuda_set_array_elements(1, d_workarea+EPAIREXT,  0.0, 1); // set this contribution to 0.0
#endif
    
#ifdef ENABLE_VELOCITY_EXT
    kernel_compute_energy_velocity_ext<<<nblocks, nthreads>>>(it,
                                                    rho_a, rho_b, nu, 
                                                    j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                    tau_a, tau_b, 
                                                    delta,
                                                    E_ext
                                                    );
    ierr = local_reductionR(d_workarea +  4*NX, NX, d_workarea +  4*NX, nthreads, 0);
    if(ierr!=0) return ierr;
    if( cudaMemcpy( d_workarea+EVELEXT , d_workarea +  4*NX , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -105;
#else
    cuda_set_array_elements(1, d_workarea+EVELEXT,  0.0, 1); // set this contribution to 0.0
#endif

    
    // Step 5: compute particle number
    ierr = local_reductionR(rho_a, NX, d_workarea +  NPARTA, nthreads, 0);
    if(ierr!=0) return ierr;  
    ierr = local_reductionR(rho_b, NX, d_workarea +  NPARTB, nthreads, 0);
    if(ierr!=0) return ierr;
    cuda_scale_array_elements(2, d_workarea +  NPARTA,  DXYZ*NY*NZ, 1); // add missing volume element 
    
    // Step 6: commpute angular momentum
    // no angular momentum for 1D
    
    return 0;
}

// =======================================================================================
// ================================== apply_hamiltonian ==================================
// =======================================================================================
extern "C" void *pca_cufft_work_area;
extern "C" int compute_gradient_real_f(double *f, double *df_dx, double *df_dy, double *df_dz, int nthreads);
extern "C" int compute_derivative_real_vector_f(double *fx, double *fy, double *fz, double *dfx_dx, double *dfy_dy, double *dfz_dz,int nthreads);
extern "C" int compute_laplace_real_f(double *f, double *laplace_f, int nthreads);
extern "C" int compute_divergence_real_vector_f(double *fx, double *fy, double *fz, double *divf, int nthreads);

__global__ void kernel_form_alpha_j_corr(double *rho_a, double *rho_b,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                         double *alph_a, double *alph_b,
                                         double *j_corr_a_x, double *j_corr_a_y, double *j_corr_a_z, double *j_corr_b_x, double *j_corr_b_y, double *j_corr_b_z,
                                         double cccoeff
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double aa, ab;
    double na, nb, p;
#ifdef CURRENT_CORRECTIONS
    double ja, jb;
    double fra, frb; // regularization functions
#endif
    if(ixyz<NX)
    {
        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz]; 
        
        p=polarization(na, nb);
        aa=alpha_a(p);
        ab=alpha_b(p);
        alph_a[ixyz]=aa;
        alph_b[ixyz]=ab;
        
#ifdef CURRENT_CORRECTIONS
        // current corrections
//         p=na+nb; // total density
//         fr =p_regularization(p );
        
        // Spin a component
        fra=p_regularization(na) * cccoeff;
        if(fra==0.0)
        {
            j_corr_a_x[ixyz]=0.0;
            // j_corr_a_y[ixyz]=0.0; 
            // j_corr_a_z[ixyz]=0.0;
        }
        else
        {
            // x-coordinate
            ja=j_a_x[ixyz]; 
            j_corr_a_x[ixyz] = fra*(1.-aa)*ja/na;
            // y-coordinate - no current
            // j_corr_a_y[ixyz]=0.0; 
            // z-coordinate - no current
            // j_corr_a_z[ixyz]=0.0;
        }
        
        // Spin b component
        frb=p_regularization(nb) * cccoeff;
        if(frb==0.0)
        {
            j_corr_b_x[ixyz]=0.0;
            // j_corr_b_y[ixyz]=0.0;
            // j_corr_b_z[ixyz]=0.0;
        }
        else
        {
            // x-coordinate
            jb=j_b_x[ixyz]; 
            j_corr_b_x[ixyz] = frb*(1.-ab)*jb/nb;
            // y-coordinate - no current
            // j_corr_b_y[ixyz]=0.0;
            // z-coordinate - no current
            // j_corr_b_z[ixyz]=0.0;
        }      
#endif
    }
}


__global__ void kernel_apply_hamiltonian(int it, double *rho_a, double *rho_b,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z,
                                         double *alph_a_x, double *alph_a_y, double *alph_a_z, double *alph_b_x, double *alph_b_y, double *alph_b_z, double * laplace_alpha_a, double *laplace_alpha_b,
                                         double *j_corr_a_x, double *j_corr_a_y, double *j_corr_a_z, double *j_corr_b_x, double *j_corr_b_y, double *j_corr_b_z,
                                         double *V_a, double *V_b, Complex *delta, 
                                         size_t n, Complex *wf_in, Complex *wf_out, 
                                         Complex *wf_d_dx, double *d_kky, double *d_kkz, Complex *wf_laplace, Complex *alphawf_laplace,
                                         double cccoeff,
                                         double *vx_a, double *divv_a, double *vx_b, double *divv_b
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double aa, ab;
    double na, nb, p;
    double Va, Vb;
    Complex D;
    double cja=0.0, cjb=0.0;  
    double ja, jb/*, jp*/;
#ifdef CURRENT_CORRECTIONS
    Complex gax, gbx;
    double /*fr,*/ fra, frb; // regularization functions
#endif
    
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
    int ix;
#endif
    size_t iwf;
    Complex u, v, tu, tv;
    
    if(ixyz<NX)
    {
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
        ix=ixyz;
#endif
        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        
        p=polarization(na, nb);
        aa=alpha_a(p); // will be multiplied by -1/2 later 
        ab=alpha_b(p); // will be multiplied by -1/2 later 
        
        // read potentials
        Va=V_a[ixyz];
        Vb=V_b[ixyz];
        D=delta[ixyz];
        
#ifdef ENABLE_DELTA_EXT
        D = D + macro_delta_ext(ix, 0, 0, it, D);
#endif
      
#ifdef CURRENT_CORRECTIONS
        // reset
        gax=Complex(0.0, 0.0);
        gbx=Complex(0.0, 0.0);
#endif
        
        // read gradient corrections
#if FUNCTIONAL==ASLDA
        cja+=-0.5*(j_corr_a_x[ixyz]);
        cjb+=-0.5*(j_corr_b_x[ixyz]); 
        
//         p=na+nb;
//         fr =p_regularization(p );
        
        // Spin a component
        fra=p_regularization(na) * cccoeff;
        if(fra==0.0)
        {
            gax+=Complex(0.0, 0.0);
        }
        else
        {
            // x-coordinate
            ja=j_a_x[ixyz];
            gax+=Complex(0.0,  -1.*fra*(1.-aa)*ja/na); 
            // y-coordinate - no current
            // z-coordinate - no current
        }
        
        // Spin b component
        frb=p_regularization(nb) * cccoeff;
        if(frb==0.0)
        {
            gbx+=Complex(0.0, 0.0);
        }
        else
        {
            // x-coordinate
            jb=j_b_x[ixyz];
            gbx+=Complex(0.0,   1.*frb*(1.-ab)*jb/nb); // note conjugate of complex number (beacuse of "-h*" operator)
            // y-coordinate - no current
            // z-coordinate - no current
        }
#endif

#ifdef FAST_CONST_EFFECTIVE_MASS_MODE    
        // multiply mass by -1/2
        aa*=-0.5;
        ab*=-0.5;
#else
        // multiply mass by -1/4
        aa*=-0.25;
        ab*=-0.25;

        // Read laplace of effective mass - I use na and nb as working buffers
        na=0.25*laplace_alpha_a[ixyz];
        nb=0.25*laplace_alpha_b[ixyz];
#endif
        
#ifdef ENABLE_VELOCITY_EXT
        gax+=Complex(0.0, 1.0*vx_a[ixyz]);
        gbx+=Complex(0.0,-1.0*vx_b[ixyz]); // NOTE: complex conjugate included
                        
        cja+=  0.5*divv_a[ixyz];
        cjb+=  0.5*divv_b[ixyz];
#endif
    
        // apply to each wave-function
        for(iwf=0; iwf<n; iwf++)
        {
            // reset
            u=Complex(0.0, 0.0);
            v=Complex(0.0, 0.0);
            
            // read wf
            tu=wf_in[      iwf*NX+ixyz];
            tv=wf_in[n*NX +iwf*NX+ixyz];
            
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE    
            // ja = (1/2)*alpha_a*kz^2, NOTE aa has been multiplied by -0.5 already
            ja = d_kky[iwf]; // kz, jb as working buffer   
            jb = d_kkz[iwf]; // kz, jb as working buffer
            jb = jb*jb + ja*ja;      // ky^2 + kz^2, jb as working buffer
            ja=-1.0*aa*jb;
            jb=-1.0*ab*jb;
#else
            // ja = (1/2)*alpha_a*kz^2, NOTE aa has been multiplied by -0.25 already
            ja = d_kky[iwf]; // kz, jb as working buffer   
            jb = d_kkz[iwf]; // kz, jb as working buffer
            jb = jb*jb + ja*ja;      // ky^2 + kz^2, jb as working buffer
            ja=-2.0*aa*jb;
            jb=-2.0*ab*jb;
#endif
            u+=tu*Complex(Va-dc_mu_a+ja,     cja); // V_a*u
            v-=tv*Complex(Vb-dc_mu_b+jb,-1.0*cjb); // V_b*v, note conjugate of complex number (beacuse of "-h*" operator)
            
            u+=tv*D;                // delta   * v
            v+=tu*thrust::conj(D);  // delta^* * u

#ifndef FAST_CONST_EFFECTIVE_MASS_MODE            
            // kinetic part - contribution: (1/4)*laplace(alpha)*u, see further
            u+=tu*na; // note: na keeps (1/4)*laplace(alpha_a)
            v-=tv*nb; // note: nb keeps (1/4)*laplace(alpha_b)
#endif

#ifdef CURRENT_CORRECTIONS            
            // read gradients of wf: x-coordinate
            tu=wf_d_dx[      iwf*NX+ixyz];
            tv=wf_d_dx[n*NX +iwf*NX+ixyz];
            u+=tu*gax; 
            v-=tv*gbx;
#endif

            
            // kinetic part
            // note: (-1/2)\nabla(alpha * nabla u) = -(1/4)*laplace(alpha*u) - (1/4)*alpha*laplace(u) + (1/4)*laplace(alpha)*u
            // (1/4)*laplace(alpha)*u already done
            // read laplace of wf
            tu=wf_laplace[      iwf*NX+ixyz];
            tv=wf_laplace[n*NX +iwf*NX+ixyz];          
            u+=tu*aa; // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_a is given by aa
            v-=tv*ab; // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_b is given by ab    
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
            // read laplace(alpha*u)
            u+=alphawf_laplace[      iwf*NX+ixyz]*(-0.25);
            v-=alphawf_laplace[n*NX +iwf*NX+ixyz]*(-0.25);
#endif      
            // save to global memory wf 
            wf_out[      iwf*NX+ixyz]=u;
            wf_out[n*NX +iwf*NX+ixyz]=v;
        }
    }
}

__global__ void kernel_apply_hamiltonian_bdg(int it,
                                         double *V_a, double *V_b, Complex *delta, 
                                         size_t n, Complex *wf_in, Complex *wf_out, 
                                         Complex *wf_d_dx, double *d_kky, double *d_kkz, Complex *wf_laplace,
                                         double *vx_a, double *divv_a, double *vx_b, double *divv_b
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double p,k;
    double Va, Vb;
    Complex D;

    size_t iwf;
    Complex u, v, tu, tv;
    
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
    int ix;
#endif
    
    double cja=0.0, cjb=0.0;
#ifdef CURRENT_CORRECTIONS
    Complex gax;
    Complex gbx;
#endif
    
    if(ixyz<NX)
    {
#if defined(ENABLE_DELTA_EXT) || defined(ENABLE_VELOCITY_EXT)
        ix=ixyz;
#endif       
        // read potentials
        Va=V_a[ixyz];
        Vb=V_b[ixyz];
        D=delta[ixyz];
        
#ifdef ENABLE_DELTA_EXT
        D = D + macro_delta_ext(ix, 0, 0, it, D);
#endif
        
#ifdef ENABLE_VELOCITY_EXT
        gax=Complex(0.0, 1.0*vx_a[ixyz]);
        gbx=Complex(0.0,-1.0*vx_b[ixyz]); // NOTE: complex conjugate included
                        
        cja=  0.5*divv_a[ixyz];
        cjb=  0.5*divv_b[ixyz];
#endif
                
        // apply to each wave-function
        for(iwf=0; iwf<n; iwf++)
        {
            // reset
            u=Complex(0.0, 0.0);
            v=Complex(0.0, 0.0);
            
            // read wf
            tu=wf_in[      iwf*NX+ixyz];
            tv=wf_in[n*NX +iwf*NX+ixyz];
               
            // p = (1/2)*kz^2
            k = d_kky[iwf]; // ky
            p = d_kkz[iwf]; // kz
            p = 0.5*(p*p+k*k);    // kz^2/2 + ky^2/2

            u+=tu*Complex(Va-dc_mu_a+p,       cja); // V_a*u
            v-=tv*Complex(Vb-dc_mu_b+p,  -1.0*cjb); // V_b*v, note conjugate of complex number (beacuse of "-h*" operator)
            
            u+=tv*D;                // delta   * v
            v+=tu*thrust::conj(D);  // delta^* * u

#ifdef CURRENT_CORRECTIONS            
            // read gradients of wf: x-coordinate
            tu=wf_d_dx[      iwf*NX+ixyz];
            tv=wf_d_dx[n*NX +iwf*NX+ixyz];
            u+=tu*gax; 
            v-=tv*gbx;
#endif
            
            // kinetic part
            
            // read laplace of wf
            tu=wf_laplace[      iwf*NX+ixyz];
            tv=wf_laplace[n*NX +iwf*NX+ixyz];          
            u+=tu*(-0.5); // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_a is given by aa
            v-=tv*(-0.5); // (-1/4)*alpha*laplace(u), note (-1/4)*alpha_b is given by ab    
     
            // save to global memory wf 
            wf_out[      iwf*NX+ixyz]=u;
            wf_out[n*NX +iwf*NX+ixyz]=v;
        }
    }
}

__global__ void kernel_add_quantum_friction(double *rho_a, double *rho_b, 
                                            double *djax_dx, double *djay_dy, double *djaz_dz, double *djbx_dx, double *djby_dy, double *djbz_dz,
                                            double *V_a, double *V_b, double qfalpha)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    
    int ix; 
    double coeff=0.0, r;
    if(ixyz<NX)
    {
        // see Eq.(3) in paper https://arxiv.org/abs/1305.6891
        V_a[ixyz]-=qfalpha*(djax_dx[ixyz])/dc_nF; // here I divide be reference density, to avoid problems of division by zero
        V_b[ixyz]-=qfalpha*(djbx_dx[ixyz])/dc_nF; // here I divide be reference density, to avoid problems of division by zero
    }
}

__global__ void kernel_compute_qpe(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<NX)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DX;
    }
}


__global__ void kernel_compute_qpe_norm(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *qpe_re, double *norm_re)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    Complex _wf1_u, _wf1_v;
    if(ixyz<NX)
    {
        // qpe
        _wf1_u=wf1_u[ixyz];
        _wf1_v=wf1_v[ixyz];
        p=thrust::conj(_wf1_u)*wf2_u[ixyz] + thrust::conj(_wf1_v)*wf2_v[ixyz];
        qpe_re[ixyz]=p.real()*DX;
        // norm
        norm_re[ixyz]=(thrust::norm(_wf1_u)+thrust::norm(_wf1_v))*DX;
    }
}

__global__ void kernel_subtruct_qpe(int n, Complex *wf, Complex *Hwf, double *qpe)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double e;
    int iwf;
    if(ixyz<NX)
    {
        for(iwf=0; iwf<n; iwf++)
        {       
            e=qpe[iwf];
            Hwf[ixyz+iwf*NX      ]-=wf[ixyz+iwf*NX      ]*e;
            Hwf[ixyz+iwf*NX +n*NX]-=wf[ixyz+iwf*NX +n*NX]*e;
        }
    }
}

__global__ void kernel_subtruct_qpe_norm(int n, Complex *wf, Complex *Hwf, double *qpe, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double e, c;
    int iwf;
    Complex _u, _v;
    if(ixyz<NX)
    {
        for(iwf=0; iwf<n; iwf++)
        {       
            e=qpe[iwf];
            c=1.0/sqrt(norm[iwf]);
            
            // normalize wf
            _u=wf[ixyz+iwf*NX      ]; _u*=c;
            _v=wf[ixyz+iwf*NX +n*NX]; _v*=c;
            
            //subtruct <H>, note: <c*psi|H|c*psi>=c*c*<H>
            Hwf[ixyz+iwf*NX      ]=(Hwf[ixyz+iwf*NX      ]*c - _u*e*c*c);
            Hwf[ixyz+iwf*NX +n*NX]=(Hwf[ixyz+iwf*NX +n*NX]*c - _v*e*c*c);
            
            // save normalized wave-functions
            wf[ixyz+iwf*NX      ]=_u;
            wf[ixyz+iwf*NX +n*NX]=_v;
        }
    }
}

#ifdef ENABLE_VELOCITY_EXT
__global__ void kernel_get_vector_vext(int it, int spin, double *vx)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    if(ixyz<NX)
    {
        ix=ixyz;
      
        vx[ixyz]=velocity_ext(ix, 0, 0, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        
//         if(ixyz==0) printf("TEST1: %12.6f %12.6f %12.6f %12.6f\n", 0.1*it, vx[ixyz], vy[ixyz], vz[ixyz]);
    }
}
#endif

/**
 * Function applies hamiltonian (H-<H>)*Psi.
 * NOTE: this function executes cuFFT
 * @param it  iteration number
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf_in array with wave-functions (INPUT)
 * @param wf_out array with wave-functions (OUTPUT),
 *               It can be the same as wf_d_dx, ..., wf_d_laplace, but it CANNOT be wf_in
 * @param wf_d_dx derivative with respect to dx (INPUT)
 * @param wf_d_dy derivative with respect to dy (INPUT)
 * @param wf_d_dz derivative with respect to dz (INPUT)
 * @param wf_d_laplace laplace of wave-functions (INPUT)
 * @param d_densites (INPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NX,
 *                          nu - double complex array of size NX,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NX
 *                   In total size of d_densites is 12*NX
 * @param d_potentials (INPUT)
 *                     collective array with potentials [V_a, V_b, delta]
 *                     where: V_a, V_b - double arrays of size NX
 *                            delta - double complex array of size NX
 *                     In total size of d_potentials is 4*NX
 * @param qfalpha coefficient for quantum friction term, beta coefficient in Eq.(3) in paper https://arxiv.org/abs/1305.6891,
 *                if qfalpha=0.0 then quantum friction is NOT active
 * @param useqpe array of size [n]
 *               if NULL then quasiparticle energies will be computed from wf_in,
 *               otherwise given array will be used, 
 * @param cccoeff the current corrections coefficient
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int apply_hamiltonian(int it, int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out, 
                            cufftDoubleComplex *wf_d_dx, double *d_kkyz, cufftDoubleComplex *wf_laplace, cufftDoubleComplex *alphawf_laplace,
                            double *d_densities, double *d_potentials, double qfalpha, double *useqpe, double cccoeff, 
                            int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    int ierr;
    
    // Set pointers for to simplify notation
    // densities 
//     Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
//     double *tau_a = (double *)(d_densities +  3*NX);
    double *j_a_x = (double *)(d_densities +  4*NX);
    double *j_a_y = (double *)(d_densities +  5*NX);
    double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
//     double *tau_b = (double *)(d_densities +  8*NX);
    double *j_b_x = (double *)(d_densities +  9*NX);
    double *j_b_y = (double *)(d_densities + 10*NX);
    double *j_b_z = (double *)(d_densities + 11*NX);
    
    // pontentials
    double *V_a = (double *)(d_potentials +  0*NX);
    double *V_b = (double *)(d_potentials +  1*NX);
    Complex *delta = (Complex *)(d_potentials +  2*NX);  
    
    double *d_kky = d_kkyz;
    double *d_kkz = d_kkyz+n;
    
    double * grad_alpha_a  = (double *)pca_cufft_work_area; // Use here work area of cuFFT
    double * grad_alpha_b  = grad_alpha_a + NX*3;
    double * grad_j_corr_a = grad_alpha_a + NX*6; // (j+/n+ - alpha_a*ja/na), see GW notes
    double * grad_j_corr_b = grad_alpha_a + NX*9; // (j+/n+ - alpha_b*jb/nb), see GW notes
    double * laplace_alpha_a  = grad_alpha_a + NX*12;
    double * laplace_alpha_b  = grad_alpha_b + NX*13;
    
    // Step 1: if quantum friction is active, update mean-field potentials
    if(qfalpha>0.0)
    {
//         printf("QUANTUM FRICTION ACTIVE! qfalpha=%f\n", qfalpha);
        // compute nabla*j, use grad_j_corr_a and grad_j_corr_b as temporary buffers
        ierr=compute_derivative_real_vector_f(j_a_x, NULL, NULL, grad_j_corr_a, NULL, NULL, nthreads);
        if(ierr!=0) return ierr;
        ierr=compute_derivative_real_vector_f(j_b_x, NULL, NULL, grad_j_corr_b, NULL, NULL, nthreads);
        if(ierr!=0) return ierr;  
        
        // update mean field potential by friction term
        kernel_add_quantum_friction<<<nblocks, nthreads>>>(rho_a, rho_b,
                                                    grad_j_corr_a, NULL, NULL, grad_j_corr_b, NULL, NULL,
                                                    V_a, V_b, qfalpha);        
    }
    
    double *vecvext_a    = NULL;
    double *divvext_a    = NULL;
    double *vecvext_b    = NULL;
    double *divvext_b    = NULL;   
#ifdef ENABLE_VELOCITY_EXT
    // set pointers
    vecvext_a    = (double *)(grad_alpha_a + 14*NX); // storage for keeping vext=[vx(r),vy(r),vz(r)]
    divvext_a    = (double *)(grad_alpha_a + 17*NX); // storage for keeping div(vext) 
    vecvext_b    = (double *)(grad_alpha_a + 18*NX); // storage for keeping vext=[vx(r),vy(r),vz(r)]
    divvext_b    = (double *)(grad_alpha_a + 21*NX); // storage for keeping div(vext)
    
    // fill arrays with data
    kernel_get_vector_vext<<<nblocks, nthreads>>>(it, SPINA, vecvext_a); // NOTE - only SPINA
    ierr=compute_divergence_real_vector_f(vecvext_a, NULL, NULL, divvext_a, nthreads);
    if(ierr!=0) return ierr+300;
    
    kernel_get_vector_vext<<<nblocks, nthreads>>>(it, SPINB, vecvext_b); // NOTE - only SPINB
    ierr=compute_divergence_real_vector_f(vecvext_b, NULL, NULL, divvext_b, nthreads);
    if(ierr!=0) return ierr+400;
#endif
    
#ifdef BDG_MODE
    // Step 3: apply hamiltonian
    kernel_apply_hamiltonian_bdg<<<nblocks, nthreads>>>(it,  
                                            V_a, V_b, delta, 
                                            n, (Complex *)wf_in, (Complex *)wf_out, 
                                            (Complex *)wf_d_dx, d_kky, d_kkz, (Complex *)wf_laplace,
                                            vecvext_a, divvext_a, vecvext_b, divvext_b
                                                   ); 
#else
    // Step 2: Prepare data neded for current corrections and effective mass handling
    kernel_form_alpha_j_corr<<<nblocks, nthreads>>>(rho_a, rho_b,
                                                    j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z,
                                                    grad_alpha_a, grad_alpha_b,
                                                    grad_j_corr_a, NULL, NULL, grad_j_corr_b, NULL, NULL,
                                                    cccoeff
                                                   );
    
#ifndef FAST_CONST_EFFECTIVE_MASS_MODE
    // and compute gradient and laplace of effective mass
    ierr=compute_laplace_real_f(grad_alpha_a, laplace_alpha_a, nthreads);
    if(ierr!=0) return ierr;
    ierr=compute_laplace_real_f(grad_alpha_b, laplace_alpha_b, nthreads);
    if(ierr!=0) return ierr;
#endif
    
#ifdef CURRENT_CORRECTIONS
    ierr=compute_derivative_real_vector_f(grad_j_corr_a, NULL, NULL, grad_j_corr_a, NULL, NULL, nthreads);
    if(ierr!=0) return ierr;
    ierr=compute_derivative_real_vector_f(grad_j_corr_b, NULL, NULL, grad_j_corr_b, NULL, NULL, nthreads);
    if(ierr!=0) return ierr;  
#endif
    
    // Step 3: apply hamiltonian
    kernel_apply_hamiltonian<<<nblocks, nthreads>>>(it, rho_a, rho_b,
                                            j_a_x, NULL, NULL, j_b_x, NULL, NULL,
                                            grad_alpha_a, NULL, NULL, grad_alpha_b, NULL, NULL, laplace_alpha_a, laplace_alpha_b,
                                            grad_j_corr_a, NULL, NULL, grad_j_corr_b, NULL, NULL, 
                                            V_a, V_b, delta, 
                                            n, (Complex *)wf_in, (Complex *)wf_out, 
                                            (Complex *)wf_d_dx, d_kky, d_kkz, (Complex *)wf_laplace, (Complex *)alphawf_laplace,
                                            cccoeff,
                                            vecvext_a, divvext_a, vecvext_b, divvext_b
                                                   ); 
#endif
    
    // Step 4: to increas stability - subtruct <H>*wf, where <H> is quasi particle energy
    // and normalize
    // use grad_alpha_a as working buffer
    double *gpe;
    
    if(useqpe==NULL) // compute quasiparticle energies
    {
        gpe=grad_alpha_a; // for easier notation
        size_t shift=0;
        int iwf;
        
        // compute quasi-particle energy for each wave-function
        for(iwf=0; iwf<n; iwf++) // for each wave-function
        {
            kernel_compute_qpe<<<nblocks, nthreads>>>((Complex *)wf_in+shift,        (Complex *)wf_out+shift, 
                                                      (Complex *)wf_in+shift+n*NX,   (Complex *)wf_out+shift+n*NX, 
                                                       gpe+iwf);
            
            ierr = local_reductionR(gpe+iwf, NX, gpe+iwf, nthreads, 0);
            if(ierr!=0) return ierr;    
            
            shift+=NX; // move pointer to next wf
        }

    }
    else // use given values
    {
        gpe=useqpe;
    }
    
    // subtruct <H>*wf
    kernel_subtruct_qpe<<<nblocks, nthreads>>>(n, (Complex *)wf_in, (Complex *)wf_out, gpe);  
 
    return 0;
}

// ================================================================================================
// ========================================= compute_ovelap =======================================
// ================================================================================================
__global__ void kernel_compute_ovelap(Complex *wf1_u, Complex *wf2_u, Complex *wf1_v, Complex *wf2_v, double *re, double *im)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex p;
    if(ixyz<NX)
    {
        p=thrust::conj(wf1_u[ixyz])*wf2_u[ixyz] + thrust::conj(wf1_v[ixyz])*wf2_v[ixyz];
        re[ixyz]=p.real()*DX;
        im[ixyz]=p.imag()*DX;
    }
}
/**
 * Function computes overlaps between two wave-functions (wf1,wf2)
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf1 array with wave-functions (INPUT)
 * @param wf2 array with wave-functions (INPUT)
 * @param overlap_re computed values of overlaps, real parts, array of size n (INPUT/OUTPUT)
 *                   If pointer is set as NULL on input then computation of real part is skipped.
 * @param overlap_im computed values of overlaps, real parts, array of size n (INPUT/OUTPUT) 
 *                   If pointer is set as NULL on input then computation of imaginary part is skipped.
 * @param workarea work space of size 2*NX, it can be the same workspace as used by cuFFT
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int compute_ovelap(int n, cufftDoubleComplex *wf1, cufftDoubleComplex *wf2, double *overlap_re, double *overlap_im, 
                              double *workarea, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    
    double *re = workarea;
    double *im = workarea + NX;
    int iwf;
    size_t shift=0;
    int ierr;
    
    for(iwf=0; iwf<n; iwf++) // for each wave-function
    {
        kernel_compute_ovelap<<<nblocks, nthreads>>>((Complex *)wf1+shift, (Complex *)wf2+shift, (Complex *)wf1+shift+n*NX, (Complex *)wf2+shift+n*NX, re, im);
        shift+=NX; // move pointer to next wf
        
        if(overlap_re!=NULL)
        {
            ierr = local_reductionR(re, NX, re, nthreads, 0);
            if(ierr!=0) return ierr;
            if( cudaMemcpy( overlap_re+iwf , re , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -333;
        }
        
        if(overlap_im!=NULL)
        {
            ierr = local_reductionR(im, NX, im, nthreads, 0);
            if(ierr!=0) return ierr;
            if( cudaMemcpy( overlap_im+iwf , im , sizeof(double), cudaMemcpyDeviceToDevice )!= cudaSuccess ) return -334;            
        }        
    }
    return 0;
}

// =======================================================================================
// ====================================== amb_step1 ======================================
// =======================================================================================
__global__ void kernel_amb_step1(size_t n, Complex *ykm1, 
                         Complex *fkm1, Complex *fkm2, Complex *fkm3)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm1, _fkm2, _fkm3;
    Complex dti;
    
    if(ixyz<n*2*NX)
    {    
        // read data
        _ykm1=ykm1[ixyz];
        _fkm1=fkm1[ixyz];
        _fkm2=fkm2[ixyz];
        _fkm3=fkm3[ixyz];
        
        dti = Complex(0.0, -1.0*dc_dt);
        
        // save data
        ykm1[ixyz] = _ykm1 + dti*(_fkm1*(23./12.) - _fkm2*(16./12.) + _fkm3*(5./12.));
        fkm3[ixyz] = _ykm1 + dti*(_fkm1*(19./24.) - _fkm2*( 5./24.) + _fkm3*(1./24.));
        
    }
}

/**
 * Functions perform step 1 from intgration.pdf
 * @param n  number of wave-functions (u,v pairs) to process
 * @param ykm1 array y_{k-1} of size 2*n*NX
 * @param fkm1 array f_{k-1} of size 2*n*NX
 * @param fkm2 array f_{k-2} of size 2*n*NX
 * @param fkm3 array f_{k-3} of size 2*n*NX
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int amb_step1(int n, cufftDoubleComplex *ykm1, 
                         cufftDoubleComplex *fkm1, cufftDoubleComplex *fkm2, cufftDoubleComplex *fkm3, 
                         int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)2*NX*n/nthreads);
    kernel_amb_step1<<<nblocks, nthreads>>>(n, (Complex *)ykm1, (Complex *)fkm1, (Complex *)fkm2, (Complex *)fkm3);
    
    return 0;
}

// =======================================================================================
// ====================================== amb45_step1 ======================================
// =======================================================================================
__global__ void kernel_amb45_step1(size_t n, Complex *ykm1, 
                         Complex *fkm1, Complex *fkm2, Complex *fkm3, Complex *fkm4)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm1, _fkm2, _fkm3, _fkm4;
    Complex dti;
    
    if(ixyz<n*2*NX)
    {    
        // read data
        _ykm1=ykm1[ixyz];
        _fkm1=fkm1[ixyz];
        _fkm2=fkm2[ixyz];
        _fkm3=fkm3[ixyz];
        _fkm4=fkm4[ixyz];
        
        dti = Complex(0.0, -1.0*dc_dt);
        
        // save data
        ykm1[ixyz] = _ykm1 + dti*(_fkm1*(55./24.  ) - _fkm2*(59./24.  ) + _fkm3*(37./24.  ) - _fkm4*(9./24.  ));
        fkm4[ixyz] = _ykm1 + dti*(_fkm1*(646./720.) - _fkm2*(264./720.) + _fkm3*(106./720.) - _fkm4*(19./720.));
    }
}

/**
 * Functions perform step 1 from intgration.pdf for algorithm AB4AM5
 * @param n  number of wave-functions (u,v pairs) to process
 * @param ykm1 array y_{k-1} of size 2*n*NX
 * @param fkm1 array f_{k-1} of size 2*n*NX
 * @param fkm2 array f_{k-2} of size 2*n*NX
 * @param fkm3 array f_{k-3} of size 2*n*NX
 * @param fkm4 array f_{k-4} of size 2*n*NX
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int amb45_step1(int n, cufftDoubleComplex *ykm1, 
                         cufftDoubleComplex *fkm1, cufftDoubleComplex *fkm2, cufftDoubleComplex *fkm3, cufftDoubleComplex *fkm4, 
                         int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)2*NX*n/nthreads);
    kernel_amb45_step1<<<nblocks, nthreads>>>(n, (Complex *)ykm1, (Complex *)fkm1, (Complex *)fkm2, (Complex *)fkm3, (Complex *)fkm4);
    
    return 0;
}

// =======================================================================================
// ====================================== amb_step4 ======================================
// =======================================================================================
__global__ void kernel_amb_step4(size_t n, Complex *ykm1_in, Complex *ykm1_out, Complex *fkm3)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm3;
    Complex dti;
    
    if(ixyz<n*2*NX)
    {    
        // read data
        _ykm1=ykm1_in[ixyz];
        _fkm3=fkm3[ixyz];
        
        dti = Complex(0.0, -1.0*dc_dt);
        
        // save data
        ykm1_out[ixyz] = dti*_ykm1*(9./24.) + _fkm3;
    }
}

/**
 * Functions perform step 4 from intgration.pdf
 * @param n  number of wave-functions (u,v pairs) to process
 * @param ykm1_in array y_{k-1} of size 2*n*NX (INPUT)
 * @param ykm1_out array y_{k-1} of size 2*n*NX (OUTPUT)
 * @param fkm3 array f_{k-3} of size 2*n*NX
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int amb_step4(int n, cufftDoubleComplex *ykm1_in, cufftDoubleComplex *ykm1_out, 
                         cufftDoubleComplex *fkm3, 
                         int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)2*NX*n/nthreads);
    kernel_amb_step4<<<nblocks, nthreads>>>(n, (Complex *)ykm1_in, (Complex *)ykm1_out, (Complex *)fkm3);
    
    return 0;
}

// =======================================================================================
// ====================================== amb45_step4 ======================================
// =======================================================================================
__global__ void kernel_amb45_step4(size_t n, Complex *ykm1_in, Complex *ykm1_out, Complex *fkm4)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    Complex _ykm1, _fkm4;
    Complex dti;
    
    if(ixyz<n*2*NX)
    {    
        // read data
        _ykm1=ykm1_in[ixyz];
        _fkm4=fkm4[ixyz];
        
        dti = Complex(0.0, -1.0*dc_dt);
        
        // save data
        ykm1_out[ixyz] = dti*_ykm1*(251./720.) + _fkm4;
    }
}

/**
 * Functions perform step 4 from intgration.pdf for algorithm AB4AM5
 * @param n  number of wave-functions (u,v pairs) to process
 * @param ykm1_in array y_{k-1} of size 2*n*NX (INPUT)
 * @param ykm1_out array y_{k-1} of size 2*n*NX (OUTPUT)
 * @param fkm4 array f_{k-3} of size 2*n*NX
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int amb45_step4(int n, cufftDoubleComplex *ykm1_in, cufftDoubleComplex *ykm1_out, 
                         cufftDoubleComplex *fkm4, 
                         int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)2*NX*n/nthreads);
    kernel_amb45_step4<<<nblocks, nthreads>>>(n, (Complex *)ykm1_in, (Complex *)ykm1_out, (Complex *)fkm4);
    
    return 0;
}

// =======================================================================================
// ==================================== normalize_wf =====================================
// =======================================================================================
__global__ void kernel_compute_norm(Complex *wf_u, Complex *wf_v, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<NX)
    {
        norm[ixyz]=(thrust::norm(wf_u[ixyz])+thrust::norm(wf_v[ixyz]))*DX;
    }
}

__global__ void kernel_normalize_wf(int n, Complex *wf, double *norm)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    double c;
    int iwf;
    if(ixyz<NX)
    {
        for(iwf=0; iwf<n; iwf++)
        {       
            c=1.0/sqrt(norm[iwf]);
            
            // normalize wf
            wf[ixyz+iwf*NX     ]*=c;
            wf[ixyz+iwf*NX+n*NX]*=c;
        }
    }
}


/**
 * Function normalizes wave-functions
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf pointer to wave-functions
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int normalize_wf(int n, cufftDoubleComplex *wf, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    int ierr;
    
    double * norm  = (double *)pca_cufft_work_area; // Use here work area of cuFFT as working buffer

    size_t shift=0;
    int iwf;
    
    // compute norm for each wave-function
    for(iwf=0; iwf<n; iwf++) // for each wave-function
    {
        kernel_compute_norm<<<nblocks, nthreads>>>((Complex *)wf+shift, (Complex *)wf+shift+n*NX, norm+iwf);
        
        ierr = local_reductionR(norm+iwf, NX, norm+iwf, nthreads, 0);
        if(ierr!=0) return ierr;    
                
        shift+=NX; // move pointer to next wf
        
//         // TEST
//         double tt;
//         if( cudaMemcpy( &tt , norm+iwf , sizeof(double), cudaMemcpyDeviceToHost )!= cudaSuccess ) return -99;
//         printf("iwf=%d norm=%f\n", iwf, tt);
        
    }
    
    // normalize wf
    kernel_normalize_wf<<<nblocks, nthreads>>>(n, (Complex *)wf, norm); 

    return 0;
    
}

// =======================================================================================
// ============================== multiply_wf_by_alpha ===================================
// =======================================================================================
__global__ void kernel_multiply_wf_by_alpha(int n, double * rho_a, double * rho_b, Complex *wf_in, Complex *wf_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    
    // registers
    int iwf;
    double na, nb, p;
    Complex u, v;
        
    if(ixyz<NX)
    {
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        p = polarization(na, nb);
        na=alpha_a(p); // na as working buffer
        nb=alpha_b(p); // nb as working buffer
        
        for(iwf=0; iwf<n; iwf++)
        {    
            // read u and v
            u=wf_in[       iwf*NX+ixyz];
            v=wf_in[n*NX +iwf*NX+ixyz];
            
            wf_out[       iwf*NX+ixyz]=u*na; // multiply by effective mass and save
            wf_out[n*NX +iwf*NX+ixyz]=v*nb; // multiply by effective mass and save
        }
    }
}

/**
 * Function applies hamiltonian.
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf_in array with wave-functions (INPUT)
 * @param wf_out array with wave-functions (OUTPUT),
 * @param d_densites (INPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NX,
 *                          nu - double complex array of size NX,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NX
 *                   In total size of d_densites is 12*NX
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/ 
extern "C" int multiply_wf_by_alpha(int n, cufftDoubleComplex *wf_in, cufftDoubleComplex *wf_out, double *d_densities, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    Complex *nu   =(Complex *)(d_densities +  0*NX);
    double *rho_a = (double *)(d_densities +  2*NX);
//     double *tau_a = (double *)(d_densities +  3*NX);
//     double *j_a_x = (double *)(d_densities +  4*NX);
//     double *j_a_y = (double *)(d_densities +  5*NX);
//     double *j_a_z = (double *)(d_densities +  6*NX);
    double *rho_b = (double *)(d_densities +  7*NX);
//     double *tau_b = (double *)(d_densities +  8*NX);
//     double *j_b_x = (double *)(d_densities +  9*NX);
//     double *j_b_y = (double *)(d_densities + 10*NX);
//     double *j_b_z = (double *)(d_densities + 11*NX);
    
    kernel_multiply_wf_by_alpha<<<nblocks, nthreads>>>(n, rho_a, rho_b, (Complex *)wf_in, (Complex *)wf_out);
    
    return 0;
}


// =======================================================================================
// ======================== taylor_expansion_contribution ================================
// =======================================================================================
__global__ void kernel_taylor_expansion_contribution(int n, Complex *wf_hpsi, Complex *wf_update, Complex *wf_contr, double t_coeff)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    
    // registers
    int iwf;

    Complex u, v;
        
    if(ixyz<NX)
    {

        for(iwf=0; iwf<n; iwf++)
        {    
            // read u and v
            u=wf_hpsi[       iwf*NX+ixyz];
            v=wf_hpsi[n*NX +iwf*NX+ixyz];
            
            // add Taylor coefficients
            u*=Complex(0.0, t_coeff);
            v*=Complex(0.0, t_coeff);
            
            
            wf_contr[       iwf*NX+ixyz]=u; // save taylor contribution
            wf_contr[n*NX +iwf*NX+ixyz]=v; // save taylor contribution
            
            wf_update[       iwf*NX+ixyz]+=u; // add to rest of Taylor expansion
            wf_update[n*NX +iwf*NX+ixyz]+=v; // add to rest of Taylor expansion            
        }
    }
}

/**
 * Function add it-order conytribution from Taylor expansion: (1/it!)(-i*(H-<H>)dt)*psi
 * @param it Taylor order contribution
 * @param dt time step
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf_hpsi result of apply_hamiltonian (INPUT)
 * @param wf_update array with wave-functions contructed by Taylor expansion (OUTPUT),
 * @param wf_contr array with wave-functions contructed by Taylor expansion (OUTPUT),
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR
 **/ 
extern "C" int taylor_expansion_contribution(int it, double dt, int n, cufftDoubleComplex *wf_hpsi, 
                                             cufftDoubleComplex *wf_update, cufftDoubleComplex *wf_contr, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NX/nthreads);
    
    kernel_taylor_expansion_contribution<<<nblocks, nthreads>>>(n, (Complex *)wf_hpsi, (Complex *)wf_update, (Complex *)wf_contr, 
                                                                -1.0*dt/(double)(it));
    
    return 0;
}


// =======================================================================================
// ======================== get_ext_potentials - handlers ================================
// =======================================================================================

__global__ void kernel_get_v_ext(int it, int spin, double *data)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
        
    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        data[ixyz]=u_ext(ix,0,0,it,spin);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_v_ext(int datadim, int spin, int it, double *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NX/nthreads);
    
    // allocate cuda memory
    double *wrkspace = (double *)pca_cufft_work_area; // reuse workspace
    kernel_get_v_ext<<<nblocks, nthreads>>>(it, spin, wrkspace);
    
    if( cudaMemcpy( data , wrkspace, sizeof(double)*NX, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;
    
    return 0;
}

__global__ void kernel_get_delta_ext(int it, Complex *deltain, Complex *data)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
        
    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        data[ixyz]=macro_delta_ext(ix, 0, 0, it, deltain[ixyz]);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_delta_ext(int datadim, int it, void *deltain, void *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NX/nthreads);
    
    // allocate cuda memory
    Complex *wrkspace = (Complex *)pca_cufft_work_area; // reuse workspace
    kernel_get_delta_ext<<<nblocks, nthreads>>>(it, (Complex *)deltain, (Complex *)wrkspace);
    
    if( cudaMemcpy( data , wrkspace, sizeof(Complex)*NX, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;
    
    return 0;
}

__global__ void kernel_get_velocity_ext(int it, int spin, double *datax, double *datay, double *dataz)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
        
    if(ixyz<NX)
    {
        ix = ixyz; // decode cartesian coordinates
        datax[ixyz]=velocity_ext(ix, 0, 0, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        datay[ixyz]=velocity_ext(ix, 0, 0, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        dataz[ixyz]=velocity_ext(ix, 0, 0, it, spin, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
    }
}

/**
 * @return array with external potential
 * */
extern "C" int get_velocity_ext(int datadim, int spin, int it, double *data)
{
    int nthreads = 256;
    int nblocks = (int)ceil((float)NX/nthreads);
    
    double *wrkspace = (double *)pca_cufft_work_area; // reuse workspace
    double *vx = wrkspace + 0*NX;
    double *vy = wrkspace + 1*NX;
    double *vz = wrkspace + 2*NX;
    
    kernel_get_velocity_ext<<<nblocks, nthreads>>>(it, spin, vx, vy, vz);
    
    if( cudaMemcpy( data , wrkspace, sizeof(double)*NX*3, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;    
    
    return 0;
}
