/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * The file implements varius utils function in CUDA
 *  
 * @author Gabriel Wlazlowski
 * @date 25.10.2020
 * */  
#ifndef __WSLDA_CUDA_UTILS__
#define __WSLDA_CUDA_UTILS__

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
__constant__ double dc_sclgth;

__constant__ void *dc_extra_data;
__constant__ size_t dc_extra_data_size;

__constant__ double dc_params[MAX_USER_PARAMS]; // array with params from input file

extern "C" int memcopy_extra_data(size_t extra_data_size, void *extra_data)
{
    if( cudaMemcpyToSymbol(dc_extra_data,      &extra_data,      sizeof(void *))!= cudaSuccess ) return 1;
    if( cudaMemcpyToSymbol(dc_extra_data_size, &extra_data_size, sizeof(size_t))!= cudaSuccess ) return 2;
    return 0;
}

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

/**
 * This function copies BdG functional data
 * */
extern "C" int memcopy_const_BdG(double aBdG)
{
    double gBdG = aBdG;
    if( cudaMemcpyToSymbol(dc_sclgth, &gBdG, sizeof(double))!= cudaSuccess ) return 1;

    return 0;
}

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


// ===========================================================================
// ============================ AUXLIARY FUNCTIONS ===========================
// ===========================================================================
__global__ void kernel_cuda_set_array_elements(int array_dim, double *array,  double value)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<array_dim) array[ixyz]=value;
}

/**
 * Function sets array elements to desired value.
 * @param array_dim size of array
 * @param array pointer to array (device)
 * @param value elements will be set to this value
 * @param nthreads number of threads to be used for kernel invocation
 * @return 0-ok
 * */
int cuda_set_array_elements(int array_dim, double *array,  double value, int nthreads)
{
 
    // number of blocks
    int nblocks = (int)ceil((float)array_dim/nthreads);
    if(nblocks<1) nblocks=1;
    kernel_cuda_set_array_elements<<<nblocks, nthreads>>>(array_dim, array, value);
    return 0;
}


__global__ void kernel_cuda_scale_array_elements(int array_dim, double *array,  double value)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    if(ixyz<array_dim) array[ixyz]*=value;
}

/**
 * Function sets array elements to desired value.
 * @param array_dim size of array
 * @param array pointer to array (device)
 * @param value elements will be scaled by this value
 * @param nthreads number of threads to be used for kernel invocation
 * @return 0-ok
 * */
int cuda_scale_array_elements(int array_dim, double *array,  double value, int nthreads)
{
 
    // number of blocks
    int nblocks = (int)ceil((float)array_dim/nthreads);
    if(nblocks<1) nblocks=1;
    kernel_cuda_scale_array_elements<<<nblocks, nthreads>>>(array_dim, array, value);
    return 0;
}

#endif
