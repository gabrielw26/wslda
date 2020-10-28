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

