// Author: Stanisław Tabisz
// Date: 20-12-2021
// Version: 1.9

#include <iostream>
#include <cstdio>

template <unsigned int blockSize>
__device__ void warpReduce(volatile double *shm, int tid) {
    if (blockSize >= 64) shm[tid] += shm[tid + 32];
    if (blockSize >= 32) shm[tid] += shm[tid + 16];
    if (blockSize >= 16) shm[tid] += shm[tid +  8];
    if (blockSize >=  8) shm[tid] += shm[tid +  4];
    if (blockSize >=  4) shm[tid] += shm[tid +  2];
    if (blockSize >=  2) shm[tid] += shm[tid +  1]; 
}

template <unsigned int blockSize>
__global__ void reduce_many(double *g_idata, double *g_odata, int distBetweenSums, int nElementsInVector, int nBlocksPerArray, int endVal) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int currArrayIdx = blockIdx.x / nBlocksPerArray;
    int currBlockIdx = blockIdx.x - currArrayIdx*nBlocksPerArray;
    int idxInArray = currBlockIdx*blockDim.x*distBetweenSums + tid*distBetweenSums;
    int idxInGlobalMemory = currArrayIdx*nElementsInVector + idxInArray;
    int gridStride = blockDim.x*nBlocksPerArray*distBetweenSums;    
    sdata[tid] = 0;

    // multiple adds per thread
    for(int positionInArray=0; positionInArray < endVal; positionInArray+=gridStride) {
        if(idxInArray + positionInArray < endVal) {
            sdata[tid] += g_idata[idxInGlobalMemory + positionInArray];
        }
    }
    __syncthreads();

    // reduce
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }  
    if (tid < 32) 
        warpReduce<blockSize>(sdata, tid);

    // write result for this block to global mem
    if (tid == 0) 
       g_odata[idxInGlobalMemory] = sdata[0]; 
}

__global__ void reduce_many_32(double *g_idata, double *g_odata, int distBetweenSums, int nElementsInVector, int nBlocksPerArray, int endVal) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int currArrayIdx = blockIdx.x / nBlocksPerArray;
    int currBlockIdx = blockIdx.x - currArrayIdx*nBlocksPerArray;
    int idxInArray = currBlockIdx*blockDim.x*distBetweenSums + tid*distBetweenSums;
    int idxInGlobalMemory = currArrayIdx*nElementsInVector + idxInArray;
    int gridStride = blockDim.x*nBlocksPerArray*distBetweenSums;    
    sdata[tid] = 0;

    // multiple adds per thread
    for(int positionInArray=0; positionInArray < endVal; positionInArray+=gridStride) {
        if(idxInArray + positionInArray < endVal) {
            sdata[tid] += g_idata[idxInGlobalMemory + positionInArray];
        }
    }
    __syncthreads();

    // reduce for block == warp
    if (tid < 16) 
        warpReduce<32>(sdata, tid);

    // write result for this block to global mem
    if (tid == 0) 
       g_odata[idxInGlobalMemory] = sdata[0]; 
}

void call_template_reduce_many_kernels(int dimGrid, int dimBlock, double *d_idata, double *d_odata, 
    int distBetweenSums, int nElementsInVector, int nBlocksPerArray, int endVal) {
    int smemSize=dimBlock*sizeof(double);
    switch (dimBlock)
    {
        case 1024:
            reduce_many<1024><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal); break;
        case 512:
            reduce_many< 512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal); break;
        case 256:
            reduce_many< 256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal); break;
        case 128:
            reduce_many< 128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal); break;
        case 64:
            reduce_many<  64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal); break;
        case 32:
            reduce_many_32<<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal); break;
    }   
}

__global__ void write_at_begining(double *g_odata, int nVectors, int nElementsInVector, int idxFirstArrayToRewrite) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int currArrayIdx = idxFirstArrayToRewrite + tid;
    
    if(currArrayIdx < nVectors) {
        sdata[tid] = g_odata[currArrayIdx*nElementsInVector];
    }

    __syncthreads();

    if(currArrayIdx < nVectors) {
        g_odata[currArrayIdx] = sdata[tid];
    }
}

int getNumberOfBlockDim(int nElementsToSum, int optNumOperationPerThread) {
    int nThreads = (int)ceil((float)nElementsToSum/optNumOperationPerThread);
    
    if(nThreads <= 32) 
        return 32;
    else if(nThreads <= 64)
        return 64;
    else if(nThreads <= 128) 
        return 128;
    else if(nThreads <= 256) 
        return 256;
    else 
        return 512;
}

/**
 * Function reduces many arrays.
 * Mathematical operation: sum (+)
 * @param n number of vectors to be reduced
 * @param m length of the vector
 * @param d_in = A pointer to memory (device) of size n*m*sizeof(double), (INPUT)
 *          i-th element of k-th vector is located A_i^k=A[k*m + I],
 *          NOTE: the array can be overwritten by the computation process, 
 *                i.e. it can be used as working space
 * @param d_out = r pointer to memory (device) of size n*sizeof(double) (OUTPUT)
 *          where reductions will be stored, namely 
 *           r[k] = sum_i^m  A_i^{k}
 *          NOTE: it is allowed to provide as pointer r address of input array A.  
 * @return 0 – ok, -1 - error before first calling kernel, >0 - otherwise error code.
 **/
int local_reductions_many(int n, int m, double *d_in, double *d_out) {
    const int optNumOperationPerThread = 1024;
    const int nElementsInVector = m;
    const int nVectors = n;

    cudaError_t error = cudaPeekAtLastError();
    if(error != cudaSuccess) 
        return -1;

    // Do first reduction
    int nElementsToSum = m;
    int dimBlock= getNumberOfBlockDim(nElementsToSum, optNumOperationPerThread);
    int nMemoryBlocksPerArray = (int)ceil((float)nElementsToSum/(dimBlock));
    int nBlocksPerArray = (int)ceil((float)nMemoryBlocksPerArray/optNumOperationPerThread);
    int nBlocksTotal = nBlocksPerArray * nVectors;
    int endVal = nElementsInVector;
    int distBetweenSums = 1;    
    call_template_reduce_many_kernels(nBlocksTotal, dimBlock, d_in, d_out, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal);
    error = cudaGetLastError();
    if(error != cudaSuccess) 
        return error;
    
    // Do iteratively reduction of partial sums
    while(nBlocksPerArray>1) 
    {
        nElementsToSum = nBlocksPerArray;
        distBetweenSums *= dimBlock; 
        endVal = nElementsToSum*distBetweenSums;               
        dimBlock = getNumberOfBlockDim(nElementsToSum, optNumOperationPerThread);
        nMemoryBlocksPerArray = (int)ceil((float)nElementsToSum/(dimBlock));
        nBlocksPerArray = (int)ceil((float)nMemoryBlocksPerArray/optNumOperationPerThread);
        nBlocksTotal = nBlocksPerArray * nVectors;        
        call_template_reduce_many_kernels(nBlocksTotal, dimBlock, d_out, d_out, distBetweenSums, nElementsInVector, nBlocksPerArray, endVal);
        error = cudaGetLastError();
        if(error != cudaSuccess) 
            return error;
    }
    
    // Write at the begining of output
    int maxThreadsPerBlock = 512;
    for(int idxFirstArrayToRewrite=0; idxFirstArrayToRewrite<nVectors; idxFirstArrayToRewrite+=maxThreadsPerBlock) {
        int smemSize=maxThreadsPerBlock*sizeof(double);
        write_at_begining<<< 1, maxThreadsPerBlock, smemSize>>>(d_out, nVectors, nElementsInVector, idxFirstArrayToRewrite);
    }
    
    return 0;
}