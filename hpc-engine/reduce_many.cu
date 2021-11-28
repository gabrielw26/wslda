// Author: Stanisław Tabisz
// Date: 25-11-2021
// Version: 1.6

#include <iostream>

__device__ int getEndValue(int currArrayIdx, int nElementsInVector, int nMemoryBlocksPerArray, int currMemoryBlockIdx) {
    int idxLastMemoryBlockInCurrArray = (currArrayIdx+1)*nMemoryBlocksPerArray - 1;
    int nthreadsUsedInLastMemoryBlock = nElementsInVector - (nMemoryBlocksPerArray-1)*blockDim.x;
    return currMemoryBlockIdx < idxLastMemoryBlockInCurrArray ? blockDim.x : nthreadsUsedInLastMemoryBlock;
}

template <unsigned int blockSize>
__device__ void warpReduce(volatile double *shm, int tid, int nEndValue) {
    if (blockSize >= 64 && tid + 32 < nEndValue) shm[tid] += shm[tid + 32];
    if (blockSize >= 32 && tid + 16 < nEndValue) shm[tid] += shm[tid + 16];
    if (blockSize >= 16 && tid + 8  < nEndValue) shm[tid] += shm[tid +  8];
    if (blockSize >=  8 && tid + 4  < nEndValue) shm[tid] += shm[tid +  4];
    if (blockSize >=  4 && tid + 2  < nEndValue) shm[tid] += shm[tid +  2];
    if (blockSize >=  2 && tid + 1  < nEndValue) shm[tid] += shm[tid +  1]; 
}

template <unsigned int blockSize>
__global__ void reduce_many(double *g_idata, double *g_odata, 
    int nElementsInVector, int nblocksPerArray, int nMemoryBlocksPerArray) {
    extern __shared__ double sdata[];

    int tid = threadIdx.x;
    int currArrayIdx = blockIdx.x / nblocksPerArray;
    int currBlockIdx = blockIdx.x - currArrayIdx*nblocksPerArray;
    int currMemoryBlockIdx = currArrayIdx*nMemoryBlocksPerArray + currBlockIdx;
    int currMemoryThreadIdx = currMemoryBlockIdx*blockDim.x + threadIdx.x;
    int shift = nMemoryBlocksPerArray*blockDim.x - nElementsInVector;
	int idxInGlobalMemory = currMemoryThreadIdx - (shift*currArrayIdx);
    int nEndValue = getEndValue(currArrayIdx, nElementsInVector, nMemoryBlocksPerArray, currMemoryBlockIdx);
    int numThreadsPerGrid = blockSize*nblocksPerArray;
    int i = currBlockIdx*blockDim.x + tid;
    sdata[tid] = 0;

    // multiple adds per thread
    for(int processedNum=0; processedNum < nElementsInVector; processedNum+=numThreadsPerGrid) {
        if(i + processedNum < nElementsInVector) {
            sdata[tid] += g_idata[idxInGlobalMemory + processedNum];
        }
    }
    __syncthreads();

    // reduce
    if (blockSize >= 512) { if (tid < 256 && tid + 256 < nEndValue) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128 && tid + 128 < nEndValue) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid < 64  && tid + 64  < nEndValue) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }  
    if (tid < 32)
        warpReduce<blockSize>(sdata, tid, nEndValue);

    // write result for this block to global mem
    if (tid == 0) 
        g_odata[currArrayIdx*nblocksPerArray + currBlockIdx] = sdata[0];
}

void call_template_kernels(int dimGrid, int dimBlock, double *d_idata, double *d_odata, 
    int nElementsInVector, int nblocksPerArray, int nMemoryBlocksPerArray) {
    int smemSize=dimBlock*sizeof(double);
    switch (dimBlock)
    {
        case 1024:
            reduce_many<1024><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 512:
            reduce_many< 512><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 256:
            reduce_many< 256><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 128:
            reduce_many< 128><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 64:
            reduce_many<  64><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 32:
            reduce_many<  32><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 16:
            reduce_many<  16><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 8:
            reduce_many<   8><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 4:
            reduce_many<   4><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 2:
            reduce_many<   2><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;
        case 1:
            reduce_many<   1><<< dimGrid, dimBlock, smemSize >>>(d_idata, d_odata, nElementsInVector, nblocksPerArray, nMemoryBlocksPerArray); break;  
    }   
}

int getNumberOfGridDim(int m, int optNumOperation, int blockDim) {
    return (int)ceil((float)m/(2*blockDim*optNumOperation));
}

int getNumberOfBlockDim(int m, int optNumOperationPerThread) {
    int nThreads = (int)ceil((float)m/optNumOperationPerThread);
    
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
 * @param n number of vectors to be reduced -> vectorNo
 * @param m length of the vector            -> indexNo
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

    cudaError_t error = cudaPeekAtLastError();
    if(error != cudaSuccess) 
        return -1;

    int dimBlock=getNumberOfBlockDim(m, optNumOperationPerThread);
    int nMemoryBlocksPerArray = (int)ceil((float)m/(dimBlock));
    int nblocksPerArray = (int)ceil((float)nMemoryBlocksPerArray/optNumOperationPerThread);
    int nblocksTotal = nblocksPerArray * n;

    // Do first reduction
    call_template_kernels(nblocksTotal, dimBlock, d_in, d_out, m, nblocksPerArray, nMemoryBlocksPerArray);
    error = cudaGetLastError();
    if(error != cudaSuccess) 
        return error;

    // Do iteratively reduction of partial_sums
    while(nblocksPerArray>1) {
        m = nblocksPerArray;
        dimBlock = getNumberOfBlockDim(m, optNumOperationPerThread);
        nMemoryBlocksPerArray = (int)ceil((float)m/(dimBlock));
        nblocksPerArray = (int)ceil((float)nMemoryBlocksPerArray/optNumOperationPerThread);
        nblocksTotal = nblocksPerArray * n;
        call_template_kernels(nblocksTotal, dimBlock, d_out, d_out, m, nblocksPerArray, nMemoryBlocksPerArray);
        error = cudaGetLastError();
        if(error != cudaSuccess) 
            return error;
    }     
    return 0;
}