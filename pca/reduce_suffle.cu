#include <cuda.h> 
#include <stdio.h>

__inline__ __device__
double warpReduceSum(double val) {
  for (int offset = warpSize/2; offset > 0; offset /= 2) 
    val += __shfl_down(val, offset);
  return val;
}

__inline__ __device__
double blockReduceSum(double val) {

  static __shared__ double shared[32]; // Shared mem for 32 partial sums
  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;

  val = warpReduceSum(val);     // Each warp performs partial reduction

  if (lane==0) shared[wid]=val; // Write reduced value to shared memory

  __syncthreads();              // Wait for all partial reductions

  //read from shared memory only if that warp existed
  val = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  
  if (wid==0) val = warpReduceSum(val); //Final reduce within first warp
  
  return val;
}

__global__ void deviceReduceKernel(double* in, double* out, int N) {
  double sum = 0;
  //reduce multiple elements per thread
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
       		i < N; 
       		i += blockDim.x * gridDim.x) {
       sum += in[i];
  }
  sum = blockReduceSum(sum);

  if (threadIdx.x==0)
    out[blockIdx.x]=sum;
}



void deviceReduce_suffle(double *in, double* out, int N, int threads) {
//  int threads = 512;
  int blocks = min((N + threads - 1) / threads, 512);

  deviceReduceKernel<<<blocks, threads>>>(in, out, N);
if(blocks>1);
  deviceReduceKernel<<<1, 512>>>(out, out, blocks);
}

void deviceReduce_suffle_st(double *in, double* out, int N, int threads, cudaStream_t stream) {
  int blocks = min((N + threads - 1) / threads, 512);

  deviceReduceKernel<<<blocks, threads, 0, stream>>>(in, out, N);
if(blocks>1);
  deviceReduceKernel<<<1, 512, 0, stream>>>(out, out, blocks);
}

__global__ void deviceReduceKernel_MANY(double* in, double* out, int N, int M) {
  double sum ;
  for (int m = 0; m<M; m++){
  sum = 0;
//reduce multiple elements per thread
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
      	i < N;
      	i += blockDim.x * gridDim.x) {
      sum += in[i+N*m];
      }
  sum = blockReduceSum(sum);

  if (threadIdx.x==0)
  out[blockIdx.x+m*gridDim.x]=sum;

//if (threadIdx.x==0) printf("%d   BId=%d  %d  m=%d N=%d  %24.22f\n",blockIdx.x+m*gridDim.x, blockIdx.x, gridDim.x, m, N, sum);

  }//m
}

void deviceReduce_suffle_MANY_st(double *in, double* out, int N, int threads, cudaStream_t stream, int M, double* work) {
  int blocks = min((N + threads - 1) / threads, 512);

  deviceReduceKernel_MANY<<<blocks, threads, 0, stream>>>(in, work, N, M);
  deviceReduceKernel_MANY<<<1, 512, 0, stream>>>(work, out, blocks, M);
}

void deviceReduce_suffle_MANY(double *in, double* out, int N, int threads, int M, double* work) {
  int blocks = min((N + threads - 1) / threads, 512);

  deviceReduceKernel_MANY<<<blocks, threads>>>(in, work, N, M);
  deviceReduceKernel_MANY<<<1, 512>>>(work, out, blocks, M);
}

      
