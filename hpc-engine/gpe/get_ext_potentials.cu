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
#define double_complex Complex
#define __externc
#include "wslda_potdens.h"

#include "get_ext_potentials.h"

// =======================================================================================
// ======================== get_ext_potentials - handlers ================================
// =======================================================================================

// __global__ void kernel_get_v_ext(int it, int spin, double *data)
// {
//     size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
//     int ix;

//     if(ixyz<NX)
//     {
//         ix = ixyz; // decode cartesian coordinates
//         data[ixyz]=u_ext(ix,0,0,it,spin);
//     }
// }

/**
 * @return array with external potential
 * */
extern "C" int get_v_ext(int datadim, int spin, int it, double *data)
{
    // int nthreads = 256;
    // int nblocks = (int)ceil((float)NX/nthreads);

    // // allocate cuda memory
    // double *wrkspace = (double *)pca_cufft_work_area; // reuse workspace
    // kernel_get_v_ext<<<nblocks, nthreads>>>(it, spin, wrkspace);

    // if( cudaMemcpy( data , wrkspace, sizeof(double)*NX, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}

// __global__ void kernel_get_delta_ext(int it, Complex *deltain, Complex *data)
// {
//     size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
//     int ix;

//     if(ixyz<NX)
//     {
//         ix = ixyz; // decode cartesian coordinates
//         data[ixyz]=macro_delta_ext(ix, 0, 0, it, deltain[ixyz]);
//     }
// }

/**
 * @return array with external potential
 * */
extern "C" int get_delta_ext(int datadim, int it, void *deltain, void *data)
{
    // int nthreads = 256;
    // int nblocks = (int)ceil((float)NX/nthreads);

    // // allocate cuda memory
    // Complex *wrkspace = (Complex *)pca_cufft_work_area; // reuse workspace
    // kernel_get_delta_ext<<<nblocks, nthreads>>>(it, (Complex *)deltain, (Complex *)wrkspace);

    // if( cudaMemcpy( data , wrkspace, sizeof(Complex)*NX, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}

// __global__ void kernel_get_velocity_ext(int it, int spin, double *datax, double *datay, double *dataz)
// {
//     size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
//     int ix;

//     if(ixyz<NX)
//     {
//         ix = ixyz; // decode cartesian coordinates
//         datax[ixyz]=velocity_ext(ix, 0, 0, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
//         datay[ixyz]=velocity_ext(ix, 0, 0, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
//         dataz[ixyz]=velocity_ext(ix, 0, 0, it, spin, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
//     }
// }

/**
 * @return array with external potential
 * */
extern "C" int get_velocity_ext(int datadim, int spin, int it, double *data)
{
    // int nthreads = 256;
    // int nblocks = (int)ceil((float)NX/nthreads);

    // double *wrkspace = (double *)pca_cufft_work_area; // reuse workspace
    // double *vx = wrkspace + 0*NX;
    // double *vy = wrkspace + 1*NX;
    // double *vz = wrkspace + 2*NX;

    // kernel_get_velocity_ext<<<nblocks, nthreads>>>(it, spin, vx, vy, vz);

    // if( cudaMemcpy( data , wrkspace, sizeof(double)*NX*3, cudaMemcpyDeviceToHost )!= cudaSuccess ) return 1;

    return 0;
}