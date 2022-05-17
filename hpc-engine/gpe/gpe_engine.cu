/***************************************************************************
 *   Copyright (C) 2015 by                                                 *
 *   WARSAW UNIVERSITY OF TECHNOLOGY                                       *
 *   FACULTY OF PHYSICS                                                    *
 *   NUCLEAR THEORY GROUP                                                  *
 *   See also AUTHORS file                                                 *
 *                                                                         *
 *   This file is a part of GPE for GPU project.                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/ 
#include <math.h>
#include <complex>      // std::complex

#include "predefines.h"
#include "gpe_engine.h"

typedef std::complex<double> cplx;

    
/***************************************************************************/ 
/******************************** GLOBALS **********************************/
/***************************************************************************/
typedef struct
{
    double *kkx;
    double *kky;
    double *kkz;
    cufftHandle plan; // cufft plan
    Complex * d_wrk; // workspace for cufft - device
    Complex * d_wrk2; // additional work space - device
    Complex * d_psi; // wave function - device
    Complex * d_psi2; // copy of wave function - device
    double * d_wrk2R; // d_wrk2R = (double *) d_wrk2; - for convinience
    double * d_wrk3R; // for quantum friction computation
    Complex * d_wrk3C; // for quantum friction computation
    Complex * d_psi_ref; // reference psi on device
    Complex * d_psi_Lz; // adhoc
    Complex * d_psi_copy; // adhoc
    
    // other variables
    uint it; // number of steps performed from last call gpe_set_psi
    double dt;
    double t0;
    double alpha;
    double beta;
    double npart;
    double qfcoeff;
    int threads;
    int blocks;
    
} gpe_mem_t;

gpe_mem_t gpe_mem;


// CONST MEMORY - device
__constant__ double d_alpha;
__constant__ double d_beta;
__constant__ double d_qfcoeff; // quantum friction coeff

__constant__ double d_kkx[NX];
__constant__ double d_kky[NY];
__constant__ double d_kkz[NZ];
__constant__ Complex d_step_coeff; // 0.5*dt/(i*alpha-beta)
__constant__ Complex d_exp_kkx2[NX]; // exp( (dt/(i*alpha-beta)) * 1/(2gamma) * kx^2 )
__constant__ Complex d_exp_kky2[NY]; // exp( (dt/(i*alpha-beta)) * 1/(2gamma) * ky^2 )
__constant__ Complex d_exp_kkz2_over_nxyz[NZ]; // exp( (dt/(i*alpha-beta)) * 1/(2gamma) * kz^2 ) / nxyz

__constant__ double d_user_param[MAX_USER_PARAMS];
__constant__ uint d_nx; // lattice size in x direction
__constant__ uint d_ny; // lattice size in y direction
__constant__ uint d_nz; // lattice size in z direction
__constant__ double d_dt;
__constant__ double d_t0;
__constant__ double d_npart;
__constant__ Complex *d_psi_ref; // pointer to reference psi on device - use gpe_set_psi_ref() function to set it

__constant__ void *dc_extra_data;
__constant__ size_t dc_extra_data_size;

#define PARTICLES 1
#define DIMERS 2

// TODO -> enable wslda_density in problem-definition.h file
#define double_complex Complex
#define __externc
#include "wslda_potdens.h"
#include "problem-definition.h"

#if GPE_FOR == PARTICLES
#define GAMMA 1.0
#elif GPE_FOR == DIMERS
#define GAMMA 2.0
#endif

/***************************************************************************/ 
/******************************** MACROS ***********************************/
/***************************************************************************/ 
#define nx NX
#define ny NY
#define nz NZ
#define nxyz (nx*ny*nz)

#define ixyz2ixiyiz(ixyz,_ix,_iy,_iz,i)     \
    i=ixyz;                                 \
    _ix=i/(ny*nz);                          \
    i=i-_ix * ny * nz;                      \
    _iy=i/nz;                               \
    _iz=i-_iy * nz;

#define constgpu(c) (double)(c)
    
#define myerrchecksimple(cmd,err)                                                                       \
    {                                                                                                   \
        err=cmd ;                                                                                       \
        if(err != 0)                                                                                    \
        {                                                                                               \
            fprintf( stderr, "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ;                      \
            return (int)(err);                                                                          \
        }                                                                                               \
    } 
    
#define myerrcheck(cmd)                                                                                 \
    {                                                                                                   \
        err=cmd ;                                                                                       \
        if(err != cudaSuccess)                                                                          \
        {                                                                                               \
            fprintf( stderr, "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ;                      \
            fprintf( stderr, "CUDA ERROR %d: %s\n", err, cudaGetErrorString((cudaError_t)(err)));       \
            return (int)(err);                                                                          \
        }                                                                                               \
    } 
    
#define GPE_QF_EPSILON 1.0e-12
    
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
__device__ void warpReduce(volatile double *sdata, unsigned int tid)
{
    if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
    if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
    if (blockSize >= 16) sdata[tid] += sdata[tid +  8];
    if (blockSize >=  8) sdata[tid] += sdata[tid +  4];
    if (blockSize >=  4) sdata[tid] += sdata[tid +  2];
    if (blockSize >=  2) sdata[tid] += sdata[tid +  1];
}

template <unsigned int blockSize>
__global__ void __reduce_kernel__(double *g_idata, double *g_odata, int n)
{
    extern __shared__ double sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    unsigned int ishift=i+blockDim.x;

    if(ishift<n) sdata[tid] = g_idata[i] + g_idata[ishift];
    else if(i<n) sdata[tid] = g_idata[i];
    else         sdata[tid] = 0.0;
    __syncthreads();

    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    if (tid < 32) warpReduce<blockSize>(sdata, tid);

    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

void call_reduction_kernel(int dimGrid, int dimBlock, int size, double *d_idata, double *d_odata, cudaStream_t stream)
{
    int smemSize=dimBlock*sizeof(double);
    switch (dimBlock)
    {
        case 1024:
            __reduce_kernel__<1024><<< dimGrid, dimBlock, smemSize , stream >>>(d_idata, d_odata, size); break;
        case 512:
            __reduce_kernel__< 512><<< dimGrid, dimBlock, smemSize , stream >>>(d_idata, d_odata, size); break;
        case 256:
            __reduce_kernel__< 256><<< dimGrid, dimBlock, smemSize , stream >>>(d_idata, d_odata, size); break;
        case 128:
            __reduce_kernel__< 128><<< dimGrid, dimBlock, smemSize , stream >>>(d_idata, d_odata, size); break;
        case 64:
            __reduce_kernel__<  64><<< dimGrid, dimBlock, smemSize , stream >>>(d_idata, d_odata, size); break;
    }
}
/**
 * Function does fast reduction (sum of elements) of array.
 * Result is located in partial_sums[0] element
 * If partial_sums==array then array will be destroyed
 * */
int local_reduction(double *array, int size, double *partial_sums, int threads, cudaStream_t stream)
{
    int blocks=(int)ceil((float)size/threads);
    unsigned int lthreads=threads/2; // Threads is always power of 2
    if(lthreads<64) lthreads=64; // at least 2*warp_size
    unsigned int new_blocks, current_size;

    // First reduction of the array
    call_reduction_kernel(blocks, lthreads, size, array, partial_sums,stream);

    // Do iteratively reduction of partial_sums
    current_size=blocks;
    while(current_size>1)
    {
        new_blocks=(int)ceil((float)current_size/threads);
        lthreads=opt_threads(new_blocks,threads, current_size)/2;
        if(lthreads<64) lthreads=64; // at least 2*warp_size
        call_reduction_kernel(new_blocks, lthreads, current_size, partial_sums, partial_sums,stream);
        current_size=new_blocks;
    }

    return 0;
}


////////////////////////////////////////////////////////////////////////////////
// Complex operations
////////////////////////////////////////////////////////////////////////////////

// Complex addition
static __device__ __host__ inline Complex cplxAdd(Complex a, Complex b)
{
    Complex c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}

static __device__ __host__ inline Complex cplxSub(Complex a, Complex b)
{
    Complex c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    return c;
}

// Complex scale
static __device__ __host__ inline Complex cplxScale(Complex a, double s)
{
    Complex c;
    c.x = s * a.x;
    c.y = s * a.y;
    return c;
}

// Complex scalei = multiplication  by purely imaginay number, i.e: c = a*(i*s)
static __device__ __host__ inline Complex cplxScalei(Complex a, double s)
{
    Complex c;
    c.x = s * a.y * constgpu(-1.0);
    c.y = s * a.x;
    return c;
}

static __device__ __host__ inline Complex cplxConj(Complex a)
{
    a.y=-1.0*a.y;
    return a;
}

// norm2=|a|^2
static __device__ __host__ inline double cplxNorm2(Complex a)
{
    return (double)(a.x*a.x + a.y*a.y);
}

// abs=|a|
static __device__ __host__ inline double cplxAbs(Complex a)
{
    return (double)(sqrt(cplxNorm2(a)));
}

static __device__ __host__ inline double cplxArg(Complex a)
{
    return (double)(atan2(a.y,a.x));
}


// Complex multiplication
static __device__ __host__ inline Complex cplxMul(Complex a, Complex b)
{
    Complex c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

// Complex multiplication - and return only real part
static __device__ __host__ inline double cplxMulR(Complex a, Complex b)
{
    double c;
    c = a.x * b.x - a.y * b.y;
    return c;
}

// Complex multiplication - and return only imag part
static __device__ __host__ inline double cplxMulI(Complex a, Complex b)
{
    double c;
    c = a.x * b.y + a.y * b.x;
    return c;
}

static __device__ __host__ inline Complex cplxDiv(Complex a, Complex b)
{
    Complex c;
    double t=b.x*b.x + b.y*b.y;
    c.x = (a.x * b.x + a.y * b.y)/t;
    c.y = (a.y * b.x - a.x * b.y)/t;
    return c;
}

static __device__ __host__ inline Complex cplxSqrt(Complex z)
{
    Complex r;
    double norm=sqrt(z.x*z.x + z.y*z.y);

    r.x=sqrt((norm+z.x)/2.);
    if(z.y>0.0) r.y=sqrt((norm-z.x)/2.);
    else r.y=constgpu(-1.)*sqrt((norm-z.x)/2.);
    return r;
}

static __device__ __host__ inline Complex cplxLog(Complex z)
{
    Complex r;
    double norm=z.x*z.x + z.y*z.y;
    double arg=atan2(z.y,z.x);
    r.x=constgpu(0.5)*log(norm);
    r.y=arg;
    return r;
}

// z=exp(i*x)=cos(x)+i*sin(x), x - real number
static __device__ __host__ inline Complex cplxExpi(double x)
{
    Complex r;
    r.x=cos(x);
    r.y=sin(x);
    return r;
}

// z=exp(x), x - complex number
static __device__ __host__ inline Complex cplxExp(Complex x)
{
    double t;
    Complex r;
    t=exp(x.x);
    r.x=t*cos(x.y);
    r.y=t*sin(x.y);
    return r;
}

// returns 1/c
static __device__ __host__ inline Complex cplxInv(Complex c)
{
    Complex r;
    double n=c.x*c.x + c.y*c.y;
    r.x=c.x/n;
    r.y=constgpu(-1.0)*c.y/n;
    return r;
}

/***************************************************************************/ 
/****************************** FUNCTIONS **********************************/
/***************************************************************************/
int gpe_set_extra_data(void* extra_data, size_t extra_data_size)
{
    if( cudaMemcpyToSymbol(dc_extra_data,      &extra_data,      sizeof(void *))!= cudaSuccess ) return 1;
    if( cudaMemcpyToSymbol(dc_extra_data_size, &extra_data_size, sizeof(size_t))!= cudaSuccess ) return 2;
    return 0;
}

void gpe_get_lattice(int *_nx, int *_ny, int *_nz)
{
    *_nx = nx;
    *_ny = ny;
    *_nz = nz;
}

int gpe_create_engine(double alpha, double beta, double dt, double npart, int nthreads)
{
    cudaError err;
    uint ui;
    int i,j;
    double r;
    
    #ifndef GAMMA
        return -99; // not supported mode
    #endif
    
    // Set number of blocks, if number of threads is given
    gpe_mem.threads=nthreads;
    gpe_mem.blocks=(int)ceil((float)nxyz/nthreads);

//     printf("GPU SETTING: THREADS=%d, BLOCKS=%d, THREADS*BLOCKS=%d, nxyz=%d\n",gpe_mem.threads,gpe_mem.blocks,gpe_mem.threads*gpe_mem.blocks,nxyz);
    
    // Fill const memory
    ui=nx;
    myerrcheck( cudaMemcpyToSymbol(d_nx, &ui, sizeof(uint)) ) ;
    ui=ny;
    myerrcheck( cudaMemcpyToSymbol(d_ny, &ui, sizeof(uint)) ) ;
    ui=nz;
    myerrcheck( cudaMemcpyToSymbol(d_nz, &ui, sizeof(uint)) ) ;   
    myerrcheck( cudaMemcpyToSymbol(d_alpha, &alpha, sizeof(double)) ) ;
    gpe_mem.alpha=alpha;
    myerrcheck( cudaMemcpyToSymbol(d_beta, &beta, sizeof(double)) ) ;
    gpe_mem.beta=beta;
    myerrcheck( cudaMemcpyToSymbol(d_dt, &dt, sizeof(double)) ) ;
    gpe_mem.dt=dt;
    r=0.0;
    myerrcheck( cudaMemcpyToSymbol(d_t0, &r, sizeof(double)) ) ;
    myerrcheck( cudaMemcpyToSymbol(d_qfcoeff, &r, sizeof(double)) ) ;
    gpe_mem.t0=0.0;
    gpe_mem.it=0;
    gpe_mem.qfcoeff=0.0;
    myerrcheck( cudaMemcpyToSymbol(d_npart, &npart, sizeof(double)) ) ;
    gpe_mem.npart=npart;
    
    // Generate arrays
    gpemalloc(gpe_mem.kkx,nx,double);
    gpemalloc(gpe_mem.kky,ny,double);
    gpemalloc(gpe_mem.kkz,nz,double);  
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */
//     printf("CREATING LATTICE: %dx%dx%d...\n", nx,ny,nz);
    
    for ( i = 0 ; i <= nx / 2 - 1 ; i++ ) {
        gpe_mem.kkx[ i ] = 2. * ( double ) M_PI / nx * ( double ) i ;  }
    j = - i ;
    for ( i = nx / 2 ; i < nx ; i++ ) 
    {
        gpe_mem.kkx[ i ] = 2. * ( double ) M_PI / nx * ( double ) j ; 
        j++ ;
    }

    for ( i = 0 ; i <= ny / 2 - 1 ; i++ ) {
        gpe_mem.kky[ i ] = 2. * ( double ) M_PI / ny * ( double ) i ;  }
    j = - i ;
    for ( i = ny / 2 ; i < ny ; i++ ) 
    {
        gpe_mem.kky[ i ] = 2. * ( double ) M_PI / ny * ( double ) j ; 
        j++ ;
    }

    for ( i = 0 ; i <= nz / 2 - 1 ; i++ ) {
        gpe_mem.kkz[ i ] = 2. * ( double ) M_PI / nz * ( double ) i ;  }
    j = - i ;
    for ( i = nz / 2 ; i < nz ; i++ ) 
    {
        gpe_mem.kkz[ i ] = 2. * ( double ) M_PI / nz * ( double ) j ; 
        j++ ;
    }    
    myerrcheck( cudaMemcpyToSymbol(d_kkx, gpe_mem.kkx, nx*sizeof(double)) ) ;
    myerrcheck( cudaMemcpyToSymbol(d_kky, gpe_mem.kky, ny*sizeof(double)) ) ;
    myerrcheck( cudaMemcpyToSymbol(d_kkz, gpe_mem.kkz, nz*sizeof(double)) ) ;
        
    // 0.5*dt/(i*alpha-beta)*GAMMA
    cplx c1(GAMMA*0.5*dt,0.0);
    cplx c2(-1.0*beta, alpha);
    cplx c3=c1/c2;
    Complex c4 = {c3.real(), c3.imag()};
    myerrcheck( cudaMemcpyToSymbol(d_step_coeff, &c4, sizeof(Complex)) ) ;
    
    // kinetic operator mulipliers
    Complex *carr;
    
    // nx direction
    gpemalloc(carr,nx,Complex);
    for(ui=0; ui<nx; ui++)
    {
        c1=exp(c3*gpe_mem.kkx[ui]*gpe_mem.kkx[ui]/(GAMMA*GAMMA));
        carr[ui].x=c1.real(); carr[ui].y=c1.imag();
//         printf("%d %f %f\n", ui, carr[ui].x, carr[ui].y);
    }
    myerrcheck( cudaMemcpyToSymbol(d_exp_kkx2, carr, nx*sizeof(Complex)) ) ;
    free(carr);
    
    // ny direction
    gpemalloc(carr,ny,Complex);
    for(ui=0; ui<ny; ui++)
    {
        c1=exp(c3*gpe_mem.kky[ui]*gpe_mem.kky[ui]/(GAMMA*GAMMA));
        carr[ui].x=c1.real(); carr[ui].y=c1.imag();
//         printf("%d %f %f\n", ui, carr[ui].x, carr[ui].y);
    }
    myerrcheck( cudaMemcpyToSymbol(d_exp_kky2, carr, ny*sizeof(Complex)) ) ;
    free(carr);
    
    // nz direction
    gpemalloc(carr,nz,Complex);
    for(ui=0; ui<nz; ui++)
    {
        c1=exp(c3*gpe_mem.kkz[ui]*gpe_mem.kkz[ui]/(GAMMA*GAMMA)) / (double)(nxyz);
        carr[ui].x=c1.real(); carr[ui].y=c1.imag();
//         printf("%d %f %f\n", ui, carr[ui].x, carr[ui].y);
    }
    myerrcheck( cudaMemcpyToSymbol(d_exp_kkz2_over_nxyz, carr, nz*sizeof(Complex)) ) ;
    free(carr);
    
    // create cufft plans
    cufftResult cufft_result;
    cufft_result=cufftCreate(&gpe_mem.plan); if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    cufft_result=cufftSetAutoAllocation(gpe_mem.plan, 0); if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    size_t workSize;
    cufft_result=cufftMakePlan3d(gpe_mem.plan, nx, ny, nz, CUFFT_Z2Z, &workSize);
    if(workSize<sizeof(Complex)*nxyz) workSize=sizeof(Complex)*nxyz;
    myerrcheck( cudaMalloc( &gpe_mem.d_wrk , workSize ) );
    cufft_result=cufftSetWorkArea(gpe_mem.plan, gpe_mem.d_wrk); if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    
    // create buffers
    myerrcheck( cudaMalloc( &gpe_mem.d_wrk2 , sizeof(Complex)*nxyz ) );
    myerrcheck( cudaMalloc( &gpe_mem.d_psi , sizeof(Complex)*nxyz ) );
    myerrcheck( cudaMalloc( &gpe_mem.d_psi2 , sizeof(Complex)*nxyz ) );
    gpe_mem.d_wrk2R = (double *) gpe_mem.d_wrk2; 
    
    gpe_mem.d_wrk3R = NULL;
    gpe_mem.d_wrk3C = NULL;
    
    // reference psi
    gpe_mem.d_psi_ref = NULL;
    myerrcheck( cudaMemcpyToSymbol(d_psi_ref, &gpe_mem.d_psi_ref, sizeof(Complex *)) ) ;
    
    gpe_mem.d_psi_Lz = NULL;
    gpe_mem.d_psi_copy = NULL;
    
    return 0; // success
}

int gpe_destroy_engine()
{
    cudaError err;
    cufftResult cufft_result;
    free(gpe_mem.kkx);
    free(gpe_mem.kky);
    free(gpe_mem.kkz);
    cufft_result=cufftDestroy(gpe_mem.plan); if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    myerrcheck( cudaFree(gpe_mem.d_wrk) );
    myerrcheck( cudaFree(gpe_mem.d_wrk2) );
    myerrcheck( cudaFree(gpe_mem.d_psi) );
    myerrcheck( cudaFree(gpe_mem.d_psi2) );
    if(gpe_mem.d_wrk3R != NULL) myerrcheck( cudaFree(gpe_mem.d_wrk3R) );
    if(gpe_mem.d_wrk3C != NULL) myerrcheck( cudaFree(gpe_mem.d_wrk3C) );
    if(gpe_mem.d_psi_ref != NULL) myerrcheck( cudaFree(gpe_mem.d_psi_ref) );
    if(gpe_mem.d_psi_Lz != NULL) myerrcheck( cudaFree(gpe_mem.d_psi_Lz) );
    if(gpe_mem.d_psi_copy != NULL) myerrcheck( cudaFree(gpe_mem.d_psi_copy) );
    
    return 0; // success
}


int gpe_change_alpha_beta(double alpha, double beta)
{
    cudaError err;
    uint ui;
    
    myerrcheck( cudaMemcpyToSymbol(d_alpha, &alpha, sizeof(double)) ) ;
    gpe_mem.alpha=alpha;
    myerrcheck( cudaMemcpyToSymbol(d_beta, &beta, sizeof(double)) ) ;
    gpe_mem.beta=beta;
    
    double dt = gpe_mem.dt;
    
    // 0.5*dt/(i*alpha-beta)*GAMMA
    cplx c1(GAMMA*0.5*dt,0.0);
    cplx c2(-1.0*beta, alpha);
    cplx c3=c1/c2;
    Complex c4 = {c3.real(), c3.imag()};
    myerrcheck( cudaMemcpyToSymbol(d_step_coeff, &c4, sizeof(Complex)) ) ;
    
    // kinetic operator mulipliers
    Complex *carr;
    
    // nx direction
    gpemalloc(carr,nx,Complex);
    for(ui=0; ui<nx; ui++)
    {
        c1=exp(c3*gpe_mem.kkx[ui]*gpe_mem.kkx[ui]/(GAMMA*GAMMA));
        carr[ui].x=c1.real(); carr[ui].y=c1.imag();
    }
    myerrcheck( cudaMemcpyToSymbol(d_exp_kkx2, carr, nx*sizeof(Complex)) ) ;
    free(carr);
    
    // ny direction
    gpemalloc(carr,ny,Complex);
    for(ui=0; ui<ny; ui++)
    {
        c1=exp(c3*gpe_mem.kky[ui]*gpe_mem.kky[ui]/(GAMMA*GAMMA));
        carr[ui].x=c1.real(); carr[ui].y=c1.imag();
    }
    myerrcheck( cudaMemcpyToSymbol(d_exp_kky2, carr, ny*sizeof(Complex)) ) ;
    free(carr);
    
    // nz direction
    gpemalloc(carr,nz,Complex);
    for(ui=0; ui<nz; ui++)
    {
        c1=exp(c3*gpe_mem.kkz[ui]*gpe_mem.kkz[ui]/(GAMMA*GAMMA)) / (double)(nxyz);
        carr[ui].x=c1.real(); carr[ui].y=c1.imag();
    }
    myerrcheck( cudaMemcpyToSymbol(d_exp_kkz2_over_nxyz, carr, nz*sizeof(Complex)) ) ;
    free(carr);
    
    return 0;
}

int gpe_set_time(double t0)
{
    cudaError err;
    
    myerrcheck( cudaMemcpyToSymbol(d_t0, &t0, sizeof(double)) ) ;
    gpe_mem.t0=t0;
    gpe_mem.it=0;    
    
    return 0;
}

int gpe_set_user_params(int size, double *params)
{
    cudaError err;
    
    if(size>MAX_USER_PARAMS) return -9;
    
    myerrcheck( cudaMemcpyToSymbol(d_user_param, params, MAX_USER_PARAMS*sizeof(double)) ) ;
    
    return 0;
}

int gpe_set_quantum_friction_coeff(double qfcoeff)
{
    cudaError err;
    if(qfcoeff!=0.0)
    {
        qfcoeff=qfcoeff/( GAMMA*(double)(nxyz) );
        myerrcheck( cudaMemcpyToSymbol(d_qfcoeff, &qfcoeff, sizeof(double)) ) ;
        gpe_mem.qfcoeff=qfcoeff;
        
        if(gpe_mem.d_wrk3R==NULL) myerrcheck( cudaMalloc( &gpe_mem.d_wrk3R , sizeof(double)*nxyz ) );
        if(gpe_mem.d_wrk3C==NULL) myerrcheck( cudaMalloc( &gpe_mem.d_wrk3C , sizeof(Complex)*nxyz ) );
    }
    else
    {
        myerrcheck( cudaMemcpyToSymbol(d_qfcoeff, &qfcoeff, sizeof(double)) ) ;
        gpe_mem.qfcoeff=qfcoeff;
        
        if(gpe_mem.d_wrk3R != NULL) myerrcheck( cudaFree(gpe_mem.d_wrk3R) );
        if(gpe_mem.d_wrk3C != NULL) myerrcheck( cudaFree(gpe_mem.d_wrk3C) );   
        
        gpe_mem.d_wrk3R = NULL;
        gpe_mem.d_wrk3C = NULL;
    }
    
    return 0;
}


int gpe_set_psi_ref(Complex * psi)
{
    cudaError err;
    
    // allocate memory
    if(gpe_mem.d_psi_ref==NULL) 
    {
        myerrcheck( cudaMalloc( &gpe_mem.d_psi_ref , sizeof(Complex)*nxyz ) );
        myerrcheck( cudaMemcpyToSymbol(d_psi_ref, &gpe_mem.d_psi_ref, sizeof(Complex *)) ) ;
    }
    // copy buffer and pointer
    myerrcheck( cudaMemcpy( gpe_mem.d_psi_ref , psi , sizeof(Complex)*nxyz, cudaMemcpyHostToDevice ) );  
    
    return 0;
}

int gpe_set_psi_ref_from_present_state()
{
    cudaError err;
    
    // allocate memory
    if(gpe_mem.d_psi_ref==NULL) 
    {
        myerrcheck( cudaMalloc( &gpe_mem.d_psi_ref , sizeof(Complex)*nxyz ) );
        myerrcheck( cudaMemcpyToSymbol(d_psi_ref, &gpe_mem.d_psi_ref, sizeof(Complex *)) ) ;
    }
    // copy buffer and pointer
    myerrcheck( cudaMemcpy( gpe_mem.d_psi_ref , gpe_mem.d_psi , sizeof(Complex)*nxyz, cudaMemcpyDeviceToDevice ) );  
    
    return 0;
}

int gpe_clear_psi_ref()
{
    cudaError err;
    
    if(gpe_mem.d_psi_ref != NULL) 
    {
        myerrcheck( cudaFree(gpe_mem.d_psi_ref) );
        gpe_mem.d_psi_ref=NULL;
        myerrcheck( cudaMemcpyToSymbol(d_psi_ref, &gpe_mem.d_psi_ref, sizeof(Complex *)) ) ;
    }
    
    
    return 0;
}

__device__  double gpe_density(Complex psi)
{
    return GAMMA * (psi.x*psi.x + psi.y*psi.y); // |psi|^2 * GAMMA, where GAMMA=1 for particles, GAMMA=2 for dimers
}

__global__ void __gpe_compute_density__(Complex *psi_in, double *rho_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    if(ixyz<nxyz)
    {
        rho_out[ixyz] = gpe_density(psi_in[ixyz]);
    }
}

__global__ void __gpe_normalize__(Complex *psi_inout, double *sumrho)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    if(ixyz<nxyz)
    {
//         if(ixyz==0) printf("sumrho[0]=%f\n", sumrho[0]);
        psi_inout[ixyz] = cplxScale(psi_inout[ixyz], sqrt(d_npart/sumrho[0]));
    }
}

// Normalizes wave function
int gpe_normalize(Complex *psi, double *wrk)
{
    __gpe_compute_density__<<<gpe_mem.blocks, gpe_mem.threads>>>(psi, wrk);
    int ierr;
    myerrchecksimple( local_reduction(wrk, nxyz, wrk, gpe_mem.threads, 0), ierr );
    __gpe_normalize__<<<gpe_mem.blocks, gpe_mem.threads>>>(psi, wrk);
    
    return 0;
}

int gpe_normalize_psi()
{
    return gpe_normalize(gpe_mem.d_psi, gpe_mem.d_wrk2R);
}

int gpe_set_psi(double t, Complex * psi)
{
    cudaError err;
    myerrcheck( cudaMemcpyToSymbol(d_t0, &t, sizeof(double)) ) ;
    gpe_mem.it=0;
    gpe_mem.t0=t;
    myerrcheck( cudaMemcpy( gpe_mem.d_psi , psi , sizeof(Complex)*nxyz, cudaMemcpyHostToDevice ) );
    
    return 0;
}

int gpe_get_psi(double *t, Complex * psi)
{
    cudaError err;
    myerrcheck( cudaMemcpy( psi , gpe_mem.d_psi , sizeof(Complex)*nxyz, cudaMemcpyDeviceToHost ) );
    *t = gpe_mem.t0 + gpe_mem.dt*gpe_mem.it;
    
    return 0;
}


/**
 * construct  exp(-i*dt*V/2) and apply exp(-i*dt*V/2) * psi 
 * */
__global__ void __gpe_exp_Vstep1_(uint it, Complex *psi_in, Complex *psi_out, double * wrkR)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    // registers
    Complex lpsi, exp_lv;
    double lrho, lv, dEDFdn;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        
        lpsi = psi_in[ixyz]; // psi to register
        gpe_modify_psi(ix, iy, iz, it, &lpsi, d_user_param, 0, NULL); // modify psi
        lrho = gpe_density(lpsi); // compute density
        compute_potentials_gpe(it, lrho, &dEDFdn, d_user_param, NULL, NULL); //compute potentials defining Hamiltonian
        lv=v_ext(ix, iy, iz, it, NULL, d_user_param, NULL, NULL) + dEDFdn; // external potential + mean field
        
        wrkR[ixyz]=lv; // it will be use later
        exp_lv = cplxExp( cplxScale(d_step_coeff,lv) );
        
        psi_out[ixyz] = cplxMul(lpsi, exp_lv); // apply and save
    }    
}

__global__ void __gpe_exp_Vstep1_qf_(uint it, Complex *psi_in, Complex *psi_out, double * wrkR, double *qfpotential)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    // registers
    Complex lpsi, exp_lv;
    double lrho, lv, dEDFdn;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        
        lpsi = psi_in[ixyz]; // psi to register
        gpe_modify_psi(ix, iy, iz, it, &lpsi, d_user_param, 0, NULL); // modify psi
        lrho = gpe_density(lpsi); // compute density
        compute_potentials_gpe(it, lrho, &dEDFdn, d_user_param, NULL, NULL); //compute potentials defining Hamiltonian
        lv=v_ext(ix, iy, iz, it, NULL, d_user_param, NULL, NULL) + dEDFdn + qfpotential[ixyz]; 
           // external potential + mean field + quantum friction potential
        
        wrkR[ixyz]=lv; // it will be use later
        exp_lv = cplxExp( cplxScale(d_step_coeff,lv) );
        
        psi_out[ixyz] = cplxMul(lpsi, exp_lv); // apply and save
    }    
}

__global__ void __gpe_exp_Vstep2_(uint it, Complex *psi_in, Complex *psi_out, double * wrkR, Complex * wrkC)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    // registers
    Complex lpsi, exp_lv;
    double lrho, lv, dEDFdn;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        
        lpsi = psi_in[ixyz]; // psi to register
        lv = wrkR[ixyz]; // potentials to register
        exp_lv = cplxExp( cplxScale(d_step_coeff,lv) );
        lpsi=cplxMul(lpsi, exp_lv); // finalize step from predictor
        
        lrho = gpe_density(lpsi); // compute density
        compute_potentials_gpe(it+1, lrho, &dEDFdn, d_user_param, NULL, NULL); //compute potentials defining Hamiltonian
        lv=0.5*(lv + v_ext(ix, iy, iz, it+1, NULL, d_user_param, NULL, NULL) + dEDFdn); // external potential + mean field - take average
        exp_lv = cplxExp( cplxScale(d_step_coeff,lv) );
        wrkC[ixyz]=exp_lv; // it will be used later
        
        lpsi = psi_out[ixyz]; // psi to register
        gpe_modify_psi(ix, iy, iz, it, &lpsi, d_user_param, 0, NULL); // modify psi
        psi_out[ixyz] = cplxMul(lpsi, exp_lv); // apply and save      
    }    
}

__global__ void __gpe_exp_Vstep2_part1_(Complex *psi_in, Complex *psi_out, double * wrkR)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    
    // registers
    Complex exp_lv;

    if(ixyz<nxyz)
    {
        exp_lv = cplxExp( cplxScale(d_step_coeff, wrkR[ixyz]) );
        psi_out[ixyz]=cplxMul(psi_in[ixyz], exp_lv); // finalize step from predictor     
    }    
}

__global__ void __gpe_exp_Vstep2_part2_(uint it, Complex *psi_in, Complex *psi_out, double * wrkR, Complex * wrkC)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    // registers
    Complex lpsi, exp_lv;
    double lrho, lv, dEDFdn;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        
        lpsi = psi_in[ixyz]; // psi to register
        lv = wrkR[ixyz]; // potentials to register
        lrho = gpe_density(lpsi); // compute density
        compute_potentials_gpe(it+1, lrho, &dEDFdn, d_user_param, NULL, NULL); //compute potentials defining Hamiltonian
        lv=0.5*(lv + v_ext(ix, iy, iz, it+1, NULL, d_user_param, NULL, NULL) + dEDFdn); // external potential + mean field - take average
        exp_lv = cplxExp( cplxScale(d_step_coeff,lv) );
        wrkC[ixyz]=exp_lv; // it will be used later
        
        lpsi = psi_out[ixyz]; // psi to register
        gpe_modify_psi(ix, iy, iz, it, &lpsi, d_user_param, 0, NULL); // modify psi
        psi_out[ixyz] = cplxMul(lpsi, exp_lv); // apply and save      
    }    
}

__global__ void __gpe_exp_Vstep2_part2_qf_(uint it, Complex *psi_in, Complex *psi_out, double * wrkR, Complex * wrkC, double *qfpotential)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    // registers
    Complex lpsi, exp_lv;
    double lrho, lv, dEDFdn;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        
        lpsi = psi_in[ixyz]; // psi to register
        lv = wrkR[ixyz]; // potentials to register
        lrho = gpe_density(lpsi); // compute density
        compute_potentials_gpe(it+1, lrho, &dEDFdn, d_user_param, NULL, NULL); //compute potentials defining Hamiltonian
        lv=0.5*(lv + v_ext(ix, iy, iz, it+1, NULL, d_user_param, NULL, NULL) + dEDFdn + qfpotential[ixyz]) ; 
          // external potential + mean field + quantum friction potential - take average 
        exp_lv = cplxExp( cplxScale(d_step_coeff,lv) );
        wrkC[ixyz]=exp_lv; // it will be used later
        
        lpsi = psi_out[ixyz]; // psi to register
        gpe_modify_psi(ix, iy, iz, it, &lpsi, d_user_param, 0, NULL); // modify psi
        psi_out[ixyz] = cplxMul(lpsi, exp_lv); // apply and save      
    }    
}

__global__ void __gpe_exp_Vstep3_(Complex *psi_inout, Complex * wrkC)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    
    if(ixyz<nxyz)
    {
        psi_inout[ixyz] = cplxMul(psi_inout[ixyz], wrkC[ixyz]); // apply and save      
    }    
}

__global__ void __gpe_multiply_by_expT__(Complex *psi_in, Complex *psi_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    Complex _wavef;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        _wavef=psi_in[ixyz]; // bring to register
        _wavef=cplxMul(_wavef,d_exp_kkx2[ix]);
        _wavef=cplxMul(_wavef,d_exp_kky2[iy]);
        _wavef=cplxMul(_wavef,d_exp_kkz2_over_nxyz[iz]); // note - normalization factor is included here
        psi_out[ixyz]=_wavef; // send to global memory
    }    
}

__global__ void __gpe_multiply_by_k2_qf__(Complex *psi_in, Complex *psi_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 

        psi_out[ixyz]=cplxScale(psi_in[ixyz], d_qfcoeff*( d_kkx[ix]*d_kkx[ix] + d_kky[iy]*d_kky[iy] + d_kkz[iz]*d_kkz[iz] ) ); 
            
    }    
}

__global__ void __gpe_overlap_imag_qf__(Complex *psi1, Complex *psi2, double *overlap)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    double lrho;
    Complex lpsi;
    if(ixyz<nxyz)
    {
        lpsi = psi1[ixyz]; // psi to register
        lrho = gpe_density(lpsi); // compute density
        overlap[ixyz]= cplxMulI( cplxConj(lpsi),  psi2[ixyz] )/(lrho+GPE_QF_EPSILON);
    }
}

int gpe_compute_qf_potential(Complex *psi, Complex *wrk, double *qfpotential)
{
    cufftResult cufft_result;
    
        cufft_result=cufftExecZ2Z(gpe_mem.plan, psi, wrk, CUFFT_FORWARD);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        __gpe_multiply_by_k2_qf__<<<gpe_mem.blocks, gpe_mem.threads>>>(wrk, wrk);
        cufft_result=cufftExecZ2Z(gpe_mem.plan, wrk, wrk, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        __gpe_overlap_imag_qf__<<<gpe_mem.blocks, gpe_mem.threads>>>(psi, wrk, qfpotential);

        
    return 0;
}

/**
 * Function evolves wave funcion nt steps 
 * */
int gpe_evolve(int nt)
{
    cufftResult cufft_result;
    int i;
    int ierr;
        
    for(i=0; i<nt; i++)
    {
                
        if(gpe_mem.qfcoeff!=0.0) // quantum friction is active
        {
            myerrchecksimple( gpe_compute_qf_potential(gpe_mem.d_psi, gpe_mem.d_wrk3C, gpe_mem.d_wrk3R), ierr);
            __gpe_exp_Vstep1_qf_<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_psi, gpe_mem.d_psi2, gpe_mem.d_wrk2R, gpe_mem.d_wrk3R);
        }
        else
        {
            // potential part exp(V/2)
            __gpe_exp_Vstep1_<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_psi, gpe_mem.d_psi2, gpe_mem.d_wrk2R);
        }
        
        // kinetic part exp(T)
        cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi2, gpe_mem.d_psi2, CUFFT_FORWARD);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        __gpe_multiply_by_expT__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2, gpe_mem.d_psi2);
        cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi2, gpe_mem.d_psi2, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        
        // potential part exp(V/2)
        if(gpe_mem.beta==0.0 && gpe_mem.qfcoeff==0.0)
        {
            // without normalization
            __gpe_exp_Vstep2_<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_psi2, gpe_mem.d_psi, gpe_mem.d_wrk2R, gpe_mem.d_psi2);
        }
        else
        {
            __gpe_exp_Vstep2_part1_<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2, gpe_mem.d_psi2, gpe_mem.d_wrk2R);
            
            if(gpe_mem.beta!=0.0)
            {
                // with normalization between
                myerrchecksimple( gpe_normalize(gpe_mem.d_psi2, gpe_mem.d_wrk2R+nxyz), ierr);
            }
            
            if(gpe_mem.qfcoeff!=0.0) // quantum friction is active
            {
                myerrchecksimple( gpe_compute_qf_potential(gpe_mem.d_psi, gpe_mem.d_wrk3C, gpe_mem.d_wrk3R), ierr);
                __gpe_exp_Vstep2_part2_qf_<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_psi2, gpe_mem.d_psi, gpe_mem.d_wrk2R, gpe_mem.d_psi2, gpe_mem.d_wrk3R);
            }
            else
            {
                __gpe_exp_Vstep2_part2_<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_psi2, gpe_mem.d_psi, gpe_mem.d_wrk2R, gpe_mem.d_psi2);
            }
        }
        
        // kinetic part exp(T)
        cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi, gpe_mem.d_psi, CUFFT_FORWARD);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        __gpe_multiply_by_expT__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_psi);
        cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi, gpe_mem.d_psi, CUFFT_INVERSE);
        if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
        
        // potential part exp(V/2)
        __gpe_exp_Vstep3_<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_psi2);
        
        if(gpe_mem.beta!=0.0)
        {
            // normalize
            myerrchecksimple( gpe_normalize(gpe_mem.d_psi, gpe_mem.d_wrk2R+nxyz), ierr);
        }
        
        gpe_mem.it = gpe_mem.it + 1;
    }
    
    return 0;
}


__global__ void __gpe_multiply_by_kinetic_kernel__(Complex *psi)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 

        psi[ixyz]=cplxScale(psi[ixyz], ( d_kkx[ix]*d_kkx[ix] + d_kky[iy]*d_kky[iy] + d_kkz[iz]*d_kkz[iz] )/(2.*GAMMA*nxyz) ); 
            
    }    
}

__global__ void __gpe_add_potential_part__(uint it, Complex *psi_in, Complex *psi_add)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    // registers
    Complex lpsi, exp_lv;
    double lrho, lv, dEDFdn;
    Complex Hpsi;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        
        lpsi = psi_in[ixyz]; // psi to register
        lrho = gpe_density(lpsi); // compute density
        compute_potentials_gpe(it, lrho, &dEDFdn, d_user_param, NULL, NULL); //compute potentials defining Hamiltonian
        lv=v_ext(ix, iy, iz, it, NULL, d_user_param, NULL, NULL) + dEDFdn; // external potential + mean field
        
        // K*psi -> K*psi + V*psi
        psi_add[ixyz] = cplxAdd(psi_add[ixyz] , cplxScale(lpsi, lv));
    }    
}

__global__ void __gpe_add_potential_part2__(Complex *psi, Complex *psi_add, double *U)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    
    if(ixyz<nxyz)
    {        
        // K*psi -> K*psi + V*psi
        psi_add[ixyz] = cplxAdd(psi_add[ixyz] , cplxScale(psi[ixyz], U[ixyz]*GAMMA));
    }    
}

__global__ void __gpe_update_psi_imag__(double dt, Complex *psi, Complex *Hpsi)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    
    if(ixyz<nxyz)
    {      
        // K*psi -> K*psi + V*psi
        psi[ixyz] = cplxSub(psi[ixyz] , cplxScale(Hpsi[ixyz], dt));
    }    
}

__global__ void __gpe_imprint_psi__(int it, Complex *psi)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    if(ixyz<nxyz)
    {      
        ixyz2ixiyiz(ixyz,ix,iy,iz,i);
        gpe_modify_psi(ix, iy, iz, it, &(psi[ixyz]), d_user_param, 0, NULL); // modify psi
    }    
}

__global__ void __gpe_imprint_psi_and_compute_potential__(int it, Complex *psi, double *U, Complex *psi_copy)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    Complex lpsi;
    double lrho, dEDFdn;
    
    if(ixyz<nxyz)
    {      
        ixyz2ixiyiz(ixyz,ix,iy,iz,i);
        gpe_modify_psi(ix, iy, iz, it, &(psi[ixyz]), d_user_param, 0, NULL); // modify psi
        lrho = gpe_density(psi[ixyz]); // compute density
        compute_potentials_gpe(it, lrho, &dEDFdn, d_user_param, NULL, NULL); //compute potentials defining Hamiltonian
        U[ixyz]=v_ext(ix, iy, iz, it, NULL, d_user_param, NULL, NULL) + dEDFdn; // external potential + mean field        
    }    
}

__global__ void __gpe_update_psi_by_next_order__(double dt, int k, Complex *psi, Complex *Hpsi, Complex *Hpsi_copy)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
//     double coeff = -1.*dt/(double)k;
    Complex lHpsi;
    Complex coeffC = cplxScale(d_step_coeff, 2.0 / (double)k / (double)GAMMA);
    
    if(ixyz<nxyz)
    {      
//         lHpsi = cplxScale(Hpsi[ixyz], coeff);
        lHpsi = cplxMul(Hpsi[ixyz], coeffC);

        psi[ixyz] = cplxAdd(psi[ixyz] , lHpsi);
        Hpsi_copy[ixyz] = lHpsi;
    }    
}

__global__ void __gpe_multiply_by_kx__(Complex *psi_in, Complex *psi_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 

        psi_out[ixyz]=cplxScale(psi_in[ixyz], d_kkx[ix]/(constgpu(GAMMA*nxyz)) ); 
    }    
}

__global__ void __gpe_multiply_by_ky__(Complex *psi_in, Complex *psi_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 

        psi_out[ixyz]=cplxScale(psi_in[ixyz], d_kky[iy]/(constgpu(GAMMA*nxyz)) ); 
    }    
}

__global__ void __gpe_multiply_by_kz__(Complex *psi_in, Complex *psi_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 

        psi_out[ixyz]=cplxScale(psi_in[ixyz], d_kkz[iz]/(constgpu(GAMMA*nxyz)) ); 
    }    
}

__global__ void __gpe_y_d_dx__(Complex *d_dx_psi, Complex *psi_add, double Omega, int it)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    double _iy;
//     double Omega_z, s, time;
    double Omega_z, _x, _y, _r;
    
    if(ixyz<nxyz)
    {      
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        _iy = (double)(iy) - 1.0*(NY/2);
//         time=d_t0 + d_dt*it;
// //         s = smooth_step(time, 0.0, d_user_param[17], d_user_param[18], 1.0);
// //         Omega_z = Omega + s*d_user_param[19]*0.02875*(exp(-0.005*(iz-45)*(iz-45)) + 0.5*exp(-0.007*(iz-55)*(iz-55))) ;
//         
//         s = sin(time * (2.*M_PI/d_user_param[18])) * smooth_step(time, -1000., 270.0/0.5, 20/0.5, 1.0);
//         
// //         Omega_z = Omega + s*d_user_param[19]*0.02875*(exp(-0.005*(iz-45)*(iz-45)) + 0.5*exp(-0.007*(iz-55)*(iz-55))) ;
//         Omega_z = Omega*(1.0 + s*d_user_param[19]*(exp(-0.009*(iz-45)*(iz-45)) + exp(-0.007*(iz-55)*(iz-55)))) ;
        
//         Omega_z = Omega;
//         _r = 8.0*(exp(-0.005*(iz-45)*(iz-45)) + 0.5*exp(-0.007*(iz-55)*(iz-55)));
//         _x = _r*sin(2.0*M_PI/(1.*NZ/3.)*iz);
//         _y = _r*cos(2.0*M_PI/(1.*NZ/3.)*iz);
//         _iy += _y;
        
        Omega_z = Omega; // no modification
        psi_add[ixyz] = cplxAdd(psi_add[ixyz] , cplxScale(d_dx_psi[ixyz], Omega_z*_iy*GAMMA));
    }    
}

__global__ void __gpe_x_d_dy__(Complex *d_dy_psi, Complex *psi_add, double Omega, int it)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    double _ix;
//     double Omega_z, s, time;
    double Omega_z, _x, _y, _r;
    
    if(ixyz<nxyz)
    {      
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        _ix = (double)(ix) - 1.0*(NX/2);   
//         time=d_t0 + d_dt*it;
// //         s = smooth_step(time, 0.0, d_user_param[17], d_user_param[18], 1.0);
// //         Omega_z = Omega + s*d_user_param[19]*0.02875*(exp(-0.005*(iz-45)*(iz-45)) + 0.5*exp(-0.007*(iz-55)*(iz-55))) ;
//         
//         s = sin(time * (2.*M_PI/d_user_param[18])) * smooth_step(time, -1000., 270.0/0.5, 20/0.5, 1.0);
//         
// //         Omega_z = Omega + s*d_user_param[19]*0.02875*(exp(-0.005*(iz-45)*(iz-45)) + 0.5*exp(-0.007*(iz-55)*(iz-55))) ;
//         Omega_z = Omega*(1.0 + s*d_user_param[19]*(exp(-0.009*(iz-45)*(iz-45)) + exp(-0.007*(iz-55)*(iz-55)))) ;
        
        
//         Omega_z = Omega;
//         _r = 8.0*(exp(-0.005*(iz-45)*(iz-45)) + 0.5*exp(-0.007*(iz-55)*(iz-55)));
//         _x = _r*sin(2.0*M_PI/(1.*NZ/3.)*iz);
//         _y = _r*cos(2.0*M_PI/(1.*NZ/3.)*iz);
//         _ix += _x;  
        
        Omega_z = Omega; // no modification
        psi_add[ixyz] = cplxAdd(psi_add[ixyz] , cplxScale(d_dy_psi[ixyz], -1.0*Omega_z*_ix*GAMMA));
    }    
}

/**
 * Function applies hamiltonian to wave function, i.e. psi_out = H*psi_in
 * */
int gpe_apply_hamiltonian(Complex *psi_in, Complex *psi_out)
{
    cufftResult cufft_result;
    
    __gpe_imprint_psi__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, psi_in);
    
    // apply kinetic part, psi_out as working buffer
    cufft_result=cufftExecZ2Z(gpe_mem.plan, psi_in, psi_out, CUFFT_FORWARD);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    __gpe_multiply_by_kinetic_kernel__<<<gpe_mem.blocks, gpe_mem.threads>>>(psi_out);
    cufft_result=cufftExecZ2Z(gpe_mem.plan, psi_out, psi_out, CUFFT_INVERSE);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;    
    
    // add to kinetic part potential part
    __gpe_add_potential_part__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, psi_in, psi_out);    
    
    return 0;
}

/**
 * Function evolves in imaginary time using Taylor expansion up to n-th order
 * */
int gpe_evolve_taylor(int nt, int order, double Omega)
{
    cufftResult cufft_result;
    int i, n;
    int ierr, err;
    
    // I need extra memory here
    if(gpe_mem.d_wrk3C==NULL) myerrcheck( cudaMalloc( &gpe_mem.d_wrk3C , sizeof(Complex)*nxyz ) );
    if(gpe_mem.d_psi_Lz==NULL) myerrcheck( cudaMalloc( &gpe_mem.d_psi_Lz , sizeof(Complex)*nxyz ) );
    if(gpe_mem.d_psi_copy==NULL) myerrcheck( cudaMalloc( &gpe_mem.d_psi_copy , sizeof(Complex)*nxyz ) );
    
    for(i=0; i<nt; i++)
    {
        // -------------------- PREDICTOR -----------------------
        // make copy of psi
        myerrcheck( cudaMemcpy( gpe_mem.d_psi_copy , gpe_mem.d_psi , sizeof(Complex)*nxyz, cudaMemcpyDeviceToDevice ) );
        
        // imprint the wave-function and compute the potential
        __gpe_imprint_psi_and_compute_potential__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_psi, gpe_mem.d_wrk2R, gpe_mem.d_wrk3C);
        
        for(n=0; n<order; n++)
        {
            // apply hamitonian to psi 
            // apply kinetic part, gpe_mem.d_psi2 as working buffer
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk3C, gpe_mem.d_psi2, CUFFT_FORWARD);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
            __gpe_multiply_by_kinetic_kernel__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2);
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi2, gpe_mem.d_psi2, CUFFT_INVERSE);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;   
            
            // add to kinetic part potential part
            __gpe_add_potential_part2__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_wrk3C, gpe_mem.d_psi2, gpe_mem.d_wrk2R); 
            // gpe_mem.d_psi2 - keeps H*psi
            
            // add contribution from Omega*Lz
            // derivative along d/dx
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk3C, gpe_mem.d_psi_Lz, CUFFT_FORWARD);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
            __gpe_multiply_by_kx__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz);
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz, CUFFT_INVERSE);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;     
            __gpe_y_d_dx__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi2, Omega, gpe_mem.it);
            // derivative along d/dy
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk3C, gpe_mem.d_psi_Lz, CUFFT_FORWARD);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
            __gpe_multiply_by_ky__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz);
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz, CUFFT_INVERSE);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;     
            __gpe_x_d_dy__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi2, Omega, gpe_mem.it);
            
            // update psi by next order - half step
            __gpe_update_psi_by_next_order__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.dt*0.5, n+1, gpe_mem.d_psi, gpe_mem.d_psi2, gpe_mem.d_wrk3C);
        }
        
        // normalize
        myerrchecksimple( gpe_normalize(gpe_mem.d_psi, gpe_mem.d_wrk2R+nxyz), ierr);
        
        // -------------------- CORRECTOR -----------------------
        
        // imprint the wave-function and compute the potential
        __gpe_imprint_psi_and_compute_potential__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_psi, gpe_mem.d_wrk2R, gpe_mem.d_wrk3C);
        
        // revert original psi
        myerrcheck( cudaMemcpy( gpe_mem.d_psi,    gpe_mem.d_psi_copy , sizeof(Complex)*nxyz, cudaMemcpyDeviceToDevice ) );
        myerrcheck( cudaMemcpy( gpe_mem.d_wrk3C , gpe_mem.d_psi_copy , sizeof(Complex)*nxyz, cudaMemcpyDeviceToDevice ) );
        
        for(n=0; n<order; n++)
        {
            // apply hamitonian to psi 
            // apply kinetic part, gpe_mem.d_psi2 as working buffer
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk3C, gpe_mem.d_psi2, CUFFT_FORWARD);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
            __gpe_multiply_by_kinetic_kernel__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2);
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi2, gpe_mem.d_psi2, CUFFT_INVERSE);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;   
            
            // add to kinetic part potential part
            __gpe_add_potential_part2__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_wrk3C, gpe_mem.d_psi2, gpe_mem.d_wrk2R); 
            // gpe_mem.d_psi2 - keeps H*psi
            
            // add contribution from Omega*Lz
            // derivative along d/dx
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk3C, gpe_mem.d_psi_Lz, CUFFT_FORWARD);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
            __gpe_multiply_by_kx__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz);
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz, CUFFT_INVERSE);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;     
            __gpe_y_d_dx__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi2, Omega, gpe_mem.it);
            // derivative along d/dy
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk3C, gpe_mem.d_psi_Lz, CUFFT_FORWARD);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
            __gpe_multiply_by_ky__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz);
            cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi_Lz, gpe_mem.d_psi_Lz, CUFFT_INVERSE);
            if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;     
            __gpe_x_d_dy__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi_Lz, gpe_mem.d_psi2, Omega, gpe_mem.it);
            
            // update psi by next order - full step
            __gpe_update_psi_by_next_order__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.dt, n+1, gpe_mem.d_psi, gpe_mem.d_psi2, gpe_mem.d_wrk3C);
        }
        
        // normalize
        myerrchecksimple( gpe_normalize(gpe_mem.d_psi, gpe_mem.d_wrk2R+nxyz), ierr);
        
        // increase iterator
        gpe_mem.it = gpe_mem.it + 1;
    }
    
    return 0;
}

__global__ void __gpe_compute_vext__(uint it, double *rho, double *wrk)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        wrk[ixyz]=rho[ixyz]*v_ext(ix, iy, iz, it, NULL, d_user_param, NULL, NULL);
    }
}

__global__ void __gpe_compute_vint__(uint it, double *rho, double *wrk)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    if(ixyz<nxyz)
    {
        compute_energy_gpe(it, rho[ixyz], &(wrk[ixyz]), d_user_param, NULL, NULL);
    }
}

__global__ void __gpe_compute_vext_vint__(uint it, double *rho, double *wrk1, double *wrk2)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    double lrho;
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 
        lrho=rho[ixyz];
        wrk1[ixyz]=lrho*v_ext(ix, iy, iz, it, NULL, d_user_param, NULL, NULL);
        compute_energy_gpe(it, lrho, &(wrk2[ixyz]), d_user_param, NULL, NULL);
    }
}

__global__ void __gpe_multiply_by_k2__(Complex *psi_in, Complex *psi_out)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    uint ix, iy, iz, i;
    
    if(ixyz<nxyz)
    {
        ixyz2ixiyiz(ixyz,ix,iy,iz,i); 

        psi_out[ixyz]=cplxScale(psi_in[ixyz], ( d_kkx[ix]*d_kkx[ix] + d_kky[iy]*d_kky[iy] + d_kkz[iz]*d_kkz[iz] )/(constgpu(2.0*GAMMA*nxyz)) ); 
            
    }    
}

__global__ void __gpe_overlap_real__(Complex *psi1, Complex *psi2, double *overlap)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x;
    if(ixyz<nxyz)
    {
        overlap[ixyz]= cplxMulR( cplxConj(psi1[ixyz]),  psi2[ixyz] );
    }
}

int gpe_energy(double *t, double *ekin, double *eint, double *eext)
{
    int ierr;
    cudaError err;
    cufftResult cufft_result;
   
    *t = gpe_mem.t0 + gpe_mem.dt*gpe_mem.it;
    
    if(gpe_mem.beta==0.0) // normalize - otherwise is normalized every time step
    {
        ierr=gpe_normalize_psi();
        if(ierr!=0) return ierr;
    }
    
    // Compute density
    __gpe_compute_density__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_wrk2R);
    
    // Compute <V_ext> and <V_int>
    double * wrk1 = (double *)gpe_mem.d_wrk;
    double * wrk2 = wrk1 + nxyz;
    __gpe_compute_vext_vint__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.it, gpe_mem.d_wrk2R, wrk1, wrk2);
    myerrchecksimple( local_reduction(wrk1, nxyz, wrk1, gpe_mem.threads, 0), ierr );
    myerrchecksimple( local_reduction(wrk2, nxyz, wrk2, gpe_mem.threads, 0), ierr );
    myerrcheck( cudaMemcpy( eext , wrk1 , sizeof(double), cudaMemcpyDeviceToHost ) );
    myerrcheck( cudaMemcpy( eint , wrk2 , sizeof(double), cudaMemcpyDeviceToHost ) );
        
    // Compute <T> - kinetic energy
    cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi, gpe_mem.d_psi2, CUFFT_FORWARD);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    __gpe_multiply_by_k2__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2, gpe_mem.d_psi2);
    cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi2, gpe_mem.d_psi2, CUFFT_INVERSE);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;    
    __gpe_overlap_real__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_psi2, gpe_mem.d_wrk2R);
    myerrchecksimple( local_reduction(gpe_mem.d_wrk2R, nxyz, gpe_mem.d_wrk2R, gpe_mem.threads, 0), ierr );
    myerrcheck( cudaMemcpy( ekin , gpe_mem.d_wrk2R , sizeof(double), cudaMemcpyDeviceToHost ) );
    
    return 0;
}

int gpe_density(double *t, double * density)
{
    __gpe_compute_density__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_wrk2R);
    cudaError err;
    myerrcheck( cudaMemcpy( density , gpe_mem.d_wrk2R , sizeof(double)*nxyz, cudaMemcpyDeviceToHost ) );
    
    *t = gpe_mem.t0 + gpe_mem.dt*gpe_mem.it;
    return 0;
}

int gpe_currents(double *t, double * currents)
{
    int ierr;
    cudaError err;
    cufftResult cufft_result;
   
    *t = gpe_mem.t0 + gpe_mem.dt*gpe_mem.it;
    
    int alloc=0;
    if(gpe_mem.d_wrk3R==NULL) // I need extra memory
    {
        alloc=1;
        myerrcheck( cudaMalloc( &gpe_mem.d_wrk3R , sizeof(double)*nxyz ) );
    }
    
    // move to momentum space
    cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_psi, gpe_mem.d_psi2, CUFFT_FORWARD);  
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;
    
    // Compute d / dx and jx
    __gpe_multiply_by_kx__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2, gpe_mem.d_wrk2);
    cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk2, gpe_mem.d_wrk2, CUFFT_INVERSE);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;    
    __gpe_overlap_real__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_wrk2, gpe_mem.d_wrk3R);    
    myerrcheck( cudaMemcpy( currents , gpe_mem.d_wrk3R , sizeof(double)*nxyz, cudaMemcpyDeviceToHost ) );    
    
    // Compute d / dy and jy
    __gpe_multiply_by_ky__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2, gpe_mem.d_wrk2);
    cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk2, gpe_mem.d_wrk2, CUFFT_INVERSE);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;    
    __gpe_overlap_real__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_wrk2, gpe_mem.d_wrk3R);    
    myerrcheck( cudaMemcpy( currents + nxyz , gpe_mem.d_wrk3R , sizeof(double)*nxyz, cudaMemcpyDeviceToHost ) );
    
    // Compute d / dz and jz
    __gpe_multiply_by_kz__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi2, gpe_mem.d_wrk2);
    cufft_result=cufftExecZ2Z(gpe_mem.plan, gpe_mem.d_wrk2, gpe_mem.d_wrk2, CUFFT_INVERSE);
    if(cufft_result!= CUFFT_SUCCESS) return (int)cufft_result;    
    __gpe_overlap_real__<<<gpe_mem.blocks, gpe_mem.threads>>>(gpe_mem.d_psi, gpe_mem.d_wrk2, gpe_mem.d_wrk3R);    
    myerrcheck( cudaMemcpy( currents + 2 * nxyz , gpe_mem.d_wrk3R , sizeof(double)*nxyz, cudaMemcpyDeviceToHost ) );

    // Free memory
    if(alloc) 
    {
        myerrcheck( cudaFree(gpe_mem.d_wrk3R) );  
        gpe_mem.d_wrk3R = NULL;      
    }


    return 0;
}