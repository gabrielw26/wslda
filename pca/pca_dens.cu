#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>
typedef thrust::complex<double> Complex;
#include "pca_settings.h"


// GW - it is so compilcated...

__global__ void density_one_GPU (int N, int batch, int Nbatch,
			thrust::complex<double>* y, 
			thrust::complex<double>* dx, thrust::complex<double>* dy, thrust::complex<double>* dz,	//=Psi
	                double* nab, double* tauab, thrust::complex<double>* v,
			double* j_a_x, double* j_a_y, double* j_a_z,
                        double* j_b_x, double* j_b_y, double* j_b_z,
			double* betaE
			)
/*
 * INPUT:
 * 	N - liczba wezlow siatki
 *	batch - liczba funkcji falowych
 *	Nband = N * band
 *	y - funkcja falowa
 *	dy - pochodna funkcji falowej
 */
{
 int s = threadIdx.x + blockIdx.x * blockDim.x;
// int krok = gridDim.x * blockDim.x;
 int tt = threadIdx.x;
 int t = 2*tt;
 int i;
 thrust::complex<double>* A;
 __shared__ double A_p[blockSize_d*2];
 A = (thrust::complex<double>*) A_p;

 thrust::complex<double>* psi;
 __shared__ double psi_p[blockSize_d*2 * 2];
 psi = (thrust::complex<double>*) psi_p;

 __shared__ double bE[blockSize_d]; 

 const thrust::complex<double> I2 = thrust::complex<double>(0., 0.5);

//Jakas redukcja?
 if (s<N){
//*
//	psi[t  ] = y[s+(batch-1)*N];
//	psi[t+1] = y[s+Nbatch];
//	bE[tt] = betaE[0];
	nab[s] = 0;//psi[t].real();//bE[tt] * thrust::norm(psi[t]);
	v[s] = 0;//.5*(1-2*bE[tt]) * psi[t] * thrust::conj(psi[t+1]);
        nab[s+N] = 0;//(1-bE[tt]) * thrust::norm(psi[t+1]);
//*/
	for (i=0; i<batch; i++){
	        psi[t  ] = y[s+i*N];
        	psi[t+1] = y[s+i*N+Nbatch];
	        bE[tt] = betaE[i];
                nab[s] += bE[tt] * thrust::norm(psi[t]) * 10000.;
                v[s] += .5*(1.-2.*bE[tt]) * psi[t] * thrust::conj(psi[t+1]) * 10000.;
                nab[s+N] += (1.-bE[tt]) * thrust::norm(psi[t+1]) * 10000.;
	}
nab[s] = nab[s] / 10000.;
v[s] = v[s]  / 10000.;
nab[s+N] = nab[s+N] / 10000.;

 }

 if (s < N){
//	for (j=0; j<3; j++){
//*
//        psi[t  ] = thrust::norm(dx[s      ]) + thrust::norm(dy[s      ]) + thrust::norm(dz[s      ]);
//        psi[t+1] = thrust::norm(dx[s+Nbatch]) + thrust::norm(dy[s+Nbatch]) + thrust::norm(dz[s+Nbatch]);
//        bE[tt] = betaE[0];
	tauab[s    ] = 0;//bE[tt] * psi[t].real() * 10000;
	tauab[s + N] = 0;//(1-bE[tt]) * psi[t+1].real() * 10000;

//*/
 	for (i=0; i<batch; i++){
		psi[t  ] = thrust::norm(dx[s       +i*N]) + thrust::norm(dy[s       +i*N]) + thrust::norm(dz[s       +i*N]);
		psi[t+1] = thrust::norm(dx[s+Nbatch+i*N]) + thrust::norm(dy[s+Nbatch+i*N]) + thrust::norm(dz[s+Nbatch+i*N]);
	        bE[tt] = betaE[i];
		tauab[s    ] += bE[tt] * psi[t].real() * 10000.;
	        tauab[s + N] += (1.-bE[tt]) * psi[t+1].real() * 10000.;
	}
tauab[s  ] = tauab[s  ] / 10000.;
tauab[s+N] = tauab[s+N] / 10000.;
 }
 if (s<N){
//	bE[tt] = betaE[0];
//	psi[t  ] = y[s];
	j_a_x[s] = 0;
	j_a_y[s] = 0;
	j_a_z[s] = 0;
	for (i=0; i<batch; i++){
		bE[tt] = betaE[i]     ;
		psi[tt  ] =  y[s+i*N];
//	        A[tt] = thrust::conj(psi[tt]) * dx[s+i*N];
//	        j_a_x[s] += bE[tt] * (I2*(A[tt] - thrust::conj(A[tt]))).real() * 10000;
j_a_x[s] -= bE[tt] * (thrust::conj(psi[tt]) * dx[s+i*N]).imag() * 10000;
//	        A[tt] = thrust::conj(psi[tt]) * dy[s+i*N];
//	        j_a_y[s] += bE[tt] * (I2*(A[tt] - thrust::conj(A[tt]))).real() * 10000;
j_a_y[s] -= bE[tt] * (thrust::conj(psi[tt]) * dy[s+i*N]).imag() * 10000;
//        	A[tt] = thrust::conj(psi[tt]) * dz[s+i*N];
//	        j_a_z[s] += bE[tt] * (I2*(A[tt] - thrust::conj(A[tt]))).real() * 10000;
j_a_z[s] -= bE[tt] * (thrust::conj(psi[tt]) * dz[s+i*N]).imag() * 10000;

	}
        j_a_x[s] = j_a_x[s] / 10000.;
        j_a_y[s] = j_a_y[s] / 10000.;
        j_a_z[s] = j_a_z[s] / 10000.;
        j_b_x[s] = 0;
        j_b_y[s] = 0;
        j_b_z[s] = 0;
        for (i=0; i<batch; i++){
                bE[tt] = (1.-betaE[i])     ;
                psi[tt  ] = y[s+i*N+Nbatch];
//                A[tt] = thrust::conj(psi[tt]) * dx[s+i*N+Nbatch];
//                j_b_x[s+N] += bE[tt] * (I2*(A[tt] - thrust::conj(A[tt]))).real() * 10000;
j_a_x[s] -= bE[tt] * (thrust::conj(psi[tt]) * dx[s+i*N+Nbatch]).imag() * 10000;
//                A[tt] = thrust::conj(psi[tt]) * dy[s+i*N+Nbatch];
//                j_b_y[s+N] += bE[tt] * (I2*(A[tt] - thrust::conj(A[tt]))).real() * 10000;
j_a_y[s] -= bE[tt] * (thrust::conj(psi[tt]) * dy[s+i*N+Nbatch]).imag() * 10000;
//                A[tt] = thrust::conj(psi[tt]) * dz[s+i*N+Nbatch];
//                j_b_z[s+N] += bE[tt] * (I2*(A[tt] - thrust::conj(A[tt]))).real() * 10000;
j_a_z[s] -= bE[tt] * (thrust::conj(psi[tt]) * dz[s+i*N+Nbatch]).imag() * 10000;
        }
        j_b_x[s] = j_b_x[s] / 10000.;
        j_b_y[s] = j_b_y[s] / 10000.;
        j_b_z[s] = j_b_z[s] / 10000.;
 }

}

int calc_dens (int N, int batch, int Nbatch,
                        double* psi, 
                        double* dpsiX, double* dpsiY, double* dpsiZ,
                        double* nab, double* tauab, double* v,
                        double* jaX, double* jaY, double* jaZ,
                        double* jbX, double* jbY, double* jbZ,
                        double* betaE)
{
//dPSI
    
 int nblocks = (int)ceil((float)Nbatch/blockSize_d);
 density_one_GPU <<< nblocks, blockSize_d >>> 
                        (N, batch, Nbatch, 
                                (thrust::complex<double>*) psi, 
                                (thrust::complex<double>*) dpsiX, 
                                (thrust::complex<double>*) dpsiY, 
                                (thrust::complex<double>*) dpsiZ, 
                                nab, tauab, (thrust::complex<double>*) v, 
                                jaX, jaY, jaZ,
                                jbX, jbY, jbZ,
                                betaE);

//  density_all (N, nab, nab_l); 
    return 0;
}

// ================================================================================================
// ========================================= calculate_densities ==================================
// ================================================================================================
__global__ void kernel_calculate_densities(size_t n, Complex *wf, 
                                         Complex *wf_d_dx, Complex *wf_d_dy, Complex *wf_d_dz, 
                                         double *fbetaEn,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;
    double taua=0.0, taub=0.0;
    Complex _nu=Complex(0.0, 0.0);
    double jax=0.0, jay=0.0, jaz=0.0;
    double jbx=0.0, jby=0.0, jbz=0.0;

    Complex u, v, wfdx, wfdy, wfdz;
    double fbEn, fbmEn;
    #define DENS_FACTOR_M 10000.
    
    size_t iwf;
    
    if(ixyz<NXYZ)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            
            // read u and v
            u=wf[       iwf*NXYZ+ixyz];
            v=wf[n*NXYZ+iwf*NXYZ+ixyz];
            
            // form na, nb, nu
            na+=thrust::norm(u)*fbEn;
            nb+=thrust::norm(v)*fbmEn;
#ifdef SPINSYMMETRY_MODE
            _nu+=u*thrust::conj(v)*fbmEn;
#else
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)/2.0;
#endif
            
            // read derivatives from u 
            wfdx=wf_d_dx[       iwf*NXYZ+ixyz];
            wfdy=wf_d_dy[       iwf*NXYZ+ixyz];
            wfdz=wf_d_dz[       iwf*NXYZ+ixyz];
            
            // form taua and j_a
            taua+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbEn;
            jax+=(thrust::conj(u)*wfdx).imag()*fbEn;
            jay+=(thrust::conj(u)*wfdy).imag()*fbEn;
            jaz+=(thrust::conj(u)*wfdz).imag()*fbEn;
            
            // read derivatives from v 
            wfdx=wf_d_dx[n*NXYZ+iwf*NXYZ+ixyz];
            wfdy=wf_d_dy[n*NXYZ+iwf*NXYZ+ixyz];
            wfdz=wf_d_dz[n*NXYZ+iwf*NXYZ+ixyz];
            
            // form taua and j_a
            taub+=(thrust::norm(wfdx)+thrust::norm(wfdy)+thrust::norm(wfdz))*fbmEn;
            jbx-=(thrust::conj(v)*wfdx).imag()*fbmEn;
            jby-=(thrust::conj(v)*wfdy).imag()*fbmEn;
            jbz-=(thrust::conj(v)*wfdz).imag()*fbmEn;
        }
        
        // save result to global memory
        rho_a[ixyz]=na/DENS_FACTOR_M;
        rho_b[ixyz]=nb/DENS_FACTOR_M;
        tau_a[ixyz]=taua/DENS_FACTOR_M;
        tau_b[ixyz]=taub/DENS_FACTOR_M;
        nu[ixyz]=_nu/DENS_FACTOR_M;
        j_a_x[ixyz]=jax/DENS_FACTOR_M;
        j_a_y[ixyz]=jay/DENS_FACTOR_M;
        j_a_z[ixyz]=jaz/DENS_FACTOR_M;
        j_b_x[ixyz]=jbx/DENS_FACTOR_M;
        j_b_y[ixyz]=jby/DENS_FACTOR_M;
        j_b_z[ixyz]=jbz/DENS_FACTOR_M;
    }
}

__global__ void kernel_calculate_densities_limited(size_t n, Complex *wf,
                                         double *fbetaEn,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;
    Complex _nu=Complex(0.0, 0.0);

    Complex u, v;
    double fbEn, fbmEn;
    #define DENS_FACTOR_M 10000.
    
    size_t iwf;
    
    if(ixyz<NXYZ)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            
            // read u and v
            u=wf[       iwf*NXYZ+ixyz];
            v=wf[n*NXYZ+iwf*NXYZ+ixyz];
            
            // form na, nb, nu
            na+=thrust::norm(u)*fbEn;
            nb+=thrust::norm(v)*fbmEn;
#ifdef SPINSYMMETRY_MODE
            _nu+=u*thrust::conj(v)*fbmEn;
#else
            _nu+=u*thrust::conj(v)*(fbmEn-fbEn)/2.0;
#endif
        }
        
        // save result to global memory
        rho_a[ixyz]=na/DENS_FACTOR_M;
        rho_b[ixyz]=nb/DENS_FACTOR_M;
        tau_a[ixyz]=0.0;
        tau_b[ixyz]=0.0;
        nu[ixyz]=_nu/DENS_FACTOR_M;
        j_a_x[ixyz]=0.0;
        j_a_y[ixyz]=0.0;
        j_a_z[ixyz]=0.0;
        j_b_x[ixyz]=0.0;
        j_b_y[ixyz]=0.0;
        j_b_z[ixyz]=0.0;
    }
}
    
/**
 * Function computes energy.
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf array with wave-functions (INPUT)
 * @param wf_d_dx derivative with respect to dx (INPUT)
 * @param wf_d_dy derivative with respect to dy (INPUT)
 * @param wf_d_dz derivative with respect to dz (INPUT)
 * @param fbetaEn weight of wave-function (INPUT)
 * @param d_densites (OUTPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NXYZ,
 *                          nu - double complex array of size NXYZ,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NXYZ
 *                   In total size of d_densites is 12*NXYZ
 * @param gradients_computed if 1 then gradients are computed, otherwise only normal and anomalus density will be computed
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int calculate_densities(int n, cufftDoubleComplex *wf,
                            cufftDoubleComplex *wf_d_dx, cufftDoubleComplex *wf_d_dy, cufftDoubleComplex *wf_d_dz, 
                            double *d_fbetaEn, 
                            double *d_densities,
                            int gradients_computed, int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXYZ/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    double *rho_a = (double *)(d_densities +  0*NXYZ);
    double *rho_b = (double *)(d_densities +  1*NXYZ);
    double *tau_a = (double *)(d_densities +  2*NXYZ);
    double *tau_b = (double *)(d_densities +  3*NXYZ);
    Complex *nu   =(Complex *)(d_densities +  4*NXYZ);
    double *j_a_x = (double *)(d_densities +  6*NXYZ);
    double *j_a_y = (double *)(d_densities +  7*NXYZ);
    double *j_a_z = (double *)(d_densities +  8*NXYZ);
    double *j_b_x = (double *)(d_densities +  9*NXYZ);
    double *j_b_y = (double *)(d_densities + 10*NXYZ);
    double *j_b_z = (double *)(d_densities + 11*NXYZ);
    
    
    if(gradients_computed) // computation of all densities
    {
        kernel_calculate_densities<<<nblocks, nthreads>>>(n, (Complex *)wf, 
                                            (Complex *)wf_d_dx, (Complex *)wf_d_dy, (Complex *)wf_d_dz, 
                                            d_fbetaEn,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);
    }
    else // only normal and anomalus density will be computed, other are set to zero
    {
        kernel_calculate_densities_limited<<<nblocks, nthreads>>>(n, (Complex *)wf, 
                                            d_fbetaEn,
                                            rho_a, rho_b,
                                            tau_a, tau_b,
                                            nu,
                                            j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);        
    }

    return 0;
}

// ----------------------------------------------------------------------------------------
// ----------------------------------calculate_densities_weighted -------------------------
// ----------------------------------------------------------------------------------------
__global__ void kernel_calculate_densities_weighted(size_t n, Complex *wf,
                                         double *fbetaEn, double *weights,
                                         double *rho_a, double *rho_b,
                                         double *tau_a, double *tau_b,
                                         Complex *nu,
                                         double *j_a_x, double *j_a_y, double *j_a_z, double *j_b_x, double *j_b_y, double *j_b_z
                                        )
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point

    double na=0.0, nb=0.0;

    Complex u, v;
    double fbEn, fbmEn, wcnt;
    #define DENS_FACTOR_M 10000.
    
    size_t iwf;
    
    if(ixyz<NXYZ)
    {
        // reduce over each wave-function
        for(iwf=0; iwf<n; iwf++)
        {    
            // weight
            fbEn=fbetaEn[iwf]*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;

            
            wcnt=weights[iwf]; // add external weight 
            
            // read u and v
            u=wf[       iwf*NXYZ+ixyz];
            v=wf[n*NXYZ+iwf*NXYZ+ixyz];
            
            // form na, nb, nu
            na+=thrust::norm(u)*fbEn *wcnt;
            nb+=thrust::norm(v)*fbmEn*wcnt;
        }
        
        // save result to global memory
        rho_a[ixyz]=na/DENS_FACTOR_M;
        rho_b[ixyz]=nb/DENS_FACTOR_M;
        tau_a[ixyz]=0.0;
        tau_b[ixyz]=0.0;
        nu[ixyz]=0.0;
        j_a_x[ixyz]=0.0;
        j_a_y[ixyz]=0.0;
        j_a_z[ixyz]=0.0;
        j_b_x[ixyz]=0.0;
        j_b_y[ixyz]=0.0;
        j_b_z[ixyz]=0.0;
    }
}


/**
 * Function computes density.
 * @param n  number of wave-functions (u,v pairs) to process
 * @param wf array with wave-functions (INPUT)
 * @param fbetaEn weight of wave-function (INPUT)
 * @param weights for density computation (INPUT)
 * @param d_densites (OUTPUT)
 *                   collective array with densities [rho_a, rho_b, tau_a, tau_b, nu, j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z]  
 *                   where: rho_a, rho_b, tau_a, tau_b - double arrays of size NXYZ,
 *                          nu - double complex array of size NXYZ,
 *                          j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z - double arrays of size NXYZ
 *                   In total size of d_densites is 12*NXYZ
 * @param nthreads number of threads per block
 * @return 0 - OK, otherwise ERROR 
 * */
extern "C" int calculate_densities_weighted(int n, cufftDoubleComplex *wf,
                            double *d_fbetaEn, 
                            double *d_weights, 
                            double *d_densities,
                            int nthreads)
{
    // number of blocks
    int nblocks = (int)ceil((float)NXYZ/nthreads);
    
    // Set pointers for to simplify notation
    // densities 
    double *rho_a = (double *)(d_densities +  0*NXYZ);
    double *rho_b = (double *)(d_densities +  1*NXYZ);
    double *tau_a = (double *)(d_densities +  2*NXYZ);
    double *tau_b = (double *)(d_densities +  3*NXYZ);
    Complex *nu   =(Complex *)(d_densities +  4*NXYZ);
    double *j_a_x = (double *)(d_densities +  6*NXYZ);
    double *j_a_y = (double *)(d_densities +  7*NXYZ);
    double *j_a_z = (double *)(d_densities +  8*NXYZ);
    double *j_b_x = (double *)(d_densities +  9*NXYZ);
    double *j_b_y = (double *)(d_densities + 10*NXYZ);
    double *j_b_z = (double *)(d_densities + 11*NXYZ);
    

    kernel_calculate_densities_weighted<<<nblocks, nthreads>>>(n, (Complex *)wf,
                                        d_fbetaEn, d_weights, 
                                        rho_a, rho_b,
                                        tau_a, tau_b,
                                        nu,
                                        j_a_x, j_a_y, j_a_z, j_b_x, j_b_y, j_b_z);        


    return 0;
}

