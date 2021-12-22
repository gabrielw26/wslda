// Author: Gabriel Wlazlowski

// Code for computing matrix elements of hamiltonian
#include <stdlib.h>
#include <stddef.h>

// global variables
extern double *dc_params; /* Declaration of the variable */
extern size_t dc_extra_data_size;
extern void *dc_extra_data;
extern double dc_mu_a;
extern double dc_mu_b;
extern double dc_ec;

extern int wsldapid; // process id - global variable

#include "pca_settings.h"
#include "pca_macro.h"
#include "s3dpca_me.h"

#include "pca_utils.h"

#define BLOCKSIZE (NX*NY*NZ)

// List of functions from "problem-definition.h"
double v_ext(int ix, int iy, int iz, int it, int spin, double *params, size_t extra_data_size, void *extra_data);
double complex delta_ext(int ix, int iy, int iz, int it, double complex delta, double *params, size_t extra_data_size, void *extra_data);
double velocity_ext(int ix, int iy, int iz, int it, int spin, int coordinate, double *params, size_t extra_data_size, void *extra_data);

#include "aslda_edf.h"

#define Complex(a,b) (a + I*b)
#define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a)) 
// Number of self-consistent iterations for U and delta computation
#undef UD_SCITERS
#define UD_SCITERS 1000
// Mixing parameter for self-consistent algorithm - fraction of new solution used for mixing
#undef UD_MIX_COEFF
#define UD_MIX_COEFF 0.25
#define UD_EPSILON 1.0e-12

int indxl2g_(int*, int*, int*, int*, int*);

// ==========================================================================
// ============================== FUNCTIONS =================================
// ==========================================================================
/** 
 * 1D kinetic energy operator matrix elements
 * see: arxiv.org/abs/1301.7354, Eq.(23)
 * @param i coordinate index i=0,...,N-1
 * @param j coordinate index j=0,...,N-1
 * @param N lattice size
 * @param a lattice spacing
 * @return value of matrix element un units fm^-1
 * */
double k_1D(int k, int l, int N, double a)
{
    if(k==l)
        return M_PI*M_PI*( 1.0+2.0/(N*N) ) / (6.0*a*a);
    else
    {
        double pm_one=1.0;
        if(abs(k-l)%2 == 1) pm_one=-1.0;
        double _sinkl = sin(M_PI*(k-l)/N);
        
        return M_PI*M_PI*pm_one / ( a*a*N*N*_sinkl*_sinkl );
    }
}

/**
 * Function computes matrix elements of BdG type hamiltonian
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials recomputed potentials (INPUT)
 * @param mdfft metadata for ffts plans execution
 * @param h hamitonian matrix of size [2*NXYZ,2*NXYZ] (OUTPUT, only local part)
 * @param me_d_dx matrix elements of (-i*d/dx) operator, matrix of size [NX x NX] (INPUT)
 * @param me_d_dy matrix elements of (-i*d/dy) operator, matrix of size [NY x NY] (INPUT)
 * @param me_d_dz matrix elements of (-i*d/dz) operator, matrix of size [NZ x NZ] (INPUT)
 * */
int compute_matrix_elements_3d(metadata_s3dpca_grid *bgrid, int it, wslda_density h_densities, wslda_potential h_potentials, metadata_s3dpca_fft *mdfft, double complex *h, double complex * me_d_dx, double complex * me_d_dy, double complex * me_d_dz)
{
    // densities - decode 
    double *rho_a = h_densities.rho_a;
    double *rho_b = h_densities.rho_b;
    double *tau_a = h_densities.tau_a;
    double *tau_b = h_densities.tau_b;
    double complex *nu = h_densities.nu;
    double *j_a_x = h_densities.j_a_x;
    double *j_a_y = h_densities.j_a_y;
    double *j_a_z = h_densities.j_a_z;
    double *j_b_x = h_densities.j_b_x;
    double *j_b_y = h_densities.j_b_y;
    double *j_b_z = h_densities.j_b_z;
    
    // potentials - decode
    double *V_a = h_potentials.V_a;
    double *V_b = h_potentials.V_b;
    double complex *delta = h_potentials.delta;
    
    double p, alph_1, alph_2;
    
    double Fx1, Fx2, Fy1, Fy2, Fz1, Fz2, t7, vext;          
    
    // iterate over all matrix elemnts
    
    int ix1, iy1, iz1, ixyz1; // row iterator
    int ix2, iy2, iz2, ixyz2; // column iterator
    int ix, iy, iz, ixyz; // global iterator
    int ci, ri; // column and row iterator - global index
    int li, lj, ij; // local indices
    int ZERO = 0, ti;
    
    // compute gradient of effective mass
    double *laplace_alpha_a;
    double *laplace_alpha_b;
    cppmallocl(laplace_alpha_a,NXYZ,double);
    cppmallocl(laplace_alpha_b,NXYZ,double);
    int ierr;
    if(bgrid->nip*bgrid->niq<NXYZ/2) return -199; // check if enough memory
    double *wrk_dble = (double *)(h); // only for temporary calculations 
    
    ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++) 
    {
        wrk_dble[ixyz] = h_potentials.alpha_a[ixyz]; 
        ixyz++;
    }
    ierr = compute_laplace_real_f(wrk_dble, laplace_alpha_a, mdfft);
    if(ierr!=0) return (100+ierr);
    
    ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++) 
    {
        wrk_dble[ixyz] = h_potentials.alpha_b[ixyz];
        ixyz++;
    }
    ierr = compute_laplace_real_f(wrk_dble, laplace_alpha_b, mdfft);
    if(ierr!=0) return (200+ierr);
    
    
    // Reset matrix elements
    for(ixyz=0; ixyz<bgrid->nip*bgrid->niq; ixyz++) h[ixyz] = 0.0 + I*0.0;
    
    // MATRIX structure
    // |   ha    |   Delta  |
    // ----------------------
    // | Delta^* |   -hb^*  |
    
    // kinetic part
    // K(ixiyiz,jxjyjz) =  k(ix,jx)*delta(iy,jy)*delta(iz,jz)
    //                     +delta(ix,jx)*k(iy,jy)*delta(iz,jz)
    //                     +delta(ix,jx)*delta(iy,jy)*k(iz,jz)
    // where k(ix,jx) is matrix element of 1D kinetic energy (see: arxiv.org/abs/1301.7354, Eq.(23))
    
    for(lj=0; lj<bgrid->niq; lj++) // column-major iteration: over local index (column)
    {
        for(li=0; li<bgrid->nip; li++) // over local index (row)
        {
            ij=li + bgrid->nip * lj; // local index of elemnt
            ixyz1 = li+1; ixyz2=lj+1; // conversion to fortran standard
            // find indices in global matrix: row and colummn
            ri = indxl2g_( &ixyz1, &bgrid->mb, &bgrid->ip, &ZERO, &bgrid->p )-1; // back to C standard
            ci = indxl2g_( &ixyz2, &bgrid->nb, &bgrid->iq, &ZERO, &bgrid->q )-1; // back to C standard            
            
            // just in case: consitency check
            if(ri>=2*NXYZ) {wprintf("ERROR: ri=%d >= 2*NXYZ\n", ri); return -99;}
            if(ci>=2*NXYZ) {wprintf("ERROR: ci=%d >= 2*NXYZ\n", ci); return -98;}
            
            if(ci<NXYZ && ri<NXYZ) // part: |   ha    |
            {
                ixyz1 = ri;
                ixyz2 = ci;
                
                // decode to cartesian coordinates
                ixyz2ixiyiz(ixyz1,ix1,iy1,iz1,ti);
                ixyz2ixiyiz(ixyz2,ix2,iy2,iz2,ti);
                
                vext=v_ext(ix1,iy1,iz1,it,SPINA,dc_params,dc_extra_data_size,dc_extra_data); // mean field does not contain external potential 
                
                alph_1 = h_potentials.alpha_a[ixyz1]; 
                alph_2 = h_potentials.alpha_a[ixyz2]; 
                
                // diagonal part: V_a-mu_a + 0.5*alpha_a*kz^2
                if(ixyz1==ixyz2) h[ij] += V_a[ixyz1] + vext - dc_mu_a + 0.25*laplace_alpha_a[ixyz1];
                
                // Kinetic term K_a                
                // k(ix,jx)*delta(iy,jy)*delta(iz,jz)
                if(iy1==iy2 && iz1==iz2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX);
                // delta(ix,jx)*k(iy,jy)*delta(iz,jz)
                if(ix1==ix2 && iz1==iz2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY);
                // delta(ix,jx)*delta(iy,jy)*k(iz,jz)
                if(ix1==ix2 && iy1==iy2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(iz1, iz2, NZ, DZ);
                
                // vector field contribution
                Fx1 = h_potentials.A_a_x[ixyz1];
                Fy1 = h_potentials.A_a_y[ixyz1];
                Fz1 = h_potentials.A_a_z[ixyz1];
                
                Fx2 = h_potentials.A_a_x[ixyz2];
                Fy2 = h_potentials.A_a_y[ixyz2];
                Fz2 = h_potentials.A_a_z[ixyz2];
                
                // k(ix,jx)*delta(iy,jy)*delta(iz,jz)
                if(iy1==iy2 && iz1==iz2) h[ij] += me_d_dx[ix1 + ix2*NX] * (Fx1 + Fx2) * ( 0.5); // (+1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                // delta(ix,jx)*k(iy,jy)*delta(iz,jz)
                if(ix1==ix2 && iz1==iz2) h[ij] += me_d_dy[iy1 + iy2*NY] * (Fy1 + Fy2) * ( 0.5); // (+1.0) because me_d_dy keeps matrix elements of (-i d/dy)
                // delta(ix,jx)*delta(iy,jy)*k(iz,jz)
                if(ix1==ix2 && iy1==iy2) h[ij] += me_d_dz[iz1 + iz2*NZ] * (Fz1 + Fz2) * ( 0.5); // (+1.0) because me_d_dz keeps matrix elements of (-i d/dy)

                
                // contribution from velocity_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = velocity_ext(ix1, iy1, iz1, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = velocity_ext(ix2, iy2, iz2, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2 && iz1==iz2) h[ij] += me_d_dx[ix1 + ix2*NX] * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = velocity_ext(ix1, iy1, iz1, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = velocity_ext(ix2, iy2, iz2, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2 && iz1==iz2) h[ij] += me_d_dy[iy1 + iy2*NY] * (Fy1 + Fy2) * ( -0.5); // note me_d_dy keeps matrix elements of (-i d/dy)   
                
                Fz1 = velocity_ext(ix1, iy1, iz1, it, SPINA, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fz2 = velocity_ext(ix2, iy2, iz2, it, SPINA, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2 && iy1==iy2) h[ij] += me_d_dz[iz1 + iz2*NZ] * (Fz1 + Fz2) * ( -0.5); // note me_d_dz keeps matrix elements of (-i d/dz)   
                
            }
            else if(ci>=NXYZ && ri<NXYZ) // part: |   Delta  |
            {
                ixyz1 = ri;
                ixyz2 = ci-NXYZ;             
                
                // decode to cartesian coordinates
                ixyz2ixiyiz(ixyz1,ix1,iy1,iz1,ti);
                ixyz2ixiyiz(ixyz2,ix2,iy2,iz2,ti);
                
                // diagonal part: Delta(r)
                if(ixyz1==ixyz2) h[ij] += delta[ixyz1] + delta_ext(ix1, iy1, iz1, it, delta[ixyz1], dc_params,dc_extra_data_size,dc_extra_data);
            }
            else if(ci<NXYZ && ri>=NXYZ) // part: | Delta^* |
            {
                ixyz1 = ri-NXYZ;
                ixyz2 = ci;          
                
                // decode to cartesian coordinates
                ixyz2ixiyiz(ixyz1,ix1,iy1,iz1,ti);
                ixyz2ixiyiz(ixyz2,ix2,iy2,iz2,ti);
                
                // diagonal part: (Delta(r))^*
                if(ixyz1==ixyz2) h[ij] += conj(delta[ixyz1] + delta_ext(ix1, iy1, iz1, it, delta[ixyz1], dc_params,dc_extra_data_size,dc_extra_data));
            }
            else // part: |   -hb^*  |
            {
                ixyz1 = ri-NXYZ;
                ixyz2 = ci-NXYZ;
                
                // decode to cartesian coordinates
                ixyz2ixiyiz(ixyz1,ix1,iy1,iz1,ti);
                ixyz2ixiyiz(ixyz2,ix2,iy2,iz2,ti);
                
                vext=v_ext(ix1,iy1,iz1,it,SPINB,dc_params,dc_extra_data_size,dc_extra_data); // mean field does not contain external potential
                
                alph_1 = h_potentials.alpha_b[ixyz1]; 
                alph_2 = h_potentials.alpha_b[ixyz2]; 
                
                // diagonal part: V_b-mu_b + 0.5*alpha_b*kz^2
                if(ixyz1==ixyz2) h[ij] -= V_b[ixyz1] + vext - dc_mu_b + 0.25*laplace_alpha_b[ixyz1]; // NOTE -= operator has minus!
                
                // Kinetic term -K_b
                // k(ix,jx)*delta(iy,jy)*delta(iz,jz)
                if(iy1==iy2 && iz1==iz2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX); // NOTE -= operator has minus!
                // delta(ix,jx)*k(iy,jy)*delta(iz,jz)
                if(ix1==ix2 && iz1==iz2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY); // NOTE -= operator has minus!
                // delta(ix,jx)*delta(iy,jy)*k(iz,jz)
                if(ix1==ix2 && iy1==iy2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(iz1, iz2, NZ, DZ); // NOTE -= operator has minus!
                
                
                // vector field contribution
                Fx1 = h_potentials.A_b_x[ixyz1];
                Fy1 = h_potentials.A_b_y[ixyz1];
                Fz1 = h_potentials.A_b_z[ixyz1];
                
                Fx2 = h_potentials.A_b_x[ixyz2];
                Fy2 = h_potentials.A_b_y[ixyz2];
                Fz2 = h_potentials.A_b_z[ixyz2];
                
                // k(ix,jx)*delta(iy,jy)*delta(iz,jz)
                if(iy1==iy2 && iz1==iz2) h[ij] -= conj(me_d_dx[ix1 + ix2*NX]) * (Fx1 + Fx2) * ( 0.5); // (+1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                // delta(ix,jx)*k(iy,jy)*delta(iz,jz)
                if(ix1==ix2 && iz1==iz2) h[ij] -= conj(me_d_dy[iy1 + iy2*NY]) * (Fy1 + Fy2) * ( 0.5); // (+1.0) because me_d_dy keeps matrix elements of (-i d/dy)
                // delta(ix,jx)*delta(iy,jy)*k(iz,jz)
                if(ix1==ix2 && iy1==iy2) h[ij] -= conj(me_d_dz[iz1 + iz2*NZ]) * (Fz1 + Fz2) * ( 0.5); // (+1.0) because me_d_dz keeps matrix elements of (-i d/dy)

                
                // contribution from v_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = velocity_ext(ix1, iy1, iz1, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = velocity_ext(ix2, iy2, iz2, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2 && iz1==iz2) h[ij] -= conj(me_d_dx[ix1 + ix2*NX]) * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = velocity_ext(ix1, iy1, iz1, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = velocity_ext(ix2, iy2, iz2, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2 && iz1==iz2) h[ij] -= conj(me_d_dy[iy1 + iy2*NY]) * (Fy1 + Fy2) * ( -0.5); // note me_d_dy keeps matrix elements of (-i d/dy)   
                
                Fz1 = velocity_ext(ix1, iy1, iz1, it, SPINB, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fz2 = velocity_ext(ix2, iy2, iz2, it, SPINB, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2 && iy1==iy2) h[ij] -= conj(me_d_dz[iz1 + iz2*NZ]) * (Fz1 + Fz2) * ( -0.5); // note me_d_dz keeps matrix elements of (-i d/dz)   
            }
            
        } // for(li=0; li<bgrid->nip; li++)
//         wprintf("ci=%d\n", ci); fflush(stdout);
    } // for(lj=0; lj<bgrid->niq; lj++)
    
    
    // clear memory
    free(laplace_alpha_a);
    free(laplace_alpha_b);
    
    return 0;
}


/**
 * Function that computes matrix elements of -i*d/dx operator
 * @param nx lattice size (INPUT)
 * @param h matrix of size [nx*nx] with computed matrix elements (OUTPUT)
 *          <x1| -id/dx |x2> is given in matrix element h[x1 + nx*x2]
 * */
int compute_matrix_elements_of_momentum_operator(int nx, double dx, double complex *me)
{
    int ix, ci, ri; // iterator
    double kx;
    
    // allocate memory 
    double complex *fft1;
    cppmallocl(fft1, nx, double complex); // for working area of plan
    
    // fftw plan
    fftw_plan plan_f_1d;
    fftw_plan plan_b_1d;
    
    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    plan_f_1d = fftw_plan_dft_1d(nx, fft1, fft1, FFTW_FORWARD, USE_FFTW_PLANNER);
    plan_b_1d = fftw_plan_dft_1d(nx, fft1, fft1, FFTW_BACKWARD, USE_FFTW_PLANNER);     
    
    // reset matrix elements
    for(ix=0; ix<nx*nx; ix++) me[ix] = 0.0 + I*0.0;
    
    for(ci=0; ci<nx; ci++) // column-major iteration fashion, for each column do:
    {
        for(ix=0; ix<nx; ix++) fft1[ix]=0.0+I*0.0;
        fft1[ci] = 1.0/(dx*nx); // normalization factor already included 
        fftw_execute(plan_f_1d); // to momentum space
        
        // multiply by momentum
        for(ix=0; ix<nx; ix++)
        {
            // extract momentum
            if(ix<nx/2) kx=2.*M_PI/( ( double )nx * dx ) * ( double )(ix   );
            else        kx=2.*M_PI/( ( double )nx * dx ) * ( double )(ix-nx);
            
            if(ix==nx/2) kx=0.0;
            
            fft1[ix]*=kx; // -i*(i*kx) = kx
        }
        
        fftw_execute(plan_b_1d); // to coordinate space
        
        // copy as matrix element
        for(ri=0; ri<nx; ri++) // for each row in column
            me[ri + nx*ci] = fft1[ri];
    }
    
    // just in case - check if hermitian
    double complex x1x2, x2x1, diff;
    int hermitian_violated=0;
    for(ri=0; ri<nx; ri++) for(ci=0; ci<nx; ci++)
    {
        x1x2 =      me[ri + nx*ci] ;
        x2x1 = conj(me[ci + nx*ri]);
        diff = x1x2 - x2x1;
        
        if(fabs(creal(diff))>1.0e-14 || fabs(cimag(diff))>1.0e-14) 
        {
            wprintf("# ERROR[compute_matrix_elements_of_momentum_operator]: hermitian problem: %6d %6d (%f,%f) <=> (%f,%f)\n",
                ri, ci, creal(x1x2), cimag(x1x2), creal(x2x1), cimag(x2x1)
            );
            hermitian_violated=1;
        }
    }
    
    
    // clear
    free(fft1);
    fftw_destroy_plan(plan_f_1d);
    fftw_destroy_plan(plan_b_1d);
    
    return hermitian_violated;
}


/**
 * Function computes angular momentum along z-direction
 * @param jx current, x-component (INPUT)
 * @param jy current, y-component (INPUT)
 * @param Lz angular momentum (OUTPUT)
 * */
int compute_angular_momentum_Lz(double *jx, double *jy, double *Lz)
{
    double _x, _y;
    int ix, iy, iz, ixyz=0;
    
    Lz[0] = 0.0; // reset
    
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        _x = DX*(double)(ix-NX/2);
        _y = DY*(double)(iy-NY/2);
        
        Lz[0] += (_x*jy[ixyz] - _y*jx[ixyz])*DXYZ;
        
        ixyz++; // go to next point
    }
    
    
    return 0;
}

/**
 * Function computes angular momentum along z-direction
 * @param it iteration number
 * @param spin - spin idicator (INPUT)
 * @param jx current, x-component (INPUT)
 * @param jy current, y-component (INPUT)
 * @param jz current, z-component (INPUT)
 * @param vext_dot_j integrated v_ext(r)*j(r) (OUTPUT)
 * */
int compute_vext_dot_j(int it, int spin, double *jx, double *jy, double *jz, double *vext_dot_j)
{
    int ix, iy, iz, ixyz=0;
    
    vext_dot_j[0] = 0.0; // reset
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
            
        vext_dot_j[0] += (
                             velocity_ext(ix, iy, iz, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)*jx[ixyz] 
                           + velocity_ext(ix, iy, iz, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data)*jy[ixyz]
                           + velocity_ext(ix, iy, iz, it, spin, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data)*jz[ixyz]
                         )*DX*DY*DZ;
        
        ixyz++; // gp to next point
    }
    
    
    return 0;
}
