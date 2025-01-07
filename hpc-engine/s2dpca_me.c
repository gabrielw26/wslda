// Author: Gabriel Wlazlowski

// Code for computing matrix elements of hamiltonian
#include <stdlib.h>
#include <stddef.h>

// globa variables
extern double *dc_params; /* Declaration of the variable */
extern size_t dc_extra_data_size;
extern void *dc_extra_data;
extern double dc_mu_a;
extern double dc_mu_b;
extern double dc_ec;

extern int wsldapid; // process id - global variable

#include "pca_settings.h"
#include "pca_macro.h"
#include "s2dpca_fft.h"
// #include "s2dpca_uext.h"
#include "s3dpca_grid.h"

#include "pca_utils.h"
#include "s2dpca_me.h"
#define BLOCKSIZE (NX*NY)

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

#define ixy2ixiy(ixy,_ix,_iy) \
    _ix=ixy/NY;               \
    _iy=ixy-_ix * NY;                   

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
 * @param bgrid stores information about the bc matrix distribution
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials recomputed potentials (INPUT)
 * @param mdfft metadata for ffts plans execution
 * @param h hamitonian matrix of size [2*NX*NY,2*NX*NY] (OUTPUT)
 * @param kz value of kz vector (INPUT)
 * @param me_d_dx matrix elements of (-i*d/dx) operator, matrix of size [NX x NX] (INPUT)
 * @param me_d_dy matrix elements of (-i*d/dy) operator, matrix of size [NY x NY] (INPUT)
 * */
int compute_matrix_elements_2d(metadata_s3dpca_grid *bgrid, int it, wslda_density h_densities, wslda_potential h_potentials, metadata_s2dpca_fft *mdfft, double complex *h, double kz, double complex * me_d_dx, double complex * me_d_dy)
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
    
    double Fx1, Fx2, Fy1, Fy2, t7, vext;          
    
    // iterate over all matrix elemnts
    
    int ix1, iy1, ixyz1; // row iterator
    int ix2, iy2, ixyz2; // column iterator
    int ix, iy, ixyz; // global iterator
    int ci, ri; // column and row iterator
    int li, lj, ij; // local indices
    int ZERO = 0;
    
    // compute gradient of effective mass
    double *laplace_alpha_a;
    double *laplace_alpha_b;
    cppmallocl(laplace_alpha_a,NX*NY,double);
    cppmallocl(laplace_alpha_b,NX*NY,double);
    int ierr;
    if(bgrid->nip*bgrid->niq<NX*NY/2) return -199; // check if enough memory
    double *wrk_dble = (double *)(h); // only for temporary calculations 
    
    ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) 
    {
        wrk_dble[ixyz] = h_potentials.alpha_a[ixyz]; 
        ixyz++;
    }
    ierr = compute_laplace_real_f(wrk_dble, laplace_alpha_a, mdfft);
    if(ierr!=0) return (100+ierr);
    
    ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) 
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
    // K(ixiy,jxjy) =  k(ix,jx)*delta(iy,jy)
    //                +delta(ix,jx)*k(iy,jy)
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
            
            if(ci<NX*NY && ri<NX*NY) // part: |   ha    |
            {
                ixyz1 = ri;
                ixyz2 = ci;
                
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                vext=v_ext(ix1,iy1,0,it,SPINA,dc_params,dc_extra_data_size,dc_extra_data); // mean field does not contain external potential 
                
                alph_1 = h_potentials.alpha_a[ixyz1]; 
                alph_2 = h_potentials.alpha_a[ixyz2]; 
                
                // diagonal part: V_a-mu_a + 0.5*alpha_a*kz^2
                if(ixyz1==ixyz2) h[ij] += V_a[ixyz1] + vext - dc_mu_a + 0.5*alph_1*kz*kz + 0.25*laplace_alpha_a[ixyz1];
                
                // Kinetic term K_a
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX);
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY);
                
                // vector field contribution
                Fx1 = h_potentials.A_a_x[ixyz1];
                Fy1 = h_potentials.A_a_y[ixyz1];
                
                Fx2 = h_potentials.A_a_x[ixyz2];
                Fy2 = h_potentials.A_a_y[ixyz2];
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] += me_d_dx[ix1 + ix2*NX] * (Fx1 + Fx2) * ( 0.5); // ( 1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] += me_d_dy[iy1 + iy2*NY] * (Fy1 + Fy2) * ( 0.5); // ( 1.0) because me_d_dx keeps matrix elements of (-i d/dx)


                // contribution from v_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = velocity_ext(ix1, iy1, 0, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = velocity_ext(ix2, iy2, 0, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2) h[ij] += me_d_dx[ix1 + ix2*NX] * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = velocity_ext(ix1, iy1, 0, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = velocity_ext(ix2, iy2, 0, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2) h[ij] += me_d_dy[iy1 + iy2*NY] * (Fy1 + Fy2) * ( -0.5); // note because me_d_dx keeps matrix elements of (-i d/dx)                
                
            }
            else if(ci>=NX*NY && ri<NX*NY) // part: |   Delta  |
            {
                ixyz1 = ri;
                ixyz2 = ci-NX*NY;             
                
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                // diagonal part: Delta(r)
                if(ixyz1==ixyz2) h[ij] += delta[ixyz1] + delta_ext(ix1, iy1, 0, it, delta[ixyz1], dc_params,dc_extra_data_size,dc_extra_data);
            }
            else if(ci<NX*NY && ri>=NX*NY) // part: | Delta^* |
            {
                ixyz1 = ri-NX*NY;
                ixyz2 = ci;          
                
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                // diagonal part: (Delta(r))^*
                if(ixyz1==ixyz2) h[ij] += conj(delta[ixyz1] + delta_ext(ix1, iy1, 0, it, delta[ixyz1], dc_params,dc_extra_data_size,dc_extra_data));
            }
            else // part: |   -hb^*  |
            {
                ixyz1 = ri-NX*NY;
                ixyz2 = ci-NX*NY;
                
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                vext=v_ext(ix1,iy1,0,it,SPINB,dc_params,dc_extra_data_size,dc_extra_data); // mean field does not contain external potential
                
                alph_1 = h_potentials.alpha_b[ixyz1]; 
                alph_2 = h_potentials.alpha_b[ixyz2];
                
                // diagonal part: V_b-mu_b + 0.5*alpha_b*kz^2
                if(ixyz1==ixyz2) h[ij] -= V_b[ixyz1] + vext - dc_mu_b + 0.5*alph_1*kz*kz + 0.25*laplace_alpha_b[ixyz1]; // NOTE -= operator has minus!
                
                // Kinetic term -K_b                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX); // NOTE -= operator has minus!
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY); // NOTE -= operator has minus!
                
                // vector field contribution
                Fx1 = h_potentials.A_b_x[ixyz1];
                Fy1 = h_potentials.A_b_y[ixyz1];
                
                Fx2 = h_potentials.A_b_x[ixyz2];
                Fy2 = h_potentials.A_b_y[ixyz2];
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] -= conj(me_d_dx[ix1 + ix2*NX]) * (Fx1 + Fx2) * ( 0.5); // (+1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] -= conj(me_d_dy[iy1 + iy2*NY]) * (Fy1 + Fy2) * ( 0.5); // (+1.0) because me_d_dx keeps matrix elements of (-i d/dx)

                
                // contribution from v_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = velocity_ext(ix1, iy1, 0, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = velocity_ext(ix2, iy2, 0, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2) h[ij] -= conj(me_d_dx[ix1 + ix2*NX]) * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = velocity_ext(ix1, iy1, 0, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = velocity_ext(ix2, iy2, 0, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2) h[ij] -= conj(me_d_dy[iy1 + iy2*NY]) * (Fy1 + Fy2) * ( -0.5); // note because me_d_dx keeps matrix elements of (-i d/dx)                
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
 * @param dx lattice constant
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
        fft1[ci] = 1.0/(nx); // normalization factor already included
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
        
        if(fabs(creal(diff))>1.0e-9 || fabs(cimag(diff))>1.0e-9)
        {
            wprintf("# ERROR[compute_matrix_elements_of_momentum_operator]: hermitian problem: %6d %6d (%f,%f) <=> (%f,%f)\n",
                ri, ci, creal(x1x2), cimag(x1x2), creal(x2x1), cimag(x2x1)
            );
            hermitian_violated=1;
        }
    }

    // // TEST check with the formula
    // // Both the formula and numerical approch should give EXACTLY the same result
    // for(ri=0; ri<nx; ri++) for(ci=0; ci<nx; ci++)
    // {
    //     double complex f=0.0;
    //     if(ri!=ci) f= M_PI/(nx*dx) * pow(-1.,ri-ci)*cos(M_PI*(ri-ci)/nx)/sin(M_PI*(ri-ci)/nx);
    //     f*=-1.0*I;
    //
    //     if(cabs(me[ri + nx*ci]-f)>1.0e-6)
    //     {
    //         wprintf("ERROR: %d %d (%f,%f) == (%f,%f)\n", ri, ci, creal(me[ri + nx*ci]),cimag(me[ri + nx*ci]), creal(f), cimag(f));
    //         hermitian_violated=2;
    //     }
    // }
    
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
    int ix, iy, ixyz=0;
    
    Lz[0] = 0.0; // reset
    
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        _x = DX*(double)(ix-NX/2);
        _y = DY*(double)(iy-NY/2);
        
        Lz[0] += (_x*jy[ixyz] - _y*jx[ixyz])*LZ*DX*DY;
        
        ixyz++; // gp to next point
    }
    
    
    return 0;
}

/**
 * Function computes angular momentum along z-direction
 * @param it iteration number
 * @param spin - spin idicator (INPUT)
 * @param jx current, x-component (INPUT)
 * @param jy current, y-component (INPUT)
 * @param vext_dot_j integrated v_ext(r)*j(r) (OUTPUT)
 * */
int compute_vext_dot_j(int it, int spin, double *jx, double *jy, double *vext_dot_j)
{
    int ix, iy, ixyz=0;
    
    vext_dot_j[0] = 0.0; // reset
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
            
        vext_dot_j[0] += (
                             velocity_ext(ix, iy, 0, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)*jx[ixyz] 
                           + velocity_ext(ix, iy, 0, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data)*jy[ixyz]
                         )*LZ*DX*DY;
        
        ixyz++; // gp to next point
    }
    
    
    return 0;
}



double fbeta(double E, double beta);
int test_density(metadata_s3dpca_grid *bgrid, double *En, int ne, double beta, double complex *U, double *rho_a, double *rho_b)
{
    int lj, li, ij;
    int ixyz1, ixyz2;
    double fbEn, fbmEn;
    int ri, ci;
    int ZERO = 0;
    for(lj=0; lj<bgrid->niq; lj++) // column-major iteration: over local index (column)
    {
        for(li=0; li<bgrid->nip; li++) // over local index (row)
        {
            ij=li + bgrid->nip * lj; // local index of elemnt
            ixyz1 = li+1; ixyz2=lj+1; // conversion to fortran standard
            // find indices in global matrix: row and colummn
            ri = indxl2g_( &ixyz1, &bgrid->mb, &bgrid->ip, &ZERO, &bgrid->p )-1; // back to C standard
            ci = indxl2g_( &ixyz2, &bgrid->nb, &bgrid->iq, &ZERO, &bgrid->q )-1; // back to C standard  
            if(ci>=ne) continue;
            
            // weight
            #define DENS_FACTOR_M 1.
            #define Complex(a,b) (a + I*b)
            #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
            
            fbEn=fbeta(En[ci], beta)*DENS_FACTOR_M;
            fbmEn = DENS_FACTOR_M - fbEn;
            
            if(ri<NX*NY) *rho_a+=cnorm(U[ij])*fbEn;
            else         *rho_b+=cnorm(U[ij])*fbmEn;
        }
    }
    return 0;
}

