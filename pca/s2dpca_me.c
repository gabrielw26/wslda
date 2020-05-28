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

// rotating frame
extern double dc_Omega_a;
extern double dc_Omega_b;

#include "pca_settings.h"
#include "pca_macro.h"
#include "s2dpca_fft.h"
#include "s2dpca_uext.h"
#include "s3dpca_grid.h"

// EDF functions
double polarization_h(double n_a, double n_b);
double der_polarization__der_na_h(double n_a, double n_b);
double der_polarization__der_nb_h(double n_a, double n_b);
double alpha_h(double p);
double alpha_a_h(double p);
double alpha_b_h(double p);
double alpha_plus_h(double p);
double alpha_minus_h(double p);
double der_alpha_plus__der_na_h(double n_a, double n_b);
double der_alpha_plus__der_nb_h(double n_a, double n_b);
double der_alpha_minus__der_na_h(double n_a, double n_b);
double der_alpha_minus__der_nb_h(double n_a, double n_b);
double funG_h(double p);
double funD_h(double n_a, double n_b);
double der_funD__der_na_h(double n_a, double n_b);
double der_funD__der_nb_h(double n_a, double n_b);
double tildeC_h(double n_a, double n_b);
double der_tildeC__der_na_h(double n_a, double n_b);
double der_tildeC__der_nb_h(double n_a, double n_b);

double polarization(double n_a, double n_b);
double der_polarization__der_na(double n_a, double n_b);
double der_polarization__der_nb(double n_a, double n_b);
double alpha(double p);
double alpha_a(double p);
double alpha_b(double p);
double alpha_plus(double p);
double alpha_minus(double p);
double der_alpha_plus__der_na(double n_a, double n_b);
double der_alpha_plus__der_nb(double n_a, double n_b);
double der_alpha_minus__der_na(double n_a, double n_b);
double der_alpha_minus__der_nb(double n_a, double n_b);
double funG(double p);
double funD(double n_a, double n_b);
double der_funD__der_na(double n_a, double n_b);
double der_funD__der_nb(double n_a, double n_b);
double tildeC(double n_a, double n_b);
double der_tildeC__der_na(double n_a, double n_b);
double der_tildeC__der_nb(double n_a, double n_b);
double p_regularization(double n);
double der_p_regularization(double n);

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
 * Function recomputes potentials
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials from PREVIOUS iteration, (INPUT) 
 *                     they are used as initial point for computation of new potentials 
 * @param h_potentials_new recomputed potentials (OUTPUT)
 * */
int recompute_potentials_aslda(int it, double *h_densities, double *h_potentials, double *h_potentials_new)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    // pontentials - decode
    double *V_a_new = (double *)(h_potentials_new +  0*NX*NY);
    double *V_b_new = (double *)(h_potentials_new +  1*NX*NY);
    double complex *delta_new = (double complex *)(h_potentials_new +  2*NX*NY);    
    
    // Code is equivivalent to the code implemented in pca_kernels.cu
    
    int ix, iy, ixyz;
    int i;
    int isconverged;
    
    // registers
    double na, nb;
    double t1, t2, t3, t4, t5, t6, t7; // working buffers
    
    double alph_plus;
#ifdef CURRENT_CORRECTIONS
    double alph_minus;
#endif
    double dalphm_dna, dalphm_dnb, dalphp_dna, dalphp_dnb;
    double Va, Vb, Vanew, Vbnew, Va_const, Vb_const;
    double complex p0, kc, wz_0, Zone, lnu, ldelta;
    
    ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) // for all points
    {
        Va_const=u_ext(ix,iy,it,SPINA,dc_params,dc_extra_data_size,dc_extra_data);
        Vb_const=u_ext(ix,iy,it,SPINB,dc_params,dc_extra_data_size,dc_extra_data);     
        
        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        t1 = polarization(na, nb);
        alph_plus = alpha_plus(t1);
#ifdef CURRENT_CORRECTIONS
        alph_minus = alpha_minus(t1);
#endif
        
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
        dalphm_dna=0.0;
        dalphm_dnb=0.0;
        dalphp_dna=0.0;
        dalphp_dnb=0.0;
#else
        // term dalphm_dna*tau_m/2.0
        t1=tau_a[ixyz]; // tau_a
        t2=tau_b[ixyz]; // tau_b       
        t3=t1 - t2; // tau_m
        dalphm_dna=der_alpha_minus__der_na(na, nb);
        dalphm_dnb=der_alpha_minus__der_nb(na, nb);
        Va_const+=dalphm_dna*t3/2.0; // dalphm_dna*tau_m/2.0
        Vb_const+=dalphm_dnb*t3/2.0; // dalphm_dnb*tau_m/2.0
        
        // term dalphp_dna*tau_p/2.0
        t3=t1 + t2; // tau_p
        dalphp_dna=der_alpha_plus__der_na(na, nb);
        dalphp_dnb=der_alpha_plus__der_nb(na, nb);
        Va_const+=dalphp_dna*t3/2.0; // dalphp_dna*tau_p/2.0
        Vb_const+=dalphp_dnb*t3/2.0; // dalphp_dnb*tau_p/2.0
        // no other terms with tau, now I can resue t1 and t2
#endif

        // term dD_dna and dD_dnb
        Va_const += der_funD__der_na(na, nb);
        Vb_const += der_funD__der_nb(na, nb);
        
#ifdef CURRENT_CORRECTIONS
        // current terms
        t1=j_a_x[ixyz];
        t2=j_a_y[ixyz];
        t3=j_a_z[ixyz];
        t4=j_b_x[ixyz];
        t5=j_b_y[ixyz];
        t6=j_b_z[ixyz];
        // terms with ja^2
        t7 = p_regularization(na);
        if(t7!=0.0)
        {
            t7 = t7*(t1*t1 + t2*t2 + t3*t3)/(2.0*na); // fr(na)*ja^2/2na
            Va_const+=( (alph_plus+alph_minus-1.0)/na - (dalphp_dna+dalphm_dna) ) *t7; // fr(na)*(alpha_a-1)*ja^2/2na^2 - fr(na)*dalpha_dna*ja^2/2na
            Va_const-=der_p_regularization(na)*(alph_plus+alph_minus-1.0)*(t1*t1 + t2*t2 + t3*t3)/(2.0*na); // derivative of regularization function
            Vb_const-=(dalphp_dnb+dalphm_dnb) *t7; // -dalpha_dnb*fr(na)*ja^2/2na            
        }

        // terms with jb^2
        t7 = p_regularization(nb);
        if(t7!=0.0)
        {
            t7 = t7*(t4*t4 + t5*t5 + t6*t6)/(2.0*nb); // fr(nb)*jb^2/2nb
            Va_const-=(dalphp_dna-dalphm_dna) *t7; // -dalphb_dna*fr(nb)*jb^2/2nb
            Vb_const+=( (alph_plus-alph_minus-1.0)/nb - (dalphp_dnb-dalphm_dnb-1.0) ) *t7; // fr(nb)*(alpha_b-1)*jb^2/2nb^2 - fr(nb)*dalphb_dnb*jb^2/2nb     
            Vb_const-=der_p_regularization(nb)*(alph_plus-alph_minus-1.0)*(t4*t4 + t5*t5 + t6*t6)/(2.0*nb); // derivative of regularization function
        }
        
// //         // terms with j+^2
// //         t7 = p_regularization(na+nb)*((t1+t4)*(t1+t4) + (t2+t5)*(t2+t5) + (t3+t6)*(t3+t6))/(2.0*(na+nb)*(na+nb)); // j+^2/2n+^2
// //         Va_const-=t7;
// //         Vb_const-=t7;
// //         t7 = der_p_regularization(na+nb)*((t1+t4)*(t1+t4) + (t2+t5)*(t2+t5) + (t3+t6)*(t3+t6))/(2.0*(na+nb)); // derivative of regularization function
// //         Va_const+=t7;
// //         Vb_const+=t7;
        // NOTE: divergence terms include when applied hamiltonian - here not needed
#endif
        
        // prepare other variables for self-consistent process
        t1=dalphp_dna/alph_plus; 
        t2=dalphp_dnb/alph_plus;
        t3=der_tildeC__der_na(na, nb) / alph_plus; // dtildeC_dna / alph_plus
        t4=der_tildeC__der_nb(na, nb) / alph_plus; // dtildeC_dnb / alph_plus
        t5 = tildeC(na, nb); // tC
        Va = V_a[ixyz]; // initial values
        Vb = V_b[ixyz]; // initial values
        lnu = nu[ixyz];
        Zone = Complex(1.0, 0.0);
        
        // computation of Va and Vb and delta
        for(i=0; i<UD_SCITERS; i++) // self-consistent loop
        {
            // pairing
            t7=(dc_mu_a-Va+dc_mu_b-Vb)/2.0;
            p0 = csqrt( Complex(2.0*t7/ alph_plus, 0.0) );
            if(cimag(p0)<0.) p0 *= -1. ;
            kc = csqrt( Complex(2.0*(dc_ec+t7)/ alph_plus, 0.0) );
            if(cimag(kc)<0.) kc *= -1. ;
            
            wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( cimag(wz_0) < 0. ) wz_0 += Complex(0.0, 2. * M_PI) ;    
            wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
            wz_0 = Zone*alph_plus / (Zone*t5 - wz_0);
            // g_eff = wz_0.real(); 
            ldelta = lnu*(-1.0*creal(wz_0));
            
            // potential
            t6=creal(conj(ldelta)*lnu); // delta^+ * nu 
            t7=cnorm(ldelta);
            Vanew = Va_const - t1*t6 - t3*t7;
            Vbnew = Vb_const - t2*t6 - t4*t7;
            
            // check convergence
            isconverged=1;
            if(fabs(Vanew-Va)>UD_EPSILON) isconverged=0; // Va not converged
            if(fabs(Vbnew-Vb)>UD_EPSILON) isconverged=0; // Vb not converged
            if(isconverged) break; 
             
            // mixing of potentials
            Va = UD_MIX_COEFF*Vanew+(1.0-UD_MIX_COEFF)*Va;
            Vb = UD_MIX_COEFF*Vbnew+(1.0-UD_MIX_COEFF)*Vb;
        }
        
//         if(i==UD_SCITERS) return -1; // error
        
        // save results to global memory
        V_a_new[ixyz]=Va;
        V_b_new[ixyz]=Vb;
#ifdef ENABLE_DELTA_EXT
        ldelta += delta_ext(ix, iy, it, ldelta, dc_params, dc_extra_data_size, dc_extra_data); // add external field
#endif
        delta_new[ixyz]=ldelta;       
        
        ixyz++; // go to next point
        
    } // for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    
    
    return 0;
}


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
 * Function computes matrix elemnts of BdG hamiltonian
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
int compute_matrix_elements_aslda(metadata_s3dpca_grid *bgrid, int it, double *h_densities, double *h_potentials, metadata_s2dpca_fft *mdfft, double complex *h, double kz, double complex * me_d_dx, double complex * me_d_dy)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    double p, alph_1, alph_2;
    
#if defined (CURRENT_CORRECTIONS) || defined (ENABLE_VELOCITY_EXT)
      double Fx1, Fx2, Fy1, Fy2, t7;          
#endif
    
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
        wrk_dble[ixyz] = alpha_a(polarization(rho_a[ixyz], rho_b[ixyz]));
        ixyz++;
    }
    ierr = compute_laplace_real_f(wrk_dble, laplace_alpha_a, mdfft);
    if(ierr!=0) return (100+ierr);
    
    ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) 
    {
        wrk_dble[ixyz] = alpha_b(polarization(rho_a[ixyz], rho_b[ixyz]));
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
                
                alph_1 = alpha_a(polarization(rho_a[ixyz1], rho_b[ixyz1]));
                alph_2 = alpha_a(polarization(rho_a[ixyz2], rho_b[ixyz2]));
                
                // diagonal part: V_a-mu_a + 0.5*alpha_a*kz^2
                if(ixyz1==ixyz2) h[ij] += V_a[ixyz1] - dc_mu_a + 0.5*alph_1*kz*kz + 0.25*laplace_alpha_a[ixyz1];
                
                // Kinetic termn K_a
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX);
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY);
                
#ifdef CURRENT_CORRECTIONS
                t7 = p_regularization(rho_a[ixyz1]);
                if(t7!=0.0)
                {
                    Fx1 = (1.0 - alph_1) * t7 *j_a_x[ixyz1] / (2.0 * rho_a[ixyz1]);
                    Fy1 = (1.0 - alph_1) * t7 *j_a_y[ixyz1] / (2.0 * rho_a[ixyz1]);
                }
                else
                {
                    Fx1 = 0.0;
                    Fy1 = 0.0;
                }

                t7 = p_regularization(rho_a[ixyz2]);
                if(t7!=0.0)
                {
                    Fx2 = (1.0 - alph_2) * t7 *j_a_x[ixyz2] / (2.0 * rho_a[ixyz2]);
                    Fy2 = (1.0 - alph_2) * t7 *j_a_y[ixyz2] / (2.0 * rho_a[ixyz2]);
                }
                else
                {
                    Fx2 = 0.0;
                    Fy2 = 0.0;
                }
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] += me_d_dx[ix1 + ix2*NX] * (Fx1 + Fx2) * ( 1.0); // ( 1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] += me_d_dy[iy1 + iy2*NY] * (Fy1 + Fy2) * ( 1.0); // ( 1.0) because me_d_dx keeps matrix elements of (-i d/dx)
#endif

                // rotating frame
                if(dc_Omega_a!=0.0)
                {
                    if(ix1==ix2) h[ij] -= me_d_dy[iy1 + iy2*NY]*dc_Omega_a*DX*(double)(ix1-NX/2)         ; 
                    if(iy1==iy2) h[ij] -= me_d_dx[ix1 + ix2*NX]*dc_Omega_a*DY*(double)(iy1-NY/2)* (-1.0) ; // (-1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                }
                
#ifdef ENABLE_VELOCITY_EXT
                // contribution from v_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = vector_vext(ix1, iy1, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = vector_vext(ix2, iy2, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2) h[ij] += me_d_dx[ix1 + ix2*NX] * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = vector_vext(ix1, iy1, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = vector_vext(ix2, iy2, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2) h[ij] += me_d_dy[iy1 + iy2*NY] * (Fy1 + Fy2) * ( -0.5); // note because me_d_dx keeps matrix elements of (-i d/dx)                
#endif
                
            }
            else if(ci>=NX*NY && ri<NX*NY) // part: |   Delta  |
            {
                ixyz1 = ri;
                ixyz2 = ci-NX*NY;             
                
                // diagonal part: Delta(r)
                if(ixyz1==ixyz2) h[ij] += delta[ixyz1];
            }
            else if(ci<NX*NY && ri>=NX*NY) // part: | Delta^* |
            {
                ixyz1 = ri-NX*NY;
                ixyz2 = ci;          
                
                // diagonal part: (Delta(r))^*
                if(ixyz1==ixyz2) h[ij] += conj(delta[ixyz1]);
            }
            else // part: |   -hb^*  |
            {
                ixyz1 = ri-NX*NY;
                ixyz2 = ci-NX*NY;
                
                alph_1 = alpha_b(polarization(rho_a[ixyz1], rho_b[ixyz1]));
                alph_2 = alpha_b(polarization(rho_a[ixyz2], rho_b[ixyz2]));
                
                // diagonal part: V_b-mu_b + 0.5*alpha_b*kz^2
                if(ixyz1==ixyz2) h[ij] -= V_b[ixyz1] - dc_mu_b + 0.5*alph_1*kz*kz + 0.25*laplace_alpha_b[ixyz1]; // NOTE -= operator has minus!
                
                // Kinetic termn -K_b
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX); // NOTE -= operator has minus!
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY); // NOTE -= operator has minus!
                
#ifdef CURRENT_CORRECTIONS
                t7 = p_regularization(rho_b[ixyz1]);
                if(t7!=0.0)
                {
                    Fx1 = (1.0 - alph_1) * t7 *j_b_x[ixyz1] / (2.0 * rho_b[ixyz1]);
                    Fy1 = (1.0 - alph_1) * t7 *j_b_y[ixyz1] / (2.0 * rho_b[ixyz1]);
                }
                else
                {
                    Fx1 = 0.0;
                    Fy1 = 0.0;
                }

                t7 = p_regularization(rho_b[ixyz2]);
                if(t7!=0.0)
                {
                    Fx2 = (1.0 - alph_2) * t7 *j_b_x[ixyz2] / (2.0 * rho_b[ixyz2]);
                    Fy2 = (1.0 - alph_2) * t7 *j_b_y[ixyz2] / (2.0 * rho_b[ixyz2]);
                }
                else
                {
                    Fx2 = 0.0;
                    Fy2 = 0.0;
                } 
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] -= conj(me_d_dx[ix1 + ix2*NX]) * (Fx1 + Fx2) * ( 1.0); // (+1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] -= conj(me_d_dy[iy1 + iy2*NY]) * (Fy1 + Fy2) * ( 1.0); // (+1.0) because me_d_dx keeps matrix elements of (-i d/dx)
#endif

                // rotating frame
                if(dc_Omega_b!=0.0)
                {
                    if(ix1==ix2) h[ij] += conj(me_d_dy[iy1 + iy2*NY])*dc_Omega_b*DX*(double)(ix1-NX/2)         ;
                    if(iy1==iy2) h[ij] += conj(me_d_dx[ix1 + ix2*NX])*dc_Omega_b*DY*(double)(iy1-NY/2)* (-1.0) ;  // (-1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                }
                
#ifdef ENABLE_VELOCITY_EXT
                // contribution from v_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = vector_vext(ix1, iy1, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = vector_vext(ix2, iy2, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2) h[ij] -= conj(me_d_dx[ix1 + ix2*NX]) * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = vector_vext(ix1, iy1, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = vector_vext(ix2, iy2, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2) h[ij] -= conj(me_d_dy[iy1 + iy2*NY]) * (Fy1 + Fy2) * ( -0.5); // note because me_d_dx keeps matrix elements of (-i d/dx)                
#endif
            }
            
        } // for(li=0; li<bgrid->nip; li++)
//         printf("ci=%d\n", ci); fflush(stdout);
    } // for(lj=0; lj<bgrid->niq; lj++)
    
    // clear memory
    free(laplace_alpha_a);
    free(laplace_alpha_b);
    
    return 0;
}

/**
 * Function that computes energy of the system
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials corresponding to the densities (INPUT) 
 * @param energy array with contributions to the energy (OUTPUT)
 * @param npart array with contributions to the particle number (OUTPUT)
 * */
int compute_energy_aslda(int it, double *h_densities, double *h_potentials, double *energy, double *npart)
{
    // Set pointers for to simplify notation
    // densities 
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    // pontentials
    // pontentials - decode
//     double *V_a = (double *)(h_potentials +  0*NX*NY);
//     double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    // buffers for energies
    double *E_kin = (double *)(energy +  0);
    double *E_pot = (double *)(energy +  1);
    double *E_pair= (double *)(energy +  2);
    double *E_CM  = (double *)(energy +  3);
    double *E_ext = (double *)(energy +  4);
    
    // reset buffers
    int i;
    for(i=0; i<5; i++) energy[i]=0.0;
    npart[SPINA]=0.0; npart[SPINB]=0.0;
    
    double na, nb;
    double p;
    double taua, taub;
#ifdef CURRENT_CORRECTIONS
    double tx1, ty1, tz1, tx2, ty2, tz2;
#endif
    
    int ix, iy;
    int ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {

        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        
        // particle number
        npart[SPINA]+=na;
        npart[SPINB]+=nb;
        
        // External potential energy
        E_ext[0]+=na*u_ext(ix,iy,it,SPINA,dc_params,dc_extra_data_size,dc_extra_data) + 
                  nb*u_ext(ix,iy,it,SPINB,dc_params,dc_extra_data_size,dc_extra_data);
        
        // kinetic energy
        p=polarization(na, nb);
        taua=tau_a[ixyz]; // tau_a
        taub=tau_b[ixyz]; // tau_b    
        
        // current corrections
#ifdef CURRENT_CORRECTIONS
        tx1=j_a_x[ixyz];
        ty1=j_a_y[ixyz];
        tz1=j_a_z[ixyz];
        tx2=j_b_x[ixyz];
        ty2=j_b_y[ixyz];
        tz2=j_b_z[ixyz];
        taua-=p_regularization(na)*(tx1*tx1+ty1*ty1+tz1*tz1)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(tx2*tx2+ty2*ty2+tz2*tz2)/nb; // -jb^2/nb: correction for tilde{tau}_b
#endif

        E_kin[0]+=0.5*(alpha_a(p)*taua + alpha_b(p)*taub);
        
        // potential energy
        E_pot[0]+=funD(na, nb);
        
        // pairing energy
        E_pair[0]+=creal(delta[ixyz]*conj(nu[ixyz]))*(-1.0);
        
        // center of mass motion energy
#ifdef CURRENT_CORRECTIONS
//         E_CM[ixyz]=p_regularization(na+nb)*((tx1+tx2)*(tx1+tx2) + (ty1+ty2)*(ty1+ty2) + (tz1+tz2)*(tz1+tz2))/(2.0*(na+nb));
        E_CM[0]+=   p_regularization(na)*(tx1*tx1 + ty1*ty1 + tz1*tz1)/(2.*na)  
                    + p_regularization(nb)*(tx2*tx2 + ty2*ty2 + tz2*tz2)/(2.*nb);  
#else
        E_CM[0]+=0.0;
#endif
        ixyz++;
    }
    
    // take into account uniformity in the "z" direction
    for(i=0; i<5; i++) energy[i]*=DX*DY*DZ*NZ;
    npart[SPINA]*=DX*DY*DZ*NZ;
    npart[SPINB]*=DX*DY*DZ*NZ;
    
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
            printf("# ERROR[compute_matrix_elements_of_momentum_operator]: hermitian problem: %6d %6d (%f,%f) <=> (%f,%f)\n",
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
                             vector_vext(ix, iy, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data)*jx[ixyz] 
                           + vector_vext(ix, iy, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data)*jy[ixyz]
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



// --------------------------------------------------------------------------------------------------
// -------------------------------------- BdG variants ----------------------------------------------
// --------------------------------------------------------------------------------------------------
extern double aBdG; // scattering length
/**
 * Function recomputes potentials
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials from PREVIOUS iteration, (INPUT) 
 *                     they are used as initial point for computation of new potentials 
 * @param h_potentials_new recomputed potentials (OUTPUT)
 * */
int recompute_potentials_bdg(int it, double *h_densities, double *h_potentials, double *h_potentials_new)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    // pontentials - decode
    double *V_a_new = (double *)(h_potentials_new +  0*NX*NY);
    double *V_b_new = (double *)(h_potentials_new +  1*NX*NY);
    double complex *delta_new = (double complex *)(h_potentials_new +  2*NX*NY);    
    
    // Code is equivivalent to the code implemented in pca_kernels.cu
    
    int ix, iy, ixyz;
    int i;
    int isconverged;
    
    // registers
    double na, nb;
    double t1, t2, t3, t4, t5, t6, t7; // working buffers
    
    double Va, Vb, Vanew, Vbnew, Va_const, Vb_const;
    double complex p0, kc, wz_0, Zone, lnu, ldelta;
    
    ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) // for all points
    {
        Va_const=u_ext(ix,iy,it,SPINA,dc_params,dc_extra_data_size,dc_extra_data);
        Vb_const=u_ext(ix,iy,it,SPINB,dc_params,dc_extra_data_size,dc_extra_data);
        
        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];

        t5 = 1.0/ (4.0*M_PI*aBdG);
        Va = Va_const; // initial values
        Vb = Vb_const; // initial values
        lnu = nu[ixyz];
        Zone = Complex(1.0, 0.0);
        

        // pairing
        t7=(dc_mu_a-Va+dc_mu_b-Vb)/2.0;
        p0 = csqrt( Complex(2.0*t7, 0.0) );
        if(cimag(p0)<0.) p0 *= -1. ;
        kc = csqrt( Complex(2.0*(dc_ec+t7), 0.0) );
        if(cimag(kc)<0.) kc *= -1. ;
        
        wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
        if ( cimag(wz_0) < 0. ) wz_0 += Complex(0.0, 2. * M_PI) ;    
        wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
        wz_0 = Zone / (Zone*t5 - wz_0);
        // g_eff = wz_0.real(); 
        ldelta = lnu*(-1.0*creal(wz_0));
            
        // save results to global memory
        V_a_new[ixyz]=Va;
        V_b_new[ixyz]=Vb;
#ifdef ENABLE_DELTA_EXT
        ldelta += delta_ext(ix, iy, it, ldelta, dc_params, dc_extra_data_size, dc_extra_data); // add external field
#endif
        delta_new[ixyz]=ldelta;       
        
        ixyz++; // go to next point
        
    } // for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    
    
    return 0;
}

/**
 * Function computes matrix elemnts of BdG hamiltonian
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
int compute_matrix_elements_bdg(metadata_s3dpca_grid *bgrid, int it, double *h_densities, double *h_potentials, metadata_s2dpca_fft *mdfft, double complex *h, double kz, double complex * me_d_dx, double complex * me_d_dy)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    double p, alph_1, alph_2;
    
    // iterate over all matrix elemnts
    
    int ix1, iy1, ixyz1; // row iterator
    int ix2, iy2, ixyz2; // column iterator
    int ix, iy, ixyz; // global iterator
    int ci, ri; // column and row iterator
    int li, lj, ij; // local indices
    int ZERO = 0;
    
#ifdef ENABLE_VELOCITY_EXT
    double Fx1, Fx2, Fy1, Fy2;
#endif
    
    
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
                
                alph_1 = 1.0; // const effective mass
                alph_2 = 1.0; // const effective mass
                
                // diagonal part: V_a-mu_a + 0.5*alpha_a*kz^2
                if(ixyz1==ixyz2) h[ij] += V_a[ixyz1] - dc_mu_a + 0.5*alph_1*kz*kz;
                
                // Kinetic termn K_a
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX);
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] += 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY);

                // rotating frame
                if(dc_Omega_a!=0.0)
                {
                    if(ix1==ix2) h[ij] -= me_d_dy[iy1 + iy2*NY]*dc_Omega_a*DX*(double)(ix1-NX/2)         ; 
                    if(iy1==iy2) h[ij] -= me_d_dx[ix1 + ix2*NX]*dc_Omega_a*DY*(double)(iy1-NY/2)* (-1.0) ; // (-1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                }
                
#ifdef ENABLE_VELOCITY_EXT
                // contribution from v_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = vector_vext(ix1, iy1, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = vector_vext(ix2, iy2, it, SPINA, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2) h[ij] += me_d_dx[ix1 + ix2*NX] * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = vector_vext(ix1, iy1, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = vector_vext(ix2, iy2, it, SPINA, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2) h[ij] += me_d_dy[iy1 + iy2*NY] * (Fy1 + Fy2) * ( -0.5); // note because me_d_dx keeps matrix elements of (-i d/dx)                
#endif
                
            }
            else if(ci>=NX*NY && ri<NX*NY) // part: |   Delta  |
            {
                ixyz1 = ri;
                ixyz2 = ci-NX*NY;             
                
                // diagonal part: Delta(r)
                if(ixyz1==ixyz2) h[ij] += delta[ixyz1];
            }
            else if(ci<NX*NY && ri>=NX*NY) // part: | Delta^* |
            {
                ixyz1 = ri-NX*NY;
                ixyz2 = ci;          
                
                // diagonal part: (Delta(r))^*
                if(ixyz1==ixyz2) h[ij] += conj(delta[ixyz1]);
            }
            else // part: |   -hb^*  |
            {
                ixyz1 = ri-NX*NY;
                ixyz2 = ci-NX*NY;
                
                alph_1 = 1.0; // const effective mass
                alph_2 = 1.0; // const effective mass
                
                // diagonal part: V_b-mu_b + 0.5*alpha_b*kz^2
                if(ixyz1==ixyz2) h[ij] -= V_b[ixyz1] - dc_mu_b + 0.5*alph_1*kz*kz; // NOTE -= operator has minus!
                
                // Kinetic termn -K_b
                // docode to cartesian coordinates
                ixy2ixiy(ixyz1,ix1,iy1);
                ixy2ixiy(ixyz2,ix2,iy2);
                
                // k(ix,jx)*delta(iy,jy)
                if(iy1==iy2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(ix1, ix2, NX, DX); // NOTE -= operator has minus!
                // delta(ix,jx)*k(iy,jy)
                if(ix1==ix2) h[ij] -= 0.5*(alph_1+alph_2)*k_1D(iy1, iy2, NY, DY); // NOTE -= operator has minus!
                
                // rotating frame
                if(dc_Omega_b!=0.0)
                {
                    if(ix1==ix2) h[ij] += conj(me_d_dy[iy1 + iy2*NY])*dc_Omega_b*DX*(double)(ix1-NX/2)         ;
                    if(iy1==iy2) h[ij] += conj(me_d_dx[ix1 + ix2*NX])*dc_Omega_b*DY*(double)(iy1-NY/2)* (-1.0) ;  // (-1.0) because me_d_dx keeps matrix elements of (-i d/dx)
                }
                
#ifdef ENABLE_VELOCITY_EXT
                // contribution from v_ext: 
                // -(1/2) {v_vext, p} = -(1/2)[v_ext*p + p*v_ext]
                
                Fx1 = vector_vext(ix1, iy1, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fx2 = vector_vext(ix2, iy2, it, SPINB, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(iy1==iy2) h[ij] -= conj(me_d_dx[ix1 + ix2*NX]) * (Fx1 + Fx2) * ( -0.5); // (note me_d_dx keeps matrix elements of (-i d/dx)
                
                Fy1 = vector_vext(ix1, iy1, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                Fy2 = vector_vext(ix2, iy2, it, SPINB, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
                if(ix1==ix2) h[ij] -= conj(me_d_dy[iy1 + iy2*NY]) * (Fy1 + Fy2) * ( -0.5); // note because me_d_dx keeps matrix elements of (-i d/dx)                
#endif
            }
            
        } // for(li=0; li<bgrid->nip; li++)
//         printf("ci=%d\n", ci); fflush(stdout);
    } // for(lj=0; lj<bgrid->niq; lj++)
        
    return 0;
}

/**
 * Function that computes energy of the system
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials corresponding to the densities (INPUT) 
 * @param energy array with contributions to the energy (OUTPUT)
 * @param npart array with contributions to the particle number (OUTPUT)
 * */
int compute_energy_bdg(int it, double *h_densities, double *h_potentials, double *energy, double *npart)
{
    // Set pointers for to simplify notation
    // densities 
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    // pontentials
    // pontentials - decode
//     double *V_a = (double *)(h_potentials +  0*NX*NY);
//     double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    // buffers for energies
    double *E_kin = (double *)(energy +  0);
    double *E_pot = (double *)(energy +  1);
    double *E_pair= (double *)(energy +  2);
    double *E_CM  = (double *)(energy +  3);
    double *E_ext = (double *)(energy +  4);
    
    // reset buffers
    int i;
    for(i=0; i<5; i++) energy[i]=0.0;
    npart[SPINA]=0.0; npart[SPINB]=0.0;
    
    double na, nb;
    double p;
    double taua, taub;
    
    int ix, iy;
    int ixyz=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {

        // densities, and correct them
        na=rho_a[ixyz];
        nb=rho_b[ixyz];
        
        // particle number
        npart[SPINA]+=na;
        npart[SPINB]+=nb;
        
        // External potential energy
        E_ext[0]+=na*u_ext(ix,iy,it,SPINA,dc_params,dc_extra_data_size,dc_extra_data) + 
                  nb*u_ext(ix,iy,it,SPINB,dc_params,dc_extra_data_size,dc_extra_data);
        
        // kinetic energy
        p=polarization(na, nb);
        taua=tau_a[ixyz]; // tau_a
        taub=tau_b[ixyz]; // tau_b    
        
        // current corrections

        E_kin[0]+=0.5*(taua + taub);
        
        // potential energy
        E_pot[0]+=0.0;
        
        // pairing energy
        E_pair[0]+=creal(delta[ixyz]*conj(nu[ixyz]))*(-1.0);
        
        // center of mass motion energy
        E_CM[0]+=0.0;

        ixyz++;
    }
    
    // take into account uniformity in the "z" direction
    for(i=0; i<5; i++) energy[i]*=DX*DY*DZ*NZ;
    npart[SPINA]*=DX*DY*DZ*NZ;
    npart[SPINB]*=DX*DY*DZ*NZ;
    
    return 0;
}

// ==========================================================================
// ============================== WRAPPER ===================================
// ==========================================================================
int recompute_potentials(int it, double *h_densities, double *h_potentials, double *h_potentials_new)
{
    if(fabs(aBdG)<1.0e-12) return recompute_potentials_aslda(it, h_densities, h_potentials, h_potentials_new);
    else                   return recompute_potentials_bdg  (it, h_densities, h_potentials, h_potentials_new);
}

int compute_matrix_elements(metadata_s3dpca_grid *bgrid, int it, double *h_densities, double *h_potentials, metadata_s2dpca_fft *mdfft, double complex *h, double kz, double complex * me_d_dx, double complex * me_d_dy)
{
    if(fabs(aBdG)<1.0e-12) return compute_matrix_elements_aslda(bgrid, it, h_densities, h_potentials, mdfft, h, kz, me_d_dx, me_d_dy);
    else                   return compute_matrix_elements_bdg  (bgrid, it, h_densities, h_potentials, mdfft, h, kz, me_d_dx, me_d_dy);
}

int compute_energy(int it, double *h_densities, double *h_potentials, double *energy, double *npart)
{
    if(fabs(aBdG)<1.0e-12) return compute_energy_aslda(it, h_densities, h_potentials, energy, npart);
    else                   return compute_energy_bdg  (it, h_densities, h_potentials, energy, npart);
}
