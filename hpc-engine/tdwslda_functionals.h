/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 11-11-2021
 * 
 * This file implements CUDA kernels that implement functionals. 
 * 
 * !!! READ THIS BEFORE YOU START TO IMPLEMENT NEW FUNCTIONAL !!!
 * 
 * To define a functional, you must provide a body of two functions:
 *  1. tdwslda_compute_potentials(int it, wslda_density h_densities, wslda_potential h_potentials, double cccoeff):
 *     where:
 *     - it is time index and time=dc_t0 + dc_dt*it;
 *     - h_densities is structure storing pointers to density (INPUT), 
 *          see: hpc-engine/wslda_potdens.h 
 *     - h_potentials is structure storing pointers to potentials (INPUT/OUTPUT). 
 *          On input: potentials from previews iteration
 *          On output: you need to provide values of updated potentials
 *          NOTE:
 *          -- V_a and V_b: must contain sum of external potential and mean field 
 *             see: hpc-engine/wslda_potdens.h 
 *     - cccoeff is external coefficient controlled from input file via ccstart, ccstop, and ccswitch. 
 *           Precisely: cccoeff=h_smooth_step(t0+it*dt, md.ccstart/eF,  md.ccstop/eF,  md.ccswitch/eF, 1.0);
 *           We keep it due to legacy issues.
 * 
 *  2. tdwslda_compute_energy(int it, wslda_density h_densities, wslda_potential h_potentials, // < INPUT
 *                            double *E_kin, double *E_pot, double *E_pair, double *E_curr)    // < OUTPUT
 *     This function computes energy contributions to the total energy. Arguments are analogus to tdwslda_compute_potentials(…) function. 
 * 
 *  Other Info:
 *     - u_ext(ix,iy,iz,it,SPINA): 
 *          it is equivalent to: 
 *        v_ext(int ix, int iy, int iz, int it, int spin, double *params, size_t extra_data_size, void *extra_data)
 *     - Use CODEDIM macro-variable to identify dimensionality of the code.
 *     - Use  decode_ixyz2ixiyiz(ixyz,ix,iy,iz, i) to decompose global index into lattice coordinate. 
 *          In case of 2D code iz will be set to iz=0, and for 1D code iz and iy will be set iz=iy=0. 
 *     - Avoid multiple reading of the same variable from global memmory, load it before to registers.
 * */  

// DO NOT REMOVE!
#include "tdwslda_functionals_framework_enable.h"

// ------------------------------------------------------------------------------------
// ---------------------------------- (A)SLDA -----------------------------------------
// ------------------------------------------------------------------------------------
#if FUNCTIONAL==SLDA || FUNCTIONAL==ASLDA
__global__ void tdwslda_compute_potentials(int it, wslda_density h_densities, wslda_potential h_potentials, double cccoeff)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;
    
    // registers
    double na, nb;
    double t1, t2, t3, t4, t5, t6, t7; // working buffers
    
    double alph_plus;
#ifdef CURRENT_CORRECTIONS
    double alph_minus;
#endif
    double dalphm_dna, dalphm_dnb, dalphp_dna, dalphp_dnb;
    double Va, Vb, Vanew, Vbnew, Va_const, Vb_const;
    Complex p0, kc, wz_0, Zone, lnu, ldelta;
    
    if(ixyz<NUMBER_ELEMENT)
    {
        decode_ixyz2ixiyiz(ixyz,ix,iy,iz, i); // decode cartesian coordinates
        
        // load data to registers from global memory and form constant part of potentials
        Va_const=u_ext(ix,iy,iz,it,SPINA);
        Vb_const=u_ext(ix,iy,iz,it,SPINB);
        
        // densities, and correct them
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        t1 = polarization(na, nb);
        alph_plus = alpha_plus(t1);
#ifdef CURRENT_CORRECTIONS
        alph_minus = alpha_minus(t1);
        
        // compute effecctive masses
        t2=alpha_a(t1);
        t3=alpha_b(t1);
        h_potentials.alpha_a[ixyz]=t2; // save
        h_potentials.alpha_b[ixyz]=t3; // save 
        
        // compute vector potentials
        // Spin a component
        t4=p_regularization(na) * cccoeff;
        if(t4==0.0)
        {
#if CODEDIM>=1
            h_potentials.A_a_x[ixyz]=0.0;
#endif
#if CODEDIM>=2
            h_potentials.A_a_y[ixyz]=0.0;
#endif
#if CODEDIM>=3
            h_potentials.A_a_z[ixyz]=0.0;
#endif
        }
        else
        {
#if CODEDIM>=1
            // x-coordinate
            t5=h_densities.j_a_x[ixyz]; 
            h_potentials.A_a_x[ixyz] = t4*(1.-t2)*t5/na;
#endif
#if CODEDIM>=2
            // y-coordinate
            t5=h_densities.j_a_y[ixyz]; 
            h_potentials.A_a_y[ixyz] = t4*(1.-t2)*t5/na;
#endif
#if CODEDIM>=3
            // z-coordinate
            t5=h_densities.j_a_z[ixyz]; 
            h_potentials.A_a_z[ixyz] = t4*(1.-t2)*t5/na;
#endif
        }
        
        // Spin b component
        t4=p_regularization(nb) * cccoeff;
        if(t4==0.0)
        {
#if CODEDIM>=1
            h_potentials.A_b_x[ixyz]=0.0;
#endif
#if CODEDIM>=2
            h_potentials.A_b_y[ixyz]=0.0;
#endif
#if CODEDIM>=3
            h_potentials.A_b_z[ixyz]=0.0;
#endif
        }
        else
        {
#if CODEDIM>=1
            // x-coordinate
            t5=h_densities.j_b_x[ixyz]; 
            h_potentials.A_b_x[ixyz] = t4*(1.-t3)*t5/nb;
#endif
#if CODEDIM>=2
            // y-coordinate
            t5=h_densities.j_b_y[ixyz]; 
            h_potentials.A_b_y[ixyz] = t4*(1.-t3)*t5/nb;
#endif
#if CODEDIM>=3
            // z-coordinate
            t5=h_densities.j_b_z[ixyz]; 
            h_potentials.A_b_z[ixyz] = t4*(1.-t3)*t5/nb;
#endif
        }
#endif
        
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
        dalphm_dna=0.0;
        dalphm_dnb=0.0;
        dalphp_dna=0.0;
        dalphp_dnb=0.0;
#else
        // term dalphm_dna*tau_m/2.0
        t1=h_densities.tau_a[ixyz]; // tau_a
        t2=h_densities.tau_b[ixyz]; // tau_b       
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
#if CODEDIM>=1
        t1=h_densities.j_a_x[ixyz];
        t4=h_densities.j_b_x[ixyz];
#endif
#if CODEDIM>=2
        t2=h_densities.j_a_y[ixyz];
        t5=h_densities.j_b_y[ixyz];
#else
        t2=0.0; // to avoid empty load
        t5=0.0; // to avoid empty load
#endif
#if CODEDIM>=3
        t3=h_densities.j_a_z[ixyz];
        t6=h_densities.j_b_z[ixyz];
#else
        t3=0.0; // to avoid empty load
        t6=0.0; // to avoid empty load
#endif
        // terms with ja^2
        t7 = p_regularization(na) * cccoeff;
        if(t7!=0.0)
        {
            t7 = t7*(t1*t1 + t2*t2 + t3*t3)/(2.0*na); // fr(na)*ja^2/2na
            Va_const+=( (alph_plus+alph_minus-1.0)/na - (dalphp_dna+dalphm_dna) ) *t7; // fr(na)*(alpha_a-1)*ja^2/2na^2 - fr(na)*dalpha_dna*ja^2/2na
            Va_const-=cccoeff*der_p_regularization(na)*(alph_plus+alph_minus-1.0)*(t1*t1 + t2*t2 + t3*t3)/(2.0*na); // derivative of regularization function
            Vb_const-=(dalphp_dnb+dalphm_dnb) *t7; // -dalpha_dnb*fr(na)*ja^2/2na            
        }

        // terms with jb^2
        t7 = p_regularization(nb) * cccoeff;
        if(t7!=0.0)
        {
            t7 = t7*(t4*t4 + t5*t5 + t6*t6)/(2.0*nb); // fr(nb)*jb^2/2nb
            Va_const-=(dalphp_dna-dalphm_dna) *t7; // -dalphb_dna*fr(nb)*jb^2/2nb
            Vb_const+=( (alph_plus-alph_minus-1.0)/nb - (dalphp_dnb-dalphm_dnb-1.0) ) *t7; // fr(nb)*(alpha_b-1)*jb^2/2nb^2 - fr(nb)*dalphb_dnb*jb^2/2nb     
            Vb_const-=cccoeff*der_p_regularization(nb)*(alph_plus-alph_minus-1.0)*(t4*t4 + t5*t5 + t6*t6)/(2.0*nb); // derivative of regularization function
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
        
        t1=dalphp_dna/alph_plus; 
        t2=dalphp_dnb/alph_plus;
        t3=der_tildeC__der_na(na, nb, 1.0) / alph_plus; // dtildeC_dna / alph_plus
        t4=der_tildeC__der_nb(na, nb, 1.0) / alph_plus; // dtildeC_dnb / alph_plus
        t5 = tildeC(na, nb, 1.0); // tC
        Va = h_potentials.V_a[ixyz]; // initial values
        Vb = h_potentials.V_b[ixyz]; // initial values
        lnu = h_densities.nu[ixyz];
        Zone = Complex(1.0, 0.0);
  
        // computation of Va and Vb and delta
        for(i=0; i<UD_SCITERS; i++) // self-consistent loop
        {
            // pairing
#ifdef USE_CUBIC_CUTOFF
            wz_0=Complex(REGULARIZATION_SCHEME_K_CONST/(4.0*M_PI*DX), 0.0);
#else
            t7=(dc_mu_a-Va+dc_mu_b-Vb)/2.0;
            p0 = thrust::sqrt( Complex(2.0*t7/ alph_plus, 0.0) );
            if(p0.imag()<0.) p0 *= -1. ;
            kc = thrust::sqrt( Complex(2.0*(dc_ec+t7)/ alph_plus, 0.0) );
            if(kc.imag()<0.) kc *= -1. ;
            
            wz_0 = thrust::log( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( wz_0.imag() < 0. ) wz_0 += Complex(0.0, 2. * M_PI) ;    
            wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
#endif
            wz_0 = Zone*alph_plus / (Zone*t5 - wz_0);
            // g_eff = wz_0.real(); 
            ldelta = lnu*(-1.0*wz_0.real());
            
            // potential
            t6=(thrust::conj(ldelta)*lnu).real(); // delta^+ * nu 
            t7=thrust::norm(ldelta);
            Vanew = Va_const - t1*t6 - t3*t7;
            Vbnew = Vb_const - t2*t6 - t4*t7;
             
            // mixing of potentials
            Va = UD_MIX_COEFF*Vanew+(1.0-UD_MIX_COEFF)*Va;
            Vb = UD_MIX_COEFF*Vbnew+(1.0-UD_MIX_COEFF)*Vb;
        }
        
        // save results to global memory
        h_potentials.V_a[ixyz]=Va;
        h_potentials.V_b[ixyz]=Vb;   
        h_potentials.delta[ixyz]=ldelta;
    }

}

__global__ void tdwslda_compute_energy(int it, wslda_density h_densities, wslda_potential h_potentials, // <-- INPUT
                                       double *E_kin, double *E_pot, double *E_pair, double *E_curr)    // <- OUTPUT
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix;
    
    double na, nb;
    double p;
    double taua, taub;
    double tx1, ty1, tz1, tx2, ty2, tz2;
    
    if(ixyz<NUMBER_ELEMENT)
    {

        // densities, and correct them
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        
        // kinetic energy
        p=polarization(na, nb);
        taua=h_densities.tau_a[ixyz]; // tau_a
        taub=h_densities.tau_b[ixyz]; // tau_b    
        
        // current corrections
#if CODEDIM>=1
        tx1=h_densities.j_a_x[ixyz];
        tx2=h_densities.j_b_x[ixyz];
#endif
#if CODEDIM>=2
        ty1=h_densities.j_a_y[ixyz];
        ty2=h_densities.j_b_y[ixyz];
#else        
        ty1=0.0;
        ty2=0.0;
#endif
#if CODEDIM>=3
        tz1=h_densities.j_a_z[ixyz];
        tz2=h_densities.j_b_z[ixyz];
#else
        tz1=0.0;
        tz2=0.0;
#endif
        taua-=p_regularization(na)*(tx1*tx1+ty1*ty1+tz1*tz1)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(tx2*tx2+ty2*ty2+tz2*tz2)/nb; // -jb^2/nb: correction for tilde{tau}_b

        // galilean invariant contribution
        E_kin[ixyz]=0.5*(alpha_a(p)*taua + alpha_b(p)*taub)*VOLUME_ELEMENT;
        
        // potential energy
        E_pot[ixyz]=funD(na, nb)*VOLUME_ELEMENT;
        
        // pairing energy
        E_pair[ixyz]=(h_potentials.delta[ixyz]*thrust::conj(h_densities.nu[ixyz])).real()*(-1.0)*VOLUME_ELEMENT;
        
        // flow energy
        E_curr[ixyz]=  p_regularization(na)*(tx1*tx1 + ty1*ty1 + tz1*tz1)/(2.*na)*VOLUME_ELEMENT  
                     + p_regularization(nb)*(tx2*tx2 + ty2*ty2 + tz2*tz2)/(2.*nb)*VOLUME_ELEMENT;  
    }
}
#endif

// ------------------------------------------------------------------------------------
// ------------------------------------ BDG -------------------------------------------
// ------------------------------------------------------------------------------------
#if FUNCTIONAL==BDG
__global__ void tdwslda_compute_potentials(int it, wslda_density h_densities, wslda_potential h_potentials, double cccoeff)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;
    
    // registers
    double t5, t7; // working buffers
    
    double Va, Vb;
    Complex p0, kc, wz_0, Zone, lnu, ldelta;
    
    if(ixyz<NUMBER_ELEMENT)
    {
        decode_ixyz2ixiyiz(ixyz,ix,iy,iz, i); // decode cartesian coordinates
        
        // load data to registers from global memory and form constant part of potentials
        Va=u_ext(ix,iy,iz,it,SPINA);
        Vb=u_ext(ix,iy,iz,it,SPINB);
        
        t5 = 1.0/ (dc_sclgth);
        lnu = h_densities.nu[ixyz];
        Zone = Complex(1.0, 0.0);
  
        // pairing
#ifdef USE_CUBIC_CUTOFF
        wz_0=Complex(REGULARIZATION_SCHEME_K_CONST/(4.0*M_PI*DX), 0.0); // FIXME: account for effective mass
#else
        t7=(dc_mu_a-Va+dc_mu_b-Vb)/2.0;
        p0 = thrust::sqrt( Complex(2.0*t7, 0.0) );
        if(p0.imag()<0.) p0 *= -1. ;
        kc = thrust::sqrt( Complex(2.0*(dc_ec+t7), 0.0) );
        if(kc.imag()<0.) kc *= -1. ;
            
        wz_0 = thrust::log( ( kc + p0 ) / ( kc - p0 ) ) ;
        if ( wz_0.imag() < 0. ) wz_0 += Complex(0.0, 2. * M_PI) ;    
        wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
#endif
        wz_0 = Zone / (Zone*t5 - wz_0);
        // g_eff = wz_0.real(); 
        ldelta = lnu*(-1.0*wz_0.real());
                    
        // save results to global memory
        h_potentials.V_a[ixyz]=Va;
        h_potentials.V_b[ixyz]=Vb;
        h_potentials.delta[ixyz]=ldelta;
    }
}


__global__ void tdwslda_compute_energy(int it, wslda_density h_densities, wslda_potential h_potentials, // <-- INPUT
                                       double *E_kin, double *E_pot, double *E_pair, double *E_curr)    // <- OUTPUT
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    
    double na, nb;
    double taua, taub;
    double tx1, ty1, tz1, tx2, ty2, tz2;
    
    if(ixyz<NUMBER_ELEMENT)
    {
        // densities, and correct them
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        
        // kinetic energy
        taua=h_densities.tau_a[ixyz]; // tau_a
        taub=h_densities.tau_b[ixyz]; // tau_b  
        
        // current corrections
#if CODEDIM>=1
        tx1=h_densities.j_a_x[ixyz];
        tx2=h_densities.j_b_x[ixyz];
#endif
#if CODEDIM>=2
        ty1=h_densities.j_a_y[ixyz];
        ty2=h_densities.j_b_y[ixyz];
#else        
        ty1=0.0;
        ty2=0.0;
#endif
#if CODEDIM>=3
        tz1=h_densities.j_a_z[ixyz];
        tz2=h_densities.j_b_z[ixyz];
#else
        tz1=0.0;
        tz2=0.0;
#endif
        taua-=p_regularization(na)*(tx1*tx1+ty1*ty1+tz1*tz1)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(tx2*tx2+ty2*ty2+tz2*tz2)/nb; // -jb^2/nb: correction for tilde{tau}_b

        // galilean invariant contribution
        E_kin[ixyz]=0.5*(h_potentials.alpha_a[ixyz]*taua + h_potentials.alpha_b[ixyz]*taub)*VOLUME_ELEMENT;
        
        // potential energy
        E_pot[ixyz]=0.0;
        
        // pairing energy
        E_pair[ixyz]=(h_potentials.delta[ixyz]*thrust::conj(h_densities.nu[ixyz])).real()*(-1.0)*VOLUME_ELEMENT;
        
        // flow energy
        E_curr[ixyz]= p_regularization(na)*(tx1*tx1 + ty1*ty1 + tz1*tz1)/(2.*na)*VOLUME_ELEMENT  
                    + p_regularization(nb)*(tx2*tx2 + ty2*ty2 + tz2*tz2)/(2.*nb)*VOLUME_ELEMENT;

    }
}
#endif

// ------------------------------------------------------------------------------------
// ----------------------------------- SLDAE ------------------------------------------
// ------------------------------------------------------------------------------------
#if FUNCTIONAL==SLDAE
__global__ void tdwslda_compute_potentials(int it, wslda_density h_densities, wslda_potential h_potentials, double cccoeff)
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;
    
    
    int pccrSLDAe=0;
//     int pccrSLDAe=1;
    
//     double na,nb;
//     double Va_ext,Vb_ext;
    
    if(ixyz<NUMBER_ELEMENT)
    {
        decode_ixyz2ixiyiz(ixyz,ix,iy,iz, i); // decode cartesian coordinates
        
//         // external potential
//         Va_ext=u_ext(ix,iy,iz,it,SPINA);
//         Vb_ext=u_ext(ix,iy,iz,it,SPINB);
//         
        // read densities
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        
//         // save potential to global memory
//         h_potentials.V_a[ixyz]=...;  // <-- mean field + external potential
//         h_potentials.V_b[ixyz]=...;  // <-- mean field + external potential 
//         h_potentials.delta[ixyz]=...;  // <-- mean field + external potential
//         ...
    }
}


__global__ void tdwslda_compute_energy(int it, wslda_density h_densities, wslda_potential h_potentials, // <-- INPUT
                                       double *E_kin, double *E_pot, double *E_pair, double *E_curr)    // <- OUTPUT
{
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    
    double na, nb;
    double taua, taub;
    double jxa, jya, jza, jxb, jyb, jzb;
    
    if(ixyz<NUMBER_ELEMENT)
    {
        // densities
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        
        // kinetic energy
        taua=h_densities.tau_a[ixyz]; // tau_a
        taub=h_densities.tau_b[ixyz]; // tau_b  
        
        // currents
#if CODEDIM>=1
        jxa=h_densities.j_a_x[ixyz];
        jxb=h_densities.j_b_x[ixyz];
#endif
#if CODEDIM>=2
        jya=h_densities.j_a_y[ixyz];
        jyb=h_densities.j_b_y[ixyz];
#else        
        jya=0.0;
        jyb=0.0;
#endif
#if CODEDIM>=3
        jza=h_densities.j_a_z[ixyz];
        jzb=h_densities.j_b_z[ixyz];
#else
        jza=0.0;
        jzb=0.0;
#endif
        taua-=p_regularization(na)*(jxa*jxa+jya*jya+jza*jza)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(jxb*jxb+jyb*jyb+jzb*jzb)/nb; // -jb^2/nb: correction for tilde{tau}_b

        // galilean invariant contribution
        E_kin[ixyz]=0.5*(h_potentials.alpha_a[ixyz]*taua + h_potentials.alpha_b[ixyz]*taub)*VOLUME_ELEMENT;
        
        // potential energy
        E_pot[ixyz]=0.0;
        
        // pairing energy
        E_pair[ixyz]=(h_potentials.delta[ixyz]*thrust::conj(h_densities.nu[ixyz])).real()*(-1.0)*VOLUME_ELEMENT;
        
        // flow energy
        E_curr[ixyz]= p_regularization(na)*(jxa*jxa + jya*jya + jza*jza)/(2.*na)*VOLUME_ELEMENT  
                    + p_regularization(nb)*(jxb*jxb + jyb*jyb + jzb*jzb)/(2.*nb)*VOLUME_ELEMENT;

    }
}
#endif

// ------------------------------------------------------------------------------------
// --------------------------------- CUSTOMEDF ----------------------------------------
// ------------------------------------------------------------------------------------
#if FUNCTIONAL==CUSTOMEDF

// Functions are provided in problem-definition.h file

#endif

// DO NOT REMOVE!
#include "tdwslda_functionals_framework_disable.h"

