/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * The file implements routines for computation of potentials from densities
 * These routines are utilized by ST codes only!
 *  
 * @author Gabriel Wlazlowski
 * @date 04.09.2020
 * */ 
#include <stdio.h>
#include "pca_settings.h"
#include "pca_utils.h"
#include "wslda_functionals.h"

extern int wsldapid; // process id - global variable
#include "problem-definition.h"

#define Complex(a,b) (a + I*b)
#define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a)) 
// Number of self-consistent iterations for U and delta computation
#undef UD_SCITERS
#define UD_SCITERS 1000
// Mixing parameter for self-consistent algorithm - fraction of new solution used for mixing
#undef UD_MIX_COEFF
#define UD_MIX_COEFF 0.25
#define UD_EPSILON 1.0e-12

extern double *dc_params; /* Declaration of the variable */
extern size_t dc_extra_data_size;
extern void *dc_extra_data;
extern double dc_mu_a;
extern double dc_mu_b;
extern double dc_ec;

// --------------------------------------------------------------------------------------------------
// ----------------------------------------- Generic ------------------------------------------------
// --------------------------------------------------------------------------------------------------
/**
 * Function computes energy contributions from external fields
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials corresponding to the densities (INPUT) 
 * @param energy array with contributions to the energy (OUTPUT)
 * */
int compute_energy_ext(int it, wslda_density h_densities, wslda_potential h_potentials, double *energy)
{
    int lNX=h_densities.nx, lNY=h_densities.ny, lNZ=h_densities.nz; // local sizes
    
    // reset energy entries
    energy[EPOTEXT]=0.0;
    energy[EPAIREXT]=0.0;
    energy[EVELEXT]=0.0;
    
    int ix, iy, iz;
    int ixyz=0;
    double na, nb;
    
    for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    {
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        
        // External potential energy
        energy[EPOTEXT]+=na*v_ext(ix,iy,iz,it,SPINA,dc_params,dc_extra_data_size,dc_extra_data) + 
                         nb*v_ext(ix,iy,iz,it,SPINB,dc_params,dc_extra_data_size,dc_extra_data);
                         
        // External pairing energy -(Delta x nu^* + Delta^* x nu)=-2Re[Delta x nu^*] 
        energy[EPAIREXT]-=2.0*creal( 
                            h_densities.nu[ixyz]*conj(delta_ext(ix,iy,iz,it,h_potentials.delta[ixyz],dc_params,dc_extra_data_size,dc_extra_data)) 
                                   );
        
        // External energy related to coupling with external velocity filed: -v^{ext}(r)*j(r)
        // x component
        energy[EVELEXT]-= h_densities.j_a_x[ixyz]*velocity_ext(ix,iy,iz,it,SPINA,XAXIS,dc_params,dc_extra_data_size,dc_extra_data)
                         +h_densities.j_b_x[ixyz]*velocity_ext(ix,iy,iz,it,SPINB,XAXIS,dc_params,dc_extra_data_size,dc_extra_data);
        // y component
        energy[EVELEXT]-= h_densities.j_a_y[ixyz]*velocity_ext(ix,iy,iz,it,SPINA,YAXIS,dc_params,dc_extra_data_size,dc_extra_data)
                         +h_densities.j_b_y[ixyz]*velocity_ext(ix,iy,iz,it,SPINB,YAXIS,dc_params,dc_extra_data_size,dc_extra_data);
        // z component
        energy[EVELEXT]-= h_densities.j_a_z[ixyz]*velocity_ext(ix,iy,iz,it,SPINA,ZAXIS,dc_params,dc_extra_data_size,dc_extra_data)
                         +h_densities.j_b_z[ixyz]*velocity_ext(ix,iy,iz,it,SPINB,ZAXIS,dc_params,dc_extra_data_size,dc_extra_data);
                         
        ixyz++;
    }
    
    // take into account variants of code
    double volume_element=0.0;
    
    if(h_densities.datadim==3) volume_element=DX*DY*DZ;
    if(h_densities.datadim==2) volume_element=DX*DY*LZ;
    if(h_densities.datadim==1) volume_element=DX*LY*LZ;
    
    energy[EPOTEXT]*=volume_element;
    energy[EPAIREXT]*=volume_element;
    energy[EVELEXT]*=volume_element;
    
    return 0;
}


// --------------------------------------------------------------------------------------------------
// -------------------------------------- SLDA variant ----------------------------------------------
// --------------------------------------------------------------------------------------------------

/**
 * Function computes potentials for ALSDA functional
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials from PREVIOUS iteration as input, 
 *                     updated values as output (INPUT/OUTPUT) 
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
int compute_potentials_aslda(int it, wslda_density h_densities, wslda_potential h_potentials, double *params, size_t extra_data_size, void *extra_data)
{
    int lNX=h_densities.nx, lNY=h_densities.ny, lNZ=h_densities.nz; // local sizes
    
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
    
    
    // Code is equivivalent to the code implemented in pca_kernels.cu
    
    int ix, iy, iz, ixyz;
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
    double v_ext_a, v_ext_b;
    
    ixyz=0;
    for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    {
        // store value of potentials in separate variables, will be used later
        v_ext_a=v_ext(ix,iy,iz,it,SPINA,params,extra_data_size,extra_data);
        v_ext_b=v_ext(ix,iy,iz,it,SPINB,params,extra_data_size,extra_data);     
        
        // start computation of delta and mean-field
        Va_const=v_ext_a;
        Vb_const=v_ext_b;
        
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
        Va = V_a[ixyz]+v_ext_a; // initial values
        Vb = V_b[ixyz]+v_ext_b; // initial values
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
        
        // save potentials
        V_a[ixyz]=Va-v_ext_a; // subtract external potential to get mean-field only 
        V_b[ixyz]=Vb-v_ext_b; // subtract external potential to get mean-field only 
        delta[ixyz]=ldelta;       
        
        // potentials that do not require self-cosistent iteration
        t1 = polarization(na, nb);
        h_potentials.alpha_a[ixyz]=alpha_a(t1); 
        h_potentials.alpha_b[ixyz]=alpha_b(t1); 
        
        t7 = p_regularization(rho_a[ixyz]);
        if(t7!=0.0)
        {
            h_potentials.A_a_x[ixyz]=(1.0 - h_potentials.alpha_a[ixyz]) * t7 *j_a_x[ixyz] / rho_a[ixyz]; // derivative with respect to j_a_x 
            h_potentials.A_a_y[ixyz]=(1.0 - h_potentials.alpha_a[ixyz]) * t7 *j_a_y[ixyz] / rho_a[ixyz]; // derivative with respect to j_a_y 
            h_potentials.A_a_z[ixyz]=(1.0 - h_potentials.alpha_a[ixyz]) * t7 *j_a_z[ixyz] / rho_a[ixyz]; // derivative with respect to j_a_z 
        }
        else
        {
            h_potentials.A_a_x[ixyz]=0.0; // we are in vacum 
            h_potentials.A_a_y[ixyz]=0.0; // we are in vacum 
            h_potentials.A_a_z[ixyz]=0.0; // we are in vacum  
        }
        
        t7 = p_regularization(rho_b[ixyz]);
        if(t7!=0.0)
        {
            h_potentials.A_b_x[ixyz]=(1.0 - h_potentials.alpha_b[ixyz]) * t7 *j_b_x[ixyz] / rho_b[ixyz]; // derivative with respect to j_b_x 
            h_potentials.A_b_y[ixyz]=(1.0 - h_potentials.alpha_b[ixyz]) * t7 *j_b_y[ixyz] / rho_b[ixyz]; // derivative with respect to j_b_y 
            h_potentials.A_b_z[ixyz]=(1.0 - h_potentials.alpha_b[ixyz]) * t7 *j_b_z[ixyz] / rho_b[ixyz]; // derivative with respect to j_b_z 
        }
        else
        {
            h_potentials.A_b_x[ixyz]=0.0; // we are in vacum 
            h_potentials.A_b_y[ixyz]=0.0; // we are in vacum 
            h_potentials.A_b_z[ixyz]=0.0; // we are in vacum  
        }
        
        ixyz++; // go to next point
        
    } // for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    
    
    return 0;
}

/**
 * Function that computes energy of the system
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials corresponding to the densities (INPUT) 
 * @param energy array with contributions to the energy (OUTPUT)
 * @param npart array with contributions to the particle number (OUTPUT)
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
int compute_energy_aslda(int it, wslda_density h_densities, wslda_potential h_potentials, double *energy, double *npart, double *params, size_t extra_data_size, void *extra_data)
{    
    int lNX=h_densities.nx, lNY=h_densities.ny, lNZ=h_densities.nz; // local sizes
    
    // reset buffers
    energy[EKIN]=0.0;       //
    energy[EPOT]=0.0;       // 
    energy[EPAIR]=0.0;      // Together they give intristic energy of the system
    energy[ECURRENT]=0.0;   // 
    npart[SPINA]=0.0; npart[SPINB]=0.0;
    
    double na, nb;
    double taua, taub;
    double tx1, ty1, tz1, tx2, ty2, tz2;
    
    int ix, iy, iz;
    int ixyz=0;
    for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    {

        // densities
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        
        // particle number
        npart[SPINA]+=na;
        npart[SPINB]+=nb;
                                 
        // kinetic energy - not that I use here Galilean invariant definition of tau
        taua=h_densities.tau_a[ixyz]; // tau_a
        taub=h_densities.tau_b[ixyz]; // tau_b    
        
        // current corrections
        tx1=h_densities.j_a_x[ixyz];
        ty1=h_densities.j_a_y[ixyz];
        tz1=h_densities.j_a_z[ixyz];
        tx2=h_densities.j_b_x[ixyz];
        ty2=h_densities.j_b_y[ixyz];
        tz2=h_densities.j_b_z[ixyz];
        taua-=p_regularization(na)*(tx1*tx1+ty1*ty1+tz1*tz1)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(tx2*tx2+ty2*ty2+tz2*tz2)/nb; // -jb^2/nb: correction for tilde{tau}_b

        // kinetic energy, only aglilean invariant contribution
        energy[EKIN]+=0.5*(h_potentials.alpha_a[ixyz]*taua + h_potentials.alpha_b[ixyz]*taub);
        
        // potential energy
        energy[EPOT]+=funD(na, nb);
        
        // pairing energy
        energy[EPAIR]-=creal(h_potentials.delta[ixyz]*conj(h_densities.nu[ixyz]));
        
        // energy related to currents (j^2/2n)
        energy[ECURRENT]+=   p_regularization(na)*(tx1*tx1 + ty1*ty1 + tz1*tz1)/(2.*na)  
                           + p_regularization(nb)*(tx2*tx2 + ty2*ty2 + tz2*tz2)/(2.*nb);  

        ixyz++;
    }
    
    // take into account variants of code
    double volume_element=0.0;
    
    if(h_densities.datadim==3) volume_element=DX*DY*DZ;
    if(h_densities.datadim==2) volume_element=DX*DY*LZ;
    if(h_densities.datadim==1) volume_element=DX*LY*LZ;

    energy[EKIN]*=volume_element;
    energy[EPOT]*=volume_element;
    energy[EPAIR]*=volume_element;
    energy[ECURRENT]*=volume_element;
    npart[SPINA]*=volume_element;
    npart[SPINB]*=volume_element;
    
    return 0;
}


// --------------------------------------------------------------------------------------------------
// -------------------------------------- BdG variant -----------------------------------------------
// --------------------------------------------------------------------------------------------------
extern double aBdG; // scattering length
/**
 * Function computes potentials for BDG functional
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials from PREVIOUS iteration as input, 
 *                     updated values as output (INPUT/OUTPUT) 
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
int compute_potentials_bdg(int it, wslda_density h_densities, wslda_potential h_potentials, double *params, size_t extra_data_size, void *extra_data)
{
    int lNX=h_densities.nx, lNY=h_densities.ny, lNZ=h_densities.nz; // local sizes
    
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
    double * alpha_a = h_potentials.alpha_a;
    double * alpha_b = h_potentials.alpha_b;
    
    // Code is equivivalent to the code implemented in pca_kernels.cu
    
    int ix, iy, iz, ixyz;
    
    // registers
    double t1, t2, t3, t4, t5, t6, t7; // working buffers
    
    double Va, Vb;
    double complex p0, kc, wz_0, Zone, lnu, ldelta;
    double v_ext_a, v_ext_b;
    double alph_plus;
    
    ixyz=0;
    for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    {
        v_ext_a=v_ext(ix,iy,iz,it,SPINA,params,extra_data_size,extra_data);
        v_ext_b=v_ext(ix,iy,iz,it,SPINB,params,extra_data_size,extra_data);     
        
        // start computation of delta
        alph_plus = 0.5*(alpha_a[ixyz]+alpha_b[ixyz]);
        t5 = 1.0/ (4.0*M_PI*aBdG*alph_plus);
        Va = V_a[ixyz]+v_ext_a; // initial values
        Vb = V_b[ixyz]+v_ext_b; // initial values
        lnu = nu[ixyz];
        Zone = Complex(1.0, 0.0);
        
        // pairing
        t7=(dc_mu_a-Va+dc_mu_b-Vb)/2.0;
        p0 = csqrt( Complex(2.0*t7/ alph_plus, 0.0) );
        if(cimag(p0)<0.) p0 *= -1. ;
        kc = csqrt( Complex(2.0*(dc_ec+t7)/ alph_plus, 0.0) );
        if(cimag(kc)<0.) kc *= -1. ;
        
        wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
        if ( cimag(wz_0) < 0. ) wz_0 += Complex(0.0, 2. * M_PI) ;    
        wz_0= kc / ( 2. * M_PI * M_PI * alph_plus) *( 1. - p0 / ( 2. * kc ) * wz_0);
        wz_0 = Zone / (Zone*t5 - wz_0);
        // g_eff = wz_0.real(); 
        ldelta = lnu*(-1.0*creal(wz_0));
        
        // save potentials
        V_a[ixyz]=0.0; // no mean-filed potential in BDG
        V_b[ixyz]=0.0; // no mean-filed potential in BDG 
        
        delta[ixyz]=ldelta; // only pairing in BDG 
        
        h_potentials.alpha_a[ixyz]=1.0; // bare mass 
        h_potentials.alpha_b[ixyz]=1.0; // bare mass
        
        h_potentials.A_a_x[ixyz]=0.0; // no vector potential in BDG
        h_potentials.A_a_y[ixyz]=0.0; // no vector potential in BDG 
        h_potentials.A_a_z[ixyz]=0.0; // no vector potential in BDG 
        h_potentials.A_b_x[ixyz]=0.0; // no vector potential in BDG
        h_potentials.A_b_y[ixyz]=0.0; // no vector potential in BDG
        h_potentials.A_b_z[ixyz]=0.0; // no vector potential in BDG 
        
        ixyz++; // go to next point
        
    } // for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    
    
    return 0;
}

/**
 * Function that computes energy of the system
 * @param it iteration number
 * @param h_densities array with all densities (INPUT)
 * @param h_potentials potentials corresponding to the densities (INPUT) 
 * @param energy array with contributions to the energy (OUTPUT)
 * @param npart array with contributions to the particle number (OUTPUT)
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
int compute_energy_bdg(int it, wslda_density h_densities, wslda_potential h_potentials, double *energy, double *npart, double *params, size_t extra_data_size, void *extra_data)
{    
    int lNX=h_densities.nx, lNY=h_densities.ny, lNZ=h_densities.nz; // local sizes
    
    // reset buffers
    energy[EKIN]=0.0;       //
    energy[EPOT]=0.0;       // 
    energy[EPAIR]=0.0;      // Together they give intristic energy of the system
    energy[ECURRENT]=0.0;   // 
    npart[SPINA]=0.0; npart[SPINB]=0.0;
    
    double tx1, ty1, tz1, tx2, ty2, tz2;
    double na, nb, taua, taub;    
    int ix, iy, iz;
    int ixyz=0;
    for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    {

        // densities
        na=h_densities.rho_a[ixyz];
        nb=h_densities.rho_b[ixyz];
        
        // particle number
        npart[SPINA]+=na;
        npart[SPINB]+=nb;
                                 
        // kinetic energy
        taua=h_densities.tau_a[ixyz]; // tau_a
        taub=h_densities.tau_b[ixyz]; // tau_b
        
        // current corrections
        tx1=h_densities.j_a_x[ixyz];
        ty1=h_densities.j_a_y[ixyz];
        tz1=h_densities.j_a_z[ixyz];
        tx2=h_densities.j_b_x[ixyz];
        ty2=h_densities.j_b_y[ixyz];
        tz2=h_densities.j_b_z[ixyz];
        taua-=p_regularization(na)*(tx1*tx1+ty1*ty1+tz1*tz1)/na; // -ja^2/na: correction for tilde{tau}_a
        taub-=p_regularization(nb)*(tx2*tx2+ty2*ty2+tz2*tz2)/nb; // -jb^2/nb: correction for tilde{tau}_b
        
        // standard contribution, gallileant invariance only of alpha_a=alpha_b=1
        energy[EKIN]+=0.5*(h_potentials.alpha_a[ixyz]*taua + h_potentials.alpha_b[ixyz]*taub);
        
        // potential energy
        energy[EPOT]+=0.0;
        
        // pairing energy
        energy[EPAIR]-=creal(h_potentials.delta[ixyz]*conj(h_densities.nu[ixyz]));
        
        // energy related to currents (j^2/2n)
        energy[ECURRENT]+=   p_regularization(na)*(tx1*tx1 + ty1*ty1 + tz1*tz1)/(2.*na)  
                           + p_regularization(nb)*(tx2*tx2 + ty2*ty2 + tz2*tz2)/(2.*nb); 
        
        ixyz++;
    }
    
    // take into account variants of code
    double volume_element=0.0;
    
    if(h_densities.datadim==3) volume_element=DX*DY*DZ;
    if(h_densities.datadim==2) volume_element=DX*DY*LZ;
    if(h_densities.datadim==1) volume_element=DX*LY*LZ;

    energy[EKIN]*=volume_element;
    energy[EPOT]*=volume_element;
    energy[EPAIR]*=volume_element;
    energy[ECURRENT]*=volume_element;
    npart[SPINA]*=volume_element;
    npart[SPINB]*=volume_element;
    
    return 0;
}

// --------------------------------------------------------------------------------------------------
// -------------------------------------- Selector --------------------------------------------------
// --------------------------------------------------------------------------------------------------
int compute_potentials(int it, wslda_density h_densities, wslda_potential h_potentials)
{
#if FUNCTIONAL==BDG
    return compute_potentials_bdg(it, h_densities, h_potentials, dc_params,dc_extra_data_size,dc_extra_data);
#elif FUNCTIONAL==SLDA    
    return compute_potentials_aslda(it, h_densities, h_potentials, dc_params,dc_extra_data_size,dc_extra_data);
#elif FUNCTIONAL==ASLDA    
    return compute_potentials_aslda(it, h_densities, h_potentials, dc_params,dc_extra_data_size,dc_extra_data);    
#elif FUNCTIONAL==CUSTOMEDF    
    return compute_potentials_custom(it, h_densities, h_potentials, dc_params,dc_extra_data_size,dc_extra_data); 
#endif
}

int compute_energy(int it, wslda_density h_densities, wslda_potential h_potentials, double *energy, double *npart)
{
    compute_energy_ext(it, h_densities, h_potentials, energy);
#if FUNCTIONAL==BDG
    return compute_energy_bdg(it, h_densities, h_potentials, energy, npart, dc_params,dc_extra_data_size,dc_extra_data);
#elif FUNCTIONAL==SLDA    
    return compute_energy_aslda(it, h_densities, h_potentials, energy, npart, dc_params,dc_extra_data_size,dc_extra_data);
#elif FUNCTIONAL==ASLDA    
    return compute_energy_aslda(it, h_densities, h_potentials, energy, npart, dc_params,dc_extra_data_size,dc_extra_data);
#elif FUNCTIONAL==CUSTOMEDF    
    return compute_energy_custom(it, h_densities, h_potentials, energy, npart, dc_params,dc_extra_data_size,dc_extra_data);
#endif
}
