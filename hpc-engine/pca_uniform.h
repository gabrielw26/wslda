// Author: Gabriel Wlazlowski
// Date: 14-07-2016

#ifndef __PCA_UNIFORM__
#define __PCA_UNIFORM__

// The functions are used to solve ASLDA for uniform system

// structure that contains info needed to fast reconstruct the solution
typedef struct
{
    double n0_a;
    double n0_b;
    double mu_a;
    double mu_b;
    double V_a;
    double V_b;
    double alph_a;
    double alph_b;
    double kc;
    double ec;
    double nu;
    double delta;
    double tau_a;
    double tau_b;
    double beta;
    double ekin;
    double epot;
    double epair;
    int nwf;
} metadata_pca_uniform_t;

// Stores results;
metadata_pca_uniform_t __md_pca_uniform;


inline double fbeta(double E, double beta)
{
    double bE=beta*E;
    if(bE>50.) return 0.0; // to avoid numerical problems
    else if(bE<-50.) return 1.0; // to avoid numerical problems
    else return 1.0/(exp(bE)+1.0);
}

/**
 * @param n0_a requested density for population "a" (INPUT)
 * @param n0_b requested density for population "b" (INPUT)
 * @param nwf number of wave-functions (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int solve_uniform_problem(double n0_a, double n0_b, int *nwf, int printout)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    }    
    
    double * kk2tau ; /* kk2 = k^2 */
    cppmallocl(kk2tau,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        double _kkx=kkx[ix];
        double _kky=kky[iy];
        double _kkz=kkz[iz];
#ifdef TAU_COMPUTATION_VIA_GRADIENTS       
        if(ix==NX/2) _kkx=0.0;
        if(iy==NY/2) _kky=0.0;
        if(iz==NZ/2) _kkz=0.0;
#endif
        kk2tau[i]=_kkx*_kkx + _kky*_kky + _kkz*_kkz;
        i++;
    }
    
    
    // Set quantities that depend only on desities
    double p = polarization_h(n0_a, n0_b);
    double alph_a = alpha_a_h(p);
    double alph_b = alpha_b_h(p);
    double dalphm_dna=der_alpha_minus__der_na_h(n0_a, n0_b);
    double dalphm_dnb=der_alpha_minus__der_nb_h(n0_a, n0_b);
    double dalphp_dna=der_alpha_plus__der_na_h(n0_a, n0_b);
    double dalphp_dnb=der_alpha_plus__der_nb_h(n0_a, n0_b);
    double alph_plus = alpha_plus_h(p);
    double dtildeC_dna = der_tildeC__der_na_h(n0_a, n0_b);
    double dtildeC_dnb = der_tildeC__der_nb_h(n0_a, n0_b);
    double dD_dna = der_funD__der_na_h(n0_a, n0_b);
    double dD_dnb = der_funD__der_nb_h(n0_a, n0_b);
    double tC = tildeC_h(n0_a, n0_b);
    double D = funD_h(n0_a, n0_b);
    double eF_a=pow(6.0*M_PI*M_PI*n0_a, 2.0/3.0) / 2.0;
    double eF_b=pow(6.0*M_PI*M_PI*n0_b, 2.0/3.0) / 2.0;
    double eF_avg=pow(3.0*M_PI*M_PI*(n0_a+n0_b), 2.0/3.0) / 2.0;
    double Effg = 0.6*n0_a*eF_a*LXYZ + 0.6*n0_b*eF_b*LXYZ;
    double kc=md.kc;
    double mu_a=0.37*eF_a;
    double mu_b=0.37*eF_b;
    
    if(printout && md.init0debug>0) wprintf("# DEBUG: n_a=%f, n_b=%f\n", n0_a, n0_b);
    if(printout && md.init0debug>0) wprintf("# DEBUG: eF_a=%f, eF_b=%f, eF_avg=%f\n", eF_a, eF_b, eF_avg);
    if(printout && md.init0debug>0) wprintf("# DEBUG: N_a=%f, N_b=%f\n", n0_a*LXYZ, n0_b*LXYZ);
    if(printout && md.init0debug>0) wprintf("# DEBUG: p=%f\n", p);
    if(printout && md.init0debug>0) wprintf("# DEBUG: alph_a=%f, alph_b=%f, alph_plus=%f\n", alph_a, alph_b, alph_plus); 
    if(printout && md.init0debug>0) wprintf("# DEBUG: dalphm_dna=%f, dalphm_dnb=%f\n", dalphm_dna, dalphm_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: dalphp_dna=%f, dalphp_dnb=%f\n", dalphp_dna, dalphp_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: dtildeC_dna=%f, dtildeC_dnb=%f\n", dtildeC_dna, dtildeC_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: dD_dna=%f, dD_dnb=%f\n", dD_dna, dD_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: D=%f, tC=%f\n", D, tC);
    if(printout && md.init0debug>0 && md.spinsymmetry>0) wprintf("# SPIN SYMMETRY MODE!\n");

    
    // Set quantities updated in s-c loop
    double tau_a=pow(6.0*M_PI*M_PI*n0_a, 5.0/3.0) / (10.0*M_PI*M_PI); // initial value
    double tau_b=pow(6.0*M_PI*M_PI*n0_b, 5.0/3.0) / (10.0*M_PI*M_PI); // initial value
    double delta;
    if(n0_a<n0_b) delta = 0.5*pow(6.0*M_PI*M_PI*n0_a, 2.0/3.0) / 2.0; // initial value: 0.5*eF
    else          delta = 0.5*pow(6.0*M_PI*M_PI*n0_b, 2.0/3.0) / 2.0; // initial value: 0.5*eF
    double nu= -1.0 * delta * tC / alph_plus; // initial value
    
    if(printout && md.init0debug>0) wprintf("# DEBUG: tau_a=%f, tau_b=%f\n", tau_a, tau_b);   
    if(printout && md.init0debug>0) wprintf("# DEBUG: delta=%f, nu=%f\n", delta, nu);  
    
    // auxliary variables
    int maxiter=md.init0maxiter;
    int iter;
    double tau_m, tau_p; 
    double V_a, V_b, mu_p, g_eff, eta_a, eta_b;
    double complex p0, wz_0;
    double complex Zzero = 0.0 + I*0.0;
    double complex Zone  = 1.0 + I*0.0;
    double uk, vk, ek;
    double n_a, n_b;
    double tau_a_old, tau_b_old, delta_old, nu_old;
    double scmix=md.init0scmix, epsilon=md.init0eps;
    int is_conv;
    double beta, T;
    double energy_kin, energy_pot, energy_pair, energy_tot;
    double kc2=kc*kc;
    
    // extract number of wave-functions
    ixyz=0;
    *nwf=0;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2) // take only from sphere
            *nwf+=2; // thera are two soloutions for each momentum
        ixyz++;
    }
    
    for(T=md.init0Tstart; T>=md.init0Tstop; T-=md.init0DeltaT) 
    {
        beta=1.0/(T*eF_avg);
        for(iter=0; iter<maxiter; iter++)
        {
            // save old values of potentials
            tau_a_old=tau_a;
            tau_b_old=tau_b;
            delta_old=delta;
            nu_old=nu;
            
            // potential
            tau_p=tau_a + tau_b;
            tau_m=tau_a - tau_b;
            
            V_a = dalphm_dna*tau_m/2.0 + dalphp_dna*(tau_p/2.0 - delta*nu/alph_plus) - dtildeC_dna*delta*delta/alph_plus + dD_dna;
            V_b = dalphm_dnb*tau_m/2.0 + dalphp_dnb*(tau_p/2.0 - delta*nu/alph_plus) - dtildeC_dnb*delta*delta/alph_plus + dD_dnb;
            
            // pairing
#ifdef USE_CUBIC_CUTOFF
            wz_0=REGULARIZATION_SCHEME_K_CONST/(4.0*M_PI*DX) + I*0.0;
#else
            mu_p=(mu_a-V_a+mu_b-V_b)/2.0;
            p0 = csqrt( 2.0*mu_p/ alph_plus) ;
            if ( cimag(p0) < 0. ) p0 *= -1. ;
            
            // kc is fixed, and finally will be trananslated into ec
            wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( cimag(wz_0) < 0. ) wz_0 += I * 2. * M_PI ;    
                    
            wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
#endif
            g_eff = creal( Zone*alph_plus / (Zone*tC - wz_0) );
            delta = -1.0*g_eff*nu;
            
            // contribution from states to densities
            n_a=0.0;
            n_b=0.0;
            tau_a=0.0;
            tau_b=0.0;
            nu=0.0;
            ixyz=0;
            *nwf=0;
            for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
            {
                if(kk2[ixyz]<kc2) // take only from sphere
                {
                    eta_a = alph_a*kk2[ixyz]/2.0 + V_a - mu_a;
                    eta_b = alph_b*kk2[ixyz]/2.0 + V_b - mu_b;
                    
                    // solution 1:
                    ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
                    vk = delta*delta / ( pow(0.5*(eta_a+eta_b)+0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta);
                    uk = 1.0 - vk;
    //                 if(printout) if(vk<0) wprintf("S1v: problem\n");
    //                 if(printout) if(uk<0) wprintf("S1u: problem\n");
                    vk = sqrt(vk); uk=sqrt(uk);
                    if(eta_b+ek<0.0) vk=-1.0*vk;
                    
                    if(md.spinsymmetry>0)
                    {
                        if(ek>0.0) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b 
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        n_a+=uk*uk*fbeta(ek,beta);
                        tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                        n_b+=vk*vk*fbeta(-1.0*ek,beta);
                        tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                        nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                        (*nwf)++;
                    }
                    
                    // solution 2:
                    ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
                    vk = delta*delta / ( pow(0.5*(eta_a+eta_b)-0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta) ;
                    uk = 1.0 - vk;
    //                 if(printout) if(vk<0) wprintf("S2v: problem\n");
    //                 if(printout) if(uk<0) wprintf("S2u: problem: %f %f\n", vk, uk);
                    vk = sqrt(vk); uk=sqrt(uk);
                    if(eta_b+ek<0.0) vk=-1.0*vk;
                    
                    if(md.spinsymmetry>0)
                    {
                        if(ek>0.0) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b 
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        n_a+=uk*uk*fbeta(ek,beta);
                        tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                        n_b+=vk*vk*fbeta(-1.0*ek,beta);
                        tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                        nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                        (*nwf)++;
                    }
                }
                
                ixyz++;
            }
            
            // add normalization factors
            if(md.spinsymmetry>0)
            {
                n_b/=LXYZ; n_a=n_b;
                tau_b/=LXYZ; tau_a=tau_b; 
                nu/=2.0*LXYZ;
            }
            else
            {
                n_a/=LXYZ; n_b/=LXYZ;
                tau_a/=LXYZ; tau_b/=LXYZ;
                nu/=2.0*LXYZ;
            }
            
            if(printout && md.init0debug>1) wprintf("D: iter=%d: V_a=%f, V_b=%f, delta=%f, n_a=%f, n_b=%f, tau_a=%f, tau_b=%f, nu=%f\n", iter, V_a, V_b, delta, n_a, n_b, tau_a, tau_b, nu);
            
            // check convergence
            if(printout && md.init0debug>1) wprintf("C: iter=%d: fabs(n0_a-n_a)=%g fabs(n0_b-n_b)=%g fabs(delta-delta_old)=%g, mu_a=%f, mu_b=%f\n", iter, fabs(n0_a-n_a), fabs(n0_b-n_b), fabs(delta-delta_old), mu_a, mu_b);
            is_conv=1;
            if(fabs(n0_a-n_a)>epsilon) is_conv=0; // check for density
            if(fabs(n0_b-n_b)>epsilon) is_conv=0; // check for density
            if(fabs(delta-delta_old)>epsilon) is_conv=0; // check for delta
            if(is_conv) break; // we converged
            
            // mix potentials and go to next iteration
            tau_a = (1.-scmix)*tau_a_old + scmix*tau_a;
            tau_b = (1.-scmix)*tau_b_old + scmix*tau_b;
            delta = (1.-scmix)*delta_old + scmix*delta;
            nu = (1.-scmix)*nu_old + scmix*nu;
            mu_a += md.init0muchange*(md.init0Tstart/T)*(n0_a-n_a);
            mu_b += md.init0muchange*(md.init0Tstart/T)*(n0_b-n_b);
            if(md.spinsymmetry>0) mu_b=mu_a;
        }
        if(printout && md.init0debug>0) wprintf("# TEMPCONV: T=%f, iter=%d, delta/eF_a=%f, mu_a/eF_a=%f, delta/eF_b=%f, mu_b/eF_b=%f\n", T, iter, delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b);
        if(iter==maxiter && printout && md.init0debug>0) wprintf("# WARNING: MAXITER REACHED!\n");
            
        // Compute energy 
        energy_kin=(0.5*alph_a*tau_a + 0.5*alph_b*tau_b)*LXYZ;
        energy_pot=(D)*LXYZ;
        energy_pair=-1.0*delta*nu*LXYZ;
        energy_tot=energy_kin+energy_pot+energy_pair;
        if(printout  && md.init0debug>0) wprintf("# TEMPCONV: T=%f, energy_kin=%f, energy_pot=%f, energy_pair=%f, energy_tot=%f\n", T, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
        fflush(stdout);
    }
    
    
    double ec=0.0;
    eta_a = alph_a*kc2/2.0 + V_a - mu_a;
    eta_b = alph_b*kc2/2.0 + V_b - mu_b;
    ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC1=%f\n", ek);
    if(fabs(ek)>ec) ec=fabs(ek);
    ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC2=%f\n", ek);
    if(fabs(ek)>ec) ec=fabs(ek);
    
    mu_p=(mu_a-V_a+mu_b-V_b)/2.0;
    ec=alph_plus*kc2/2.0 - mu_p;
    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC=%f\n", ec);
    
    // print results
    if(printout) wprintf("# UNIFORM SOLUTION: delta/eF_a=%8.4f, mu_a/eF_a=%8.4f, delta/eF_b=%8.4f, mu_b/eF_b=%8.4f, ec=%8.4f\n", delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b, ec);
    if(printout) wprintf("# UNIFORM SOLUTION: energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
    if(printout) wprintf("# UNIFORM SOLUTION: nwf=%d\n", *nwf);
    
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    free(kk2tau);
    
    // write data to structure for future use
    __md_pca_uniform.n0_a=n_a;
    __md_pca_uniform.n0_b=n_b;
    __md_pca_uniform.mu_a=mu_a;
    __md_pca_uniform.mu_b=mu_b;
    __md_pca_uniform.V_a=V_a;
    __md_pca_uniform.V_b=V_b;
    __md_pca_uniform.alph_a=alph_a;
    __md_pca_uniform.alph_b=alph_b;
    __md_pca_uniform.kc=kc;
    __md_pca_uniform.ec=ec;
    __md_pca_uniform.nu=nu;
    __md_pca_uniform.delta=delta;
    __md_pca_uniform.tau_a=tau_a;
    __md_pca_uniform.tau_b=tau_b;
    __md_pca_uniform.beta=beta;    
    __md_pca_uniform.nwf=*nwf;
    __md_pca_uniform.ekin=energy_kin;
    __md_pca_uniform.epot=energy_pot;
    __md_pca_uniform.epair=energy_pair;
    
    if(iter==maxiter) return 1; // not converged!
    return 0;
}

/** 
 * Function save data needed to reconstruct uniform solution
 * @return 0 - OK, otherwise error
 * */
int save_uniform()
{
    char filename[512];
    sprintf(filename, "%s/uniform.solution", md.outprefix);
    wprintf("# UNIFORM SAVE: Creating file with solution: `%s`\n", filename);
    int fexist = exists(filename); 
    if(fexist)
    {
        if(md.overwrite) 
        {
            wprintf("# UNIFORM SAVE: File `%s` exists. Removing [overwrite=%d]\n", filename, md.overwrite);
            urm(filename);
        }
        else
        {
            wprintf("# UNIFORM SAVE: Error: File `%s` exists. [overwrite=%d] \n", filename, md.overwrite);
            return WSLDA_ERR_CANNOT_OVERWRITE;
        }
    }

    FILE * fout=fopen(filename, "wb");
    
    size_t test_ele = fwrite (&__md_pca_uniform , sizeof(metadata_pca_uniform_t), 1, fout);
    if(test_ele!=1) return 2; // data not written 
    
    fclose(fout);
    
    return 0;
    
    
}

/** 
 * Function reads data needed to construct uniform solution
 * @param nwf number of wave-functions (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int read_uniform(int *nwf, int printout)
{
    char filename[512];
    sprintf(filename, "%s/uniform.solution", md.inprefix);
    wprintf("# UNIFORM READ: Reading data from file: `%s`\n", filename);
    int fexist = exists(filename); 
    if(!fexist) return 1;

    FILE * fout=fopen(filename, "rb");
    
    size_t test_ele = fread (&__md_pca_uniform , sizeof(metadata_pca_uniform_t), 1, fout);
    if(test_ele!=1) return 2; // data not read 
    
    fclose(fout);
    
    // provide info about the solution
    // Set quantities that depend only on desities
    double n0_a=__md_pca_uniform.n0_a;
    double n0_b=__md_pca_uniform.n0_b;
    double p = polarization_h(n0_a, n0_b);
    double alph_a = alpha_a_h(p);
    double alph_b = alpha_b_h(p);
    double dalphm_dna=der_alpha_minus__der_na_h(n0_a, n0_b);
    double dalphm_dnb=der_alpha_minus__der_nb_h(n0_a, n0_b);
    double dalphp_dna=der_alpha_plus__der_na_h(n0_a, n0_b);
    double dalphp_dnb=der_alpha_plus__der_nb_h(n0_a, n0_b);
    double alph_plus = alpha_plus_h(p);
    double dtildeC_dna = der_tildeC__der_na_h(n0_a, n0_b);
    double dtildeC_dnb = der_tildeC__der_nb_h(n0_a, n0_b);
    double dD_dna = der_funD__der_na_h(n0_a, n0_b);
    double dD_dnb = der_funD__der_nb_h(n0_a, n0_b);
    double tC = tildeC_h(n0_a, n0_b);
    double D = funD_h(n0_a, n0_b);
    double eF_a=pow(6.0*M_PI*M_PI*n0_a, 2.0/3.0) / 2.0;
    double eF_b=pow(6.0*M_PI*M_PI*n0_b, 2.0/3.0) / 2.0;
    double eF_avg=pow(3.0*M_PI*M_PI*(n0_a+n0_b), 2.0/3.0) / 2.0;
    double Effg = 0.6*n0_a*eF_a*LXYZ + 0.6*n0_b*eF_b*LXYZ;
    double kc=__md_pca_uniform.kc;
    double mu_a=__md_pca_uniform.mu_a;
    double mu_b=__md_pca_uniform.mu_b;
    
    if(printout && md.init0debug>0) wprintf("# DEBUG: n_a=%f, n_b=%f\n", n0_a, n0_b);
    if(printout && md.init0debug>0) wprintf("# DEBUG: eF_a=%f, eF_b=%f, eF_avg=%f\n", eF_a, eF_b, eF_avg);
    if(printout && md.init0debug>0) wprintf("# DEBUG: N_a=%f, N_b=%f\n", n0_a*LXYZ, n0_b*LXYZ);
    if(printout && md.init0debug>0) wprintf("# DEBUG: p=%f\n", p);
    if(printout && md.init0debug>0) wprintf("# DEBUG: alph_a=%f, alph_b=%f, alph_plus=%f\n", alph_a, alph_b, alph_plus); 
    if(printout && md.init0debug>0) wprintf("# DEBUG: dalphm_dna=%f, dalphm_dnb=%f\n", dalphm_dna, dalphm_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: dalphp_dna=%f, dalphp_dnb=%f\n", dalphp_dna, dalphp_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: dtildeC_dna=%f, dtildeC_dnb=%f\n", dtildeC_dna, dtildeC_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: dD_dna=%f, dD_dnb=%f\n", dD_dna, dD_dnb);
    if(printout && md.init0debug>0) wprintf("# DEBUG: D=%f, tC=%f\n", D, tC);

    double tau_a=__md_pca_uniform.tau_a;
    double tau_b=__md_pca_uniform.tau_b;
    double delta=__md_pca_uniform.delta;
    double nu=__md_pca_uniform.nu;

    if(printout && md.init0debug>0) wprintf("# DEBUG: tau_a=%f, tau_b=%f\n", tau_a, tau_b);   
    if(printout && md.init0debug>0) wprintf("# DEBUG: delta=%f, nu=%f\n", delta, nu);  
    
    double ec=__md_pca_uniform.ec;
    *nwf=__md_pca_uniform.nwf;
    
    double energy_kin, energy_pot, energy_pair, energy_tot;
    // Compute energy 
    energy_kin=(0.5*alph_a*tau_a + 0.5*alph_b*tau_b)*LXYZ;
    energy_pot=(D)*LXYZ;
    energy_pair=-1.0*delta*nu*LXYZ;
    energy_tot=energy_kin+energy_pot+energy_pair;    
    
    // print results
    if(printout) wprintf("# UNIFORM SOLUTION: delta/eF_a=%8.4f, mu_a/eF_a=%8.4f, delta/eF_b=%8.4f, mu_b/eF_b=%8.4f, ec=%8.4f\n", delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b, ec);
    if(printout) wprintf("# UNIFORM SOLUTION: energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
    if(printout) wprintf("# UNIFORM SOLUTION: nwf=%d\n", *nwf);    
    return 0;
    
    
}

/**
 * @param idxfrom extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param idxto extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param wf pointer to array for wf of size nx*ny*nz*2*(idxto-idxfrom)*sizeof(double complex)
 * @param mu_a chemical potential for population "a" (OUTPUT)
 * @param mu_b chemical potential for population "b" (OUTPUT)
 * @param ec energy cut-off for the solution (OUTPUT)
 * @param fEn weights used for computation densities, fEn=fbeta(E_n), array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param En eigen energies, E_n, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int create_uniform_wf(int idxfrom, int idxto, double complex *wf, double *mu_a, double *mu_b, double *ec, double *fEn, double *En, int printout)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    int ix2, iy2, iz2, ixyz2; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    } 
    
    // extract number of wave-functions and check consistency
    double kc2=__md_pca_uniform.kc*__md_pca_uniform.kc;
    ixyz=0;
    int nwf=0;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2) // take only from sphere
            nwf+=2; // thera are two soloutions for each momentum
        ixyz++;
    }
    if(md.spinsymmetry==0) if(nwf!=__md_pca_uniform.nwf) return 1;
    if(idxfrom>=nwf) return 2;
    if(idxto>nwf) return 3;
    if(idxfrom>=idxto)return 4;
    
    // Create wf
    size_t shift;
    double V_a=__md_pca_uniform.V_a;
    double V_b=__md_pca_uniform.V_b;
    double eta_a, eta_b;
    double alph_a=__md_pca_uniform.alph_a;
    double alph_b=__md_pca_uniform.alph_b;
    *mu_a=__md_pca_uniform.mu_a;
    *mu_b=__md_pca_uniform.mu_b;
    double uk, vk, ek;
    double delta=__md_pca_uniform.delta;
    double sqrt_volume=sqrt((double) (LXYZ));
    *ec=__md_pca_uniform.ec;
    double beta=__md_pca_uniform.beta;
    int takeit;

    if(printout) wprintf("# UNIFORM CREATE WF: Creating wave-functions.\n");
    ixyz=0;
    nwf=-1;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2) // take only from sphere
        {
            eta_a = alph_a*kk2[ixyz]/2.0 + V_a - *mu_a;
            eta_b = alph_b*kk2[ixyz]/2.0 + V_b - *mu_b;    
            
            // solution 1
            ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            vk = delta*delta / ( pow(0.5*(eta_a+eta_b)+0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta);
            uk = 1.0 - vk;
            vk = sqrt(vk); uk=sqrt(uk);
            if(eta_b+ek<0.0) vk=-1.0*vk;    
            
            takeit=0;
            if(md.spinsymmetry>0)
            {
                if(ek>0.0) 
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                nwf++;
                if(nwf>=idxfrom && nwf<idxto) takeit=1;
            }
            
            if(takeit) // this is my wave-function
            {
                // u-components
                shift=NXYZ*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ ) for ( iz2 = 0 ; iz2 < NZ ; iz2++ ) 
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] + ( double ) iz2 * DZ * kkz[iz] ) ) * uk /  sqrt_volume;
                    ixyz2++;
                }
                
                // v-components
                shift=NXYZ*(idxto-idxfrom) + NXYZ*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ ) for ( iz2 = 0 ; iz2 < NZ ; iz2++ ) 
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] + ( double ) iz2 * DZ * kkz[iz] ) ) * vk /  sqrt_volume;
                    ixyz2++;
                }
                
                // weight
                fEn[nwf-idxfrom]=fbeta(ek,beta);
                
                // eigen energy
                En[nwf-idxfrom]=ek;
                
//                 // only for tests:
//                 if(nwf==idxfrom)
//                     wprintf("!!!!!!!!! nwf=%d: kx=%f ky=%f kz=%f kk2=%f\n", nwf, kkx[ix], kky[iy], kkz[iz], kk2[ixyz]);
            }
            
             
            
            // solution 2
            ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            vk = delta*delta / ( pow(0.5*(eta_a+eta_b)-0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta) ;
            uk = 1.0 - vk;
            vk = sqrt(vk); uk=sqrt(uk);
            if(eta_b+ek<0.0) vk=-1.0*vk;   
            
            takeit=0;
            if(md.spinsymmetry>0)
            {
                if(ek>0.0) 
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                nwf++;
                if(nwf>=idxfrom && nwf<idxto) takeit=1;
            }
            
            if(takeit) // this is my wave-function
            {
                // u-components
                shift=NXYZ*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ ) for ( iz2 = 0 ; iz2 < NZ ; iz2++ ) 
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] + ( double ) iz2 * DZ * kkz[iz] ) ) * uk /  sqrt_volume;
                    ixyz2++;
                }
                
                // v-components
                shift=NXYZ*(idxto-idxfrom) + NXYZ*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ ) for ( iz2 = 0 ; iz2 < NZ ; iz2++ ) 
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] + ( double ) iz2 * DZ * kkz[iz] ) ) * vk /  sqrt_volume;
                    ixyz2++;
                }
                
                // weight
                fEn[nwf-idxfrom]=fbeta(ek,beta);
                
                // eigen energy
                En[nwf-idxfrom]=ek;
                
//                 // only for tests:
//                 if(nwf==idxfrom)
//                     wprintf("!!!!!!!!! nwf=%d: kx=%f ky=%f kz=%f kk2=%f\n", nwf, kkx[ix], kky[iy], kkz[iz], kk2[ixyz]);
            }
            
        }
        ixyz++;
    }    
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    
    return 0;
}

/**
 * Function computes eigen-values for given kz
 * @param kz value of kz
 * @param En eigen energies, E_n, array of size 2*NX*NY*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int extract_En_for_kz(double kz, double *En, int printout)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    int ix2, iy2, iz2, ixyz2; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    } 
        
    // Create En
    double V_a=__md_pca_uniform.V_a;
    double V_b=__md_pca_uniform.V_b;
    double eta_a, eta_b;
    double alph_a=__md_pca_uniform.alph_a;
    double alph_b=__md_pca_uniform.alph_b;
    double *mu_a, *mu_b;
    mu_a=&__md_pca_uniform.mu_a;
    mu_b=&__md_pca_uniform.mu_b;
    double ek;
    double delta=__md_pca_uniform.delta;
    double sqrt_volume=sqrt((double) (LXYZ));
    double *ec;
    ec=&__md_pca_uniform.ec;
    double beta=__md_pca_uniform.beta;

    if(printout) { wprintf("# extract_En_for_kz: extracting eigen-values\n"); fflush(stdout); }
    for(ixyz=0; ixyz<2*NX*NY; ixyz++) En[ixyz]=-999999999999.;
    ixyz=0;
    int nwf=0;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kkz[iz]==kz) // take only for given kz
        {
            eta_a = alph_a*kk2[ixyz]/2.0 + V_a - *mu_a;
            eta_b = alph_b*kk2[ixyz]/2.0 + V_b - *mu_b;    
            
            // solution 1
            ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            En[nwf]=ek;            
            nwf++; 
            
            // solution 2
            ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            En[nwf]=ek;                
            nwf++; 
        }
        ixyz++;
    }    
    
    if(printout) {wprintf("# extract_En_for_kz: sorting\n"); fflush(stdout); }
    
    while(1)
    {
        i = 1;
        for(ixyz=0; ixyz<(2*NX*NY-1); ixyz++) if(En[ixyz]>En[ixyz+1])
        {
            ek=En[ixyz];
            En[ixyz]=En[ixyz+1];
            En[ixyz+1]=ek;
            i=0;
        }
        
        if(i) break;
    }
    
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    
    return 0;
}

/**
 * Function computes eigen-values for uniform system - function used for testing
 * @param En eigen energies, E_n, array of size 2*NXYZ*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int extract_En(double *En, int printout)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    int ix2, iy2, iz2, ixyz2; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    } 
        
    // Create En
    double V_a=__md_pca_uniform.V_a;
    double V_b=__md_pca_uniform.V_b;
    double eta_a, eta_b;
    double alph_a=__md_pca_uniform.alph_a;
    double alph_b=__md_pca_uniform.alph_b;
    double *mu_a, *mu_b;
    mu_a=&__md_pca_uniform.mu_a;
    mu_b=&__md_pca_uniform.mu_b;
    double ek;
    double delta=__md_pca_uniform.delta;
    double sqrt_volume=sqrt((double) (LXYZ));
    double *ec;
    ec=&__md_pca_uniform.ec;
    double beta=__md_pca_uniform.beta;

    if(printout) { wprintf("# extract_En: extracting eigen-values\n"); fflush(stdout); }
    for(ixyz=0; ixyz<2*NXYZ; ixyz++) En[ixyz]=-999999999999.;
    ixyz=0;
    int nwf=0;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        
        {
            eta_a = alph_a*kk2[ixyz]/2.0 + V_a - *mu_a;
            eta_b = alph_b*kk2[ixyz]/2.0 + V_b - *mu_b;    
            
            // solution 1
            ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            En[nwf]=ek;            
            nwf++; 
            
            // solution 2
            ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            En[nwf]=ek;                
            nwf++; 
        }
        ixyz++;
    }    
    
    if(printout) {wprintf("# extract_En: sorting\n"); fflush(stdout); }
    
    while(1)
    {
        i = 1;
        for(ixyz=0; ixyz<(2*NXYZ-1); ixyz++) if(En[ixyz]>En[ixyz+1])
        {
            ek=En[ixyz];
            En[ixyz]=En[ixyz+1];
            En[ixyz+1]=ek;
            i=0;
        }
        
        if(i) break;
    }
    
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    
    return 0;
}

// -----------------------------------------------------------------
// ------------------------ BdG version ----------------------------
// -----------------------------------------------------------------
/**
 * @param n0_a requested density for population "a" (INPUT)
 * @param n0_b requested density for population "b" (INPUT)
 * @param nwf number of wave-functions (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int solve_uniform_problem_bdg(double n0_a, double n0_b, int *nwf, int printout)
{

    int i,j;
    int ix, iy, iz, ixyz; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    }  
    
    double * kk2tau ; /* kk2 = k^2 */
    cppmallocl(kk2tau,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        double _kkx=kkx[ix];
        double _kky=kky[iy];
        double _kkz=kkz[iz];
#ifdef TAU_COMPUTATION_VIA_GRADIENTS       
        if(ix==NX/2) _kkx=0.0;
        if(iy==NY/2) _kky=0.0;
        if(iz==NZ/2) _kkz=0.0;
#endif
        kk2tau[i]=_kkx*_kkx + _kky*_kky + _kkz*_kkz;
        i++;
    }
    
    
    // Set quantities that depend only on desities
    double p = polarization_h(n0_a, n0_b);

    double eF_a=pow(6.0*M_PI*M_PI*n0_a, 2.0/3.0) / 2.0;
    double eF_b=pow(6.0*M_PI*M_PI*n0_b, 2.0/3.0) / 2.0;
    double eF_avg=pow(3.0*M_PI*M_PI*(n0_a+n0_b), 2.0/3.0) / 2.0;
    double Effg = 0.6*n0_a*eF_a*LXYZ + 0.6*n0_b*eF_b*LXYZ;
    double gbare=4.0*M_PI*md.aBdG;
    double kc=md.kc;
    double alph_a=1.0;
    double alph_b=1.0;
    double mu_a=0.37*eF_a;
    double mu_b=0.37*eF_b;
    
    if(printout && md.init0debug>0) wprintf("# DEBUG: n_a=%f, n_b=%f\n", n0_a, n0_b);
    if(printout && md.init0debug>0) wprintf("# DEBUG: eF_a=%f, eF_b=%f, eF_avg=%f\n", eF_a, eF_b, eF_avg);
    if(printout && md.init0debug>0) wprintf("# DEBUG: N_a=%f, N_b=%f\n", n0_a*LXYZ, n0_b*LXYZ);
    if(printout && md.init0debug>0) wprintf("# DEBUG: p=%f\n", p);
    
    // Set quantities updated in s-c loop
    double tau_a=pow(6.0*M_PI*M_PI*n0_a, 5.0/3.0) / (10.0*M_PI*M_PI); // initial value
    double tau_b=pow(6.0*M_PI*M_PI*n0_b, 5.0/3.0) / (10.0*M_PI*M_PI); // initial value
    double delta;
    if(n0_a<n0_b) delta = 0.5*pow(6.0*M_PI*M_PI*n0_a, 2.0/3.0) / 2.0; // initial value: 0.5*eF
    else          delta = 0.5*pow(6.0*M_PI*M_PI*n0_b, 2.0/3.0) / 2.0; // initial value: 0.5*eF
    double nu= -1.0 * delta / gbare; // initial value
    
    if(printout && md.init0debug>0) wprintf("# DEBUG: tau_a=%f, tau_b=%f\n", tau_a, tau_b);   
    if(printout && md.init0debug>0) wprintf("# DEBUG: delta=%f, nu=%f\n", delta, nu);  
    if(printout && md.init0debug>0 && md.spinsymmetry>0) wprintf("# SPIN SYMMETRY MODE!\n");
    
    // auxliary variables
    int maxiter=md.init0maxiter;
    int iter;
    double tau_m, tau_p; 
    double V_a, V_b, mu_p, g_eff, eta_a, eta_b;
    double complex p0, wz_0;
    double complex Zzero = 0.0 + I*0.0;
    double complex Zone  = 1.0 + I*0.0;
    double uk, vk, ek;
    double n_a, n_b;
    double tau_a_old, tau_b_old, delta_old, nu_old;
    double scmix=md.init0scmix, epsilon=md.init0eps;
    int is_conv;
    double beta, T;
    double energy_kin, energy_pot, energy_pair, energy_tot;
    double kc2=kc*kc;
    
    // extract number of wave-functions
    ixyz=0;
    *nwf=0;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2) // take only from sphere
            *nwf+=2; // thera are two soloutions for each momentum
        ixyz++;
    }
    
    for(T=md.init0Tstart; T>=md.init0Tstop; T-=md.init0DeltaT) 
    {
        if(T<0.0) break;
        if(T<=0.0) beta=1.0e16;
        else beta=1.0/(T*eF_avg);
        for(iter=0; iter<maxiter; iter++)
        {
            // save old values of potentials
            tau_a_old=tau_a;
            tau_b_old=tau_b;
            delta_old=delta;
            nu_old=nu;
            
            // potential
            tau_p=tau_a + tau_b;
            tau_m=tau_a - tau_b;
            
            V_a = 0.0;
            V_b = 0.0;
            
            // pairing
#ifdef USE_CUBIC_CUTOFF
            wz_0=REGULARIZATION_SCHEME_K_CONST/(4.0*M_PI*DX) + I*0.0;
#else
            mu_p=(mu_a-V_a+mu_b-V_b)/2.0;
            p0 = csqrt( 2.0*mu_p) ;
            if ( cimag(p0) < 0. ) p0 *= -1. ;
            
            // kc is fixed, and finally will be trananslated into ec
            wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( cimag(wz_0) < 0. ) wz_0 += I * 2. * M_PI ;    
                    
            wz_0= kc / ( 2. * M_PI * M_PI ) *( 1. - p0 / ( 2. * kc ) * wz_0);
#endif
            g_eff = creal( Zone / (Zone/gbare - wz_0) );
            delta = -1.0*g_eff*nu;
//             wprintf("AAA: %f %f \n", delta, g_eff);
            
            // contribution from states to densities
            n_a=0.0;
            n_b=0.0;
            tau_a=0.0;
            tau_b=0.0;
            nu=0.0;
            ixyz=0;
            *nwf=0;
            for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
            {
                if(kk2[ixyz]<kc2) // take only from sphere
                {
                    eta_a = alph_a*kk2[ixyz]/2.0 + V_a - mu_a;
                    eta_b = alph_b*kk2[ixyz]/2.0 + V_b - mu_b;
                    
                    // solution 1:
                    ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
//                     if(iter>=1947) wprintf("DDD: %f %f %f %f\n", ek, eta_a, eta_b,delta);
                    vk = delta*delta / ( pow(0.5*(eta_a+eta_b)+0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta);
                    uk = 1.0 - vk;
//                     if(printout) if(vk<0.0) wprintf("S1v: problem\n");
//                     if(printout) if(uk<0.0) wprintf("S1u: problem\n");
                    vk = sqrt(vk); uk=sqrt(uk);
                    if(eta_b+ek<0.0) vk=-1.0*vk;
                    
                    if(md.spinsymmetry>0)
                    {
                        if(ek>0.0) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b 
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        n_a+=uk*uk*fbeta(ek,beta);
                        tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                        n_b+=vk*vk*fbeta(-1.0*ek,beta);
                        tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                        nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                        (*nwf)++;
                    }
                    
                    // solution 2:
                    ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
                    vk = delta*delta / ( pow(0.5*(eta_a+eta_b)-0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta) ;
                    uk = 1.0 - vk;
//                     if(printout) if(vk<0.0) wprintf("S2v: problem\n");
//                     if(printout) if(uk<0.0) wprintf("S2u: problem: %f %f\n", vk, uk);
                    vk = sqrt(vk); uk=sqrt(uk);
                    if(eta_b+ek<0.0) vk=-1.0*vk;
                    
                    if(md.spinsymmetry>0)
                    {
                        if(ek>0.0) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b 
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        n_a+=uk*uk*fbeta(ek,beta);
                        tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                        n_b+=vk*vk*fbeta(-1.0*ek,beta);
                        tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                        nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                        (*nwf)++;
                    }
                }
                
                ixyz++;
            }
            
            // add normalization factors
            if(md.spinsymmetry>0)
            {
                n_b/=LXYZ; n_a=n_b;
                tau_b/=LXYZ; tau_a=tau_b; 
                nu/=2.0*LXYZ;
            }
            else
            {
                n_a/=LXYZ; n_b/=LXYZ;
                tau_a/=LXYZ; tau_b/=LXYZ;
                nu/=2.0*LXYZ;
            }
            
            if(printout && md.init0debug>1) wprintf("D: iter=%d: V_a=%f, V_b=%f, delta=%f, n_a=%f, n_b=%f, tau_a=%f, tau_b=%f, nu=%f\n", iter, V_a, V_b, delta, n_a, n_b, tau_a, tau_b, nu);
            
            // check convergence
            if(printout && md.init0debug>1) wprintf("C: iter=%d: fabs(n0_a-n_a)=%g fabs(n0_b-n_b)=%g fabs(delta-delta_old)=%g, mu_a=%f, mu_b=%f\n", iter, fabs(n0_a-n_a), fabs(n0_b-n_b), fabs(delta-delta_old), mu_a, mu_b);
            is_conv=1;
            if(fabs(n0_a-n_a)>epsilon) is_conv=0; // check for density
            if(fabs(n0_b-n_b)>epsilon) is_conv=0; // check for density
            if(fabs(delta-delta_old)>epsilon) is_conv=0; // check for delta
            if(is_conv) break; // we converged
            
            // mix potentials and go to next iteration
            tau_a = (1.-scmix)*tau_a_old + scmix*tau_a;
            tau_b = (1.-scmix)*tau_b_old + scmix*tau_b;
            delta = (1.-scmix)*delta_old + scmix*delta;
            nu = (1.-scmix)*nu_old + scmix*nu;
            mu_a += md.init0muchange*(n0_a-n_a);
            mu_b += md.init0muchange*(n0_b-n_b);
            if(md.spinsymmetry>0) mu_b=mu_a;
        }
        if(printout && md.init0debug>0) wprintf("# TEMPCONV: T=%f, iter=%d, delta/eF_a=%f, mu_a/eF_a=%f, delta/eF_b=%f, mu_b/eF_b=%f\n", T, iter, delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b);
        if(iter==maxiter && printout && md.init0debug>0) wprintf("# WARNING: MAXITER REACHED!\n");
            
        // Compute energy 
        energy_kin=(0.5*alph_a*tau_a + 0.5*alph_b*tau_b)*LXYZ;
        energy_pot=(0.0)*LXYZ;
        energy_pair=-1.0*delta*nu*LXYZ;
        energy_tot=energy_kin+energy_pot+energy_pair;
        if(printout  && md.init0debug>0) wprintf("# TEMPCONV: T=%f, energy_kin=%f, energy_pot=%f, energy_pair=%f, energy_tot=%f\n", T, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
        fflush(stdout);
    }
    
    
    double ec=0.0;
    eta_a = alph_a*kc2/2.0 + V_a - mu_a;
    eta_b = alph_b*kc2/2.0 + V_b - mu_b;
    ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC1=%f\n", ek);
    if(fabs(ek)>ec) ec=fabs(ek);
    ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC2=%f\n", ek);
    if(fabs(ek)>ec) ec=fabs(ek);
    
    mu_p=(mu_a-V_a+mu_b-V_b)/2.0;
    ec=1.0*kc2/2.0 - mu_p;
    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC=%f\n", ec);
    
    // print results
    if(printout) wprintf("# UNIFORM SOLUTION: delta/eF_a=%8.4f, mu_a/eF_a=%8.4f, delta/eF_b=%8.4f, mu_b/eF_b=%8.4f, ec=%8.4f\n", delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b, ec);
    if(printout) wprintf("# UNIFORM SOLUTION: energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
    if(printout) wprintf("# UNIFORM SOLUTION: nwf=%d\n", *nwf);
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    free(kk2tau);
    
    // write data to structure for future use
    __md_pca_uniform.n0_a=n_a;
    __md_pca_uniform.n0_b=n_b;
    __md_pca_uniform.mu_a=mu_a;
    __md_pca_uniform.mu_b=mu_b;
    __md_pca_uniform.V_a=V_a;
    __md_pca_uniform.V_b=V_b;
    __md_pca_uniform.alph_a=alph_a;
    __md_pca_uniform.alph_b=alph_b;
    __md_pca_uniform.kc=kc;
    __md_pca_uniform.ec=ec;
    __md_pca_uniform.nu=nu;
    __md_pca_uniform.delta=delta;
    __md_pca_uniform.tau_a=tau_a;
    __md_pca_uniform.tau_b=tau_b;
    __md_pca_uniform.beta=beta;    
    __md_pca_uniform.nwf=*nwf;
    __md_pca_uniform.ekin=energy_kin;
    __md_pca_uniform.epot=energy_pot;
    __md_pca_uniform.epair=energy_pair;
    
//     // Append results to the file
//     FILE * ff = fopen("bdg.txt", "a");
//     fprintf(ff, "%12.6f %12.6f %12.6f %16.8f %16.8f %16.8f\n", eF_a, sqrt(2.0*eF_a), md.aBdG, delta/eF_a, mu_a/eF_a, energy_tot/Effg);
//     fclose(ff);
    
    if(iter==maxiter) return 1; // not converged!
    return 0;
}

// auxliary functions for cpca code

/**
 * Function resturns number of wave-functions that need to be evolved
 * It takes into accout that hamilonian for values kz and -kz is the same. 
 * @param nwf number of wave-functions (OUTPUT)
 * @return 0 - OK, otherwise error
 * */
int get_nwf_to_evolve_2d(int *nwf)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    }   
    
    // extract number of wave-functions
    ixyz=0;
    *nwf=0;

    double kc=md.kc;
    double kc2=kc*kc;
    int total_nwf=0;
    int total_nwf2=0;
    
    int *states_selected;
    cppmallocl(states_selected,NZ, int);
    states_selected[0]=1;
    for(i=1; i<NZ/2; i++) states_selected[i]=2;
#ifdef USE_CUBIC_CUTOFF
    states_selected[NZ/2]=1;
#else
    states_selected[NZ/2]=0;
#endif
    for(i=NZ/2+1; i<NZ; i++) states_selected[i]=0;
    
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2 && states_selected[iz]>0) // take only from sphere and for non-negative kz values 
            if(md.spinsymmetry==1) *nwf+=1; // only positive energy state
            else *nwf+=2; // thera are two soloutions for each momentum
            
        // make test for correctness
        if(md.spinsymmetry==1)
        {
            if(kk2[ixyz]<kc2) 
            {
                // take only positive energy states
                if(states_selected[iz]==2) total_nwf += 1*2; // two solutions x (-kz, +kz)
                else if(states_selected[iz]==1) total_nwf += 1;
                
                total_nwf2+=1;
            }
        }
        else
        {
            if(kk2[ixyz]<kc2)
            {
                if(states_selected[iz]==2) total_nwf += 2*2; // two solutions x (-kz, +kz)
                else if(states_selected[iz]==1) total_nwf += 2;
                
                total_nwf2+=2;
            }
        }
             
//         wprintf("TTT: %9d %4d %4d %4d %9d %9d\n", ixyz, ix, iy, iz, total_nwf, nwf);
        
        ixyz++;
    }
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    free(states_selected);
    
    if(total_nwf!=__md_pca_uniform.nwf)
    {
        wprintf("# ERROR: get_nwf_to_evolve_2d: total_nwf[%d]!=__md_pca_uniform.nwf[%d], total_nwf2[%d]\n", total_nwf, __md_pca_uniform.nwf, total_nwf2);
        fflush(stdout);
        return 1;
    }
    return 0;
}

/**
 * @param idxfrom extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param idxto extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param wf pointer to array for wf of size nx*ny*nz*2*(idxto-idxfrom)*sizeof(double complex)
 * @param mu_a chemical potential for population "a" (OUTPUT)
 * @param mu_b chemical potential for population "b" (OUTPUT)
 * @param ec energy cut-off for the solution (OUTPUT)
 * @param fEn weights used for computation densities, fEn=fbeta(E_n), array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param kkzvals values of corresponding kkz values, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param En eigen energies, E_n, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int create_uniform_wf_2d(int idxfrom, int idxto, double complex *wf, double *mu_a, double *mu_b, double *ec, double *fEn, double *kkzvals, double *En, int printout)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    int ix2, iy2, iz2, ixyz2; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    } 
    
    // extract number of wave-functions and check consistency
    double kc=md.kc;
    double kc2=kc*kc;
    ixyz=0;
    int nwf=0;
    int total_nwf=0;
    int takeit;
    
    int *states_selected;
    cppmallocl(states_selected,NZ, int);
    states_selected[0]=1;
    for(i=1; i<NZ/2; i++) states_selected[i]=2;
#ifdef USE_CUBIC_CUTOFF
    states_selected[NZ/2]=1;
#else
    states_selected[NZ/2]=0;
#endif
    for(i=NZ/2+1; i<NZ; i++) states_selected[i]=0;
    
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2 && states_selected[iz]>0) // only positive energy state
            if(md.spinsymmetry==1) nwf+=1; // thera are two soloutions for each momentum
            else nwf+=2; // thera are two soloutions for each momentum
            
        // make test for correctness
        if(md.spinsymmetry==1)
        {
            if(kk2[ixyz]<kc2)
            {
                if(states_selected[iz]==2) total_nwf += 1*2; // two solutions x (-kz, +kz)
                else if(states_selected[iz]==1) total_nwf += 1;
            }
        }
        else 
        {
            if(kk2[ixyz]<kc2)
            {
                if(states_selected[iz]==2) total_nwf += 2*2; // two solutions x (-kz, +kz)
                else if(states_selected[iz]==1) total_nwf += 2;
            }
        }
        
        ixyz++;
    }
    if(total_nwf!=__md_pca_uniform.nwf) return 1;
    if(idxfrom>=nwf) return 2;
    if(idxto>nwf) return 3;
    if(idxfrom>=idxto)return 4;
    
    // Create wf
    size_t shift;
    double V_a=__md_pca_uniform.V_a;
    double V_b=__md_pca_uniform.V_b;
    double eta_a, eta_b;
    double alph_a=__md_pca_uniform.alph_a;
    double alph_b=__md_pca_uniform.alph_b;
    *mu_a=__md_pca_uniform.mu_a;
    *mu_b=__md_pca_uniform.mu_b;
    double uk, vk, ek;
    double delta=__md_pca_uniform.delta;
    double sqrt_volume=sqrt((double) (LXY));
    *ec=__md_pca_uniform.ec;
    double beta=__md_pca_uniform.beta;

    if(printout) wprintf("# UNIFORM CREATE WF: Creating wave-functions.\n");
    ixyz=0;
    nwf=-1;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2 && states_selected[iz]>0) // take only from sphere
        {
            eta_a = alph_a*kk2[ixyz]/2.0 + V_a - *mu_a;
            eta_b = alph_b*kk2[ixyz]/2.0 + V_b - *mu_b;    
            
            // solution 1
            ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            vk = delta*delta / ( pow(0.5*(eta_a+eta_b)+0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta);
            uk = 1.0 - vk;
            vk = sqrt(vk); uk=sqrt(uk);
            if(eta_b+ek<0.0) vk=-1.0*vk;  
            
            takeit=0;
            if(md.spinsymmetry>0)
            {
                if(ek>0.0) 
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                nwf++;
                if(nwf>=idxfrom && nwf<idxto) takeit=1;
            }
            
            if(takeit) // this is my wave-function
            {
                // u-components
                shift=NXY*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ )
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] ) ) * uk /  sqrt_volume;
                    ixyz2++;
                }
                
                // v-components
                shift=NXY*(idxto-idxfrom) + NXY*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ )
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] ) ) * vk /  sqrt_volume;
                    ixyz2++;
                }
                
                // weight
                fEn[nwf-idxfrom]=fbeta(ek,beta);
                
                // kkz value
                kkzvals[nwf-idxfrom]=kkz[iz];
                
                // eigen energy
                En[nwf-idxfrom]=ek;
                
//                 // only for tests:
//                 if(nwf==idxfrom)
//                     wprintf("!!!!!!!!! nwf=%d: kx=%f ky=%f kz=%f kk2=%f\n", nwf, kkx[ix], kky[iy], kkz[iz], kk2[ixyz]);
            }
            
            
            
            // solution 2
            ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            vk = delta*delta / ( pow(0.5*(eta_a+eta_b)-0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta) ;
            uk = 1.0 - vk;
            vk = sqrt(vk); uk=sqrt(uk);
            if(eta_b+ek<0.0) vk=-1.0*vk;   
            
            takeit=0;
            if(md.spinsymmetry>0)
            {
                if(ek>0.0) 
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                nwf++;
                if(nwf>=idxfrom && nwf<idxto) takeit=1;
            }
            
            if(takeit) // this is my wave-function
            {
                // u-components
                shift=NXY*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ )
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] ) ) * uk /  sqrt_volume;
                    ixyz2++;
                }
                
                // v-components
                shift=NXY*(idxto-idxfrom) + NXY*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) for ( iy2 = 0 ; iy2 < NY ; iy2++ ) 
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] + ( double ) iy2 * DY * kky[iy] ) ) * vk /  sqrt_volume;
                    ixyz2++;
                }
                
                // weight
                fEn[nwf-idxfrom]=fbeta(ek,beta);
                
                // kkz value
                kkzvals[nwf-idxfrom]=kkz[iz];
                
                // eigen energy
                En[nwf-idxfrom]=ek;
                
//                 // only for tests:
//                 if(nwf==idxfrom)
//                     wprintf("!!!!!!!!! nwf=%d: kx=%f ky=%f kz=%f kk2=%f\n", nwf, kkx[ix], kky[iy], kkz[iz], kk2[ixyz]);
            }
            
        }
        ixyz++;
    }    
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    free(states_selected);
    
    return 0;
}

/**
 * Function resturns number of wave-functions that need to be evolved
 * It takes into accout that hamilonian for values kz and -kz and ky and -ky is the same. 
 * @param nwf number of wave-functions (OUTPUT)
 * @return 0 - OK, otherwise error
 * */
int get_nwf_to_evolve_1d(int *nwf)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    }   
    
    // extract number of wave-functions
    ixyz=0;
    *nwf=0;

    double kc=md.kc;
    double kc2=kc*kc;
    int total_nwf=0;
    int total_nwf2=0;
    int deg;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
            
        if(kk2[ixyz]<kc2) // 
        {
            deg=wslda_kmodes_1d_get_weight(kky[iy], kkz[iz]);
            if(deg>0) deg=1; else deg=0;
            if(md.spinsymmetry==1) *nwf+=1*deg; // only positive energy state
            else *nwf+=2*deg; // there are two solutions for each momentum
        }
            
        // make test for correctness
        if(md.spinsymmetry==1)
        {
            if(kk2[ixyz]<kc2) 
            {
                // take only positive energy states
                total_nwf += 1*wslda_kmodes_1d_get_weight(kky[iy], kkz[iz]);
                total_nwf2+= 1;
            }
        }
        else
        {
            if(kk2[ixyz]<kc2)
            {
                // there are two solutions for each momentum
                total_nwf += 2*wslda_kmodes_1d_get_weight(kky[iy], kkz[iz]);
                total_nwf2+= 2;
            }
        }
             
//         wprintf("TTT: %9d %4d %4d %4d %9d %9d\n", ixyz, ix, iy, iz, total_nwf, nwf);
        
        ixyz++;
    }
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    
    if(total_nwf!=__md_pca_uniform.nwf)
    {
        wprintf("# ERROR: get_nwf_to_evolve_1d: total_nwf[%d]!=__md_pca_uniform.nwf[%d], total_nwf2[%d]\n", total_nwf, __md_pca_uniform.nwf, total_nwf2);
        fflush(stdout);
        return 1;
    }
    return 0;
}

/**
 * @param idxfrom extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param idxto extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param wf pointer to array for wf of size nx*ny*nz*2*(idxto-idxfrom)*sizeof(double complex)
 * @param mu_a chemical potential for population "a" (OUTPUT)
 * @param mu_b chemical potential for population "b" (OUTPUT)
 * @param ec energy cut-off for the solution (OUTPUT)
 * @param fEn weights used for computation densities, fEn=fbeta(E_n), array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param kkzvals values of corresponding kkz values, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param cnt degenerecies of states (OUTPUT)
 * @param En eigen energies, E_n, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int create_uniform_wf_1d(int idxfrom, int idxto, double complex *wf, double *mu_a, double *mu_b, double *ec, double *fEn, double *kkyzvals, int *cnt, double *En, int printout)
{
    int i,j;
    int ix, iy, iz, ixyz; 
    int ix2, iy2, iz2, ixyz2; 
    
    // Allocate memory - Double values
    // momentum on the lattice
    double * kkx , * kky , * kkz ; /* values */
    cppmallocl(kkx,NX,double);
    cppmallocl(kky,NY,double);
    cppmallocl(kkz,NZ,double);
    
    /* NOTE : nx , ny , nz = 2j forall j integers (e.g. even numbers for the lattice dimensions) */
    // Initialize lattice in momentum space (first Brullion zone)
    /* initialize the k-space lattice */    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }

    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }    
    
    double * kk2 ; /* kk2 = k^2/2m */
    cppmallocl(kk2,NXYZ,double);
    i=0;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        kk2[i]=(kkx[ix]*kkx[ix] + kky[iy]*kky[iy] + kkz[iz]*kkz[iz]);
        i++;
    } 
    
    // extract number of wave-functions and check consistency
    double kc=md.kc;
    double kc2=kc*kc;
    ixyz=0;
    int nwf=0;
    int total_nwf=0;
    int takeit;
    int deg; 
    
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        if(kk2[ixyz]<kc2) // take only from sphere and for non-negative ky and non-negative kz values 
        {
            deg=wslda_kmodes_1d_get_weight(kky[iy], kkz[iz]);
            if(deg>0) deg=1; else deg=0;
            if(md.spinsymmetry==1) nwf+=1*deg; // only positive energy state
            else nwf+=2*deg; // there are two solutions for each momentum
        }
        
        // make test for correctness
        if(md.spinsymmetry==1)
        {
            if(kk2[ixyz]<kc2)
            {
                // take only positive energy states
                total_nwf += 1*wslda_kmodes_1d_get_weight(kky[iy], kkz[iz]);
            }
        }
        else 
        {
            if(kk2[ixyz]<kc2)
            {
                // there are two solutions for each momentum
                total_nwf += 2*wslda_kmodes_1d_get_weight(kky[iy], kkz[iz]);
            }
        }
        
        ixyz++;
    }
    if(total_nwf!=__md_pca_uniform.nwf) return 1;
    if(idxfrom>=nwf) return 2;
    if(idxto>nwf) return 3;
    if(idxfrom>=idxto)return 4;
    
    // Create wf
    size_t shift;
    double V_a=__md_pca_uniform.V_a;
    double V_b=__md_pca_uniform.V_b;
    double eta_a, eta_b;
    double alph_a=__md_pca_uniform.alph_a;
    double alph_b=__md_pca_uniform.alph_b;
    *mu_a=__md_pca_uniform.mu_a;
    *mu_b=__md_pca_uniform.mu_b;
    double uk, vk, ek;
    double delta=__md_pca_uniform.delta;
    double sqrt_volume=sqrt((double) (LX));
    *ec=__md_pca_uniform.ec;
    double beta=__md_pca_uniform.beta;

    if(printout) wprintf("# UNIFORM CREATE WF: Creating wave-functions.\n");
    ixyz=0;
    nwf=-1;
    for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ ) 
    {
        deg=wslda_kmodes_1d_get_weight(kky[iy], kkz[iz]);
        if(kk2[ixyz]<kc2 && deg>0) // take only from sphere
        {
            eta_a = alph_a*kk2[ixyz]/2.0 + V_a - *mu_a;
            eta_b = alph_b*kk2[ixyz]/2.0 + V_b - *mu_b;    
            
            // solution 1
            ek = -0.5*(eta_b-eta_a) + 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            vk = delta*delta / ( pow(0.5*(eta_a+eta_b)+0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta);
            uk = 1.0 - vk;
            vk = sqrt(vk); uk=sqrt(uk);
            if(eta_b+ek<0.0) vk=-1.0*vk;  
            
            takeit=0;
            if(md.spinsymmetry>0)
            {
                if(ek>0.0) 
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                nwf++;
                if(nwf>=idxfrom && nwf<idxto) takeit=1;
            }
            
            if(takeit) // this is my wave-function
            {
                // u-components
                shift=NX*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ )
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] ) ) * uk /  sqrt_volume;
                    ixyz2++;
                }
                
                // v-components
                shift=NX*(idxto-idxfrom) + NX*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ )
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] ) ) * vk /  sqrt_volume;
                    ixyz2++;
                }
                
                // weight
                fEn[nwf-idxfrom]=fbeta(ek,beta);
                
                // kkz value
                kkyzvals[                  nwf-idxfrom]=kky[iy];
                kkyzvals[(idxto-idxfrom) + nwf-idxfrom]=kkz[iz];
                
                // degeneracy
                cnt[nwf-idxfrom]=deg;
                
                // eigen energy
                En[nwf-idxfrom]=ek;
                
//                 // only for tests:
//                 if(nwf==idxfrom)
//                     wprintf("!!!!!!!!! nwf=%d: kx=%f ky=%f kz=%f kk2=%f\n", nwf, kkx[ix], kky[iy], kkz[iz], kk2[ixyz]);
            }
            
            
            
            // solution 2
            ek = -0.5*(eta_b-eta_a) - 0.5*sqrt( (eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta));
            vk = delta*delta / ( pow(0.5*(eta_a+eta_b)-0.5*sqrt((eta_b-eta_a)*(eta_b-eta_a) + 4.0*(eta_a*eta_b + delta*delta)),2) + delta*delta) ;
            uk = 1.0 - vk;
            vk = sqrt(vk); uk=sqrt(uk);
            if(eta_b+ek<0.0) vk=-1.0*vk;   
            
            takeit=0;
            if(md.spinsymmetry>0)
            {
                if(ek>0.0) 
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                nwf++;
                if(nwf>=idxfrom && nwf<idxto) takeit=1;
            }
            
            if(takeit) // this is my wave-function
            {
                // u-components
                shift=NX*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ )
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] ) ) * uk /  sqrt_volume;
                    ixyz2++;
                }
                
                // v-components
                shift=NX*(idxto-idxfrom) + NX*(nwf-idxfrom); // local index
                ixyz2=0;
                for ( ix2 = 0 ; ix2 < NX ; ix2++ ) 
                {
                    wf[shift+ixyz2]=cexp( I * ( ( double ) ix2 * DX * kkx[ix] ) ) * vk /  sqrt_volume;
                    ixyz2++;
                }
                
                // weight
                fEn[nwf-idxfrom]=fbeta(ek,beta);
                
                // kkz value
                kkyzvals[                  nwf-idxfrom]=kky[iy];
                kkyzvals[(idxto-idxfrom) + nwf-idxfrom]=kkz[iz];
                
                // degeneracy
                cnt[nwf-idxfrom]=deg;
                
                // eigen energy
                En[nwf-idxfrom]=ek;
                
//                 // only for tests:
//                 if(nwf==idxfrom)
//                     wprintf("!!!!!!!!! nwf=%d: kx=%f ky=%f kz=%f kk2=%f\n", nwf, kkx[ix], kky[iy], kkz[iz], kk2[ixyz]);
            }
            
        }
        ixyz++;
    }    
    
    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    
    return 0;
}

#endif
