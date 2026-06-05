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
    double S;
    double ekin;
    double epot;
    double epair;
    int nwf;
} metadata_pca_uniform_t;

// Stores results;
metadata_pca_uniform_t __md_pca_uniform;

void testsuite_file_uniform(double nerr, double eerr, int wmu, double muerr, int went, double errent)
{
    char fname [512];
    sprintf(fname,"%s_uniform.ref", md.outprefix);

    FILE *f=fopen(fname, "w");

    fprintf(f,"npart[SPINA]     %20.10g     %20.10g\n", __md_pca_uniform.n0_a*LXYZ, nerr);
    fprintf(f,"npart[SPINB]     %20.10g     %20.10g\n", __md_pca_uniform.n0_b*LXYZ, nerr);

    fprintf(f,"energy[EKIN]     %20.10g     %20.10g\n", __md_pca_uniform.ekin, eerr);
    fprintf(f,"energy[EPOT]     %20.10g     %20.10g\n", __md_pca_uniform.epot, eerr);
    fprintf(f,"energy[EPAIR]    %20.10g     %20.10g\n", __md_pca_uniform.epair, eerr);
    fprintf(f,"energy[ECURRENT] %20.10g     %20.10g\n", 0.0, eerr);
    fprintf(f,"energy[EPOTEXT]  %20.10g     %20.10g\n", 0.0, eerr);
    fprintf(f,"energy[EPAIREXT] %20.10g     %20.10g\n", 0.0, eerr);
    fprintf(f,"energy[EVELEXT]  %20.10g     %20.10g\n", 0.0, eerr);

    if(wmu)
    {
        fprintf(f,"mu[SPINA]        %20.10g     %20.10g\n", __md_pca_uniform.mu_a, muerr);
        fprintf(f,"mu[SPINB]        %20.10g     %20.10g\n", __md_pca_uniform.mu_b, muerr);
    }

    if(went)
    {
        fprintf(f,"entropy          %20.10g     %20.10g\n", __md_pca_uniform.S, errent);
    }

    fclose(f);
}

double fbeta(double E, double beta)
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
    // double Effg = 0.6*n0_a*eF_a*LXYZ + 0.6*n0_b*eF_b*LXYZ;
    double Effg = 0.6*(n0_a+n0_b)*eF_avg*LXYZ;
    double kc=md.kc;
    double mu_a=0.37*eF_a;
    double mu_b=0.37*eF_b;

    if(printout && md.init0debug>0) wprintf("# DEBUG: solve_uniform_problem\n");
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
    if(printout && md.init0debug>0) wprintf("# DEBUG: alph_a=%f, alph_b=%f, beta=%f, gamma=%f\n", alph_a, alph_b, D*LXYZ/Effg, GAMMA0);
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
    double V_a=0.0, V_b=0.0, mu_p, g_eff, eta_a, eta_b;
    double V_a_old, V_b_old;
    double complex p0, wz_0;
    double complex Zzero = 0.0 + I*0.0;
    double complex Zone  = 1.0 + I*0.0;
    double uk, vk, ek;
    double n_a, n_b;
    double tau_a_old, tau_b_old, delta_old, nu_old;
    double scmix=md.init0scmix, epsilon=md.init0eps;
    int is_conv;
    double beta, T, S;
    double energy_kin, energy_pot, energy_pair, energy_tot;
    double ec=md.ec;

    // iterate over temperatures range
    double md_init0Tstart=fabs(md.init0Tstart);
    double md_init0Tstop=fabs(md.init0Tstop);
    double md_init0DeltaT=fabs(md.init0DeltaT);
    double dsT=1.0, TT;
    if(md_init0Tstart<md_init0Tstop) dsT=-1.0; // it is trick to allow increasing or decreasing temparture
    for(TT=md_init0Tstart*dsT; TT>=md_init0Tstop*dsT; TT-=md_init0DeltaT)
    {
        T=TT*dsT; // it should be always positive
        if(T<0.0) break;
        if(T<=1.0e-16) beta=1.0e16;
        else beta=1.0/(T*eF_avg);
        for(iter=0; iter<maxiter; iter++)
        {
            // save old values of potentials
            tau_a_old=tau_a;
            tau_b_old=tau_b;
            delta_old=delta;
            nu_old=nu;
            V_a_old=V_a;
            V_b_old=V_b;

            // potential
            tau_p=tau_a + tau_b;
            tau_m=tau_a - tau_b;

            V_a = dalphm_dna*tau_m/2.0 + dalphp_dna*(tau_p/2.0 - delta*nu/alph_plus) - dtildeC_dna*delta*delta/alph_plus + dD_dna;
            V_b = dalphm_dnb*tau_m/2.0 + dalphp_dnb*(tau_p/2.0 - delta*nu/alph_plus) - dtildeC_dnb*delta*delta/alph_plus + dD_dnb;

            // pairing
            mu_p=(mu_a-V_a_old+mu_b-V_b_old)/2.0;
            p0 = csqrt( 2.0*mu_p/ alph_plus) ;
            if ( cimag(p0) < 0. ) p0 *= -1. ;

            //## kc is fixed, and it will be translated into ec
            #ifdef USE_CUBIC_CUTOFF
            // no change of ec, which is set to infinity
            kc = M_PI/DX;
            #else
            ec = alph_plus*kc*kc/2.0 - mu_p;
            #endif

            // ec is fixed, and it will be translated into kc
            // kc = sqrt( 2.0*(ec+mu_p)/ alph_plus) ;

            wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( cimag(wz_0) < 0. ) wz_0 += I * 2. * M_PI ;

            wz_0= kc *REG_COEFF_R0 *( 1. - ( p0 / kc ) * REG_COEFF_R1*wz_0);

            g_eff = creal( Zone*alph_plus / (Zone*tC - wz_0) );
            delta = -1.0*g_eff*nu;

            #ifdef USE_CUBIC_CUTOFF
            // correction to the mean-field due to regularization
            double kF=pow(3.0*M_PI*M_PI*(n0_a+n0_b), 1.0/3.0);
            double bcoeff=creal(p0)/(kF+1.0e-12); // to avoid numerical problems
            double x = bcoeff*kF * DX / M_PI;
            if(g_eff<-1.0e-10 && x>1.0e-10) // to avoid numerical problems
            {
                double Lam_0 = bcoeff*REG_COEFF_R0/x *(1.0-REG_COEFF_R1*x*log((1.0+x)/(1.0-x)));
                double dLam_0_dx = (-2.*bcoeff* REG_COEFF_R0*REG_COEFF_R1)/(1.-x*x) - bcoeff*REG_COEFF_R0/x/x; // derivative of Lam_0 with respect to x

                double Lam = Lam_0*kF/alph_plus; // regularizator in codes units, with A correction
                double dLam_dn = (kF/(3.*(n0_a+n0_b)*alph_plus))*(Lam_0+dLam_0_dx*bcoeff*DX*kF/M_PI);// derivative of Lam with respect to n
            
                V_a+=dLam_dn*pow(delta,2);
                V_b+=dLam_dn*pow(delta,2);
            }
            #endif

            // contribution from states to densities
            S=0.0;
            n_a=0.0;
            n_b=0.0;
            tau_a=0.0;
            tau_b=0.0;
            nu=0.0;
            ixyz=0;
            *nwf=0;
            for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ )
            {
                #ifndef USE_CUBIC_CUTOFF
                #if CODEDIM==2 || CODEDIM==1
                if(iz==NZ/2) {ixyz++; continue;}
                #endif
                #if CODEDIM==1
                if(iy==NY/2) {ixyz++; continue;}
                #endif
                #endif
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
                        if(ek>0.0 && ek<ec) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        if(ek>-ec && ek<ec)
                        {
                            n_a+=uk*uk*fbeta(ek,beta);
                            tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                            n_b+=vk*vk*fbeta(-1.0*ek,beta);
                            tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
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
                        if(ek>0.0 && ek<ec) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        if(ek>-ec && ek<ec)
                        {
                            n_a+=uk*uk*fbeta(ek,beta);
                            tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                            n_b+=vk*vk*fbeta(-1.0*ek,beta);
                            tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
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
//             mu_a += md.init0muchange*(md.init0Tstart/T)*(n0_a-n_a);
//             mu_b += md.init0muchange*(md.init0Tstart/T)*(n0_b-n_b);
            mu_a += md.init0muchange*(n0_a-n_a)/n0_a*eF_a;
            mu_b += md.init0muchange*(n0_b-n_b)/n0_b*eF_b;
            if(md.spinsymmetry>0) mu_b=mu_a;
        }
        eF_a=pow(6.0*M_PI*M_PI*n_a, 2.0/3.0) / 2.0;
        eF_b=pow(6.0*M_PI*M_PI*n_b, 2.0/3.0) / 2.0;
        eF_avg=pow(3.0*M_PI*M_PI*(n_a+n_b), 2.0/3.0) / 2.0;
        // Effg = 0.6*n_a*eF_a*LXYZ + 0.6*n_b*eF_b*LXYZ;
        Effg = 0.6*(n_a+n_b)*eF_avg*LXYZ;
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

    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC=%f\n", ec);

    // print results
    if(printout) wprintf("# UNIFORM SOLUTION: delta/eF_a=%8.4f, mu_a/eF_a=%8.4f, delta/eF_b=%8.4f, mu_b/eF_b=%8.4f, ec=%8.4f\n", delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b, ec);
    if(printout) wprintf("# UNIFORM SOLUTION: energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
    if(printout) wprintf("# ENTROPY PER PARTICLE: S/NkB=%16.12f\n", S/((n0_a+n0_a)*LXYZ));
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
    __md_pca_uniform.S=S;
    __md_pca_uniform.nwf=*nwf;
    __md_pca_uniform.ekin=energy_kin;
    __md_pca_uniform.epot=energy_pot;
    __md_pca_uniform.epair=energy_pair;

#ifdef TESTSUITE
#ifdef TDWSLDA
    if(printout) testsuite_file_uniform(TS_NERR, TS_EERR, 0, TS_MUERR, 0, TS_SERR);
#endif
#ifdef WSLDA
    if(printout) testsuite_file_uniform(TS_NERR, TS_EERR, 1, TS_MUERR, 1, TS_SERR);
#endif
#endif

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
    // double Effg = 0.6*n0_a*eF_a*LXYZ + 0.6*n0_b*eF_b*LXYZ;
    double Effg = 0.6*(n0_a+n0_b)*eF_avg*LXYZ;
    double kc=__md_pca_uniform.kc;
    double mu_a=__md_pca_uniform.mu_a;
    double mu_b=__md_pca_uniform.mu_b;
    double S=__md_pca_uniform.S;

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
    if(printout) wprintf("# ENTROPY PER PARTICLE: S/NkB=%16.12f\n", S/((n0_a+n0_a)*LXYZ));
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
 * @param fEn weights used for computation densities, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
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
    int nwf=__md_pca_uniform.nwf;
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
                if(ek>0.0 && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                if(ek>-(*ec) && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
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
                fEn[nwf-idxfrom]=ek;

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
                if(ek>0.0 && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                if(ek>-(*ec) && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
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
                fEn[nwf-idxfrom]=ek;

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
    // double Effg = 0.6*n0_a*eF_a*LXYZ + 0.6*n0_b*eF_b*LXYZ;
    double Effg = 0.6*(n0_a+n0_b)*eF_avg*LXYZ;
    if(md.akF!=0.0) md.sclgth=md.akF/sqrt(2.*eF_avg);
    double gbare=4.0*M_PI*md.sclgth;
    double kc=md.kc;
    double alph_a=1.0;
    double alph_b=1.0;
    double mu_a=0.96087287*eF_a; // value for akF=-1
    double mu_b=0.96087287*eF_b; // value for akF=-1

    if(printout && md.init0debug>0) wprintf("# DEBUG: solve_uniform_problem_bdg\n");
    if(printout && md.init0debug>0) wprintf("# DEBUG: n_a=%f, n_b=%f\n", n0_a, n0_b);
    if(printout && md.init0debug>0) wprintf("# DEBUG: eF_a=%f, eF_b=%f, eF_avg=%f\n", eF_a, eF_b, eF_avg);
    if(printout && md.init0debug>0) wprintf("# DEBUG: N_a=%f, N_b=%f\n", n0_a*LXYZ, n0_b*LXYZ);
    if(printout && md.init0debug>0) wprintf("# DEBUG: p=%f\n", p);
    if(printout && md.init0debug>0) wprintf("# DEBUG: akF=%f\n", md.sclgth*sqrt(2.*eF_avg));

    // Set quantities updated in s-c loop
    double tau_a=pow(6.0*M_PI*M_PI*n0_a, 5.0/3.0) / (10.0*M_PI*M_PI); // initial value
    double tau_b=pow(6.0*M_PI*M_PI*n0_b, 5.0/3.0) / (10.0*M_PI*M_PI); // initial value
    double delta;
    if(n0_a<n0_b) delta = 0.21*pow(6.0*M_PI*M_PI*n0_a, 2.0/3.0) / 2.0; // value for akF=-1
    else          delta = 0.21*pow(6.0*M_PI*M_PI*n0_b, 2.0/3.0) / 2.0; // value for akF=-1
    double nu= -1.0 * delta / gbare; // initial value

    if(printout && md.init0debug>0) wprintf("# DEBUG: tau_a=%f, tau_b=%f\n", tau_a, tau_b);
    if(printout && md.init0debug>0) wprintf("# DEBUG: delta=%f, nu=%f\n", delta, nu);
    if(printout && md.init0debug>0 && md.spinsymmetry>0) wprintf("# SPIN SYMMETRY MODE!\n");

    // auxliary variables
    int maxiter=md.init0maxiter;
    int iter;
    double tau_m, tau_p;
    double V_a=0.0, V_b=0.0, mu_p, g_eff, eta_a, eta_b;
    double V_a_old, V_b_old;
    double complex p0, wz_0;
    double complex Zzero = 0.0 + I*0.0;
    double complex Zone  = 1.0 + I*0.0;
    double uk, vk, ek;
    double n_a, n_b;
    double tau_a_old, tau_b_old, delta_old, nu_old;
    double scmix=md.init0scmix, epsilon=md.init0eps;
    int is_conv;
    double beta, T, S;
    double energy_kin, energy_pot, energy_pair, energy_tot;
    double ec=md.ec;

    // iterate over temperatures range
    double md_init0Tstart=fabs(md.init0Tstart);
    double md_init0Tstop=fabs(md.init0Tstop);
    double md_init0DeltaT=fabs(md.init0DeltaT);
    double dsT=1.0, TT;
    if(md_init0Tstart<md_init0Tstop) dsT=-1.0; // it is trick to allow increasing or decreasing temparture
    for(TT=md_init0Tstart*dsT; TT>=md_init0Tstop*dsT; TT-=md_init0DeltaT)
    {
        T=TT*dsT; // it should be always positive
        if(T<0.0) break;
        if(T<=1.0e-16) beta=1.0e16;
        else beta=1.0/(T*eF_avg);
        for(iter=0; iter<maxiter; iter++)
        {
            // save old values of potentials
            tau_a_old=tau_a;
            tau_b_old=tau_b;
            delta_old=delta;
            nu_old=nu;
            V_a_old=V_a; 
            V_b_old=V_b;

            // potential
            tau_p=tau_a + tau_b;
            tau_m=tau_a - tau_b;

            V_a = 0.0;
            V_b = 0.0;

            // pairing
            mu_p=(mu_a-V_a_old+mu_b-V_b_old)/2.0;
            p0 = csqrt( 2.0*mu_p) ;
            if ( cimag(p0) < 0. ) p0 *= -1. ;

            //## kc is fixed, and it will be translated into ec
            #ifdef USE_CUBIC_CUTOFF
            // no change of ec, which is set to infinity
            kc = M_PI/DX;
            #else
            ec = kc*kc/2.0 - mu_p;
            #endif

            wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( cimag(wz_0) < 0. ) wz_0 += I * 2. * M_PI ;

            wz_0= kc *REG_COEFF_R0 *( 1. - ( p0 / kc ) * REG_COEFF_R1*wz_0);

            g_eff = creal( Zone / (Zone/gbare - wz_0) );
            delta = -1.0*g_eff*nu;
            // printf("# T=%f iter=%d: V_a=%f V_b=%f delta=%f nu=%f g_eff=%f wz_0=(%f,%f)\n", T, iter, V_a, V_b, delta, nu, g_eff, creal(wz_0), cimag(wz_0) );

            #ifdef USE_CUBIC_CUTOFF
            // correction to the mean-field due to regularization
            double kF=pow(3.0*M_PI*M_PI*(n0_a+n0_b), 1.0/3.0);
            double bcoeff=creal(p0)/(kF+1.0e-12); // to avoid numerical problems
            double x = bcoeff*kF * DX / M_PI;
            if(g_eff<-1.0e-10 && x>1.0e-10) // to avoid numerical problems
            {
                double Lam_0 = bcoeff*REG_COEFF_R0/x *(1.0-REG_COEFF_R1*x*log((1.0+x)/(1.0-x)));
                double dLam_0_dx = (-2.*bcoeff* REG_COEFF_R0*REG_COEFF_R1)/(1.-x*x) - bcoeff*REG_COEFF_R0/x/x; // derivative of Lam_0 with respect to x

                double Lam = Lam_0*kF; // regularizator in codes units, with A correction
                double dLam_dn = (kF/(3.*(n0_a+n0_b)))*(Lam_0+dLam_0_dx*bcoeff*DX*kF/M_PI);// derivative of Lam with respect to n
            
                V_a+=dLam_dn*pow(delta,2);
                V_b+=dLam_dn*pow(delta,2);
            }
            #endif

            // contribution from states to densities
            S=0.0;
            n_a=0.0;
            n_b=0.0;
            tau_a=0.0;
            tau_b=0.0;
            nu=0.0;
            ixyz=0;
            *nwf=0;
            for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ )
            {
                #ifndef USE_CUBIC_CUTOFF
                #if CODEDIM==2 || CODEDIM==1
                if(iz==NZ/2) {ixyz++; continue;}
                #endif
                #if CODEDIM==1
                if(iy==NY/2) {ixyz++; continue;}
                #endif
                #endif
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
                        if(ek>0.0 && ek<ec) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        if(ek>-ec && ek<ec)
                        {
                            n_a+=uk*uk*fbeta(ek,beta);
                            tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                            n_b+=vk*vk*fbeta(-1.0*ek,beta);
                            tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
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
                        if(ek>0.0 && ek<ec) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        if(ek>-ec && ek<ec)
                        {
                            n_a+=uk*uk*fbeta(ek,beta);
                            tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                            n_b+=vk*vk*fbeta(-1.0*ek,beta);
                            tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
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
        eF_a=pow(6.0*M_PI*M_PI*n_a, 2.0/3.0) / 2.0;
        eF_b=pow(6.0*M_PI*M_PI*n_b, 2.0/3.0) / 2.0;
        eF_avg=pow(3.0*M_PI*M_PI*(n_a+n_b), 2.0/3.0) / 2.0;
        // Effg = 0.6*n_a*eF_a*LXYZ + 0.6*n_b*eF_b*LXYZ;
        Effg = 0.6*(n_a+n_b)*eF_avg*LXYZ;
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

    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC=%f\n", ec);

    // print results
    if(printout) wprintf("# UNIFORM SOLUTION: delta/eF_a=%8.4f, mu_a/eF_a=%8.4f, delta/eF_b=%8.4f, mu_b/eF_b=%8.4f, ec=%8.4f\n", delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b, ec);
    if(printout) wprintf("# UNIFORM SOLUTION: energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
    if(printout) wprintf("# ENTROPY PER PARTICLE: S/NkB=%16.12f\n", S/((n0_a+n0_a)*LXYZ));
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
    __md_pca_uniform.S=S;
    __md_pca_uniform.nwf=*nwf;
    __md_pca_uniform.ekin=energy_kin;
    __md_pca_uniform.epot=energy_pot;
    __md_pca_uniform.epair=energy_pair;

#ifdef TESTSUITE
#ifdef TDWSLDA
    if(printout) testsuite_file_uniform(TS_NERR, TS_EERR, 0, TS_MUERR, 0, TS_SERR);
#endif
#ifdef WSLDA
    if(printout) testsuite_file_uniform(TS_NERR, TS_EERR, 1, TS_MUERR, 1, TS_SERR);
#endif
#endif

    if(iter==maxiter) return 1; // not converged!
    return 0;
}



// -----------------------------------------------------------------
// ----------------------- SLDAE version ---------------------------
// -----------------------------------------------------------------
// -----------------------------------------------------------------
//
//     ######  ##       ########     ###    ########
//    ##    ## ##       ##     ##   ## ##   ##
//    ##       ##       ##     ##  ##   ##  ##
//     ######  ##       ##     ## ##     ## ######
//          ## ##       ##     ## ######### ##
//    ##    ## ##       ##     ## ##     ## ##
//     ######  ######## ########  ##     ## ########
//
// -----------------------------------------------------------------
/**
 * SLDAe variant
 * Author: Antoine Boulet
 * Update: 2021.09.21
 * */

#include "sldae_functional.h"

/**
 * @param n0_a requested density for population "a" (INPUT)
 * @param n0_b requested density for population "b" (INPUT)
 * @param nwf number of wave-functions (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int solve_uniform_problem_sldae(double n0_a, double n0_b, int *nwf, int printout)
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
    //##
    double as_, x_, kF_, eF_; // local Fermi momentum and Fermi energy
    double alpha_, beta_, inverse_gamma_; // HFB paremeters
    double alpha_p; // HFB parameter derivative
    double af_, bf_, cf_;    // functional parameters
    double af_p, bf_p, cf_p; // functional parameters (fderiv)
    double nt_, nt_1o3, nt_2o3; // power of total local density
    double dx_dnt_;  // derivative of x_ according to nt_ = x_ / (3.*nt_)
    double deF_dnt_; // derivative of eF_ according to nt_ = kF_ * kF_ / (3.*nt_)
    double ctilde_, ctilde_p;
        // ctilde_ ~ alpha_ * nt_1o3 * inverse_gamma_
        // but depends on regularization scheme
    double g_eff, inverse_gamma_eff; // renormalized pairing coupling constants
    double lambda_, lmu_sc; // spherical cutoff integral
    //##

    // register for total density
    nt_ = n0_a + n0_b;
    nt_1o3 = pow(nt_, 1. / 3.);
    nt_2o3 = pow(nt_, 2. / 3.);

    // register for local Fermi momentum and Fermi energy
    kF_ = pow(3. * M_PI_SQ * nt_, 1. / 3.);
    eF_ = pow(kF_, 2) / 2.;
    if(md.akF!=0.0) md.sclgth=md.akF/kF_;
    as_ = md.sclgth; // s-wave scattering length
    x_ = fabs(as_ * kF_); // density-dependent coupling constant
    dx_dnt_ = x_ / (3. * nt_);
    deF_dnt_ = pow(kF_, 2) / (3. * nt_);

    //##
    /*
      - functional derivative of quantity Z_
        according to the total density is noted Z_p
      - the renormalized coupling consants due to pairing
        are ended by _eff, e.g. g_eff
      (Note that in-medium renormalization procedure does not require cf_ and cf_p)
    */
    alpha_ = alpha_parameter_d0(x_);
    beta_ = beta_parameter_d0(x_);
    inverse_gamma_ = inverse_gamma_parameter_d0(x_);
    alpha_p = dx_dnt_ * alpha_parameter_d1(x_);
    af_ = a_functional_d0(x_);
    bf_ = b_functional_d0(x_);
    cf_ = c_functional_d0(x_);
    // definition independent of the functional and pairing form used
    af_p = alpha_p;
    bf_p = 5. / 3. * (beta_ - bf_) / nt_;
    cf_p = cf_ / (3. * nt_) * (1. - cf_ * inverse_gamma_);

    // initial values
    ctilde_ = af_ * nt_1o3 / cf_;
    inverse_gamma_eff = (1 / cf_) * (1. - 3. * nt_ / cf_ * cf_p);
    ctilde_p = inverse_gamma_eff * af_ / (3. * nt_2o3) + af_p * nt_1o3 / cf_;

    if(printout && md.init0debug>0) wprintf("# DEBUG SLDAE: as = %f, kF = %f, |askF| = %f\n", as_, kF_, x_);
    if(printout && md.init0debug>0) wprintf("# DEBUG SLDAE: xi = %f, zeta = %f, eta = %f\n", ground_state_energy_d0(x_), chemical_potential_d0(x_), pairing_gap_d0(x_));
    if(printout && md.init0debug>0) wprintf("# DEBUG SLDAE: A = %f, B = %f, C = %f\n", af_, bf_, cf_);
    if(printout && md.init0debug>0) wprintf("# DEBUG SLDAE: alpha = %f, beta = %f, gamma = %f\n", alpha_, beta_, 1./inverse_gamma_);

    // Set quantities that depend only on desities (spin-symmetry)
    double p = polarization_h(n0_a, n0_b); // = 0.0
    double alph_a = af_; //alpha_a_h(p);
    double alph_b = af_; //alpha_b_h(p);
    double dalphm_dna= 0.0; //der_alpha_minus__der_na_h(n0_a, n0_b);
    double dalphm_dnb= 0.0; //der_alpha_minus__der_nb_h(n0_a, n0_b);
    double dalphp_dna= af_p; //der_alpha_plus__der_na_h(n0_a, n0_b);
    double dalphp_dnb= af_p; //der_alpha_plus__der_nb_h(n0_a, n0_b);
    double alph_plus = af_; //alpha_plus_h(p);
    double dtildeC_dna = ctilde_p; //der_tildeC__der_na_h(n0_a, n0_b);
    double dtildeC_dnb = ctilde_p; //der_tildeC__der_nb_h(n0_a, n0_b);
    double dD_dna = (3. / 5.) * (bf_p * nt_ + bf_) * eF_ +
                    (3. / 5.) * bf_ * nt_ * deF_dnt_;
                    //der_funD__der_na_h(n0_a, n0_b);
    double dD_dnb = (3. / 5.) * (bf_p * nt_ + bf_) * eF_ +
                    (3. / 5.) * bf_ * nt_ * deF_dnt_;
                    //der_funD__der_nb_h(n0_a, n0_b);
    double tC = ctilde_; //tildeC_h(n0_a, n0_b);
    double D = (3. / 5.) * bf_ * nt_ * eF_; //funD_h(n0_a, n0_b);
    double eF_a=pow(6.0*M_PI*M_PI*n0_a, 2.0/3.0) / 2.0;
    double eF_b=pow(6.0*M_PI*M_PI*n0_b, 2.0/3.0) / 2.0;
    double eF_avg=pow(3.0*M_PI*M_PI*(n0_a+n0_b), 2.0/3.0) / 2.0;
    // double Effg = 0.6*n0_a*eF_a*LXYZ + 0.6*n0_b*eF_b*LXYZ;
    double Effg = 0.6*(n0_a+n0_b)*eF_avg*LXYZ;
    double kc=md.kc;
    double mu_a=chemical_potential_d0(x_)*eF_a;
    double mu_b=chemical_potential_d0(x_)*eF_b;

    if(printout && md.init0debug>0) wprintf("# DEBUG: solve_uniform_problem_sldae\n");
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
    double nu= -1.0 * delta * tC / alpha_; // initial value

    if(printout && md.init0debug>0) wprintf("# DEBUG: tau_a=%f, tau_b=%f\n", tau_a, tau_b);
    if(printout && md.init0debug>0) wprintf("# DEBUG: delta=%f, nu=%f\n", delta, nu);

    // auxliary variables
    int maxiter=md.init0maxiter;
    int iter;
    double V_a=0.0, V_b=0.0, eta_a, eta_b;
    double V_a_old, V_b_old;
    double complex p0, wz_0;
    double complex Zone  = 1.0 + I*0.0;
    double uk, vk, ek;
    double n_a, n_b;
    double tau_a_old, tau_b_old, delta_old, nu_old;
    double scmix=md.init0scmix, epsilon=md.init0eps;
    int is_conv;
    double beta, T, S;
    double energy_kin, energy_pot, energy_pair, energy_tot;
    double ec=md.ec;


    // iterate over temperatures range
    double md_init0Tstart=fabs(md.init0Tstart);
    double md_init0Tstop=fabs(md.init0Tstop);
    double md_init0DeltaT=fabs(md.init0DeltaT);
    double dsT=1.0, TT;
    if(md_init0Tstart<md_init0Tstop) dsT=-1.0; // it is trick to allow increasing or decreasing temparture
    for(TT=md_init0Tstart*dsT; TT>=md_init0Tstop*dsT; TT-=md_init0DeltaT)
    {
        T=TT*dsT; // it should be always positive
        if(T<0.0) break;
        if(T<=1.0e-16) beta=1.0e16;
        else beta=1.0/(T*eF_avg);
        for(iter=0; iter<maxiter; iter++)
        {

            // save old values of potentials
            tau_a_old=tau_a;
            tau_b_old=tau_b;
            delta_old=delta;
            nu_old=nu;
            V_a_old=V_a;
            V_b_old=V_b;

            // potential (kinetic contribution)
            V_a = af_p * (tau_a+tau_b) / 2.;
            V_b = af_p * (tau_a+tau_b) / 2.;
            // potential (mean-field contribution)
            V_a += beta_ * eF_;
            V_b += beta_ * eF_;
            // potential (pairing contribution)
            V_a += -ctilde_p*delta*delta/af_ - af_p*delta*nu/af_;
            V_b += -ctilde_p*delta*delta/af_ - af_p*delta*nu/af_;

            // effective pairing coupling constants and pairing field
            lmu_sc = (mu_a - V_a_old + mu_b - V_b_old) / 2.;
            // p0_ = sqrt (fabs (2. * (0. + lmu_sc) / af_));

            // if (lmu_sc >= 0.) {
            //   lambda_ = (kc + p0_) / (kc - p0_);
            //   lambda_ = 1. - p0_ / (2. * kc) * log(lambda_);
            //   lambda_ *= kc / (2. * M_PI_SQ);
            // } else {
            //   lambda_ = p0_ / kc;
            //   lambda_ = 1. + p0_ / kc * atan(lambda_);
            //   lambda_ *= kc / (2. * M_PI_SQ);
            // }

            // unified expression for spherical and cubic cut-off regularization schemes
            p0 = csqrt( 2.0*lmu_sc/ af_) ;
            if ( cimag(p0) < 0. ) p0 *= -1. ;

            //## kc is fixed, and it will be translated into ec
            #ifdef USE_CUBIC_CUTOFF
            // no change of ec, which is set to infinity
            kc = M_PI/DX;
            #else
            ec = af_*kc*kc/2.0 - lmu_sc;
            #endif

            wz_0 = clog( ( kc + p0 ) / ( kc - p0 ) ) ;
            if ( cimag(wz_0) < 0. ) wz_0 += I * 2. * M_PI ;

            wz_0= kc *REG_COEFF_R0 *( 1. - ( p0 / kc ) * REG_COEFF_R1*wz_0);
            lambda_=creal(wz_0);

            // effective coupling constant
            g_eff = af_ / (ctilde_ - lambda_);
            // pairing field
            delta = -nu * g_eff;

            #ifdef USE_CUBIC_CUTOFF
            // correction to the mean-field due to regularization
            double kF=pow(3.0*M_PI*M_PI*(n0_a+n0_b), 1.0/3.0);
            double bcoeff=creal(p0)/(kF+1.0e-12); // to avoid numerical problems
            double x = bcoeff*kF * DX / M_PI;
            if(g_eff<-1.0e-10 && x>1.0e-10) // to avoid numerical problems
            {
                double Lam_0 = bcoeff*REG_COEFF_R0/x *(1.0-REG_COEFF_R1*x*log((1.0+x)/(1.0-x)));
                double dLam_0_dx = (-2.*bcoeff* REG_COEFF_R0*REG_COEFF_R1)/(1.-x*x) - bcoeff*REG_COEFF_R0/x/x; // derivative of Lam_0 with respect to x

                double dLam_dn = (kF/(3.*(n0_a+n0_b)*af_))*(Lam_0+dLam_0_dx*bcoeff*DX*kF/M_PI);// derivative of Lam with respect to n
            
                V_a+=dLam_dn*pow(delta,2);
                V_b+=dLam_dn*pow(delta,2);
            }
            #endif

            // contribution from states to densities
            S=0.0;
            n_a=0.0;
            n_b=0.0;
            tau_a=0.0;
            tau_b=0.0;
            nu=0.0;
            ixyz=0;
            *nwf=0;
            for ( ix = 0 ; ix < NX ; ix++ ) for ( iy = 0 ; iy < NY ; iy++ ) for ( iz = 0 ; iz < NZ ; iz++ )
            {
                #ifndef USE_CUBIC_CUTOFF
                #if CODEDIM==2 || CODEDIM==1
                if(iz==NZ/2) {ixyz++; continue;}
                #endif
                #if CODEDIM==1
                if(iy==NY/2) {ixyz++; continue;}
                #endif
                #endif
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
                        if(ek>0.0 && ek<ec) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        if(ek>-ec && ek<ec)
                        {
                            n_a+=uk*uk*fbeta(ek,beta);
                            tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                            n_b+=vk*vk*fbeta(-1.0*ek,beta);
                            tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
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
                        if(ek>0.0 && ek<ec) // only positive energy states
                        {
                            //n_a+=... will be taken the same as n_b
                            //tau_a+=... will be taken the same as tau_b
                            n_b+=vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta);
                            tau_b+=kk2tau[ixyz]*(vk*vk*fbeta(-1.0*ek,beta) + uk*uk*fbeta(ek,beta));
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta))*2.0; // will be divided by two later
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=2.*fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
                    }
                    else
                    {
                        if(ek>-ec && ek<ec)
                        {
                            n_a+=uk*uk*fbeta(ek,beta);
                            tau_a+=kk2tau[ixyz]*uk*uk*fbeta(ek,beta);
                            n_b+=vk*vk*fbeta(-1.0*ek,beta);
                            tau_b+=kk2tau[ixyz]*vk*vk*fbeta(-1.0*ek,beta);
                            nu+=uk*vk*(fbeta(-1.0*ek,beta)-fbeta(ek,beta));
                            if(fbeta(     ek,beta)>NUMERICAL_ZERO) S-=fbeta(     ek,beta)*log(fbeta(     ek,beta));
                            if(fbeta(-1.0*ek,beta)>NUMERICAL_ZERO) S-=fbeta(-1.0*ek,beta)*log(fbeta(-1.0*ek,beta));
                            (*nwf)++;
                        }
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
            mu_a += md.init0muchange*(n0_a-n_a)/n0_a*eF_a;
            mu_b += md.init0muchange*(n0_b-n_b)/n0_b*eF_b;
            // divided by n0_c to improve convergeance
            // by factor 4 in terms of iterations
            if(md.spinsymmetry>0) mu_b=mu_a;
        }
        eF_a=pow(6.0*M_PI*M_PI*n_a, 2.0/3.0) / 2.0;
        eF_b=pow(6.0*M_PI*M_PI*n_b, 2.0/3.0) / 2.0;
        eF_avg=pow(3.0*M_PI*M_PI*(n_a+n_b), 2.0/3.0) / 2.0;
        // Effg = 0.6*n_a*eF_a*LXYZ + 0.6*n_b*eF_b*LXYZ;
        Effg = 0.6*(n_a+n_b)*eF_avg*LXYZ;
        if(printout && md.init0debug>0) wprintf("# TEMPCONV: T=%f, iter=%d, delta/eF_a=%f, mu_a/eF_a=%f, delta/eF_b=%f, mu_b/eF_b=%f\n", T, iter, delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b);
        if(iter==maxiter && printout && md.init0debug>0) wprintf("# WARNING: MAXITER REACHED!\n");

        // Compute energy
        // energy_kin=(0.5*alph_a*tau_a + 0.5*alph_b*tau_b)*LXYZ;
        energy_kin=(af_*tau_a + af_*tau_b) / 2. *LXYZ;
        energy_pot=((3. / 5.) * bf_ * eF_ * nt_)*LXYZ;
        energy_pair=-1.0*delta*nu*LXYZ;
        energy_tot=energy_kin+energy_pot+energy_pair;
        if(printout  && md.init0debug>0) wprintf("# TEMPCONV: T=%f, energy_kin=%f, energy_pot=%f, energy_pair=%f, energy_tot=%f\n", T, energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
        if(printout  && md.init0debug>0) wprintf("# TEMPCONV: T=%f, S/NkB=%16.12f\n", T, S/((n0_a+n0_b)*LXYZ));
        fflush(stdout);
    }

    if(printout && md.init0debug>0) wprintf("# ENERGY CUT-OFF: EC=%f\n", ec);

    // print results
    if(printout) wprintf("# UNIFORM SOLUTION: delta/eF_a=%8.4f, mu_a/eF_a=%8.4f, delta/eF_b=%8.4f, mu_b/eF_b=%8.4f, ec=%8.4f\n", delta/eF_a, mu_a/eF_a, delta/eF_b, mu_b/eF_b, ec);
    if(printout) wprintf("# UNIFORM SOLUTION: energy_kin=%16.12f, energy_pot=%16.12f, energy_pair=%16.12f, energy_tot=%16.12f\n", energy_kin/Effg, energy_pot/Effg, energy_pair/Effg, energy_tot/Effg);
    if(printout) wprintf("# ENTROPY PER PARTICLE: S/NkB=%16.12f\n", S/((n0_a+n0_a)*LXYZ));
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
    //##
    __md_pca_uniform.V_a=V_a;
    __md_pca_uniform.V_b=V_b;
    __md_pca_uniform.alph_a=af_;
    __md_pca_uniform.alph_b=af_;
    //##
    __md_pca_uniform.kc=kc;
    __md_pca_uniform.ec=ec;
    __md_pca_uniform.nu=nu;
    __md_pca_uniform.delta=delta;
    __md_pca_uniform.tau_a=tau_a;
    __md_pca_uniform.tau_b=tau_b;
    __md_pca_uniform.beta=beta;
    __md_pca_uniform.S=S;
    __md_pca_uniform.nwf=*nwf;
    __md_pca_uniform.ekin=energy_kin;
    __md_pca_uniform.epot=energy_pot;
    __md_pca_uniform.epair=energy_pair;

#ifdef TESTSUITE
#ifdef TDWSLDA
    if(printout) testsuite_file_uniform(TS_NERR, TS_EERR, 0, TS_MUERR, 0, TS_SERR);
#endif
#ifdef WSLDA
    if(printout) testsuite_file_uniform(TS_NERR, TS_EERR, 1, TS_MUERR, 1, TS_SERR);
#endif
#endif

    if(iter==maxiter) return 1; // not converged!
    return 0;
}



// auxliary functions for cpca code
/**
 * @param idxfrom extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param idxto extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param wf pointer to array for wf of size nx*ny*nz*2*(idxto-idxfrom)*sizeof(double complex)
 * @param mu_a chemical potential for population "a" (OUTPUT)
 * @param mu_b chemical potential for population "b" (OUTPUT)
 * @param ec energy cut-off for the solution (OUTPUT)
 * @param fEn weights used for computation densities, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param kkzvals values of corresponding kkz values, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param En eigen energies, E_n, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @param countonly if 1 the function counts only the states to be generated
 * @param cnwf number of counted wave-functions
 * @return 0 - OK, otherwise error
 * */
int create_uniform_wf_2d_generic(int idxfrom, int idxto, double complex *wf, double *mu_a, double *mu_b, double *ec, double *fEn, double *kkzvals, double *En, int printout, int countonly, int *cnwf)
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
    ixyz=0;
    int nwf=0;
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
        if(states_selected[iz]>0) // take only from sphere
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
                if(ek>0.0 && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                if(ek>-(*ec) && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }

            if(countonly!=1 && takeit) // this is my wave-function
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
                fEn[nwf-idxfrom]=ek;

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
                if(ek>0.0 && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                if(ek>-(*ec) && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }

            if(countonly!=1 && takeit) // this is my wave-function
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
                fEn[nwf-idxfrom]=ek;

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

    if(countonly==1) *cnwf=nwf+1; // save counted wave-functions

    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);
    free(states_selected);

    return 0;
}

/**
 * @param idxfrom extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param idxto extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param wf pointer to array for wf of size nx*ny*nz*2*(idxto-idxfrom)*sizeof(double complex)
 * @param mu_a chemical potential for population "a" (OUTPUT)
 * @param mu_b chemical potential for population "b" (OUTPUT)
 * @param ec energy cut-off for the solution (OUTPUT)
 * @param fEn weights used for computation densities, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param kkzvals values of corresponding kkz values, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param En eigen energies, E_n, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int create_uniform_wf_2d(int idxfrom, int idxto, double complex *wf, double *mu_a, double *mu_b, double *ec, double *fEn, double *kkzvals, double *En, int printout)
{
    return create_uniform_wf_2d_generic(idxfrom, idxto, wf, mu_a, mu_b, ec, fEn, kkzvals, En, printout, 0, NULL);
}

/**
 * Function resturns number of wave-functions that need to be evolved
 * It takes into accout that hamilonian for values kz and -kz is the same.
 * @param nwf number of wave-functions (OUTPUT)
 * @return 0 - OK, otherwise error
 * */
int get_nwf_to_evolve_2d(int *nwf)
{
    double mu_a, mu_b, ec;
    return create_uniform_wf_2d_generic(-2, -1, NULL, &mu_a, &mu_b, &ec, NULL, NULL, NULL, 0, 1, nwf);
}


/**
 * @param idxfrom extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param idxto extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param wf pointer to array for wf of size nx*ny*nz*2*(idxto-idxfrom)*sizeof(double complex)
 * @param mu_a chemical potential for population "a" (OUTPUT)
 * @param mu_b chemical potential for population "b" (OUTPUT)
 * @param ec energy cut-off for the solution (OUTPUT)
 * @param fEn weights used for computation densities, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param kkzvals values of corresponding kkz values, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param cnt degenerecies of states (OUTPUT)
 * @param En eigen energies, E_n, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @param countonly if 1 the function counts only the states to be generated
 * @param cnwf number of counted wave-functions
 * @return 0 - OK, otherwise error
 * */
int create_uniform_wf_1d_generic(int idxfrom, int idxto, double complex *wf, double *mu_a, double *mu_b, double *ec, double *fEn, double *kkyzvals, int *cnt, double *En, int printout, int countonly, int *cnwf)
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
    ixyz=0;
    int nwf=0;
    int total_nwf=0;
    int takeit;
    int deg;

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
        if(deg>0) // take only from sphere
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
                if(ek>0.0 && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                if(ek>-(*ec) && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }

            if(countonly!=1 && takeit) // this is my wave-function
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
                fEn[nwf-idxfrom]=ek;

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
                if(ek>0.0 && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }
            else
            {
                if(ek>-(*ec) && ek<(*ec))
                {
                    nwf++;
                    if(nwf>=idxfrom && nwf<idxto) takeit=1;
                }
            }

            if(countonly!=1 && takeit) // this is my wave-function
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
                fEn[nwf-idxfrom]=ek;

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

    if(countonly==1) *cnwf=nwf+1; // save counted wave-functions

    // Clear memory
    free(kkx);
    free(kky);
    free(kkz);
    free(kk2);

    return 0;
}

/**
 * @param idxfrom extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param idxto extract only wf with indices [idxfrom, idxto) (INPUT)
 * @param wf pointer to array for wf of size nx*ny*nz*2*(idxto-idxfrom)*sizeof(double complex)
 * @param mu_a chemical potential for population "a" (OUTPUT)
 * @param mu_b chemical potential for population "b" (OUTPUT)
 * @param ec energy cut-off for the solution (OUTPUT)
 * @param fEn weights used for computation densities, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param kkzvals values of corresponding kkz values, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param cnt degenerecies of states (OUTPUT)
 * @param En eigen energies, E_n, array of size (idxto-idxfrom)*sizeof(double) (OUTPUT)
 * @param printout, function prints on output info if printout is true
 * @return 0 - OK, otherwise error
 * */
int create_uniform_wf_1d(int idxfrom, int idxto, double complex *wf, double *mu_a, double *mu_b, double *ec, double *fEn, double *kkyzvals, int *cnt, double *En, int printout)
{
    return create_uniform_wf_1d_generic(idxfrom, idxto, wf, mu_a, mu_b, ec, fEn, kkyzvals, cnt, En, printout, 0, NULL);
}

/**
 * Function resturns number of wave-functions that need to be evolved
 * It takes into accout that hamilonian for values kz and -kz and ky and -ky is the same.
 * @param nwf number of wave-functions (OUTPUT)
 * @return 0 - OK, otherwise error
 * */
int get_nwf_to_evolve_1d(int *nwf)
{
    double mu_a, mu_b, ec;
    return create_uniform_wf_1d_generic(-2, -1, NULL, &mu_a, &mu_b, &ec, NULL, NULL, NULL, NULL, 0, 1, nwf);
}

#endif
