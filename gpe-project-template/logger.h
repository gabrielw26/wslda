
#include "stdio.h"

static int lineid; // line id 

/**
 * This function computes Fermi momentum, which is used as the reference value. 
 * Other reference scales are set automatically to: eF=kF^2/2, Effg=(3/5)*N*eF (N-total number of particles)
 * For more details see: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Reference%20scales
 * NOTE units are: hbar=m=k_b=1
 * @param it iteration number
 * @param h_densities structure with densities, see (wiki) documentation for list of fields
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of Fermi momentum for your problem
 * */
double referencekF(int it, wslda_density h_densities, double *params, size_t extra_data_size, void *extra_data)
{
    if(input->referencekF>0.0) return input->referencekF; // take it from input file
    
    // define here your prescription for computing kF
    // ...
    // default: extract max density and use it for definition of kF
    double max_dens=0.0, kF;
    int ixyz;
    for(ixyz=0; ixyz<h_densities.nx*h_densities.ny*h_densities.nz; ixyz++) 
        if(h_densities.rho_a[ixyz]+h_densities.rho_b[ixyz]>max_dens) max_dens=h_densities.rho_a[ixyz]+h_densities.rho_b[ixyz];
    
    // depending on dimensionality of the problem
    if(NY==1 && NZ==1) kF = 0.5*M_PI*max_dens;                // 1D
    else if(NZ==1)     kF = pow(2.0*M_PI*max_dens,1./2.);     // 2D
    else               kF = pow(3.*M_PI*M_PI*max_dens,1./3.); // 3D
    
    return kF;
}

/**
 * This function adds new entry to `outprefix`.wlog file.
 * It is executed at the end of each iteration
 * @param log pointer to file 
 * @param it iteration number
 * @param h_densities structure with densities, see (wiki) documentation for list of fields.
 * @param h_potentials struture with potentials, see (wiki) documentation for list of fields.
 * @param kF typical Fermi momentum scale of the problem. 
 * @param observable array with observables: 
 *                      contributions to the energy: EKIN, EPOT, EPAIR, ECURRENT, EPOTEXT, EPAIREXT, EVELEXT
 *                      entropy: ENTROPY
 * @param npart array with computed particle numbers: npart[SPINA] and  npart[SPINB]
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if entry has been added successfuly, otherwise return error code. If nonzero value is returned the main code will terminate.
 * 
 * NOTES: 
 *   - in order to access fields from INPUT file use `md` global structure, ie.: md.inittype, md.outprefix, etc.
 * */
int logger(FILE *log, 
           int it, 
           wslda_density h_densities, wslda_potential h_potentials, 
           double kF, double *mu, //<-- IGNORE
           double *observable, double *npart, 
           double *params, size_t extra_data_size, void *extra_data)
{
    
    // Time stamp
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [20];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,20,"%x-%X",timeinfo);
    
    double eF = 0.5 * kF*kF;
    double Effg = 0.6 * (npart[SPINA]+npart[SPINB]) * eF;
    double E_tot = observable[EKIN]+observable[EPOT]+observable[EPOTEXT];  
    
    if(lineid==0) // HEADER
    {
        fprintf(log,"#\n");
        fprintf(log,"# ========================= LATTICE =========================\n");
        fprintf(log,"# LATTICE    : %d x %d x %d\n", NX, NY, NZ);
        fprintf(log,"# SPACING    : %.2f x %.2f x %.2f\n", DX, DY, DZ);
        fprintf(log,"# VOLUME     : %.2f x %.2f x %.2f\n", LX, LY, LZ);
        fprintf(log,"#\n");
        fprintf(log,"# ========================= COLUMNS ========================\n");
        fprintf(log,"#  1: iteration number\n");
        fprintf(log,"#  2: E_tot/Effg\n");
        fprintf(log,"#  3: observable[EKIN]/Effg\n");
        fprintf(log,"#  4: observable[EPOT]/Effg\n");
        fprintf(log,"#  5: observable[EPOTEXT]/Effg\n");
        fprintf(log,"#  6: kF\n");
        fprintf(log,"#  7: eF = 0.5 * kF*kF\n");
        fprintf(log,"#  8: Effg = 0.6 * (npart[SPINA]+npart[SPINB]) * eF\n");
        fprintf(log,"#  9: time per iteration (sec)\n");
        fprintf(log,"# 10: time & date of entry\n");
    }
    
    
    // add entry
    fprintf(log, "%6d %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %10.2f %20s\n",
        it, // 1
        E_tot/Effg, // 2
        observable[EKIN]/Effg, // 3
        observable[EPOT]/Effg, // 4
        observable[EPOTEXT]/Effg, //5
        kF, // 6
        eF, // 7
        Effg, // 8
        logger_get_time_from_last_entry(), //9
        buffer // 10
    );
    
    lineid++; // new line 
    return 0;
}
