/**
 * W-SLDA Toolkit
 * Engine version: 2021.11.12
 * */


static int lineid; // line id 

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
           double kF, double *mu,
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
    double E_tot = observable[EKIN]+observable[EPOT]+observable[EPAIR]+observable[ECURRENT]+observable[EPOTEXT]+observable[EPAIREXT]+observable[EVELEXT];     
    double Emax =  M_PI*M_PI/(2.*DX*DX);
    
    if(lineid==0) // HEADER
    {
        fprintf(log,"#\n");
        fprintf(log,"# ========================= LATTICE =========================\n");
        fprintf(log,"# LATTICE    : %d x %d x %d\n", NX, NY, NZ);
        fprintf(log,"# SPACING    : %.2f x %.2f x %.2f\n", DX, DY, DZ);
        fprintf(log,"# VOLUME     : %.2f x %.2f x %.2f\n", LX, LY, LZ);
        
        fprintf(log,"#\n");
        fprintf(log,"# ==================== PHYSICAL SETTINGS ====================\n");
        fprintf(log,"# kF                 =%14.6g\n", kF);
        fprintf(log,"# eF                 =%14.6g\n", eF);
        fprintf(log,"# Effg               =%14.6g\n", Effg);
        fprintf(log,"# mu_a/eF            =%14.6g\n", mu[SPINA]/eF);
        fprintf(log,"# mu_b/eF            =%14.6g\n", mu[SPINB]/eF);   
        fprintf(log,"# E_cut              =%14.6g\n", dc_ec);
        fprintf(log,"# E_cut/eF           =%14.6g\n", dc_ec/eF);
        fprintf(log,"#\n");
        fprintf(log,"# ==================== ALGORITHM SETTINGS ===================\n");
        fprintf(log,"# np                 =%14d\n", dc_np);
        fprintf(log,"# nwf                =%14d\n", dc_nwf);
        fprintf(log,"# nwfip              =%14d\n", dc_nwfip);
        fprintf(log,"# dt*emax            =%14.6g\n", md.dt/eF*Emax);
        fprintf(log,"# dt*eF              =%14.6g\n", md.dt);
        fprintf(log,"# dt                 =%14.6g\n", md.dt/eF);
        fprintf(log,"# t0*eF              =%14.6g\n", dc_t0*eF);
        fprintf(log,"# t0                 =%14.6g\n", dc_t0);
        fprintf(log,"# ========================= COLUMNS ========================\n");
        fprintf(log,"#  1: measurement id\n");
        fprintf(log,"#  2: time*eF\n");
        fprintf(log,"#  3: npart[SPINA]\n");
        fprintf(log,"#  4: npart[SPINB]\n");
        fprintf(log,"#  5: npart[SPINA]+npart[SPINB]\n");
        fprintf(log,"#  6: E_tot/Effg\n");
        fprintf(log,"#  7: observable[EKIN]/Effg\n");
        fprintf(log,"#  8: observable[EPOT]/Effg\n");
        fprintf(log,"#  9: observable[EPAIR]/Effg\n");
        fprintf(log,"# 10: observable[ECURRENT]/Effg\n");
        fprintf(log,"# 11: observable[EPOTEXT]/Effg\n");
        fprintf(log,"# 12: observable[EPAIREXT]/Effg\n");
        fprintf(log,"# 13: observable[EVELEXT]/Effg\n");
        fprintf(log,"# 14: time per iteration (sec)\n");
        fprintf(log,"# 15: time & date of entry\n");
    }
    
    
    // add entry
    fprintf(log, "%6d %12.4f %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %10.2f %20s\n",
        lineid, // 1
        (dc_t0+it*md.dt/eF*md.timesteps) * eF, // 2
        npart[SPINA], // 3
        npart[SPINB], // 4
        npart[SPINA]+npart[SPINB], // 5
        E_tot/Effg, // 6
        observable[EKIN]/Effg, // 7
        observable[EPOT]/Effg, // 8
        observable[EPAIR]/Effg, // 9
        observable[ECURRENT]/Effg, // 10
        observable[EPOTEXT]/Effg, //11
        observable[EPAIREXT]/Effg, //12
        observable[EVELEXT]/Effg, //13
        logger_get_time_from_last_entry(), //14
        buffer // 15
    );
    
    lineid++; // new line 
    return 0;
}
