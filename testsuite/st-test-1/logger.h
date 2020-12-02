
static int lineid; // line id 

/**
 * This function adds new entry to `outprefix`.wlog file.
 * It is executed at the end of each iteration
 * @param log pointer to file 
 * @param it iteration number
 * @param h_densities structure with densities, see (wiki) documentation for list of fields.
 * @param h_potentials struture with potentials, see (wiki) documentation for list of fields.
 * @param kF typical Fermi momentum scale of the problem. 
 * @param energy array with contributions to energy: EKIN, EPOT, EPAIR, ECURRENT, EPOTEXT, EPAIREXT, EVELEXT
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
           double *energy, double *npart, 
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
    double E_tot = energy[EKIN]+energy[EPOT]+energy[EPAIR]+energy[ECURRENT]+energy[EPOTEXT]+energy[EPAIREXT]+energy[EVELEXT];    
    
    if(lineid==0) // HEADER
    {
        fprintf(log,"#\n");
        fprintf(log,"# ========================= LATTICE =========================\n");
        fprintf(log,"# LATTICE    : %d x %d x %d\n", NX, NY, NZ);
        fprintf(log,"# SPACING    : %.2f x %.2f x %.2f\n", DX, DY, DZ);
        fprintf(log,"# VOLUME     : %.2f x %.2f x %.2f\n", LX, LY, LZ);
        
        fprintf(log,"#\n");
        fprintf(log,"# ==================== PHYSICAL SETTINGS ====================\n");
        fprintf(log,"# E_cut              =%14.6g\n", dc_ec);
        fprintf(log,"# E_cut/eF           =%14.6g\n", dc_ec/eF);
        fprintf(log,"#\n");
        fprintf(log,"# ========================= COLUMNS ========================\n");
        fprintf(log,"#  1: iteration number\n");
        fprintf(log,"#  2: npart[SPINA]\n");
        fprintf(log,"#  3: npart[SPINB]\n");
        fprintf(log,"#  4: npart[SPINA]+npart[SPINB]\n");
        fprintf(log,"#  5: E_tot/Effg\n");
        fprintf(log,"#  6: energy[EKIN]/Effg\n");
        fprintf(log,"#  7: energy[EPOT]/Effg\n");
        fprintf(log,"#  8: energy[EPAIR]/Effg\n");
        fprintf(log,"#  9: energy[ECURRENT]/Effg\n");
        fprintf(log,"# 10: energy[EPOTEXT]/Effg\n");
        fprintf(log,"# 11: energy[EPAIREXT]/Effg\n");
        fprintf(log,"# 12: energy[EVELEXT]/Effg\n");
        fprintf(log,"# 13: mu[SPINA]/eF\n");
        fprintf(log,"# 14: mu[SPINB]/eF\n");
        fprintf(log,"# 15: kF\n");
        fprintf(log,"# 16: eF = 0.5 * kF*kF\n");
        fprintf(log,"# 17: Effg = 0.6 * (npart[SPINA]+npart[SPINB]) * eF\n");
        fprintf(log,"# 18: time per iteration (sec)\n");
        fprintf(log,"# 19: time & date of entry\n");
    }
    
    
    // add entry
    fprintf(log, "%6d %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %10.2f %20s\n",
        it, // 1
        npart[SPINA], // 2
        npart[SPINB], // 3
        npart[SPINA]+npart[SPINB], // 4
        E_tot/Effg, // 5
        energy[EKIN]/Effg, // 6
        energy[EPOT]/Effg, // 7
        energy[EPAIR]/Effg, // 8
        energy[ECURRENT]/Effg, // 9
        energy[EPOTEXT]/Effg, //10
        energy[EPAIREXT]/Effg, //11
        energy[EVELEXT]/Effg, //12
        mu[SPINA]/eF, //13
        mu[SPINB]/eF, //14
        kF, // 15
        eF, // 16
        Effg, // 17
        logger_get_time_from_last_entry(), //18
        buffer
    );
    
    lineid++; // new line 
    return 0;
}
