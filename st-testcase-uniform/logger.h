
static int lineid; // line id 

/**
 * This function defines the unit in which energies are printed in stdout.
 * @param kF typical Fermi momentum scale of the problem, value returned by referencekF() function.
 * @param mu array with chemical potentials: mu[SPINA], mu[SPINB].
 * @param npart array with computed particle numbers: npart[SPINA] and npart[SPINB].
 * @param params array of input parameters, before the call of this routine, the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
double energy_unit(double kF, double *mu, double *npart,
           double *params, size_t extra_data_size, void *extra_data)
{
    double Effg;
    double eF = kF*kF/2.0; // Fermi energy
    double N = npart[SPINA]+npart[SPINB]; // total number of particles

    // depending on the dimensionality of the problem
    if(NY==1 && NZ==1) Effg=(1./3.)*N*eF;   // 1D
    else if(NZ==1)     Effg=(1./2.)*N*eF;   // 2D
    else               Effg=(3./5.)*N*eF;   // 3D

    return Effg;
}

/**
 * This function adds a new entry to `outprefix`.wlog file.
 * It is executed at the end of each iteration
 * @param log pointer to file
 * @param it iteration number
 * @param h_densities structure with densities, see (wiki) documentation for the list of fields.
 * @param h_potentials structure with potentials, see (wiki) documentation for the list of fields.
 * @param kF typical Fermi momentum scale of the problem.
 * @param mu array with values of chemical potentials
 * @param observable array with observables:
 *                      contributions to the energy: EKIN, EPOT, EPAIR, ECURRENT, EPOTEXT, EPAIREXT, EVELEXT
 *                      entropy: ENTROPY
 * @param npart array with computed particle numbers: npart[SPINA] and  npart[SPINB]
 * @param params array of input parameters, before the call of this routine, the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if the entry has been added successfully, otherwise return the error code. If a nonzero value is returned, the main code will terminate.
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
    double Effg = energy_unit(kF, mu, npart, params, extra_data_size, extra_data);
    double E_tot = observable[EKIN]+observable[EPOT]+observable[EPAIR]+observable[ECURRENT]+observable[EPOTEXT]+observable[EPAIREXT]+observable[EVELEXT];    
    
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
        fprintf(log,"#  6: observable[EKIN]/Effg\n");
        fprintf(log,"#  7: observable[EPOT]/Effg\n");
        fprintf(log,"#  8: observable[EPAIR]/Effg\n");
        fprintf(log,"#  9: observable[ECURRENT]/Effg\n");
        fprintf(log,"# 10: observable[EPOTEXT]/Effg\n");
        fprintf(log,"# 11: observable[EPAIREXT]/Effg\n");
        fprintf(log,"# 12: observable[EVELEXT]/Effg\n");
        fprintf(log,"# 13: observable[ENTROPY]/(npart[SPINA]+npart[SPINB])\n");
        fprintf(log,"# 14: mu[SPINA]/eF\n");
        fprintf(log,"# 15: mu[SPINB]/eF\n");
        fprintf(log,"# 16: kF\n");
        fprintf(log,"# 17: eF = 0.5 * kF*kF\n");
        fprintf(log,"# 18: Effg = 0.6 * (npart[SPINA]+npart[SPINB]) * eF\n");
        fprintf(log,"# 19: time per iteration (sec)\n");
        fprintf(log,"# 20: time & date of entry\n");
    }
    
    
    // add entry
    fprintf(log, "%6d %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %18.10g %10.2f %20s\n",
        it, // 1
        npart[SPINA], // 2
        npart[SPINB], // 3
        npart[SPINA]+npart[SPINB], // 4
        E_tot/Effg, // 5
        observable[EKIN]/Effg, // 6
        observable[EPOT]/Effg, // 7
        observable[EPAIR]/Effg, // 8
        observable[ECURRENT]/Effg, // 9
        observable[EPOTEXT]/Effg, //10
        observable[EPAIREXT]/Effg, //11
        observable[EVELEXT]/Effg, //12
        observable[ENTROPY]/(npart[SPINA]+npart[SPINB]), //13
        mu[SPINA]/eF, //14
        mu[SPINB]/eF, //15
        kF, // 16
        eF, // 17
        Effg, // 18
        logger_get_time_from_last_entry(), //19
        buffer // 20
    );
    
    lineid++; // new line 
    return 0;
}

/**
 * Use this routine to customize the metadata file (wtxt) of data sets.
 * This function is executed once at the beginning of the code.
 * @param wdmd pointer wdata_metadata structure, see Wiki->W-data format for more info.
 * @param params array of input parameters, before the call of this routine, the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if the entry has been added successfully, otherwise return the error code. If a nonzero value is returned, the main code will terminate.
 *
 * NOTES:
 *   - in order to access fields from INPUT file use `md` global structure, ie.: md.inittype, md.outprefix, etc.
 * */
int add_custom_variable_to_wdata_metadata(wdata_metadata *wdmd,
           double *params, size_t extra_data_size, void *extra_data)
{
    // // To add a variable, use this template
    // //                       var_name       type    unit
    // wdata_variable var1 = {"real_var_name", "real", "none", "wdat"}; // scalar variable
    // wdata_add_variable(wdmd, &var1);
    // wdata_variable var2 = {"complex_var_name", "complex", "none", "wdat"}; // complex variable
    // wdata_add_variable(wdmd, &var2);
    // wdata_variable var3 = {"vector_var_name", "vector", "none", "wdat"}; // vector variable
    // wdata_add_variable(wdmd, &var3);

    return 0;
}

/**
 * Use this routine to write a custom variable to the wdata set.
 * This function is executed for each writing event of the observables.
 * @param wdmd pointer wdata_metadata structure, see Wiki->W-data format for more info.
 * @param it iteration number.
 * @param h_densities structure with densities, see (wiki) documentation for the list of fields.
 * @param h_potentials structure with potentials, see (wiki) documentation for the list of fields.
 * @param kF typical Fermi momentum scale of the problem.
 * @param mu array with values of chemical potentials
 * @param params array of input parameters, before the call of this routine, the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if the entry has been added successfully, otherwise return the error code. If a nonzero value is returned, the main code will terminate.
 *
 * NOTES:
 *   - in order to access fields from INPUT file use `md` global structure, ie.: md.inittype, md.outprefix, etc.
 * */
int write_custom_variable_to_wdata_set(wdata_metadata *wdmd,
           int it,
           wslda_density h_densities, wslda_potential h_potentials,
           double kF, double *mu,
           double *params, size_t extra_data_size, void *extra_data)
{
    // // DETERMINE LOCAL SIZES OF ARRAYS (CODE DIMENSIONALITY DEPENDENT)
    // int lNX=h_densities.nx, lNY=h_densities.ny, lNZ=h_densities.nz; // local sizes
    // int ix, iy, iz, ixyz;
    //
    // // ITERATE OVER ALL POINTS
    // ixyz=0;
    // for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    // {
    //     double x = DX*(ix-lNX/2);
    //     double y = DY*(iy-lNY/2); // for 1d code y will be always 0
    //     double z = DZ*(iz-lNZ/2); // for 1d and 2d codes z will be always 0
    //
    //     ixyz++; // go to the next point, it should be the last line of the triple loop
    // }

    // // to add a variable to a binary file, use this function
    // wdata_write_cycle(wdmd, "var_name", pointer_to_data);

    return 0;
}
