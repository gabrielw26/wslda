
static int lineid; // line id 

/**
 * This function defines the unit in which energies are printed in stdout.
 * @param kF typical Fermi momentum scale of the problem, value returned by referencekF() function.
 * @param mu array with chemical potentials: mu[SPINA], mu[SPINB].
 * @param npart array with computed particle numbers: npart[SPINA] and npart[SPINB].
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
double energy_unit(double kF, double *mu, double *npart,
           double *params, size_t extra_data_size, void *extra_data)
{
    double Effg;
    double eF = kF*kF/2.0; // Fermi energy
    double N = npart[SPINA]+npart[SPINB]; // total number of particles

    // depending on dimensionality of the problem
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
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
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

/**
 * Use this routine to customize the metadata file (wtxt) of data sets.
 * This function is executed once at the beginning of the code.
 * @param wdmd pointer wdata_metadata structure, see Wiki->W-data format for more info.
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
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
    // // To add variable use this template
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
 * Use this routine to write custom variable to wdata set.
 * This function is executed once at the beginning of the code.
 * @param wdmd pointer wdata_metadata structure, see Wiki->W-data format for more info.
 * @param it iteration number.
 * @param h_densities structure with densities, see (wiki) documentation for the list of fields.
 * @param h_potentials structure with potentials, see (wiki) documentation for the list of fields.
 * @param kF typical Fermi momentum scale of the problem.
 * @param mu array with values of chemical potentials
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
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

    // // to add variable to binary file use this function
    // wdata_write_cycle(wdmd, "var_name", pointer_to_data);

    return 0;
}


/**
 * Use this routine to write wave functions to file.
 * NOTE: to use this routine, you need to be familiar with the MPI and data layout that WSLDA exploits.
 * @param it iteration number.
 * @param nxyz number of lattice points
 * @param nwfip number of wave functions per MPI process
 * @param nwf total number of wave functions
 * @param beta inverse of temperature beta=1/T
 * @param comm MPI communicator
 * @param h_wf buffer for wave functions of size 2*nwfip*nxyz*sizeof(double complex).
 *             Initially, it is empty. You need to copy the wave functions from the device first.
 *             Decoding: u_n(ixyz):=h_wf[n*nxyz+ixyz]; v_n(ixyz):=h_wf[nwfip*nxyz+n*nxyz+ixyz]
 *                       where n is the index of the wave function, and ixyz is the coordinate index
 * @param h_qpe instantaneous values of quasiparticle energies, the buffer of size nwfip*sizeof(double).
 *              It is computed as expectation values of the single particle hamiltonian over quasiparticle wave functions.
 * @param h_kky wave vectors along y coordinate for plane wave representation, size nwfip*sizeof(double)
 *              Meaningful only in the case of 1d. Otherwise it is set to NULL.
 * @param h_kkz wave vectors along z coordinate for plane wave representation, size nwfip*sizeof(double)
 *              Meaningful only in the case of 1d and 2d. Otherwise it is set to NULL.
 * @param h_cnt degeneracy of states, size nwfip*sizeof(int).
 * @param d_wf device buffer with wave functions. It has the same structure as h_wf.
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if the entry has been added successfully, otherwise return the error code. If a nonzero value is returned, the main code will terminate.
 *
 * NOTES:
 *   - in order to access fields from INPUT file use `md` global structure, ie.: md.inittype, md.outprefix, etc.
 *   - inspect example code (commented out) to learn how to use this routine.
 * */
int write_wave_functions(int it, int nxyz, int nwfip, int nwf, double beta, MPI_Comm comm,
                         double complex *h_wf, double *h_qpe, double *h_kky, double *h_kkz, int *h_cnt,  // Host pointers
                         double complex *d_wf,                                                           // Device pointers
                         double *params, size_t extra_data_size, void *extra_data
                        )
{
    // // example code demonstrating writing wave functions to data format (readable by VisIt).
    // int ip, np;
    // MPI_Comm_size(comm, &np); // np = total number of processes
    // MPI_Comm_rank(comm, &ip); // id of process st 0 <= ip < np
    //
    // // Copy wave functions from device to host
    // memcopy_gpu2host(d_wf, h_wf, (size_t)2*nxyz*nwfip*sizeof(double complex));
    // // set pointers to u and v components
    // double complex *wf_u=h_wf+(size_t)0*nxyz*nwfip; // pointer to nwfip wave functions
    // double complex *wf_v=h_wf+(size_t)1*nxyz*nwfip; // pointer to nwfip wave functions
    //
    // // get info about the number of wave functions kept by each MPI process
    // int * wf_tbl = (int *)malloc(np*sizeof(int));
    // MPI_Gather( &nwfip , 1 , MPI_INT , wf_tbl , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    // MPI_Bcast( wf_tbl , np , MPI_INT , 0 , MPI_COMM_WORLD ) ;
    //
    // // use MPI I/O to write a file
    // char file_name[128];
    // MPI_File fh;
    // MPI_Offset my_offset=0;
    // MPI_Status status;
    // int i, cnt ;
    // for(i=0; i<ip; i++) my_offset+=sizeof(double complex)*nxyz*wf_tbl[i];
    //
    // // u component
    // sprintf(file_name, "%s_%06d_wf_u.wdat", md.outprefix, it);
    // if(ip==0) printf("# WRITING U COMPONENTS TO FILE `%s`\n", file_name);
    // MPI_File_open(comm, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    // MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
    // MPI_File_write_all(fh, wf_u , nxyz*nwfip, MPI_DOUBLE_COMPLEX, &status);
    // MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &cnt);
    // if(cnt!=nxyz*nwfip) { printf("# PROBLEM WITH WRITING OF U COMPONENTS BY PROCESS %s!\n", ip); return 1; }
    // MPI_File_close(&fh);
    //
    // // v component
    // sprintf(file_name, "%s_%06d_wf_v.wdat", md.outprefix, it);
    // if(ip==0) printf("# WRITING V COMPONENTS TO FILE `%s`\n", file_name);
    // MPI_File_open(comm, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    // MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
    // MPI_File_write_all(fh, wf_v , nxyz*nwfip, MPI_DOUBLE_COMPLEX, &status);
    // MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &cnt);
    // if(cnt!=nxyz*nwfip) { printf("# PROBLEM WITH WRITING OF V COMPONENTS BY PROCESS %s!\n", ip); return 2; }
    // MPI_File_close(&fh);
    //
    // // to be able to see the wave functions in VisIt add wtxt file
    // // see documentation of WDATA for more details
    // wdata_metadata wdmd;
    // wdmd.datadim = CODEDIM;
    // wdmd.nx = NX; wdmd.ny = NY; wdmd.nz = NZ;
    // wdmd.dx = DX; wdmd.dy = DY; wdmd.dz = DZ;
    // sprintf(wdmd.prefix, "%s_%06d", md.outprefix, it);
    // wdmd.t0 = 0.0; wdmd.dt = 1.0; wdmd.cycles=nwf;
    // wdata_variable vwf_u = {"wf_u", "complex", "none", "wdat"}; wdata_add_variable(&wdmd, &vwf_u);
    // wdata_variable vwf_v = {"wf_v", "complex", "none", "wdat"}; wdata_add_variable(&wdmd, &vwf_v);
    // sprintf(file_name, "%s_%06d.wtxt", md.outprefix, it);
    // if(ip==0) wdata_write_metadata_to_file(&wdmd,file_name);
    //
    // // clear memory
    // free(wf_tbl);

    return 0; // return OK status
}


