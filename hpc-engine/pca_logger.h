// Author: Gabriel Wlazlowski
// Date: 08-02-2017

// Here I implement functions needed for reporting results

int create_header_of_runlog(const char *execcmd, double kF, double Effg, double *mu, double ec,
                            int nwf, int np, int nwfip
                           )
{
    FILE * log;
    char file_name[256];
    int i;
    sprintf(file_name, "%s_run.log", md.outprefix);
        
    // open file
    log = fopen (file_name,"w");
    if(log==NULL) // error - cannot create the file
        return 1;    

    time_t rawtime;
    time ( &rawtime );
    fprintf(log,"# CREATION TIME OF THE LOG: %s", ctime (&rawtime));
    fprintf(log,"# EXECUTION COMMAND       : %s\n", execcmd);
    fprintf(log,"# CODE NAME               : %s\n",STRINGIZE(CODE));
    fprintf(log,"# VERSION OF THE CODE     : %s\n", VERSION);
    fprintf(log,"# COMPILATION DATE & TIME : %s, %s\n", __DATE__,__TIME__);
    fprintf(log,"# \n");
    fprintf(log,"# ============== COMPILATION MACRO-VARIABLES ================\n");
    fprintf(log,"# LATTICE    : %d x %d x %d\n", NX, NY, NZ);
    fprintf(log,"# NXYZ       : %d\n", NXYZ);
    fprintf(log,"# SPACING    : %.2f x %.2f x %.2f\n", DX, DY, DZ);
    fprintf(log,"# VOLUME     : %.2f x %.2f x %.2f\n", LX, LY, LZ);
    fprintf(log,"# INTEGRATION_SCHEME=%d\n", INTEGRATION_SCHEME);
#ifdef TARGET_MACHINE
    fprintf(log,"# TARGET_MACHINE=%s\n",STRINGIZE(TARGET_MACHINE));
#endif
#ifdef VERBOSE
    fprintf(log,"# VERBOSE\n");
#endif
#ifdef DEBUG
    fprintf(log,"# DEBUG\n");
#endif
#ifdef TIMING
    fprintf(log,"# TIMING\n");
#endif
#ifdef EPSILON
    fprintf(log,"# EPSILON=%g\n",EPSILON);  
#endif
#ifdef MAX_USER_PARAMS
    fprintf(log,"# MAX_USER_PARAMS=%d\n", MAX_USER_PARAMS);
#endif
    fprintf(log,"# UD_SCITERS=%d\n", UD_SCITERS);
    fprintf(log,"# UD_MIX_COEFF=%f\n", UD_MIX_COEFF);
    fprintf(log,"# DENSEPSILON=%g\n", DENSEPSILON);
    fprintf(log,"# N_STABILITY_CRITERIA=%f\n", N_STABILITY_CRITERIA);
#ifdef WORK_IN_ROTATING_FRAME
    fprintf(log,"# WORK_IN_ROTATING_FRAME: YES\n");
#else
    fprintf(log,"# WORK_IN_ROTATING_FRAME: NO\n");
#endif 
#ifdef STORE_QPE
    fprintf(log,"# STORE_QPE: YES\n");
#else
    fprintf(log,"# STORE_QPE: NO\n");
#endif
#ifdef ENABLE_DELTA_EXT
    fprintf(log,"# ENABLE_DELTA_EXT: YES\n");
#else
    fprintf(log,"# ENABLE_DELTA_EXT: NO\n");
#endif
#ifdef ENABLE_VELOCITY_EXT
    fprintf(log,"# ENABLE_VELOCITY_EXT: YES\n");
#else
    fprintf(log,"# ENABLE_VELOCITY_EXT: NO\n");
#endif
#if FUNCTIONAL==BDG
    fprintf(log,"# FUNCTIONAL: BDG\n");
#elif FUNCTIONAL==SLDA
    fprintf(log,"# FUNCTIONAL: SLDA\n");
#elif FUNCTIONAL==ASLDA
    fprintf(log,"# FUNCTIONAL: ASLDA\n");
#endif
    fprintf(log,"# \n");
    fprintf(log,"# ================== EDF MACRO-VARIABLES ====================\n");
    if(fabs(md.aBdG)<1.0e-12) fprintf(log, "# ENERGY DENSITY FUNCTIONAL: (A)SLDA [UNITARITY]\n");
    else                      fprintf(log, "# ENERGY DENSITY FUNCTIONAL: BdG [a=%16.8f]\n", md.aBdG);    
#ifdef BDG_MODE
    fprintf(log,"# BDG_MODE\n");
    fprintf(log,"# a=%f\n", md.aBdG);
    fprintf(log,"# a*kF=%f\n", md.aBdG*kF);
#endif
#ifdef CURRENT_CORRECTIONS
    fprintf(log,"# CURRENT_CORRECTIONS: YES\n");
#else
    fprintf(log,"# CURRENT_CORRECTIONS: NO\n");
#endif
#ifdef SPINSYMMETRY_MODE
    fprintf(log,"# SPINSYMMETRY_MODE: YES\n");
#else    
    fprintf(log,"# SPINSYMMETRY_MODE: NO\n");
#endif
#ifdef FAST_CONST_EFFECTIVE_MASS_MODE
    fprintf(log,"# FAST_CONST_EFFECTIVE_MASS_MODE: YES\n");
#else
    fprintf(log,"# FAST_CONST_EFFECTIVE_MASS_MODE: NO\n");
#endif  
#ifdef USE_CUBIC_CUTOFF
    fprintf(log,"# USE_CUBIC_CUTOFF: YES\n");
#else
    fprintf(log,"# USE_CUBIC_CUTOFF: NO\n");
#endif
#ifdef UNIFORM_TEST_MODE
    fprintf(log,"# UNIFORM_TEST_MODE: YES\n");
#else
    fprintf(log,"# UNIFORM_TEST_MODE: NO\n");
#endif
#ifdef TAU_COMPUTATION_VIA_GRADIENTS
    fprintf(log,"# TAU_COMPUTATION_VIA_GRADIENTS: YES\n");
#else
    fprintf(log,"# TAU_COMPUTATION_VIA_GRADIENTS: NO\n");
#endif
    fprintf(log,"# A0=%f\n", A0);
    fprintf(log,"# A1=%f\n", A1);
    fprintf(log,"# A2=%f\n", A2);
    fprintf(log,"# G0=%f\n", G0);
    fprintf(log,"# G1=%f\n", G1);
    fprintf(log,"# GAMMA0=%f\n", GAMMA0);
    fprintf(log,"# P_NMIN=%g\n", P_NMIN);
    fprintf(log,"# P_NMAX=%g\n", P_NMAX);
    fprintf(log,"# P_ALPHA=%f\n", P_ALPHA);
    fprintf(log,"# \n");
    fprintf(log,"# ==================== PHYSICAL SETTINGS ====================\n");
    fprintf(log,"# kF                 =%14.6g\n", kF);
    double eF = 0.5 * kF*kF;
    fprintf(log,"# eF                 =%14.6g\n", eF);
    fprintf(log,"# Effg               =%14.6g\n", Effg);
    fprintf(log,"# mu_a/eF            =%14.6g\n", mu[SPINA]/eF);
    fprintf(log,"# mu_b/eF            =%14.6g\n", mu[SPINB]/eF);   
    fprintf(log,"# E_cut              =%14.6g\n", ec);
    fprintf(log,"# E_cut/eF           =%14.6g\n", ec/eF);    
#ifdef TDWSLDA
    double Emax =  M_PI*M_PI/2.;
    fprintf(log,"# dt*emax            =%14.6g\n", md.dt/eF*Emax);
    fprintf(log,"# dt*eF              =%14.6g\n", md.dt);
    fprintf(log,"# dt                 =%14.6g\n", md.dt/eF);
#else
    fprintf(log,"# dt*emax            =%14.6g\n", md.dt);
    double Emax =  M_PI*M_PI/2.;
    fprintf(log,"# dt*eF              =%14.6g\n", md.dt/Emax*eF);
    fprintf(log,"# dt                 =%14.6g\n", md.dt/Emax);
#endif
    fprintf(log,"# nwf                =%14d\n", nwf);
    fprintf(log,"# \n");
    fprintf(log,"# ==================== QUANTUM FRICTION =====================\n");
    fprintf(log,"# Quant.fr. alpha    =%14.6g\n", md.qfalpha);
    fprintf(log,"# Quant.fr. beta     =%14.6g\n", md.qfbeta);
    fprintf(log,"# Quant.fr. gamma    =%14.6g\n", md.qfgamma);
    fprintf(log,"# Quant.fr. start    =%14.6g\n", md.qfstart);
    fprintf(log,"# Quant.fr. stop     =%14.6g\n", md.qfstop);
    fprintf(log,"# Quant.fr. switch   =%14.6g\n", md.qfswitch);
    fprintf(log,"# \n");
    fprintf(log,"# =================== CURRENT CORRECTIONS ===================\n");
#ifdef CURRENT_CORRECTIONS
    fprintf(log,"# CURRENT_CORRECTIONS: YES\n");
#else
    fprintf(log,"# CURRENT_CORRECTIONS: NO\n");
#endif
    fprintf(log,"# Current.cor. start =%14.6g\n", md.ccstart);
    fprintf(log,"# Current.cor. stop  =%14.6g\n", md.ccstop);
    fprintf(log,"# Current.cor. switch=%14.6g\n", md.ccswitch);
    fprintf(log,"# \n");
    fprintf(log,"# ==================== ALGORITHM SETTINGS ===================\n");
    fprintf(log,"# np                 =%14d\n", np);
    fprintf(log,"# nwfip              =%14d\n", nwfip);
    fprintf(log,"# nthreads           =%14d\n", md.nthreads);
    fprintf(log,"# batch              =%14d\n", md.batch);
    fprintf(log,"# inittype           =%14d\n", md.inittype);
    fprintf(log,"# measurements       =%14d\n", md.measurements);
    fprintf(log,"# timesteps          =%14d\n", md.timesteps);
    fprintf(log,"# Prefix of reports  =  %s\n", md.outprefix);
    fprintf(log,"# Prefix of input    =  %s\n", md.inprefix); 
    fprintf(log,"# overwrite          =%14d\n", md.overwrite);
    fprintf(log,"# selfstart          =%14d\n", md.selfstart);
    fprintf(log,"# checkpoint         =%14d\n", md.checkpoint);
    fprintf(log,"# walltime           =%14.2f\n", md.walltime);
    fprintf(log,"# spinsymmetry       =%14d\n", md.spinsymmetry);
    fprintf(log,"\n");
    fprintf(log,"# ======================= PARAMS ARRAY ======================\n");
    for(i=0; i<MAX_USER_PARAMS; i++)
        fprintf(log,"# params[%2d]        =%14.6g\n", i, md.params[i]);
    fprintf(log,"\n");
    
    // close file
    fclose(log);
    
    return 0;
}

int add_line_to_file(int imeasurement, double rt, int nr_of_items, double *items)
{
    FILE * log;
    char file_name[256];
    sprintf(file_name, "%s_run.log", md.outprefix);
        
    log = fopen (file_name,"a");
    if(log==NULL) // error - cannot create the file
        return 1;
    
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [20];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,20,"%x-%X",timeinfo);
    
    char line[1024] = "";
    char item[32];
    int i;
    sprintf(item, "%8d", imeasurement);
    strcat(line,item);
    for(i=0; i< nr_of_items; i++) 
    {
        sprintf(item, " %18.10g ", items[i]);
        strcat (line,item);
    }
    sprintf(item, " %8.3f ", rt);
    strcat(line,item);
    sprintf(item, " %20s", buffer);
    strcat(line,item);
    
    fprintf(log, "%s\n", line);
    
    fclose(log);
    
    return 0;   
}
// ------------------------------- CUSTOM LOGGER -------------------------------------
static double __logger_t_gettimeofday ;
static struct timeval __logger_s ;
int logger_create_header(const char *execcmd)
{
    FILE * log;
    char file_name[512];
    int i;
    sprintf(file_name, "%s.wlog", md.outprefix);
    
    if(md.overwrite==0) // break if file exists
    {
        log = fopen (file_name,"r"); // try to open
        if(log!=NULL)
        {
            fclose(log);
            return 1;
        }
    }
    
    // create file
    log = fopen (file_name,"w");
    if(log==NULL) // error - cannot create the file
        return 2;
    
    time_t rawtime;
    time ( &rawtime );
    fprintf(log,"# CREATION TIME OF THE LOG: %s", ctime (&rawtime));
    fprintf(log,"# EXECUTION COMMAND       : %s\n", execcmd);
    fprintf(log,"# CODE NAME               : %s\n",STRINGIZE(CODE));
    fprintf(log,"# VERSION OF THE CODE     : %s\n", VERSION);
    fprintf(log,"# COMPILATION DATE & TIME : %s, %s\n", __DATE__,__TIME__);
    
    fclose(log);
    
    // reset timers
    gettimeofday( &__logger_s , NULL ) ;
    __logger_t_gettimeofday = __logger_s.tv_sec + 1e-6 * __logger_s.tv_usec ;
    
    return 0;
}

double logger_get_time_from_last_entry()
{
    gettimeofday( &__logger_s , NULL ) ;
    double rt = __logger_s.tv_sec + 1e-6 * __logger_s.tv_usec - __logger_t_gettimeofday ;
    __logger_t_gettimeofday = __logger_s.tv_sec + 1e-6 * __logger_s.tv_usec ;
    return rt;
}

// implemented in logger.h
int logger(FILE *log, 
           int it, 
           wslda_density h_densities, wslda_potential h_potentials, 
           double kF, double *mu,
           double *energy, double *npart, 
           double *params, size_t extra_data_size, void *extra_data);


int logger_add_entry(int it, 
           wslda_density h_densities, wslda_potential h_potentials, 
           double kF, double *mu,
           double *energy, double *npart, 
           double *params, size_t extra_data_size, void *extra_data)
{
    FILE * log;
    char file_name[512];
    int i;
    sprintf(file_name, "%s.wlog", md.outprefix);
    
    // append
    log = fopen (file_name,"a");
    if(log==NULL) // error - cannot create the file
        return -1;
    
    int ierr = logger(log, it, h_densities, h_potentials, kF, mu, energy, npart, params, extra_data_size, extra_data);
    
    fclose(log);
    
    return ierr;
}
