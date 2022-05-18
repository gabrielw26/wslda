
#include "stdio.h"

static int lineid; // line id 

/**
 * This function provides size of extra_data array, in bytes.
 * The extra_data of specified size will be allocated by the main process.
 * This function is thread-safe.
 * @param params with input file parameters. 
 *              NOTE: the array contains bare input file values, not processed by process_params()!
 * @return size of the extra_data array that needs to be allocated, if 0 then extra_data will not be allocated.
 * */
size_t get_extra_data_size(double *params)
{
    return 0;
}

/**
 * This function loads data into extra_data array.
 * This function is thread-safe.
 * @param size size of array computed using function get_extra_data_size()
 * @param extra_data pointer to array that should be filled with data
 * @param params with input file parameters. 
 *               NOTE: the array contains bare input file values, not processed by process_params()!
 * @return 0 if load is successful, otherwise return error code. If nonzero value is returned the main code will terminate.
 * */
int load_extra_data(size_t size, void *extra_data, double *params)
{
    return 0;
}

/** 
 * THIS FUNCTION IS CALLED DURING THE SELF-CONSISTENT PROCESS.
 * After loading params array from input file, the parameters are processed by this routine.
 * The routine is executed at beginning of each iteration.
 * @param params array of size MAX_USER_PARAMS with parameters from input file. 
 * @param kF typical Fermi momentum scale of the problem. 
 *           kF=referencekF if the referencekF tag is indicated in the input file, 
 *           otherwise to kF value is assigned according formula kF=(3*pi^2*n)^{1/3}, where n corresponds to maximal density.
 *           You can also set kF at request in this function using (*kF)=myvalue;
 * @param mu array with chemical potentials: mu[SPINA] and mu[SPINB]. 
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
void process_params(double *params, double *kF, double *mu, size_t extra_data_size, void *extra_data)
{
    // PROCESS INPUT FILE PARAMETERS
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
