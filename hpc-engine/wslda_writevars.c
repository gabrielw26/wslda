/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 31-07-2020
 * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "wdata.h"

#include "pca_settings.h"
#include "pca_utils.h"
#include <mpi.h>
#include "wslda_writevars.h"

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        wfprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        wfprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }

#ifdef WSLDA
extern double *dc_params; /* Declaration of the variable */
extern size_t dc_extra_data_size;
extern void *dc_extra_data;
double v_ext(int ix, int iy, int iz, int it, int spin, double *params, size_t extra_data_size, void *extra_data);
double complex delta_ext(int ix, int iy, int iz, int it, double complex delta, double *params, size_t extra_data_size, void *extra_data);
double velocity_ext(int ix, int iy, int iz, int it, int spin, int coordinate, double *params, size_t extra_data_size, void *extra_data);

int get_v_ext(int datadim, int spin, int it, double *data)
{
    int ix, iy, iz, ixyz=0;
    int nx=NX, ny=NY, nz=NZ;
    if(datadim==1) {ny=1; nz=1;}
    if(datadim==2) {nz=1;}
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        data[ixyz] = v_ext(ix, iy, iz, it, spin, dc_params, dc_extra_data_size, dc_extra_data);
        ixyz++;
    }
    return 0;
}

int get_delta_ext(int datadim, int it, void *deltain, void *data)
{
    double complex *_deltain = (double complex *)deltain;
    double complex *_data = (double complex *)data;
    
    int ix, iy, iz, ixyz=0;
    int nx=NX, ny=NY, nz=NZ;
    if(datadim==1) {ny=1; nz=1;}
    if(datadim==2) {nz=1;}
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        _data[ixyz] = delta_ext(ix, iy, iz, it, _deltain[ixyz], dc_params, dc_extra_data_size, dc_extra_data);
        ixyz++;
    }
    return 0;
}

int get_velocity_ext(int datadim, int spin, int it, double *data)
{
    int ix, iy, iz, ixyz=0;
    int nx=NX, ny=NY, nz=NZ;
    if(datadim==1) {ny=1; nz=1;}
    if(datadim==2) {nz=1;}
    
    double *datax = data;
    double *datay = datax+nx*ny*nz;
    double *dataz = datay+nx*ny*nz;
    
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
    {
        datax[ixyz] = velocity_ext(ix, iy, iz, it, spin, XAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        datay[ixyz] = velocity_ext(ix, iy, iz, it, spin, YAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        dataz[ixyz] = velocity_ext(ix, iy, iz, it, spin, ZAXIS, dc_params, dc_extra_data_size, dc_extra_data);
        ixyz++;
    }
    return 0;
}
#endif

#ifdef TDWSLDA
int get_v_ext(int datadim, int spin, int it, double *data); // implemented on cuda

void * ptr_d_delta; //pointer to delta on device 
int get_delta_ext(int datadim, int it, void *deltain, void *data); // implemented on cuda
void set_ptr_d_delta(void * ptr) {ptr_d_delta=ptr;}

int get_velocity_ext(int datadim, int spin, int it, double *data); // implemented on cuda

#endif

void wdata_variable_process_var(metadata_t *input, wdata_variable *va)
{
    strcpy(va->format, input->dataformat); // set format of output results

    if(input->writeprec=='d') // double precision
    {
        if     (va->type[0]=='r') strcpy(va->type, "real");   // FIXME: change to real8 when setting default writeprec=f
        else if(va->type[0]=='c') strcpy(va->type, "complex");// FIXME: change to complex16 when setting default writeprec=f
        else if(va->type[0]=='v') strcpy(va->type, "vector"); // FIXME: change to vector8(3) when setting default writeprec=f
    }

    if(input->writeprec=='f') // float precision
    {
        if     (va->type[0]=='r') strcpy(va->type, "real4");
        else if(va->type[0]=='c') strcpy(va->type, "complex8");
        else if(va->type[0]=='v') strcpy(va->type, "vector4"); // FIXME: change to vector4(3) when setting default writeprec=f
    }
}

/**
 * @param input input structure (INPUT)
 * @param datadim data dimensonality (INPUT)
 * @param t0 initial  value of time
 * @param dt time increment between measurments
 * @param spinsymmetry flag inicating if system is spin-symmetric
 * @param wdmd metadata for wdata format (OUTPUT)
 * @return exit status 0: ok, 
 * */
int create_wdata_metadata(metadata_t *input, int datadim, double t0, double dt, int spinsymmetry, wdata_metadata *wdmd)
{
    // create artificial data for visulisation in visit
    wdata_metadata tmd = {NX, NY, NZ, DX, DY, DZ, 0, 0.0, 0.0, 0.0, "none", 0, 0.0, 0.0, 0, 0, 0, 0};
    tmd.datadim=datadim;
    sprintf(tmd.prefix, "%s", md.outprefix);
    tmd.t0=t0;
    tmd.dt=dt;
    
    // add variables
    int lnvars;
    char lvars[MAX_WRITEVARS][MAX_VARNAME_LGTH]; // and their names
    int i;
    
    // simple copy of variables to list
    lnvars = md.nwritevar;
    for(i=0; i<lnvars; i++) strcpy(lvars[i],md.writevar[i]);
    
    for(i=0; i<lnvars; i++)
    {        
        if(strcmp (lvars[i],"density") == 0 || strcmp (lvars[i],"rho") == 0)
        {
            wdata_variable va = {"rho_a", "real", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"rho_b", "real", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"rho_b", "rho_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
            
            // subset densities
            if(input->subsetMinEn!=input->subsetMaxEn)
            {
                wdata_variable sva = {"subset_rho_a", "real", "none", "wdat"}; wdata_variable_process_var(input,&sva);
                wdata_variable svb = {"subset_rho_b", "real", "none", "wdat"}; wdata_variable_process_var(input,&svb);
                wdata_link sl = {"subset_rho_b", "subset_rho_a"};
                wdata_add_variable(&tmd, &sva);
                if(spinsymmetry==0) wdata_add_variable(&tmd, &svb);
                else                wdata_add_link(&tmd, &sl);
            }

        }
        else if(strcmp (lvars[i],"current") == 0 || strcmp (lvars[i],"j") == 0)
        {
            wdata_variable va = {"j_a", "vector", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"j_b", "vector", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"j_b", "j_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
            
            // subset densities
            if(input->subsetMinEn!=input->subsetMaxEn)
            {
                wdata_variable sva = {"subset_j_a", "vector", "none", "wdat"}; wdata_variable_process_var(input,&sva);
                wdata_variable svb = {"subset_j_b", "vector", "none", "wdat"}; wdata_variable_process_var(input,&svb);
                wdata_link sl = {"subset_j_b", "subset_j_a"};
                wdata_add_variable(&tmd, &sva);
                if(spinsymmetry==0) wdata_add_variable(&tmd, &svb);
                else                wdata_add_link(&tmd, &sl);
            }
            
        }
        else if(strcmp (lvars[i],"tau") == 0)
        {
            wdata_variable va = {"tau_a", "real", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"tau_b", "real", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"tau_b", "tau_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"V") == 0)
        {
            wdata_variable va = {"V_a", "real", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"V_b", "real", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"V_b", "V_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"V_ext") == 0)
        {
            wdata_variable va = {"V_ext_a", "real", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"V_ext_b", "real", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"V_ext_b", "V_ext_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"velocity_ext") == 0)
        {
            wdata_variable va = {"velocity_ext_a", "vector", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"velocity_ext_b", "vector", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"velocity_ext_b", "velocity_ext_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"delta") == 0)
        {
            wdata_variable va = {"delta", "complex", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_add_variable(&tmd, &va);
        }
        else if(strcmp (lvars[i],"delta_ext") == 0)
        {
            wdata_variable va = {"delta_ext", "complex", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_add_variable(&tmd, &va);
        }
        else if(strcmp (lvars[i],"nu") == 0)
        {
            wdata_variable va = {"nu", "complex", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_add_variable(&tmd, &va);
        }
        else if(strcmp (lvars[i],"alpha") == 0)
        {
            wdata_variable va = {"alpha_a", "real", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"alpha_b", "real", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"alpha_b", "alpha_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"A") == 0)
        {
            wdata_variable va = {"A_a", "vector", "none", "wdat"}; wdata_variable_process_var(input,&va);
            wdata_variable vb = {"A_b", "vector", "none", "wdat"}; wdata_variable_process_var(input,&vb);
            wdata_link l = {"A_b", "A_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        
    }
    

    
    *wdmd = tmd; // copy to output buffers
    return 0;
}

/**
 * @return 1 in case if files cannot be deleted (mainly due to overwrite=0)
 * */
int clear_files(metadata_t *input, wdata_metadata *wdmd)
{
    int i;
    if(md.overwrite==0) for(i=0; i<wdmd->nvar; i++) if(wdata_file_exists(wdmd,wdmd->var[i].name)) return 1;
    
    wdata_clear_database(wdmd);
    
    return 0;
}

int write_wdata_metadata_file(metadata_t *input, wdata_metadata *wdmd, char *codename)
{
    char file_name[256];
    sprintf(file_name, "%s.wtxt", md.outprefix);
    FILE *f = fopen(file_name, "w");
    if(f==NULL) return 1;
    
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [20];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,20,"%x-%X",timeinfo);
    
    fprintf(f,"# Generated by %s, engine version %s [%s]\n\n", codename, VERSION, buffer);
    wdata_print_metadata(wdmd, f);
    
    fclose(f);
    
    return 0;
}

/**
 * Writes measurements
 * @param wdmd metadata for wdata format
 * @param mpi_comm mpi communicator: function utilizes parrallel method of writing
 * @param codetype "st" - static code, "td" - timependent code 
 * @param h_densities array of densities
 * @param h_potentials array of potentials
 * @return 0: ok, otherwise error
 * */
int write_measurments(wdata_metadata *wdmd, MPI_Comm mpi_comm, char *codetype, int it, wslda_density h_densities, wslda_potential h_potentials)
{
    int iam, np;
    MPI_Comm_size( mpi_comm , &np ) ; /* total number of processes */
    MPI_Comm_rank( mpi_comm , &iam ) ; /* id of process st 0 <= iam < np */    
    
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
    
    // pontentials
    double *V_a = h_potentials.V_a;
    double *V_b = h_potentials.V_b;
    double complex *delta = h_potentials.delta;
    
    int bs; // block size
    if(wdmd->datadim==1) bs = NX;
    else if(wdmd->datadim==2) bs = NX*NY;
    else bs = NX*NY*NZ;

    // write variables
    int ivar, ierr;
    for(ivar=0; ivar<wdmd->nvar; ivar++) if(ivar%np == iam) // each process handles different variable
    {

        ierr=0;
        if      (strcmp (wdmd->var[ivar].name,"rho_a") == 0) ierr = wdata_write_cycle_d(wdmd, "rho_a", rho_a);
        else if (strcmp (wdmd->var[ivar].name,"rho_b") == 0) ierr = wdata_write_cycle_d(wdmd, "rho_b", rho_b);
        else if (strcmp (wdmd->var[ivar].name,"delta") == 0) ierr = wdata_write_cycle_d(wdmd, "delta", (double *)delta);
        else if (strcmp (wdmd->var[ivar].name,"j_a") == 0) ierr = wdata_write_cycle_d(wdmd, "j_a", j_a_x);
        else if (strcmp (wdmd->var[ivar].name,"j_b") == 0) ierr = wdata_write_cycle_d(wdmd, "j_b", j_b_x);
        else if (strcmp (wdmd->var[ivar].name,"nu") == 0) ierr = wdata_write_cycle_d(wdmd, "nu", (double *)nu);
        else if (strcmp (wdmd->var[ivar].name,"tau_a") == 0) ierr = wdata_write_cycle_d(wdmd, "tau_a", tau_a);
        else if (strcmp (wdmd->var[ivar].name,"tau_b") == 0) ierr = wdata_write_cycle_d(wdmd, "tau_b", tau_b);
        else if (strcmp (wdmd->var[ivar].name,"V_a") == 0) 
        {
#ifdef WSLDA
            ierr = wdata_write_cycle_d(wdmd, "V_a", V_a);
#else
            double *towrt;
            cppmallocl(towrt,bs,double);  
            get_v_ext(wdmd->datadim, SPINA, it, towrt);
            int ixyz=0;
            if(it==0 && input->inittype>=1 && input->inittype<=3) // special case to maintain integrity with st codes
            {
                // the subtraction was already done in st code
                // there is no need to subtract again just after loading data from file
                for(ixyz=0; ixyz<bs; ixyz++) towrt[ixyz]=V_a[ixyz]; 
            }
            else
            {
                for(ixyz=0; ixyz<bs; ixyz++) towrt[ixyz]=V_a[ixyz]-towrt[ixyz]; // subtruct from mean-field contribution the external potential
            }
            ierr = wdata_write_cycle_d(wdmd, "V_a", towrt);
            free(towrt);
#endif
        }
        else if (strcmp (wdmd->var[ivar].name,"V_b") == 0)
        {
#ifdef WSLDA
            ierr = wdata_write_cycle_d(wdmd, "V_b", V_b);
#else
            double *towrt;
            cppmallocl(towrt,bs,double);  
            get_v_ext(wdmd->datadim, SPINB, it, towrt);
            int ixyz=0;
            if(it==0 && input->inittype>=1 && input->inittype<=3) // special case to maintain integrity with st codes
            {
                // the subtraction was already done in st code
                // there is no need to subtract again just after loading data from file
                for(ixyz=0; ixyz<bs; ixyz++) towrt[ixyz]=V_b[ixyz]; 
            }
            else
            {
                for(ixyz=0; ixyz<bs; ixyz++) towrt[ixyz]=V_b[ixyz]-towrt[ixyz]; // subtruct from mean-field contribution the external potential
            }
            ierr = wdata_write_cycle_d(wdmd, "V_b", towrt);
            free(towrt);
#endif
        }
        else if (strcmp (wdmd->var[ivar].name,"V_ext_a") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs,double);  
            get_v_ext(wdmd->datadim, SPINA, it, towrt);
            ierr = wdata_write_cycle_d(wdmd, "V_ext_a", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->var[ivar].name,"V_ext_b") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs,double);  
            get_v_ext(wdmd->datadim, SPINB, it, towrt);
            ierr = wdata_write_cycle_d(wdmd, "V_ext_b", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->var[ivar].name,"delta_ext") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs*2,double); 
#ifdef TDWSLDA
            get_delta_ext(wdmd->datadim, it, ptr_d_delta, towrt);
#else
            get_delta_ext(wdmd->datadim, it, delta, towrt);
#endif
            ierr = wdata_write_cycle_d(wdmd, "delta_ext", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->var[ivar].name,"velocity_ext_a") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs*3,double);  
            get_velocity_ext(wdmd->datadim, SPINA, it, towrt);
            ierr = wdata_write_cycle_d(wdmd, "velocity_ext_a", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->var[ivar].name,"velocity_ext_b") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs*3,double);  
            get_velocity_ext(wdmd->datadim, SPINB, it, towrt);
            ierr = wdata_write_cycle_d(wdmd, "velocity_ext_b", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->var[ivar].name,"alpha_a") == 0) ierr = wdata_write_cycle_d(wdmd, "alpha_a", h_potentials.alpha_a);
        else if (strcmp (wdmd->var[ivar].name,"alpha_b") == 0) ierr = wdata_write_cycle_d(wdmd, "alpha_b", h_potentials.alpha_b);
        else if (strcmp (wdmd->var[ivar].name,"A_a") == 0) ierr = wdata_write_cycle_d(wdmd, "A_a", h_potentials.A_a_x);
        else if (strcmp (wdmd->var[ivar].name,"A_b") == 0) ierr = wdata_write_cycle_d(wdmd, "A_b", h_potentials.A_b_x);
        
        if(ierr>0) return 100*iam+10*ivar+ierr;
    }
    
    // cycle is complete
    wdata_add_cycle(wdmd);
    
    return 0;
}

/**
 * Writes measurements for subset
 * @param wdmd metadata for wdata format
 * @param mpi_comm mpi communicator: function utilizes parrallel method of writing
 * @param codetype "st" - static code, "td" - timependent code 
 * @param h_densities array of densities
 * @return 0: ok, otherwise error
 * */
int write_measurments_subset(wdata_metadata *wdmd, MPI_Comm mpi_comm, char *codetype, int it, wslda_density h_densities)
{
    int iam, np;
    MPI_Comm_size( mpi_comm , &np ) ; /* total number of processes */
    MPI_Comm_rank( mpi_comm , &iam ) ; /* id of process st 0 <= iam < np */    
    
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
    
    int bs; // block size
    if(wdmd->datadim==1) bs = NX;
    else if(wdmd->datadim==2) bs = NX*NY;
    else bs = NX*NY*NZ;

    // write variables
    int ivar, ierr;
    for(ivar=0; ivar<wdmd->nvar; ivar++) if(ivar%np == iam) // each process handles different variable
    {

        ierr=0;
        if      (strcmp (wdmd->var[ivar].name,"subset_rho_a") == 0) ierr = wdata_write_cycle_d(wdmd, "subset_rho_a", rho_a);
        else if (strcmp (wdmd->var[ivar].name,"subset_rho_b") == 0) ierr = wdata_write_cycle_d(wdmd, "subset_rho_b", rho_b);
//         else if (strcmp (wdmd->var[ivar].name,"delta") == 0) ierr = wdata_write_cycle_d(wdmd, "delta", delta);
        else if (strcmp (wdmd->var[ivar].name,"subset_j_a") == 0) ierr = wdata_write_cycle_d(wdmd, "subset_j_a", j_a_x);
        else if (strcmp (wdmd->var[ivar].name,"subset_j_b") == 0) ierr = wdata_write_cycle_d(wdmd, "subset_j_b", j_b_x);
//         else if (strcmp (wdmd->var[ivar].name,"nu") == 0) ierr = wdata_write_cycle_d(wdmd, "nu", nu);
//         else if (strcmp (wdmd->var[ivar].name,"tau_a") == 0) ierr = wdata_write_cycle_d(wdmd, "tau_a", tau_a);
//         else if (strcmp (wdmd->var[ivar].name,"tau_b") == 0) ierr = wdata_write_cycle_d(wdmd, "tau_b", tau_b);

        if(ierr>0) return 100*iam+10*ivar+ierr;
    }
    
    return 0;
}

int exists(const char *filename);
int urm(const char *filename);
/**
 * Creates wtxt file for wavefunctions written by static codes
 * */
int create_wtxt_file_for_wf(const char * prefix, int iogroup,
                           int nwf, int nx, int ny, int nz, double dx, double dy, double dz,
                           double kF, double *mu, double ec, double beta, int codedim
                           )
{ 
    char file_name[512];
    sprintf(file_name, "%s/wf.%04d.wtxt", prefix, iogroup);

    if(exists(file_name))
    {
        if(md.overwrite==0) return -1; // We do not overwrite!  
        else urm(file_name);
    }
    
     // create metadata handler
    wdata_metadata md;
    wdata_reset_metadata(&md);

    // Lattice
    md.datadim = codedim;
    md.nx = nx;
    md.ny = ny;
    md.nz = nz;
    md.dx = dx;
    md.dy = dy;
    md.dz = dz;

    sprintf(md.prefix, "wf.%04d", iogroup);
    md.t0 = 0.0;
    md.dt = 1.0;
    md.cycles = nwf; // number of measurements

    // add variables to data set
    // for each variable binary file of name `prefix_`varname`.wdat will be created
    wdata_variable wfu = {"un", "complex", "none", "wdat"};
    wdata_add_variable(&md, &wfu);

    wdata_variable wfv = {"vn", "complex", "none", "wdat"};
    wdata_add_variable(&md, &wfv);

    wdata_setconst(&md, "kF", kF);
    wdata_setconst(&md, "mu_a", mu[SPINA]);
    wdata_setconst(&md, "mu_b", mu[SPINB]);
    wdata_setconst(&md, "ec", ec);
    wdata_setconst(&md, "beta", beta);

    // add txt files to data set
    wdata_txt wen = {"en.txt"};
    wdata_add_txt(&md, &wen);

    // write metadata file (with default name)
    wdata_write_metadata_to_file(&md, file_name);


    return 0;    
}