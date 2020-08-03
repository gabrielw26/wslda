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

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
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
    wdata_metadata tmd = {NX, NY, NZ, DX, DY, DZ, 0, "none", 0, 0.0, 0.0, 0, 0};
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
        if(strcmp (lvars[i],"density") == 0)
        {
            wdata_variable va = {"density_a", "real", "none"};
            wdata_variable vb = {"density_b", "real", "none"};
            wdata_link l = {"density_b", "density_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"current") == 0)
        {
            wdata_variable va = {"current_a", "vector", "none"};
            wdata_variable vb = {"current_b", "vector", "none"};
            wdata_link l = {"current_b", "current_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"tau") == 0)
        {
            wdata_variable va = {"tau_a", "real", "none"};
            wdata_variable vb = {"tau_b", "real", "none"};
            wdata_link l = {"tau_b", "tau_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"u") == 0)
        {
            wdata_variable va = {"u_a", "real", "none"};
            wdata_variable vb = {"u_b", "real", "none"};
            wdata_link l = {"u_b", "u_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"v_ext") == 0)
        {
            wdata_variable va = {"v_ext_a", "real", "none"};
            wdata_variable vb = {"v_ext_b", "real", "none"};
            wdata_link l = {"v_ext_b", "v_ext_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"velocity_ext") == 0)
        {
            wdata_variable va = {"velocity_ext_a", "vector", "none"};
            wdata_variable vb = {"velocity_ext_b", "vector", "none"};
            wdata_link l = {"velocity_ext_b", "velocity_ext_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
        else if(strcmp (lvars[i],"delta") == 0)
        {
            wdata_variable va = {"delta", "complex", "none"};
            wdata_add_variable(&tmd, &va);
        }
        else if(strcmp (lvars[i],"delta_ext") == 0)
        {
            wdata_variable va = {"delta_ext", "complex", "none"};
            wdata_add_variable(&tmd, &va);
        }
        else if(strcmp (lvars[i],"nu") == 0)
        {
            wdata_variable va = {"nu", "complex", "none"};
            wdata_add_variable(&tmd, &va);
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
    if(md.overwrite==0) for(i=0; i<wdmd->nvars; i++) if(wdata_file_exists(wdmd,wdmd->vars[i].name)) return 1;
    
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
    
    fprintf(f,"# Generated by %s [%s]\n\n", codename, buffer);
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
int write_measurments(wdata_metadata *wdmd, MPI_Comm mpi_comm, char *codetype, int it, double *h_densities, double *h_potentials)
{
    int iam, np;
    MPI_Comm_size( mpi_comm , &np ) ; /* total number of processes */
    MPI_Comm_rank( mpi_comm , &iam ) ; /* id of process st 0 <= iam < np */    
    
    double *rho_a;
    double *rho_b;
    double *tau_a;
    double *tau_b;
    double complex *nu;
    double *j_a_x;
    double *j_a_y;
    double *j_a_z;
    double *j_b_x;
    double *j_b_y;
    double *j_b_z;
    
    // pontentials
    double *V_a;
    double *V_b;
    double complex *delta;
    
    int bs; // block size
    if(wdmd->datadim==1) bs = NX;
    else if(wdmd->datadim==2) bs = NX*NY;
    else bs = NX*NY*NZ;

    if (strcmp (codetype,"st") == 0)
    {
        rho_a = (double *)(h_densities +  0*bs);
        rho_b = (double *)(h_densities +  1*bs);
        tau_a = (double *)(h_densities +  2*bs);
        tau_b = (double *)(h_densities +  3*bs);
        nu = (double complex *)(h_densities +  4*bs);
        j_a_x = (double *)(h_densities +  6*bs);
        j_a_y = (double *)(h_densities +  7*bs);
        j_a_z = (double *)(h_densities +  8*bs);
        j_b_x = (double *)(h_densities +  9*bs);
        j_b_y = (double *)(h_densities + 10*bs);
        j_b_z = (double *)(h_densities + 11*bs);
        
        // potentials
        V_a = (double *)(h_potentials +  0*bs);
        V_b = (double *)(h_potentials +  1*bs);
        delta = (double complex *)(h_potentials +  2*bs);
    }
    else if (strcmp (codetype,"td") == 0)
    {
        nu = (double complex *)(h_densities +  0*bs);
        rho_a = (double *)(h_densities +  2*bs);
        tau_a = (double *)(h_densities +  3*bs);
        j_a_x = (double *)(h_densities +  4*bs);
        j_a_y = (double *)(h_densities +  5*bs);
        j_a_z = (double *)(h_densities +  6*bs);    
        rho_b = (double *)(h_densities +  7*bs);
        tau_b = (double *)(h_densities +  8*bs);
        j_b_x = (double *)(h_densities +  9*bs);
        j_b_y = (double *)(h_densities + 10*bs);
        j_b_z = (double *)(h_densities + 11*bs);
    
        // potentials
        V_a = (double *)(h_potentials +  0*bs);
        V_b = (double *)(h_potentials +  1*bs);
        delta = (double complex *)(h_potentials +  2*bs);
    }
    else return -1;
    
    // write variables
    int ivar, ierr;
    for(ivar=0; ivar<wdmd->nvars; ivar++) if(iam==0)/*if(ivar%np == iam)*/ // each process handles different variable
    {

        ierr=0;
        if      (strcmp (wdmd->vars[ivar].name,"density_a") == 0) ierr = wdata_write_cycle(wdmd, "density_a", rho_a);
        else if (strcmp (wdmd->vars[ivar].name,"density_b") == 0) ierr = wdata_write_cycle(wdmd, "density_b", rho_b);
        else if (strcmp (wdmd->vars[ivar].name,"delta") == 0) ierr = wdata_write_cycle(wdmd, "delta", delta);
        else if (strcmp (wdmd->vars[ivar].name,"current_a") == 0) ierr = wdata_write_cycle(wdmd, "current_a", j_a_x);
        else if (strcmp (wdmd->vars[ivar].name,"current_b") == 0) ierr = wdata_write_cycle(wdmd, "current_b", j_b_x);
        else if (strcmp (wdmd->vars[ivar].name,"nu") == 0) ierr = wdata_write_cycle(wdmd, "nu", nu);
        else if (strcmp (wdmd->vars[ivar].name,"tau_a") == 0) ierr = wdata_write_cycle(wdmd, "tau_a", tau_a);
        else if (strcmp (wdmd->vars[ivar].name,"tau_b") == 0) ierr = wdata_write_cycle(wdmd, "tau_b", tau_b);
        else if (strcmp (wdmd->vars[ivar].name,"u_a") == 0) ierr = wdata_write_cycle(wdmd, "u_a", V_a);
        else if (strcmp (wdmd->vars[ivar].name,"u_b") == 0) ierr = wdata_write_cycle(wdmd, "u_b", V_b);
        else if (strcmp (wdmd->vars[ivar].name,"v_ext_a") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs,double);  
            get_v_ext(wdmd->datadim, SPINA, it, towrt);
            ierr = wdata_write_cycle(wdmd, "v_ext_a", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->vars[ivar].name,"v_ext_b") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs,double);  
            get_v_ext(wdmd->datadim, SPINB, it, towrt);
            ierr = wdata_write_cycle(wdmd, "v_ext_b", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->vars[ivar].name,"delta_ext") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs*2,double); 
#ifdef TDWSLDA
            get_delta_ext(wdmd->datadim, it, ptr_d_delta, towrt);
#else
            get_delta_ext(wdmd->datadim, it, delta, towrt);
#endif
            ierr = wdata_write_cycle(wdmd, "delta_ext", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->vars[ivar].name,"velocity_ext_a") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs*3,double);  
            get_velocity_ext(wdmd->datadim, SPINA, it, towrt);
            ierr = wdata_write_cycle(wdmd, "velocity_ext_a", towrt);
            free(towrt);
        }
        else if (strcmp (wdmd->vars[ivar].name,"velocity_ext_b") == 0) 
        {
            double *towrt;
            cppmallocl(towrt,bs*3,double);  
            get_velocity_ext(wdmd->datadim, SPINB, it, towrt);
            ierr = wdata_write_cycle(wdmd, "velocity_ext_b", towrt);
            free(towrt);
        }
        
        if(ierr>0) return 100*iam+10*ivar+ierr;
    }
    
    // cycle is complete
    wdata_add_cycle(wdmd);
    
    return 0;
}

