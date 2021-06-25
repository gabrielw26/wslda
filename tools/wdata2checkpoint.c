/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Tool for converting W-DATA set into checkpoint file for st-wslda codes
 * */ 

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h> 
#include <complex.h>

// From W-SLDA TOOLKIT
#define WSLDA_OK 0
#define WSLDA_ERR_CANNOT_CREATE_CHECKPOINT_FILE 1
#define WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE 2

int NX, NY, NZ;
double DX, DY, DZ;
char *outprefix;
int Mbroyden;

int wslda_st_write_checkpoint(int codedim, int it, 
                              int nconsts, double *consts, 
                              int npot, double *h_potentials, 
                              int ndens, double *h_densities, 
                              int nenergy, double *energy,
                              int nbroy, double **dens_in, double **dens_out
                              )
{
    int lattice1[4]={0,NX,NY,NZ};
    double lattice2[3]={DX, DY, DZ};
    char suffix[8]="dat";
    lattice1[0]=codedim;
    
    char file_name[512];
    sprintf(file_name, "%s_checkpoint.%s", outprefix, suffix);
    printf("# WRITING CHECKPOINT FILE `%s`\n", file_name);
    
    FILE * pFile = fopen(file_name, "wb");
    if(pFile==NULL) return WSLDA_ERR_CANNOT_CREATE_CHECKPOINT_FILE;
    
    // write all nescesary data to file
    int i=it+1;
    if(fwrite(lattice1     , sizeof(int)*4          , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // lattice
    if(fwrite(lattice2     , sizeof(double)*3       , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // lattice
    if(fwrite(&i           , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // iteration number
    if(fwrite(consts       , sizeof(double)*nconsts , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // constants
    if(fwrite(h_potentials , sizeof(double)*npot    , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    if(fwrite(h_densities  , sizeof(double)*ndens   , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    if(fwrite(energy       , sizeof(double)*nenergy , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    if(fwrite(&Mbroyden , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // broyden
// //     wprintf("WRITING Mbroyden=%d\n", Mbroyden);
    for (i = 0; i < (Mbroyden + 1); i++) if(fwrite(dens_in[i]   , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    for (i = 0; i < (Mbroyden + 1); i++) if(fwrite(dens_out[i]  , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
            
    fclose(pFile);
    
// //     wprintf("WW %d %d %d %d %d\n", nconsts, npot, ndens, nenergy, nbroy);
            
    return WSLDA_OK;
}




#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    } 
    
int main( int argc , char ** argv ) 
{
    printf("# CONVERTER OF WDATA TO CHECKPOINT FILE\n");
    
    if(argc!=3 && argc!=4 && argc!=5 && argc!=6)
    {
        printf("Usage: %s file.wtxt outprefix\n", argv[0]);
        printf("\tor\n");
        printf("Usage: %s file.wtxt outprefix cycleid\n", argv[0]);
        printf("\tfile.wtxt    - metadata file\n");
        printf("\toutprefix    - checkpoint will be written to outprefix\n");
        printf("\tcycleid      - use given cycleid to create checkpoint file, default: cycleid=last cycle \n");
        printf("\tec           - absolute value, default: ec=pi^2/(2*dx^2)\n");
        printf("\ttemperature  - in eF units, default: temperature=0\n");
        return 0;
    }
    
    // metadata handler
    wdata_metadata md;
    
    // read metadata from file
    int ierr;
    printf("# READING FILE: `%s`\n", argv[1]);
    ierr = wdata_parse_metadata_file(argv[1], &md);
    if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;} 
    
    NX=md.NX; NY=md.NY, NZ=md.NZ; 
    DX=md.DX; DY=md.DY; DZ=md.DZ;
    
    int cycleid = md.cycles-1;
    if(argc>=4) cycleid = atoi(argv[3]);
    printf("# USING cycleid=%d AS INPUT DATA\n", cycleid);
    double dc_ec = 3.141592653589793238462*3.141592653589793238462 / (2.*DX*DX);
    if(argc>=5) dc_ec = atof(argv[4]);
    printf("# SETTING ec=%f\n", dc_ec);
    double temperature=1.0e-9;
    if(argc>=6) temperature = atof(argv[5]);
    printf("# SETTING TEMPERATURE=%f\n", temperature);
    
    char file_name[1024];
    
    // memory
    int DENSDIM = 12*wdata_get_blocklength(&md);
    int POTDIM = 12*wdata_get_blocklength(&md);
    int ENERGYITEMS=7;
    
    double *h_densities, *h_potentials;
    cppmallocl(h_densities,DENSDIM,double);
    cppmallocl(h_potentials,POTDIM,double);
    double energy[ENERGYITEMS];
            
    int blocklength = wdata_get_blocklength(&md);
    double complex *nu = (double complex *)(h_densities +  0*blocklength);
    
    double *rho_a = (double *)(h_densities +  2*blocklength);
    double *tau_a = (double *)(h_densities +  3*blocklength);
    double *j_a_x = (double *)(h_densities +  4*blocklength);
    double *j_a_y = (double *)(h_densities +  5*blocklength);
    double *j_a_z = (double *)(h_densities +  6*blocklength);
    
    double *rho_b = (double *)(h_densities +  7*blocklength);
    double *tau_b = (double *)(h_densities +  8*blocklength);
    double *j_b_x = (double *)(h_densities +  9*blocklength);
    double *j_b_y = (double *)(h_densities + 10*blocklength);
    double *j_b_z = (double *)(h_densities + 11*blocklength);
    
    double *V_a = (double *)(h_potentials +  0*blocklength);
    double *V_b = (double *)(h_potentials +  1*blocklength);
    double complex *delta = (double complex *)(h_potentials +  2*blocklength);
    
    double *alpha_a = (double *)(h_potentials +  4*blocklength);
    double *alpha_b = (double *)(h_potentials +  5*blocklength);
    
    double *A_a_x = (double *)(h_potentials +  6*blocklength);
    double *A_a_y = (double *)(h_potentials +  7*blocklength);
    double *A_a_z = (double *)(h_potentials +  8*blocklength);
    double *A_b_x = (double *)(h_potentials +  9*blocklength);
    double *A_b_y = (double *)(h_potentials + 10*blocklength);
    double *A_b_z = (double *)(h_potentials + 11*blocklength);
    
    // reset data
    int i,j;
    for(i=0; i<DENSDIM; i++) h_densities[i]=0.0;
    for(i=0; i<POTDIM;  i++) h_potentials[i]=0.0;
    for(i=0; i< ENERGYITEMS; i++) energy[i] = 0.0;
    // alpha is special case
    for(i=0; i<blocklength; i++) alpha_a[i]=1.0;
    for(i=0; i<blocklength; i++) alpha_b[i]=1.0;
    
    // fill data
    if(wdata_has_variable(&md, "density_a"))
    {
        printf("# Loading: density\n");
        wdata_read_cycle(&md, "density_a", cycleid, rho_a);
        wdata_read_cycle(&md, "density_b", cycleid, rho_b);
    }
    if(wdata_has_variable(&md, "rho_a"))
    {
        printf("# Loading: rho\n");
        wdata_read_cycle(&md, "rho_a", cycleid, rho_a);
        wdata_read_cycle(&md, "rho_b", cycleid, rho_b);
    }
    if(wdata_has_variable(&md, "tau_a"))
    {
        printf("# Loading: tau\n");
        wdata_read_cycle(&md, "tau_a", cycleid, tau_a);
        wdata_read_cycle(&md, "tau_b", cycleid, tau_b);
    }
    if(wdata_has_variable(&md, "current_a"))
    {
        printf("# Loading: current\n");
        wdata_read_cycle(&md, "current_a", cycleid, j_a_x);
        wdata_read_cycle(&md, "current_b", cycleid, j_b_x);
    }
    if(wdata_has_variable(&md, "j_a"))
    {
        printf("# Loading: j\n");
        wdata_read_cycle(&md, "j_a", cycleid, j_a_x);
        wdata_read_cycle(&md, "j_b", cycleid, j_b_x);
    }    
    if(wdata_has_variable(&md, "nu"))
    {
        printf("# Loading: nu\n");
        wdata_read_cycle(&md, "nu", cycleid, nu);
    }
    
    if(wdata_has_variable(&md, "u_a"))
    {
        printf("# Loading: u\n");
        wdata_read_cycle(&md, "u_a", cycleid, V_a);
        wdata_read_cycle(&md, "u_b", cycleid, V_b);
    }
    if(wdata_has_variable(&md, "V_a"))
    {
        printf("# Loading: V\n");
        wdata_read_cycle(&md, "V_a", cycleid, V_a);
        wdata_read_cycle(&md, "V_b", cycleid, V_b);
    }
    if(wdata_has_variable(&md, "alpha_a"))
    {
        printf("# Loading: alpha\n");
        wdata_read_cycle(&md, "alpha_a", cycleid, alpha_a);
        wdata_read_cycle(&md, "alpha_b", cycleid, alpha_b);
    }
    if(wdata_has_variable(&md, "A_a"))
    {
        printf("# Loading: A\n");
        wdata_read_cycle(&md, "A_a", cycleid, A_a_x);
        wdata_read_cycle(&md, "A_b", cycleid, A_b_x);
    }
    if(wdata_has_variable(&md, "delta"))
    {
        printf("# Loading: delta\n");    
        wdata_read_cycle(&md, "delta", cycleid, delta);
    }
    
    double dc_mu_a = wdata_getconst_value(&md, "mu_a");
    double dc_mu_b = wdata_getconst_value(&md, "mu_b");
    double eF = wdata_getconst_value(&md, "eF");
    double kF = wdata_getconst_value(&md, "kF");
    double beta=1.0/(temperature*eF);
    double npart[2] = {0.0, 0.0};
    for(i=0; i<blocklength; i++) npart[0]+=rho_a[i]; 
    for(i=0; i<blocklength; i++) npart[1]+=rho_b[i]; 
    if     (md.datadim==1) {npart[0]*=md.DX*md.DY*md.DZ*md.NY*md.NZ; npart[1]*=md.DX*md.DY*md.DZ*md.NY*md.NZ;}
    else if(md.datadim==2) {npart[0]*=md.DX*md.DY*md.DZ*md.NZ      ; npart[1]*=md.DX*md.DY*md.DZ*md.NZ      ;}
    else if(md.datadim==3) {npart[0]*=md.DX*md.DY*md.DZ            ; npart[1]*=md.DX*md.DY*md.DZ            ;}
    double Effg = 0.6*(npart[0]+npart[1])*eF; // hard set!
    double twrt_consts[11] = {dc_mu_a, dc_mu_b, dc_mu_a, dc_mu_b, dc_ec, beta, eF, kF, Effg, npart[0], npart[1]};
    
    outprefix=argv[2];
    Mbroyden=5;
    double **dens_in;		// pointer to array of arrays of densities
    double **dens_out;		// 		---//---
    cppmallocl(dens_in, (Mbroyden + 1), double*);
    cppmallocl(dens_out, (Mbroyden + 1), double*);
    int SOLDIM = POTDIM;
    for (i = 0; i < (Mbroyden + 1); i++){
        cppmallocl(dens_in[i], SOLDIM + 2, double);
        cppmallocl(dens_out[i], SOLDIM + 2, double);
        
        // reset values
        for(j=0; j<SOLDIM + 2; j++) dens_in[i][j]=0.0;
        for(j=0; j<SOLDIM + 2; j++) dens_out[i][j]=0.0;
    }
    
    wslda_st_write_checkpoint(md.datadim, 0, 
                              11, twrt_consts, 
                              POTDIM, h_potentials, 
                              DENSDIM, h_densities,
                              ENERGYITEMS, energy,
                              SOLDIM + 2, dens_in, dens_out
                              );
    
    printf("# DONE.\n");
    
    return 0;
}
