/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Tool for converting W-DATA set into checkpoint file for ws-wslda codes
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


#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    } 
    
int main( int argc , char ** argv ) 
{
    printf("CONVERTER OF WDATA TO CHECKPOINT FILE\n");
    
    if(argc!=3)
    {
        printf("Usage: %s file.wtxt outprefix\n", argv[0]);
        printf("\tfile.wtxt    - metadata file\n");
        printf("\toutprefix    - for new files, metadata file will be written to outprefix.wtxt\n");
        return 0;
    }
    
    // metadata handler
    wdata_metadata md;
    
    // read metadata from file
    int ierr;
    printf("READING FILE: `%s`\n", argv[1]);
    ierr = wdata_parse_metadata_file(argv[1], &md);
    if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;} 
    
    // prepare forlder for checkpoint
    char cmd[1024];
    char file_name[1024];
    sprintf(cmd, "mkdir -p %s", argv[2]);
    printf("EXECUTING: `%s`\n", cmd);
    system(cmd);
    
    // set suffix
    char suffix[8];
    if     (md.datadim==1) sprintf(suffix,"s1dpca");
    else if(md.datadim==2) sprintf(suffix,"s2dpca");
    else if(md.datadim==3) sprintf(suffix,"s3dpca");
    else {printf("Incorrect value of datadim in wtxt file!\n"); return 1;}
    
    // memory
    int DENSDIM = 12*wdata_get_blocklength(&md);
    int POTDIM = 12*wdata_get_blocklength(&md);
    int ENERGYITEMS=7;
    
    double *h_densities, *h_potentials;
    cppmallocl(h_densities,DENSDIM,double);
    cppmallocl(h_potentials,POTDIM,double);
    double energy[ENERGYITEMS];
    
    // reset data
    int i;
    for(i=0; i<DENSDIM; i++) h_densities[i]=0.0;
    for(i=0; i<POTDIM;  i++) h_potentials[i]=0.0;
    for(i=0; i< ENERGYITEMS; i++) energy[i] = 0.0;
    
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
    
    // fill data
    wdata_read_cycle(&md, "density_a", md.cycles-1, rho_a);
    wdata_read_cycle(&md, "density_b", md.cycles-1, rho_b);
    wdata_read_cycle(&md, "tau_a", md.cycles-1, tau_a);
    wdata_read_cycle(&md, "tau_b", md.cycles-1, tau_b);
    wdata_read_cycle(&md, "current_a", md.cycles-1, j_a_x);
    wdata_read_cycle(&md, "current_b", md.cycles-1, j_b_x);
    wdata_read_cycle(&md, "nu", md.cycles-1, nu);
    
    wdata_read_cycle(&md, "u_a", md.cycles-1, V_a);
    wdata_read_cycle(&md, "u_b", md.cycles-1, V_b);
    wdata_read_cycle(&md, "alpha_a", md.cycles-1, alpha_a);
    wdata_read_cycle(&md, "alpha_b", md.cycles-1, alpha_b);
    wdata_read_cycle(&md, "A_a", md.cycles-1, A_a_x);
    wdata_read_cycle(&md, "A_b", md.cycles-1, A_b_x);
    wdata_read_cycle(&md, "delta", md.cycles-1, delta);
    
    sprintf(file_name, "%s/checkpoint.%s", argv[2], suffix);
    printf("CREATING CHECKPOINT FILE `%s`\n", file_name);
    FILE * pFile = fopen(file_name, "wb");
    
    // write all nescesary data to file
    i=md.cycles;
    fwrite(&i           , sizeof(int)         , 1 , pFile); // iteration number
    
    double dc_mu_a = wdata_getconst_value(&md, "mu_a");
    double dc_mu_b = wdata_getconst_value(&md, "mu_b");
    fwrite(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
    fwrite(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
    
    double dc_ec = 3.141592653589793238462*3.141592653589793238462 / (2.*md.DX*md.DX);
    fwrite(&dc_ec       , sizeof(double)      , 1 , pFile);
    
    double npart[2] = {0.0, 0.0};
    for(i=0; i<blocklength; i++) npart[0]+=rho_a[i]; 
    for(i=0; i<blocklength; i++) npart[1]+=rho_b[i]; 
    if     (md.datadim==1) {npart[0]*=md.DX*md.DY*md.DZ*md.NY*md.NZ; npart[1]*=md.DX*md.DY*md.DZ*md.NY*md.NZ;}
    else if(md.datadim==2) {npart[0]*=md.DX*md.DY*md.DZ*md.NZ      ; npart[1]*=md.DX*md.DY*md.DZ*md.NZ      ;}
    else if(md.datadim==3) {npart[0]*=md.DX*md.DY*md.DZ            ; npart[1]*=md.DX*md.DY*md.DZ            ;}

    double beta=1.0e9; // zero temperature - hard set
    fwrite(&beta        , sizeof(double)      , 1 , pFile);
    
    double eF = wdata_getconst_value(&md, "eF");
    fwrite(&eF          , sizeof(double)      , 1 , pFile); 
    
    double kF = wdata_getconst_value(&md, "kF");
    fwrite(&kF          , sizeof(double)      , 1 , pFile);
    
    double Effg = 1.0; // hard set!
    fwrite(&Effg        , sizeof(double)      , 1 , pFile);
    
    fwrite(h_potentials , sizeof(double)*POTDIM, 1 , pFile);
    fwrite(h_densities  , sizeof(double)*DENSDIM,1, pFile);
    fwrite(energy       , sizeof(double)      , ENERGYITEMS , pFile);
    fwrite(npart        , sizeof(double)      , 2 , pFile);
    fwrite(&dc_mu_a , sizeof(double)      , 1 , pFile); 
    fwrite(&dc_mu_b , sizeof(double)      , 1 , pFile);
    
    // Ignore
    printf("NOTE: Buffers for broyden are not filled!\n");
//     for (i = 0; i < (md.Mbroyden + 1); i++) fwrite(dens_in[i]   , sizeof(double) , DENSDIM + 2 , pFile);
//     for (i = 0; i < (md.Mbroyden + 1); i++) fwrite(dens_out[i]  , sizeof(double) , DENSDIM + 2 , pFile);
//             
    fclose(pFile);
    
    printf("DONE.\n");
    
    return 0;
}
