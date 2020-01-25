/**
 * This is ToolKit for pca-type codes.
 * The tool converts 2d checkpoint into 3d checkpoint
* 
 * COMPILATION:
 * GW laptop
 *      gcc tkpca_checkpoint2dto3d.c -o tkpca_checkpoint2dto3d -lm 
 * */

#define file_operationl( cmd )                                                  \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "FILE ERROR:: cannot execute: %s\n" , #cmd);          \
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        return( EXIT_FAILURE ) ;                                                \
    } }

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "pca_utils.h"
#include "pca_io.h"

typedef char * string;

int main( int argc , char ** argv ) 
{
    // settings
        
    // variables
    int i, j, k; // basic iterators
    int it;
    int ix, iy, iz, ixyz; // lattice iterators
    int ixyz2d;
    int ierr; // error flag
    char file_name[256];
    double eF, kF, Effg, beta;
    double energy[5];
    string energy_labels[5];
    energy_labels[0] = "E_kin";
    energy_labels[1] = "E_pot";
    energy_labels[2] = "E_pair";
    energy_labels[3] = "E_CM";
    energy_labels[4] = "E_ext";
    double E_tot;
    double npart[2];
    
    double dc_mu_a;
    double dc_mu_b;
    double dc_ec;
    
    md.overwrite=1;
        
    // INPUT OUTPUT
    if(argc!=9)
    {
        printf("Usage: %s filename-2d filename-3d NX NY NZ DX DY DZ\n", argv[0]);
        printf("    filename-2d - [string] 2d checkpoint filename (INPUT)\n");
        printf("    filename-3d - [string] 2d checkpoint filename (OUTPUT)\n");
    }
    
    
    // INPUT
    char *infilename = argv[1];
    char *outfilename = argv[2];
    int inNX=atoi(argv[3]);
    int inNY=atoi(argv[4]);
    int inNZ=atoi(argv[5]);
    double inDX=atof(argv[6]);
    double inDY=atof(argv[7]);
    double inDZ=atof(argv[8]);
    
    // OUTPUT
    int outNX=inNX;
    int outNY=inNY;
    int outNZ=inNZ;
    double outDX=inDX;
    double outDY=inDY;
    double outDZ=inDZ;  
    
    // arrays for input
    double *h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    cppmallocl(h_densities,12*inNX*inNY,double);
    cppmallocl(h_potentials,4*inNX*inNY,double);
    
    // For easier access to data
    // densities 
    double *rho_a = (double *)(h_densities +  0*inNX*inNY);
    double *rho_b = (double *)(h_densities +  1*inNX*inNY);
    double *tau_a = (double *)(h_densities +  2*inNX*inNY);
    double *tau_b = (double *)(h_densities +  3*inNX*inNY);
    double complex *nu = (double complex *)(h_densities +  4*inNX*inNY);
    double *j_a_x = (double *)(h_densities +  6*inNX*inNY);
    double *j_a_y = (double *)(h_densities +  7*inNX*inNY);
    double *j_a_z = (double *)(h_densities +  8*inNX*inNY);
    double *j_b_x = (double *)(h_densities +  9*inNX*inNY);
    double *j_b_y = (double *)(h_densities + 10*inNX*inNY);
    double *j_b_z = (double *)(h_densities + 11*inNX*inNY);
    
    // pontentials
    double *V_a = (double *)(h_potentials +  0*inNX*inNY);
    double *V_b = (double *)(h_potentials +  1*inNX*inNY);
    double complex *delta = (double complex *)(h_potentials +  2*inNX*inNY);
    
    
    // arrays for output
    double *out_h_densities; // pointer to array of densities [rho_a, rho_b, tau_a, tau_b, nu, \vec{j}_a, \vec{j}_b] (CPU)
    double *out_h_potentials; // pointer to array with potentials [V_a, V_b, delta] (CPU)
    cppmallocl(out_h_densities,12*outNX*outNY*outNZ,double);
    cppmallocl(out_h_potentials,4*outNX*outNY*outNZ,double);
    
    // For easier access to data
    // densities 
    double *out_rho_a = (double *)(out_h_densities +  0*outNX*outNY*outNZ);
    double *out_rho_b = (double *)(out_h_densities +  1*outNX*outNY*outNZ);
    double *out_tau_a = (double *)(out_h_densities +  2*outNX*outNY*outNZ);
    double *out_tau_b = (double *)(out_h_densities +  3*outNX*outNY*outNZ);
    double complex *out_nu = (double complex *)(out_h_densities +  4*outNX*outNY*outNZ);
    double *out_j_a_x = (double *)(out_h_densities +  6*outNX*outNY*outNZ);
    double *out_j_a_y = (double *)(out_h_densities +  7*outNX*outNY*outNZ);
    double *out_j_a_z = (double *)(out_h_densities +  8*outNX*outNY*outNZ);
    double *out_j_b_x = (double *)(out_h_densities +  9*outNX*outNY*outNZ);
    double *out_j_b_y = (double *)(out_h_densities + 10*outNX*outNY*outNZ);
    double *out_j_b_z = (double *)(out_h_densities + 11*outNX*outNY*outNZ);
    
    // pontentials
    double *out_V_a = (double *)(out_h_potentials +  0*outNX*outNY*outNZ);
    double *out_V_b = (double *)(out_h_potentials +  1*outNX*outNY*outNZ);
    double complex *out_delta = (double complex *)(out_h_potentials +  2*outNX*outNY*outNZ);


    printf("# READING CHECKPOINT FILE `%s`\n", infilename);
    printf("# INPUT LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
    printf("# INPUT SPACING: %.2f x %.2f x %.2f\n", inDX, inDY, inDZ);
    printf("# INPUT VOLUME : %.2f x %.2f x %.2f\n", inDX*inNX, inDY*inNY, inDZ*inNZ);
    
    FILE * pFile = fopen(infilename, "rb");
    
    // write all nescesary data to file
    fread(&it          , sizeof(int)         , 1 , pFile); // iteration number
    fread(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
    fread(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
    fread(&dc_ec       , sizeof(double)      , 1 , pFile); 
    fread(&beta        , sizeof(double)      , 1 , pFile); 
    fread(&eF          , sizeof(double)      , 1 , pFile); 
    fread(&kF          , sizeof(double)      , 1 , pFile);
    fread(&Effg        , sizeof(double)      , 1 , pFile);
    fread(h_potentials , sizeof(double)*inNX*inNY, 4 , pFile);
    fread(h_densities  , sizeof(double)*inNX*inNY, 12, pFile);
    fread(energy       , sizeof(double)      , 5 , pFile);
    fread(npart        , sizeof(double)      , 2 , pFile);
            
    fclose(pFile);
    
    printf("# CHECKPOINT READ: it=%d\n", it);
    printf("# CHECKPOINT READ: dc_mu_a=%16.8g  dc_mu_b=%16.8g\n", dc_mu_a, dc_mu_b);
    printf("# CHECKPOINT READ: dc_ec=%16.8g  beta=%16.8g\n", dc_ec, beta);
    printf("# CHECKPOINT READ: eF=%16.8g  kF=%16.8g  Effg=%16.8g\n", eF, kF, Effg);
    printf("# CHECKPOINT READ: ------- NPART -------\n");
    printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
        "SPINA", npart[SPINA], npart[SPINA], (npart[SPINA]-npart[SPINA]));
    printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
        "SPINB", npart[SPINB], npart[SPINB], (npart[SPINB]-npart[SPINB]));  
    printf("# CHECKPOINT READ: ------- ENERGY -------\n");
    E_tot=0.0;
    for(i=0; i<5; i++)
    {
        E_tot+=energy[i]; 
        
        printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
            energy_labels[i], energy[i]/Effg, energy[i]/Effg, (energy[i]-energy[i])/Effg);
    }
    printf("  ------------------------------------------------------------------------\n");
    printf("%8s: NEW=%16.8g OLD=%16.8g DIFF=%16.8g\n", 
        "E_tot", E_tot/Effg, E_tot/Effg, (E_tot-E_tot)/Effg);
    
    

    

    printf("#\n");
    printf("# TARGET LATTICE: %d x %d x %d\n", outNX, outNY, outNZ);
    printf("# TARGET SPACING: %.2f x %.2f x %.2f\n", outDX, outDY, outDZ);
    printf("# TARGET VOLUME : %.2f x %.2f x %.2f\n", outDX*outNX, outDY*outNY, outDZ*outNZ);
    printf("#\n");
    

    // test - particle number
    for(i=0; i<2; i++)
    {
        double _npart=0.0;
        double *rho = (h_densities +  i*inNX*inNY);
        for(ixyz=0; ixyz<(inNX*inNY); ixyz++) _npart+=rho[ixyz]*inDX*inDY;
        _npart*=inDZ*inNZ;
        printf("# CHECK IN: NPART[%d]=%f\n", i, _npart);
    }
    
    
    // ---------------------- converting 2d in 3d -----------------------------
    ixyz=0;
    for(ix=0; ix<outNX; ix++) for(iy=0; iy<outNY; iy++) for(iz=0; iz<outNZ; iz++) 
    {
        ixyz2d = iy + inNY*ix;
        
        out_rho_a[ixyz] = rho_a[ixyz2d]; 
        out_rho_b[ixyz] = rho_b[ixyz2d]; 
        out_tau_a[ixyz] = tau_a[ixyz2d]; 
        out_tau_b[ixyz] = tau_b[ixyz2d]; 
        out_nu[ixyz] = nu[ixyz2d]; 
        out_j_a_x[ixyz] = j_a_x[ixyz2d]; 
        out_j_a_y[ixyz] = j_a_y[ixyz2d]; 
        out_j_a_z[ixyz] = j_a_z[ixyz2d]; 
        out_j_b_x[ixyz] = j_b_x[ixyz2d]; 
        out_j_b_y[ixyz] = j_b_y[ixyz2d]; 
        out_j_b_z[ixyz] = j_b_z[ixyz2d]; 
    
        // pontentials
        out_V_a[ixyz] = V_a[ixyz2d]; 
        out_V_b[ixyz] = V_b[ixyz2d]; 
        out_delta[ixyz] = delta[ixyz2d]; 
    
        ixyz++; // next point
    }

    
    // checkpoint 3d
    printf("# CREATING CHECKPOINT FILE `%s`\n", outfilename);
    pFile = fopen(outfilename, "wb");
    
    // write all nescesary data to file
    i=it;
    fwrite(&i           , sizeof(int)         , 1 , pFile); // iteration number
    fwrite(&dc_mu_a     , sizeof(double)      , 1 , pFile); 
    fwrite(&dc_mu_b     , sizeof(double)      , 1 , pFile); 
    fwrite(&dc_ec       , sizeof(double)      , 1 , pFile); 
    fwrite(&beta        , sizeof(double)      , 1 , pFile); 
    fwrite(&eF          , sizeof(double)      , 1 , pFile); 
    fwrite(&kF          , sizeof(double)      , 1 , pFile);
    fwrite(&Effg        , sizeof(double)      , 1 , pFile);
    fwrite(out_h_potentials , sizeof(double)*outNX*outNY*outNZ, 4 , pFile);
    fwrite(out_h_densities  , sizeof(double)*outNX*outNY*outNZ, 12, pFile);
    fwrite(energy       , sizeof(double)      , 5 , pFile);
    fwrite(npart        , sizeof(double)      , 2 , pFile);
            
    fclose(pFile);
    

    // test - particle number
    for(i=0; i<2; i++)
    {
        double _npart=0.0;
        double *rho = (out_h_densities +  i*outNX*outNY*outNZ);
        for(ixyz=0; ixyz<(outNX*outNY*outNZ); ixyz++) _npart+=rho[ixyz]*inDX*inDY*inDZ;
        printf("# CHECK OUT: NPART[%d]=%f\n", i, _npart);
    }
    
    printf("# TKPCA: DONE.\n");
    
   
    /* messy exit here */
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
