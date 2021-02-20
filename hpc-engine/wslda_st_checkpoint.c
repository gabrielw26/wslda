/**
 * W-SLDA Toolkit
 * 
 * Warsaw University of Technology, Faculty of Physics (2021)
 * 
 * @author Gabriel Wlazlowski 
 * */ 

#include "pca_settings.h"
#include "wslda_errors.h"
#include "pca_utils.h"
#include "wslda_st_checkpoint.h"
int exists(const char *filename);

/**
 * This file implements functions for reading and writing checkpoint
 * */

int wslda_st_read_checkpoint()
{
    return WSLDA_OK;
}


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
    sprintf(file_name, "%s/checkpoint.%s", md.outprefix, suffix);
    printf("# WRITING CHECKPOINT FILE `%s`\n", file_name);
    if(md.overwrite==0) if(exists(file_name)) return WSLDA_ERR_CANNOT_OVERWRITE;
    
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
    if(fwrite(&md.Mbroyden , sizeof(double)         , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // broyden
    for (i = 0; i < (md.Mbroyden + 1); i++) if(fwrite(dens_in[i]   , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    for (i = 0; i < (md.Mbroyden + 1); i++) if(fwrite(dens_out[i]  , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
            
    fclose(pFile);
            
    return WSLDA_OK;
}


