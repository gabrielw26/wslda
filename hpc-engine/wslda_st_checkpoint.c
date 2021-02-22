/**
 * W-SLDA Toolkit
 * 
 * Warsaw University of Technology, Faculty of Physics (2021)
 * 
 * @author Gabriel Wlazlowski 
 * */ 
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "wslda_errors.h"
#include "pca_utils.h"
#include "wslda_potdens.h"
#include "wslda_st_checkpoint.h"
int exists(const char *filename);
int urm( const char * fn );

/**
 * This file implements functions for reading and writing checkpoint
 * */

int wslda_stcheckpoint_format(int codedim)
{
    char file_name[512];
    
    // Step 1: Search for present data format
    sprintf(file_name, "%s/checkpoint.dat", md.inprefix);
    if(exists(file_name)) return WSLDA_ST_CHECKPOINT_DAT;
    
    // Step 2: Search for old type of checkpoint file
    if     (codedim==1) sprintf(file_name, "%s/checkpoint.s1dpca", md.inprefix);
    else if(codedim==2) sprintf(file_name, "%s/checkpoint.s2dpca", md.inprefix);
    else                sprintf(file_name, "%s/checkpoint.s3dpca", md.inprefix);
    if(exists(file_name)) return WSLDA_ST_CHECKPOINT_OLD;
    
    // no file
    return WSLDA_ST_CHECKPOINT_NOFILE;
}

int wslda_st_read_checkpoint(int fileidx, int codedim, int *it, 
                              int nconsts, double *consts, 
                              int npot, double *h_potentials, 
                              int ndens, double *h_densities, 
                              int nenergy, double *energy,
                              int nbroy, double **dens_in, double **dens_out
                              )
{
    char file_name[512];
    if(fileidx==0) 
    {
        sprintf(file_name, "%s/checkpoint.dat", md.inprefix);
        printf("# LOADING CHECKPOINT FILE `%s`\n", file_name);
    }
    else sprintf(file_name, "%s/checkpoint.dat.%d", md.outprefix, fileidx);
    
    int lattice1[4]={0,NX,NY,NZ};
    double lattice2[3]={DX, DY, DZ};
    int mbroy, i;

    FILE * pFile = fopen(file_name, "rb");
    if(pFile==NULL) return WSLDA_ERR_CANNOT_OPEN_CHECKPOINT_FILE;
    
    // read all nescesary data to file
    if(fread(lattice1     , sizeof(int)*4          , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // lattice
    if(fread(lattice2     , sizeof(double)*3       , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // lattice
    
    // check lattice
    if(lattice1[0]!=codedim) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    if(lattice1[1]!=NX) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    if(lattice1[2]!=NY) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    if(lattice1[3]!=NZ) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    if(lattice2[0]!=DX) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    if(lattice2[1]!=DY) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    if(lattice2[2]!=DZ) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    
    if(fread(it           , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // iteration number
    if(fread(consts       , sizeof(double)*nconsts , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // constants
    if(fread(h_potentials , sizeof(double)*npot    , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(h_densities  , sizeof(double)*ndens   , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(energy       , sizeof(double)*nenergy , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(&mbroy       , sizeof(double)         , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // broyden
    if(mbroy==md.Mbroyden)
    {
        for (i = 0; i < (md.Mbroyden + 1); i++) if(fread(dens_in[i]   , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
        for (i = 0; i < (md.Mbroyden + 1); i++) if(fread(dens_out[i]  , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    }
    else report_warning(WSLDA_WRN_CHECKPOINT_NOT_CONSITENT_BROYDEN, stdout);
            
    fclose(pFile);

    
    if(fileidx>0) urm(file_name); // it is temporary file, remove it!
    return WSLDA_OK;
}

static int wslda_st_checkpoint_bs(int *lattice)
{
    if(lattice[0]==1) return lattice[1];
    else if(lattice[0]==2) return lattice[1]*lattice[2];
    else if(lattice[0]==3) return lattice[1]*lattice[2]*lattice[2];
    return -1;
}

// static void wslda_st_checkpoint_2dto3d_r(int *lattice_in, double array_in, int *lattice_out, double array_out)
// {
//     int ix, iy, iz, ixyz;
//     int NXo=lattice_out[1], NYo=lattice_out[2], NZo==lattice_out[3];
//     
//     
//     ixyz=0;
//     for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
//     {
//         ixyz++;
//     }    
// }

int wslda_st_checkpoint_convert(int operation, int fileidx, int codedim, int *it, 
                              int nconsts, double *consts, 
                              int npot, double *h_potentials, 
                              int ndens, double *h_densities, 
                              int nenergy, double *energy,
                              int nbroy, double **dens_in, double **dens_out
                              )
{
    int mylattice1[4]={0,NX,NY,NZ};
    mylattice1[0]=codedim;
    
    char file_name_in[512], file_name_out[512];
    if(fileidx==0) 
    {
        sprintf(file_name_in, "%s/checkpoint.dat", md.inprefix);
        printf("# LOADING CHECKPOINT FILE `%s`\n", file_name_in);
    }
    else sprintf(file_name_in, "%s/checkpoint.dat.%d", md.outprefix, fileidx);
    sprintf(file_name_out, "%s/checkpoint.dat.%d", md.outprefix, fileidx+1);

    int mbroy, i;
    FILE * pFile;
    
    // input buffers 
    int lattice1_in[4];
    double lattice2_in[3];
    pFile = fopen(file_name_in, "rb");
    if(pFile==NULL) return WSLDA_ERR_CANNOT_OPEN_CHECKPOINT_FILE;
    
    // read all nescesary data to file
    if(fread(lattice1_in  , sizeof(int)*4          , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // lattice
    if(fread(lattice2_in  , sizeof(double)*3       , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // lattice   
    if(fread(it           , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // iteration number
    if(fread(consts       , sizeof(double)*nconsts , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // constants
    
    double *h_densities_in, *h_potentials_in;
    int DENSDIM_in = ndens;  DENSDIM_in/=wslda_st_checkpoint_bs(mylattice1);  DENSDIM_in*=wslda_st_checkpoint_bs(lattice1_in); // rescale array lenght
    int POTDIM_in = npot;    POTDIM_in /=wslda_st_checkpoint_bs(mylattice1);  POTDIM_in *=wslda_st_checkpoint_bs(lattice1_in); // rescale array lenght
    cppmallocl(h_densities_in,DENSDIM_in,double);
    cppmallocl(h_potentials_in,POTDIM_in,double);

    
    if(fread(h_potentials , sizeof(double)*POTDIM_in, 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(h_densities  , sizeof(double)*DENSDIM_in, 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(energy       , sizeof(double)*nenergy , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(&mbroy       , sizeof(double)         , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // broyden
    
    double **dens_in_in;		// pointer to array of arrays of densities
    double **dens_out_in;		// 		---//---
    int SOLDIM_in = nbroy-2; SOLDIM_in /=wslda_st_checkpoint_bs(mylattice1);  SOLDIM_in *=wslda_st_checkpoint_bs(lattice1_in); // rescale array lenght
    SOLDIM_in+=2; // get back factor 2
    if(mbroy==md.Mbroyden)
    {
        cppmallocl(dens_in_in, (md.Mbroyden + 1), double*);
        cppmallocl(dens_out_in, (md.Mbroyden + 1), double*);
        for (i = 0; i < (md.Mbroyden + 1); i++){
            cppmallocl(dens_in_in[i], SOLDIM_in, double);
            cppmallocl(dens_out_in[i], SOLDIM_in, double);
        }
        for (i = 0; i < (md.Mbroyden + 1); i++) if(fread(dens_in_in[i]   , sizeof(double)*SOLDIM_in , 1, pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
        for (i = 0; i < (md.Mbroyden + 1); i++) if(fread(dens_out_in[i]  , sizeof(double)*SOLDIM_in , 1, pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    }
    fclose(pFile);

    // output buffers
    pFile = fopen(file_name_in, "rb");
    if(pFile==NULL) return WSLDA_ERR_CANNOT_OPEN_CHECKPOINT_FILE;
    if(fwrite(lattice1_in  , sizeof(int)*4          , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // lattice
    if(fwrite(lattice2_in  , sizeof(double)*3       , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // lattice
    if(fwrite(it           , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // iteration number
    if(fwrite(consts       , sizeof(double)*nconsts , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // constants
    
    if(operation==ST_CHECKPOINT_2D_TO_3D)
    {
        printf("# CONVERTING CHECKPOINT FILE: 2D --> 3D\n"); fflush(stdout);
        

        // TODO
        return -1;
    }
    
    fclose(pFile);
    
    // clear
    free(h_densities_in);
    free(h_potentials_in);
    for (i = 0; i < (md.Mbroyden + 1); i++){
        free(dens_in_in[i]);
        free(dens_out_in[i]);
    }
    free(dens_in_in);
    free(dens_out_in);
        
    if(fileidx>0) urm(file_name_in); // it is temporary file, remove it!
    
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


