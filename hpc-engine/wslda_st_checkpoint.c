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
#include <string.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "wslda_errors.h"
#include "pca_utils.h"
#include "wslda_potdens.h"
#include "wslda_st_checkpoint.h"

int exists(const char *filename);
int urm( const char * fn );

int wslda_resize_array_2d_to_3d(char type, int nx, int ny, void *in_2d, int nz, void *out_3d);
int wslda_resize_array_1d_to_3d(char type, int nx, void *in_1d, int ny, int nz, void *out_3d);
int wslda_resize_array_1d_to_2d(char type, int nx, void *in_1d, int ny, void *out_2d);

int wslda_interpolation_1d(char type, int nxi, void *funIn, int nxo, void *funOut);
int wslda_interpolation_2d(char type, int nxi, int nyi, void *funIn, int nxo, int nyo, void *funOut);
int wslda_interpolation_3d(char type, int nxi, int nyi, int nzi, void *funIn, int nxo, int nyo, int nzo, void *funOut);

/**
 * This file implements functions for reading and writing checkpoint
 * */

int wslda_stcheckpoint_format(int codedim)
{
    char file_name[512];
    
    // Step 1: Search for present data format
    sprintf(file_name, "%s_checkpoint.dat", md.inprefix);
    if(exists(file_name)) return WSLDA_ST_CHECKPOINT_DAT;
    
    // deprecated: remove in future
    sprintf(file_name, "%s/checkpoint.dat", md.inprefix);
    if(exists(file_name)) 
    {
        char cmd[1024];
        sprintf(cmd,"cp %s/checkpoint.dat %s_checkpoint.dat", md.inprefix,md.inprefix);
        wprintf("# OLD NAMING FOR CHECKPOINT FILE `%s`. EXECUTING: %s\n", file_name, cmd);
        system(cmd);
        return WSLDA_ST_CHECKPOINT_DAT;
    }
    
    // Step 2: Search for old type of checkpoint file
    if     (codedim==1) sprintf(file_name, "%s/checkpoint.s1dpca", md.inprefix);
    else if(codedim==2) sprintf(file_name, "%s/checkpoint.s2dpca", md.inprefix);
    else                sprintf(file_name, "%s/checkpoint.s3dpca", md.inprefix);
    if(exists(file_name)) return WSLDA_ST_CHECKPOINT_OLD;
    
    // no file
    return WSLDA_ST_CHECKPOINT_NOFILE;
}

int wslda_st_required_operations(int codedim, int *intepolation, int *resize)
{
    intepolation[0]=0; 
    resize[0]=0;
    char file_name[512];
    sprintf(file_name, "%s_checkpoint.dat", md.inprefix);
    wprintf("# INSPECTING CHECKPOINT FILE `%s`\n", file_name);
    
    FILE * pFile = fopen(file_name, "rb");
    if(pFile==NULL) return WSLDA_ERR_CANNOT_OPEN_CHECKPOINT_FILE;
    
    int lattice1[4];
    double lattice2[3];
    
    // read all nescesary data to file
    if(fread(lattice1     , sizeof(int)*4          , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // lattice
    if(fread(lattice2     , sizeof(double)*3       , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // lattice
    fclose(pFile);
    
    int incodedim=lattice1[0];
    int inNX=lattice1[1], inNY=lattice1[2], inNZ=lattice1[3];
    double inDX=lattice2[0], inDY=lattice2[1], inDZ=lattice2[2];
    double inLX=inDX*inNX, inLY=inDY*inNY, inLZ=inDZ*inNZ;  
    wprintf("# CHECKPOINT FOR %dD LATTICE: [NX,NY,NZ]=[%d,%d,%d], [DX,DY,DZ]=[%.3f,%.3f,%.3f], [LX,LY,LZ]=[%.3f,%.3f,%.3f]\n", 
    incodedim, inNX, inNY, inNZ, inDX, inDY, inDZ, inLX, inLY, inLZ);
    if(incodedim==codedim && inNX==NX && inNY==NY && inNZ==NZ && inDX==DX && inDY==DY && inDZ==DZ)
    {
        wprintf("# CHECKPOINT FULLY COMPATIBLE WITH PREDEFINES.H\n"); fflush(stdout);
        return WSLDA_OK;
    }
    
    
    if(inLX!=LX || inLY!=LY || inLZ!=LZ)
    {
        wprintf("# WARNING: [LX,LY,LZ]=[%.3f,%.3f,%.3f] FOR THE TARGET LATTICE DIFFERS FROM INPUT LATTICE!\n", LX, LY, LZ);
//         return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
        report_warning(WSLDA_WRN_CHECKPOINT_UNPREDICTED, stdout);
        intepolation[0]=incodedim; 
    }
    else if(inNX!=NX || inNY!=NY || inNZ!=NZ)
    {
        report_warning(WSLDA_WRN_CHECKPOINT_DOINTERPOLATION, stdout);
        intepolation[0]=incodedim; 
    }
    
    if(incodedim<codedim)
    {
        if     (incodedim==1 && codedim==3) resize[0]=13;
        else if(incodedim==2 && codedim==3) resize[0]=23;
        else if(incodedim==1 && codedim==2) resize[0]=12;
        else return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
        report_warning(WSLDA_WRN_CHECKPOINT_DORESIZE, stdout);
    }   
    else if(incodedim>codedim) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
    
    fflush(stdout);
    
    return WSLDA_OK;
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
        sprintf(file_name, "%s_checkpoint.dat", md.inprefix);
        wprintf("# LOADING CHECKPOINT FILE `%s`\n", file_name);
    }
    else sprintf(file_name, "%s_checkpoint.dat.%d", md.outprefix, fileidx);
    
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
    if(fread(&mbroy       , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // broyden
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
    else if(lattice[0]==3) return lattice[1]*lattice[2]*lattice[3];
    return -1;
}

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
        sprintf(file_name_in, "%s_checkpoint.dat", md.inprefix);
        wprintf("# LOADING CHECKPOINT FILE `%s`\n", file_name_in);
    }
    else sprintf(file_name_in, "%s_checkpoint.dat.%d", md.outprefix, fileidx);
    sprintf(file_name_out, "%s_checkpoint.dat.%d", md.outprefix, fileidx+1);

    int mbroy, i, j, k, l;
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
    int bs_in = wslda_st_checkpoint_bs(lattice1_in);
    int bs = wslda_st_checkpoint_bs(mylattice1);
    int POTDIM_in = npot;    POTDIM_in /=bs;  POTDIM_in *=bs_in; // rescale array lenght
    int DENSDIM_in = ndens;  DENSDIM_in/=bs;  DENSDIM_in*=bs_in; // rescale array lenght
    cppmallocl(h_potentials_in,POTDIM_in,double);
    cppmallocl(h_densities_in,DENSDIM_in,double);

    if(fread(h_potentials_in, sizeof(double)*POTDIM_in, 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(h_densities_in , sizeof(double)*DENSDIM_in, 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(energy       , sizeof(double)*nenergy , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE;
    if(fread(&mbroy       , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE; // broyden
//     wprintf("111 mbroy=%d md.Mbroyden=%d\n", mbroy, md.Mbroyden);
    double **dens_in_in;		// pointer to array of arrays of densities
    double **dens_out_in;		// 		---//---
    int SOLDIM_in = nbroy-2; SOLDIM_in /=bs;  SOLDIM_in *=bs_in; // rescale array lenght
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
    
//     wprintf("CI %d %d %d %d %d %d %d\n", nconsts, POTDIM_in, DENSDIM_in, nenergy, SOLDIM_in ,bs_in, bs);

    // output buffers  
    int interNX, interNY, interNZ;
    if(operation==ST_CHECKPOINT_2D_TO_3D) 
    { 
        if(lattice1_in[1]!=NX || lattice1_in[2]!=NY || lattice1_in[0]!=2) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
        
        lattice1_in[0]=3; // dimension go up 
        lattice1_in[3]=NZ; 
        
        wprintf("# CONVERTING CHECKPOINT FILE: 2D --> 3D\n"); fflush(stdout);
    }
    else if(operation==ST_CHECKPOINT_1D_TO_3D) 
    { 
        if(lattice1_in[1]!=NX || lattice1_in[0]!=1) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
        
        lattice1_in[0]=3; // dimension go up 
        lattice1_in[2]=NY;
        lattice1_in[3]=NZ; 
        
        wprintf("# CONVERTING CHECKPOINT FILE: 1D --> 3D\n"); fflush(stdout);
    }    
    else if(operation==ST_CHECKPOINT_1D_TO_2D) 
    { 
        if(lattice1_in[1]!=NX || lattice1_in[0]!=1) return WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE;
        
        lattice1_in[0]=2; // dimension go up 
        lattice1_in[2]=NY;
        lattice1_in[3]=NZ; 
        
        wprintf("# CONVERTING CHECKPOINT FILE: 1D --> 2D\n"); fflush(stdout);
    } 
    else if(operation==ST_CHECKPOINT_RESIZE)
    {
        wprintf("# CONVERTING CHECKPOINT WITH RESOLUTION [DX,DY,DZ]=[%.3f,%.3f,%.3f] TO NEW RESOLUTION [%.3f,%.3f,%.3f]\n",
            lattice2_in[0], lattice2_in[1], lattice2_in[2], DX, DY, DZ
        ); fflush(stdout);
        interNX=lattice1_in[1]; // make copy ...
        interNY=lattice1_in[2]; 
        interNZ=lattice1_in[3]; 
        lattice1_in[1]=NX;      // change...
        lattice1_in[2]=NY;
        lattice1_in[3]=NZ;        
        lattice2_in[0]=DX;
        lattice2_in[1]=DY;
        lattice2_in[2]=DZ;
    }
    else return WSLDA_ERR_INTRISTIC_ERROR;
    
    int bs_out = wslda_st_checkpoint_bs(lattice1_in);
    POTDIM_in = npot;    POTDIM_in /=bs;  POTDIM_in *=bs_out; // rescale array lenght
    DENSDIM_in = ndens;  DENSDIM_in/=bs;  DENSDIM_in*=bs_out; // rescale array lenght
    SOLDIM_in = nbroy-2; SOLDIM_in /=bs;  SOLDIM_in *=bs_out; // rescale array lenght
    SOLDIM_in+=2; // get back factor 2
//     wprintf("COO %d %d %d %d %d %d %d %d\n", nconsts, POTDIM_in, DENSDIM_in, nenergy, SOLDIM_in, bs_in, bs, bs_out);
    
    pFile = fopen(file_name_out, "wb");
    if(pFile==NULL) return WSLDA_ERR_CANNOT_OPEN_CHECKPOINT_FILE;
    if(fwrite(lattice1_in  , sizeof(int)*4          , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // lattice
    if(fwrite(lattice2_in  , sizeof(double)*3       , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // lattice
    if(fwrite(it           , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // iteration number
    if(fwrite(consts       , sizeof(double)*nconsts , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // constants
    
    double *ptr_in, *ptr;
    char type[32];

    // h_potentials array
    ptr_in = h_potentials_in;
    ptr = h_potentials;
    strcpy(type, POTTYPE);
    j=0; k=0;
    while(k<POTCNT)
    {
//             wprintf("POTCNT %d %d %c\n", k, j, type[j]);
        if(operation==ST_CHECKPOINT_2D_TO_3D) 
        { 
            if(wslda_resize_array_2d_to_3d(type[j], NX, NY, ptr_in, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
        }
        else if(operation==ST_CHECKPOINT_1D_TO_3D)
        {
            if(wslda_resize_array_1d_to_3d(type[j], NX, ptr_in, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
        }
        else if(operation==ST_CHECKPOINT_1D_TO_2D)
        {
            if(wslda_resize_array_1d_to_2d(type[j], NX, ptr_in, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
        }
        else if(operation==ST_CHECKPOINT_RESIZE)
        {
            if(lattice1_in[0]==1) 
            {
                if(wslda_interpolation_1d(type[j], interNX, ptr_in, NX, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
            }
            else if(lattice1_in[0]==2) 
            {
                if(wslda_interpolation_2d(type[j], interNX, interNY, ptr_in, NX, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
            }
            else
            {
                if(wslda_interpolation_3d(type[j], interNX, interNY, interNZ, ptr_in, NX, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
            }
        }
        if(type[j]=='c') { ptr_in+=bs_in*2; ptr+=bs_out*2; k+=2;}
        else { ptr_in+=bs_in; ptr+=bs_out; k+=1; }          
        j++;
    }
    if(fwrite(h_potentials , sizeof(double)*POTDIM_in , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    
    // h_densities array
    ptr_in = h_densities_in;
    ptr = h_densities;
    strcpy(type, DENSTYPE);
    j=0; k=0;
    while(k<DENSCNT)
    {
// //             wprintf("DENSCNT %d %d %c\n", k, j, type[j]);
        if(operation==ST_CHECKPOINT_2D_TO_3D) 
        { 
            if(wslda_resize_array_2d_to_3d(type[j], NX, NY, ptr_in, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
        }
        else if(operation==ST_CHECKPOINT_1D_TO_3D)
        {
            if(wslda_resize_array_1d_to_3d(type[j], NX, ptr_in, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
        }
        else if(operation==ST_CHECKPOINT_1D_TO_2D)
        {
            if(wslda_resize_array_1d_to_2d(type[j], NX, ptr_in, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
        }
        else if(operation==ST_CHECKPOINT_RESIZE)
        {
            if(lattice1_in[0]==1) 
            {
                if(wslda_interpolation_1d(type[j], interNX, ptr_in, NX, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
            }
            else if(lattice1_in[0]==2) 
            {
                if(wslda_interpolation_2d(type[j], interNX, interNY, ptr_in, NX, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
            }
            else
            {
                if(wslda_interpolation_3d(type[j], interNX, interNY, interNZ, ptr_in, NX, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
            }
        }
        if(type[j]=='c') { ptr_in+=bs_in*2; ptr+=bs_out*2; k+=2;}
        else { ptr_in+=bs_in; ptr+=bs_out; k+=1; }          
        j++;
    }
    if(fwrite(h_densities  , sizeof(double)*DENSDIM_in , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    
    if(fwrite(energy       , sizeof(double)*nenergy , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    if(fwrite(&mbroy       , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // broyden
// //         wprintf("mbroy=%d md.Mbroyden=%d\n", mbroy, md.Mbroyden);
// //         wprintf("CO %d %d %d %d %d\n", nconsts, POTDIM_in, DENSDIM_in, nenergy, SOLDIM_in);
    if(mbroy==md.Mbroyden)
    {
        if(md.mixingtype=='d') {strcpy(type, DENSTYPE); l=DENSCNT;}
        else                   {strcpy(type, POTTYPE);  l=POTCNT;}
        
        for (i = 0; i < (md.Mbroyden + 1); i++) 
        {
            ptr_in = dens_in_in[i];
            ptr = dens_in[i];
            j=0; k=0;
            while(k<l)
            {
// //                     wprintf("lin %d %d %c %d %d %d\n", k, j, type[j],bs_in,bs,bs_out);
                if(operation==ST_CHECKPOINT_2D_TO_3D) 
                { 
                    if(wslda_resize_array_2d_to_3d(type[j], NX, NY, ptr_in, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
                }
                else if(operation==ST_CHECKPOINT_1D_TO_3D)
                {
                    if(wslda_resize_array_1d_to_3d(type[j], NX, ptr_in, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
                }
                else if(operation==ST_CHECKPOINT_1D_TO_2D)
                {
                    if(wslda_resize_array_1d_to_2d(type[j], NX, ptr_in, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
                }
                else if(operation==ST_CHECKPOINT_RESIZE)
                {
                    if(lattice1_in[0]==1) 
                    {
                        if(wslda_interpolation_1d(type[j], interNX, ptr_in, NX, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
                    }
                    else if(lattice1_in[0]==2) 
                    {
                        if(wslda_interpolation_2d(type[j], interNX, interNY, ptr_in, NX, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
                    }
                    else
                    {
                        if(wslda_interpolation_3d(type[j], interNX, interNY, interNZ, ptr_in, NX, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
                    }
                }
                if(type[j]=='c') { ptr_in+=bs_in*2; ptr+=bs_out*2; k+=2;}
                else { ptr_in+=bs_in; ptr+=bs_out; k+=1; }          
                j++;
            }
            ptr[0]=ptr_in[0]; ptr[1]=ptr_in[1]; // last two elemnts 
            if(fwrite(dens_in[i]   , sizeof(double)*SOLDIM_in , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
        }
        
        for (i = 0; i < (md.Mbroyden + 1); i++) 
        {
            ptr_in = dens_out_in[i];
            ptr = dens_out[i];
            j=0; k=0;
            while(k<l)
            {
// //                     wprintf("lin %d %d %c %d %d %d\n", k, j, type[j],bs_in,bs,bs_out);
                if(operation==ST_CHECKPOINT_2D_TO_3D) 
                { 
                    if(wslda_resize_array_2d_to_3d(type[j], NX, NY, ptr_in, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
                }
                else if(operation==ST_CHECKPOINT_1D_TO_3D)
                {
                    if(wslda_resize_array_1d_to_3d(type[j], NX, ptr_in, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
                }
                else if(operation==ST_CHECKPOINT_1D_TO_2D)
                {
                    if(wslda_resize_array_1d_to_2d(type[j], NX, ptr_in, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR; 
                }
                else if(operation==ST_CHECKPOINT_RESIZE)
                {
                    if(lattice1_in[0]==1) 
                    {
                        if(wslda_interpolation_1d(type[j], interNX, ptr_in, NX, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
                    }
                    else if(lattice1_in[0]==2) 
                    {
                        if(wslda_interpolation_2d(type[j], interNX, interNY, ptr_in, NX, NY, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
                    }
                    else
                    {
                        if(wslda_interpolation_3d(type[j], interNX, interNY, interNZ, ptr_in, NX, NY, NZ, ptr)!=WSLDA_OK) return WSLDA_ERR_INTRISTIC_ERROR;
                    }
                }
                if(type[j]=='c') { ptr_in+=bs_in*2; ptr+=bs_out*2; k+=2;}
                else { ptr_in+=bs_in; ptr+=bs_out; k+=1; }          
                j++;
            }
            ptr[0]=ptr_in[0]; ptr[1]=ptr_in[1]; // last two elemnts 
            if(fwrite(dens_out[i] , sizeof(double)*SOLDIM_in , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
        }
    }

    
    fclose(pFile);
    
    // clear
    free(h_densities_in);
    free(h_potentials_in);
    if(mbroy==md.Mbroyden)
    {
        for (i = 0; i < (md.Mbroyden + 1); i++){
            free(dens_in_in[i]);
            free(dens_out_in[i]);
        }
        free(dens_in_in);
        free(dens_out_in);
    }
        
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
    sprintf(file_name, "%s_checkpoint.%s", md.outprefix, suffix);
    wprintf("# WRITING CHECKPOINT FILE `%s`\n", file_name);
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
    if(fwrite(&md.Mbroyden , sizeof(int)            , 1 , pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE; // broyden
// //     wprintf("WRITING md.Mbroyden=%d\n", md.Mbroyden);
    for (i = 0; i < (md.Mbroyden + 1); i++) if(fwrite(dens_in[i]   , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
    for (i = 0; i < (md.Mbroyden + 1); i++) if(fwrite(dens_out[i]  , sizeof(double)*nbroy , 1, pFile)!=1) return WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE;
            
    fclose(pFile);
    
// //     wprintf("WW %d %d %d %d %d\n", nconsts, npot, ndens, nenergy, nbroy);
            
    return WSLDA_OK;
}


