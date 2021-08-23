/**
 * W-SLDA Toolkit
 * 
 * This code converts existing wdata set into new one defined on lattice with different resolution. 
 * Use this code to learn how to use interpolation routines. 
 * New lattice it taken from corresponding predefines.h file
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 interpolate-dataset.c -I. -I$WSLDA/hpc-engine -I$WSLDA/lib/wdata/c -L$WSLDA/lib/wdata -lwdatac -o interpolate-dataset -lm -lfftw3
 * 
 * NOTE: you need before generate wdata lib for C compiler:
 *    cd $WSLDA/lib-wdata
 *    make libc
 * */   


// Standard libraries
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

// W-DATA Format
// Must be before wslda_toolkit.h !
#include "wdata.h"
#include "winterp.h"

// W-SLDA Toolkit API
int wsldapid;
#include "wslda_toolkit.h"

int main( int argc , char ** argv ) 
{
    wsldapid=0;
    int ierr;
    
    if(argc!=3)
    {
        printf("Usage: %s input-dataset.wtxt outprefix\n", argv[0]);
        
        return( EXIT_FAILURE ) ; 
    }
    
    // create metadata handler
    wdata_metadata wdmd;
    
    printf("# READING INPUT DATA: `%s`\n", argv[1]);
    file_operationl( wdata_parse_metadata_file(argv[1], &wdmd) );
    
    // Check
    int inNX=wdata_getNX(&wdmd), inNY=wdata_getNY(&wdmd), inNZ=wdata_getNZ(&wdmd);
    double inDX=wdata_getDX(&wdmd), inDY=wdata_getDY(&wdmd), inDZ=wdata_getDZ(&wdmd);
    double inLX=inDX*inNX, inLY=inDY*inNY, inLZ=inDZ*inNZ;
    printf("# **********************  INPUT LATTICE **********************\n");
    printf("# LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
    printf("# SPACING: %f x %f x %f\n", inDX, inDY, inDZ);
    printf("# VOLUME : %f x %f x %f\n", inLX, inLY, inLZ);
    printf("# ********************** OUTPUT LATTICE **********************\n");
    printf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
    printf("# SPACING: %f x %f x %f\n", DX, DY, DZ);
    printf("# VOLUME : %f x %f x %f\n", LX, LY, LZ);
    
    if(wdmd.datadim>=1 && fabs(LX-inLX)>0.001) {printf("# ERROR: LX FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n"); return( EXIT_FAILURE ) ;}
    if(wdmd.datadim>=2 && fabs(LY-inLY)>0.001) {printf("# ERROR: LY FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n"); return( EXIT_FAILURE ) ;}
    if(wdmd.datadim>=3 && fabs(LZ-inLZ)>0.001) {printf("# ERROR: LZ FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n"); return( EXIT_FAILURE ) ;}
    
    // storage for input variable
    int bdim = wdata_get_blocklength(&wdmd); // get block size
    double *indata;
    cppmallocl(indata, 3*bdim, double);
    
    // prepare set output set
    wdata_metadata wdmdo = wdmd;
    wdata_setNX(&wdmdo,NX); wdata_setNY(&wdmdo,NY); wdata_setNZ(&wdmdo,NZ); 
    wdata_setDX(&wdmdo,DX); wdata_setDY(&wdmdo,DY); wdata_setDZ(&wdmdo,DZ); 
    wdata_setprefix(&wdmdo,argv[2]);
    wdmdo.issetwrkdir=0;
    wdata_clear_database(&wdmdo);
    wdmdo.cycles=wdmd.cycles; // there will be same number of cycles as in dataset-1.wtxt
    wdmdo.issetwrkdir=0;
    
    int bdimo = wdata_get_blocklength(&wdmdo); // get block size
    double *outdata;
    cppmallocl(outdata, 3*bdimo, double);
    
    // interpolate for each variable
    int ivar, icycle;
    for(ivar=0; ivar<wdmdo.nvar; ivar++)
    {
        printf("# INTERPOLATING `%s`...\n", wdmdo.var[ivar].name);
        for(icycle=0; icycle<wdmdo.cycles; icycle++) // for each cycle
        {
            file_operationl( wdata_read_cycle(&wdmd, wdmdo.var[ivar].name, icycle, indata) );
            if     (wdmdo.datadim==3) winterp_interpolation_3d(wdmdo.var[ivar].type[0], inNX, inNY, inNZ, indata, NX, NY, NZ, outdata);
            else if(wdmdo.datadim==2) winterp_interpolation_2d(wdmdo.var[ivar].type[0], inNX, inNY,       indata, NX, NY,     outdata);
            else                      winterp_interpolation_1d(wdmdo.var[ivar].type[0], inNX,             indata, NX,         outdata);
            file_operationl( wdata_write_cycle(&wdmdo, wdmdo.var[ivar].name, outdata) );
        }
    }
    
    // write metadata file (with default name)
    char file_name[256];
    sprintf(file_name, "%s.wtxt", argv[2]);
    printf("# WRITING `%s`\n", file_name);
//     wdata_write_metadata_to_file(&wdmdo, file_name);
    wdata_write_metadata_to_file(&wdmdo, "");
    wdata_add_comment_to_metadata_file(file_name, "Genereted by extension:");
    char cmd[1024];
    sprintf(cmd, "\t%s %s %s", argv[0], argv[1], argv[2]);
    wdata_add_comment_to_metadata_file(file_name,cmd);
    printf("# DONE.\n");
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
