/**
 * W-SLDA Toolkit
 * 
 * This tool increases dimensionality of wdata set. 
 * Use this code to learn how to use resize routines. 
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 resize-dataset.c -I. -I$WSLDA/hpc-engine -I$WSLDA/lib-wdata -L$WSLDA/lib-wdata -lwdatac -o resize-dataset -lm -lfftw3
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

// W-SLDA Toolkit API
#include "wslda_toolkit.h"

int main( int argc , char ** argv ) 
{
    int ierr;
    
    if(argc!=4)
    {
        printf("Usage: %s input-dataset.wtxt outprefix target-dim\n", argv[0]);
        
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
    int newdim = atoi(argv[3]);
    printf("# ************************ LATTICE ***************************\n");
    printf("# LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
    printf("# SPACING: %f x %f x %f\n", inDX, inDY, inDZ);
    printf("# VOLUME : %f x %f x %f\n", inLX, inLY, inLZ);
    printf("# DIM-IN : %d\n", wdmd.datadim);
    printf("# DIM-OUT: %d\n", newdim);
    
    if(newdim>3) {printf("# ERROR: TARGET DIM >3! ERROR!\n"); return( EXIT_FAILURE ) ;}
    if(wdmd.datadim>newdim) {printf("# ERROR: TARGET DIM < PRESENT DIM! CANNOT DECREASE DIMENSONALITY!\n"); return( EXIT_FAILURE ) ;}
    if(wdmd.datadim==newdim) {printf("# ERROR: TARGET DIM == PRESENT DIM! NO NEED FOR CONVERSION!\n"); return( EXIT_FAILURE ) ;}
    
    // storage for input variable
    int bdim = wdata_get_blocklength(&wdmd); // get block size
    double *indata;
    cppmallocl(indata, 3*bdim, double);
    
    // prepare set output set
    wdata_metadata wdmdo = wdmd;
    wdata_setprefix(&wdmdo,argv[2]);
    wdmdo.issetwrkdir=0;
    wdata_clear_database(&wdmdo);
    wdmdo.cycles=wdmd.cycles; // there will be same number of cycles as in dataset-1.wtxt
    wdmdo.issetwrkdir=0;
    wdmdo.datadim=newdim;
    
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
            if     (wdmd.datadim==1 && wdmdo.datadim==3) wslda_resize_array_1d_to_3d(wdmdo.var[ivar].type[0], inNX, indata, inNY, inNZ, outdata);
            else if(wdmd.datadim==2 && wdmdo.datadim==3) wslda_resize_array_2d_to_3d(wdmdo.var[ivar].type[0], inNX, inNY, indata, inNZ, outdata);
            else                                         wslda_resize_array_1d_to_2d(wdmdo.var[ivar].type[0], inNX, indata, inNY, outdata);
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
