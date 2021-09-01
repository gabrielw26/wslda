/**
 * W-Data Tools
 * 
 * This tool increases dimensionality of wdata set. 
 * Uniformiyty alone new directions is assumed.
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
#include "wslda_resize.h"

#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "ERROR: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        return -1;                                                          \
    }
    
#define file_operationl(cmd)                                             \
    {                                                                    \
        ierr = cmd;                                                      \
        if (ierr)                                                        \
        {                                                                \
            fprintf(stderr, "FILE ERROR:: cannot execute: %s\n", #cmd);  \
            fprintf(stderr, "file=`%s`, line=%d\n", __FILE__, __LINE__); \
            fprintf(stderr, "Error=%d\nExiting!\n", ierr);               \
            return (EXIT_FAILURE);                                       \
        }                                                                \
    }

int main( int argc , char ** argv ) 
{
    int ierr;
    printf("# WDATA TOOL: DATADIM MODIFIER\n");
    if(argc!=4)
    {
        printf("Usage: %s input.wtxt outprefix target-dim\n", argv[0]);
        printf("\tinput.wtxt    -metadata file (INPUT)\n");
        printf("\toutprefix    - new metadata file will be written to outprefix.wtxt (OUTPUT)\n");
        printf("\ttarget-dim   - target datadim of output set,\n");
        printf("\t               uniformity along added dimension is assumed.\n");
        return( EXIT_FAILURE ) ; 
    }
    
    // create metadata handler
    wdata_metadata wdmd;
    
    printf("# READING INPUT DATA: `%s`\n", argv[1]);
    file_operationl( wdata_parse_metadata_file(argv[1], &wdmd) );
    
    // Check
    int inNX=wdata_getnx(&wdmd), inNY=wdata_getny(&wdmd), inNZ=wdata_getnz(&wdmd);
    double inDX=wdata_getdx(&wdmd), inDY=wdata_getdy(&wdmd), inDZ=wdata_getdz(&wdmd);
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
    
    for (ivar = 0; ivar < wdmdo.nvar; ivar++)
    {
        if( strcmp(wdmdo.var[ivar].type,"real4") == 0 || strcmp(wdmdo.var[ivar].type,"complex8") == 0 || strcmp(wdmdo.var[ivar].type,"vector4") == 0)
        {
            printf("ERROR: Unsupported type %s for variable %s!\n", wdmdo.var[ivar].type, wdmdo.var[ivar].name);
            return (EXIT_FAILURE);
        }
    }
    
    for(ivar=0; ivar<wdmdo.nvar; ivar++)
    {
        printf("# RESIZING `%s`...\n", wdmdo.var[ivar].name);
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
    wdata_add_comment_to_metadata_file(file_name, "Genereted by wdata tool:");
    char cmd[1024];
    sprintf(cmd, "\t%s %s %s", argv[0], argv[1], argv[2]);
    wdata_add_comment_to_metadata_file(file_name,cmd);
    printf("# DONE.\n");
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
