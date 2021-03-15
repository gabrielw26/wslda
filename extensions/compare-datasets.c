/**
 * W-SLDA Toolkit
 * 
 * This tool compares two different datasets
 * and creates a new dataset where differences of all variables are saved. 
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 compare-datasets.c -I. -I$WSLDA/hpc-engine -I$WSLDA/lib-wdata -L$WSLDA/lib-wdata -lwdatac -o compare-datasets -lm -lfftw3
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
        printf("Usage: %s dataset-1.wtxt dataset-2.wtxt output-dataset-prefix\n", argv[0]);
        return( EXIT_FAILURE ) ; 
    }
        
    // create metadata handler
    wdata_metadata wdmd1;
    wdata_metadata wdmd2;
    
    // read metadata from file
    printf("# READING `%s`\n", argv[1]);
    file_operationl( wdata_parse_metadata_file(argv[1], &wdmd1) );
    
    printf("# READING `%s`\n", argv[2]);
    file_operationl( wdata_parse_metadata_file(argv[2], &wdmd2) );
    
    // make a copy of the first dataset
    // NOTE: I do not check if dataset-1 and dataset-2 contain the same variables, I assume they do. 
    wdata_metadata wdmdo = wdmd1;
    sprintf(wdmdo.prefix, "%s", argv[3]);
    // just in case - clear data sets
    // it removes binary files if they alredy exists
    wdata_clear_database(&wdmdo);
    wdmdo.cycles=wdmd1.cycles; // there will be same number of cycles as in dataset-1.wtxt
    printf("# CYCLES TO PROCEED: %d\n", wdmd1.cycles);
    
    // memory for buffers
    double *buff1, *buff2, *buffo;
    int bdim = wdata_get_blocklength(&wdmdo); // get block size
    cppmallocl(buff1, bdim*3, double);
    cppmallocl(buff2, bdim*3, double);
    cppmallocl(buffo, bdim*3, double);    
    
    int ivar, ixyz;
    for(ivar=0; ivar<wdmdo.nvar; ivar++)
    {
        int bs = wdata_get_blocksize(&wdmdo, &wdmdo.var[ivar])/sizeof(double); // block size
        printf("# COMPARING `%s` [bs=%d]...\n", wdmdo.var[ivar].name, bs);
        
        int icycle, ixyz;
        for(icycle=0; icycle<wdmdo.cycles; icycle++) // for each cycle
        {
            file_operationl( wdata_read_cycle(&wdmd1, wdmdo.var[ivar].name, icycle, buff1) );
            file_operationl( wdata_read_cycle(&wdmd2, wdmdo.var[ivar].name, icycle, buff2) );    
            for(ixyz=0; ixyz<bs; ixyz++) buffo[ixyz] = buff1[ixyz]-buff2[ixyz];
            file_operationl( wdata_write_cycle(&wdmdo, wdmdo.var[ivar].name, buffo) );
        }
    }
    
    // write metadata file (with default name)
    char file_name[256];
    sprintf(file_name, "%s.wtxt", argv[3]);
    printf("# WRITING `%s`\n", file_name);
    wdata_write_metadata_to_file(&wdmdo, file_name);
    wdata_add_comment_to_metadata_file(file_name, "Genereted by extension:");
    char cmd[1024];
    sprintf(cmd, "\t%s %s %s %s", argv[0], argv[1], argv[2], argv[3]);
    wdata_add_comment_to_metadata_file(file_name,cmd);
    printf("# DONE.\n");
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
