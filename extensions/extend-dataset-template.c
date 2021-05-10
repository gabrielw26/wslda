/**
 * W-SLDA Toolkit
 * 
 * This file provides template of code that adds new variable 
 * to existing w-dataset. 
 * 
 * In this example, we extend dataset by variable
 *      w_a = j_a / sqrt(rho_a) (vector type)
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 extend-dataset-template.c -I. -I$WSLDA/hpc-engine -I$WSLDA/lib-wdata -L$WSLDA/lib-wdata -lwdatac -o extend-dataset-template -lm -lfftw3 
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
    // Read of input parameters
    char execcmd[ 256 ];
    int i,j,ierr;
    
    strcpy( execcmd , argv[ 0 ] ) ;
    for( i = 1 ; i < argc ; i++ ) 
    {
        strcat( execcmd , " " ) ; 
        strcat( execcmd , argv[ i ] ) ;
    }
    
    // parse command line
    i = readcmd( argc , argv ) ;
    if( i == -1 )
    {
        printf( "TERMINATING! NO INPUT FILE.\n" ) ;
        return( EXIT_FAILURE ) ;
    }
        
    // Read input file
    // Info from file is loaded into metadata structure
    j = parse_input_file(argv[i]);
    if ( j == 0 )
    {
        printf("PROBLEM WITH INPUT FILE: `%s`.\n" , argv[ i ] ) ;
        return( EXIT_FAILURE ) ;      
    }
        
        
    // create metadata handler
    wdata_metadata wdmd;
    
    // read metadata from file
    char file_name[256];
    sprintf(file_name, "%s.wtxt", input->outprefix);
    printf("# READING `%s`\n", file_name);
    file_operationl( wdata_parse_metadata_file(file_name, &wdmd) );
    
    if(wdata_has_variable(&wdmd, "w_a"))
    {
        printf("# DATASET `%s` ALREADY HAS w_a VARIABLE!\n");
        return( EXIT_SUCCESS ) ;
    }

    // storage for real and vector data
    int bdim = wdata_get_blocklength(&wdmd); // get block size
    double *rho; // density
    double *jv;  // current
    double *wv;  // new variable
    cppmallocl(rho, bdim  , double);
    cppmallocl(jv , bdim*3, double);
    cppmallocl(wv , bdim*3, double);
    
    // add variable to metadata file
    wdata_variable w_a = {"w_a", "vector", "none", "wdat"};
    wdata_add_variable(&wdmd, &w_a); // add to metadata structure
    wdata_add_comment_to_metadata_file(file_name, "Entries added by extend-dataset-template code");
    wdata_add_var_to_metadata_file(file_name, &w_a); // add entry to metadata file

    printf("# ADDING NEW VARIABLES TO W-DATASET...\n");
    int icycle, ixyz;
    for(icycle=0; icycle<wdmd.cycles; icycle++) // for each cycle
    {
        // load data
        file_operationl( wdata_read_cycle(&wdmd, "rho_a", icycle, rho) );
        file_operationl( wdata_read_cycle(&wdmd, "j_a", icycle, jv) );
        for(ixyz=0; ixyz<bdim; ixyz++) // for each lattice point
        {
            // new variable
            wv[ixyz+0*bdim] = jv[ixyz+0*bdim]/sqrt(rho[ixyz]+1.0e-16); // x-coordinate
            wv[ixyz+1*bdim] = jv[ixyz+1*bdim]/sqrt(rho[ixyz]+1.0e-16); // y-coordinate
            wv[ixyz+2*bdim] = jv[ixyz+2*bdim]/sqrt(rho[ixyz]+1.0e-16); // z-coordinate
        }
        
        // write to binary file
        file_operationl( wdata_write_cycle(&wdmd, "w_a", wv) );
    }
    printf("# DONE.\n");
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
