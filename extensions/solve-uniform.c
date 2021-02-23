/**
 * W-SLDA Toolkit
 * 
 * This code solves uniform problem according settings from input file.
 * It can be used for checking of correctness of solver parameters.
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 solve-uniform.c -I. -I$WSLDA/hpc-engine -o solve-uniform -lm -lfftw3
 * 
 * */  

// Standard libraries
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

// W-SLDA Toolkit API
#include "wslda_toolkit.h"

int main( int argc , char ** argv ) 
{
    // Read of input parameters
    char execcmd[ 256 ];
    int i,j;
    
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
        
    printf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
    printf("# SPACING: %f x %f x %f\n", DX, DY, DZ);
#if FUNCTIONAL==BDG
    printf("# ENERGY DENSITY FUNCTIONAL: BDG\n");
#elif FUNCTIONAL==SLDA    
    printf("# ENERGY DENSITY FUNCTIONAL: SLDA\n");
#elif FUNCTIONAL==ASLDA    
    printf("# ENERGY DENSITY FUNCTIONAL: ASLDA\n");     
#endif

    double aBdG;
#if FUNCTIONAL==BDG
    aBdG = input->aBdG; // copy to global momeory
    if ( fabs(aBdG)<1.0e-12 )
    {   
        printf("ERROR: SET aBdG IN INPUT FILE!\n");
        return( EXIT_FAILURE ) ;      
    }
#else
    aBdG = 0.0; // deactivate BdG functional
#endif

    printf("# CREATING UNIFORM SOLUTION...\n");
    
    // Generate uniform solution
    int nwf;
    if(fabs(aBdG)<1.0e-12) solve_uniform_problem    (input->init0Na/LXYZ, input->init0Nb/LXYZ, &nwf, 1);
    else                   solve_uniform_problem_bdg(input->init0Na/LXYZ, input->init0Nb/LXYZ, &nwf, 1);
            
    // Save solution
    if(input->init0save) 
    {
        create_directory(input->outprefix);
        save_uniform();
    }


    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
 
