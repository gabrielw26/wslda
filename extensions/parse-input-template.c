/**
 * W-SLDA Toolkit
 * 
 * This file demonstrate how to parse input file
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 parse-input-template.c -I. -I$WSLDA/hpc-engine -o parse-input-template -lm
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
        
        
    // ...
    // Input file tags are accessible through pointer `input`
    // ...
    printf("outprefix: %s\n", input->outprefix);
    printf("inittype:  %d\n", input->inittype);
    // ... and so on ...

    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
