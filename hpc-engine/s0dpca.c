/**
 * W-SLDA Toolkit
 *
 * This code solves uniform problem according settings from input file.
 * It can be used for checking of correctness of solver parameters.
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
    wsldapid=0; // global variable

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

#ifdef SPINSYMMETRY_MODE
    md.spinsymmetry=1; // force spin symmetry mode
#endif

    print_version("-0D");
    if(md.spinsymmetry>0)  wprintf("# SPINSYMMETRY MODE IS ACTIVE.\n");

    if(md.ec>0.0) dc_ec = md.ec;
    else          dc_ec = M_PI*M_PI/(2.*DX*DX);
    
#ifdef UNIFORM_TEST_MODE
    md.Na = ceil(1.0/(6.*M_PI*M_PI)*LXYZ);
    md.Nb = md.Na+1;
    if(md.spinsymmetry==1) md.Nb = md.Na;
    md.init0Na = md.Na;
    md.init0Nb = md.Nb;
    wprintf("# UNIFORM_TEST_MODE: Setting number of particles to be: (%f,%f)\n", md.Na,md.Nb);
#endif

    wprintf("# CREATING UNIFORM SOLUTION...\n");

    // Generate uniform solution
    int nwf;

#if FUNCTIONAL==BDG
    solve_uniform_problem_bdg(input->init0Na/LXYZ, input->init0Nb/LXYZ, &nwf, 1);
#elif FUNCTIONAL==SLDAE
    solve_uniform_problem_sldae(input->init0Na/LXYZ, input->init0Nb/LXYZ, &nwf, 1);
#else
    solve_uniform_problem(input->init0Na/LXYZ, input->init0Nb/LXYZ, &nwf, 1);
#endif

    // Save solution
    if(input->init0save)
    {
        create_directory(input->outprefix);
        save_uniform();
    }

    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
