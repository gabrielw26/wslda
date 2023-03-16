/** 
 * This code adds the phase difference with respect to the selected point in space (corner of the box).
 *
 * gcc -std=gnu99 add-relative-phasediff.c -I$WSLDA/lib/wdata/c -L$WSLDA/lib/wdata -lwdata -o add-relative-phasediff -lm
 * */

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

int main( int argc , char ** argv )
{
    int ierr;

    if(argc!=2)
    {
        printf("Usage: %s wtxt\n", argv[0]);

        return( EXIT_FAILURE ) ;
    }

    // create metadata handler
    wdata_metadata md;

    // read metadata from file
    printf("# OPENING %s.\n", argv[1]);
    ierr = wdata_parse_metadata_file(argv[1], &md);
    if (ierr != 0)
    {
        printf("Cannot read metadata file!\n");
        return 1;
    }

    if(wdata_has_variable(&md, "phase_diff")==1)
    {
        printf("# Variable already exists!\n");
        return 1;
    }

    // storage for real and vector data
    int bdim = wdata_get_blocklength(&md); // get block size

    double complex *delta;
    delta = (double complex*)malloc(sizeof(double complex) * bdim); // allocate memory
    if (delta == NULL)
    {
        printf("Cannot allocate delta!\n");
        return 1;
    }

    // add variable to metadata file
    double *phase_diff;
    phase_diff = (double *)malloc(sizeof(double) * bdim); // allocate memory
    if (phase_diff == NULL)
    {
        printf("Cannot allocate phase_diff!\n");
        return 1;
    }
    wdata_variable pd = {"phase_diff", "real", "none", "wdat"};
    wdata_add_variable(&md, &pd);                     // add to metadata structure
    wdata_add_var_to_metadata_file(argv[1], &pd);     // add entry to metadata file

    int icycle, ixyz;
    for (icycle = 0; icycle < md.cycles; icycle++) // for each cycle
    {
        // load data
        ierr = wdata_read_cycle(&md, "delta", icycle, delta);
        if (ierr != 0)
        {
            printf("ERROR: Cannot read delta [ierr=%d]!\n", ierr);
            return 1;
        }

        for (ixyz = 0; ixyz < bdim; ixyz++) // for each lattice point
        {
            phase_diff[ixyz] = carg(delta[ixyz]/delta[0]); // with respect to the corner
        }

        // write to binary file
        ierr = wdata_write_cycle(&md, "phase_diff", phase_diff);
        if (ierr != 0)
        {
            printf("ERROR: Cannot write phase_diff [ierr=%d]!\n", ierr);
            return 1;
        }
    }

    return 0;
}
