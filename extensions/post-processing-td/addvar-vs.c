 /** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2025
 * 
 * This code adds a new variable "vs" (superfluid velocity) to the wdata set.
 *
 * Compilation:
 *   gcc addvar-vs.c -o addvar-vs -I$WSLDA/lib/wdata/c -L$WSLDA/lib/wdata -lwdata -I$WSLDA/lib/wderiv/c -L$WSLDA/lib/wderiv -lwderiv -lfftw3 -lm
 *
 * */

// wdata lib
#include "wdata.h"
#include "wderiv.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

void print_help(char *pname)
{
    printf("THIS CODE COMPUTES SUPERFLUID VELOCITY AND ADDS IT TO WDATA SET\n");
    printf("Usage: %s -w wtxt -v out-varname -d order_parameter -p precision\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file [REQUIRED]\n");
    printf("\t -v, --var: name of output variable, default: vs\n");
    printf("\t -d, --delta: name of order parameter variable, default: delta\n");
    printf("\t -p, --precision: precision for output: d-double, f-float,\n");
    printf("\t                  default: the same as for the order parameter\n");
    printf("\t -h, --help: print help\n");
}

int main( int argc , char ** argv )
{
    int ierr;
    int c, i;
    char wtxt[256]; 
    char vsname[256]="vs"; // name of output variable
    char deltaname[256]="delta"; // name of order parameter variable
    char precision='x'; // precision for output variable
    int set_wtxt_name=0; // check if wtxt name provided

    while (1)
    {
        static struct option long_options[] =
        {
            /* These options don’t set a flag.
                We distinguish them by their indices. */
            {"wtxt",     required_argument,      0, 'w'},
            {"var",      required_argument,      0, 'v'},
            {"delta",    required_argument,      0, 'd'},
            {"precision",required_argument,      0, 'p'},
            {"help",     no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "w:v:d:p:h",
                        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 'h':
                print_help(argv[0]);
                exit(1);
                break;

            case 'w':
                sprintf(wtxt, "%s", optarg); // copy string
                set_wtxt_name=1;
                break;

            case 'v':
                sprintf(vsname, "%s", optarg); // copy string
                break;

            case 'd':
                sprintf(deltaname, "%s", optarg); // copy string
                break;

            case 'p':
                precision = optarg[0]; // copy char
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort ();
        }
    }

    // check provided parameters
    if(set_wtxt_name==0)
    {
        printf("Error: You must provide --wtxt filename\n");
        print_help(argv[0]);
        exit(1);
    }

    // create metadata handler
    wdata_metadata md;

    // read metadata from file
    printf("# READING METADATA FROM `%s`\n", wtxt);
    ierr = wdata_parse_metadata_file(wtxt, &md);
    if (ierr != 0)
    {
        printf("Cannot read metadata file!\n");
        exit(1);
    }

    // check if order parameter variable exists
    printf("# CHECKING IF VARIABLE `%s` EXISTS...", deltaname); fflush(stdout);
    if(wdata_has_variable(&md, deltaname)==0)
    {
        printf("\nError: Variable `%s` not present in `%s`\n", deltaname, wtxt);
        exit(1);
    }
    printf("OK\n");

    // check if output variable already exists
    printf("# CHECKING IF VARIABLE `%s` ALREADY EXISTS...", vsname); fflush(stdout);
    if(wdata_has_variable(&md, vsname)==1)
    {
        printf("\n# VARIABLE `%s` ALREADY EXISTS IN `%s`\n", vsname, wtxt);
        printf("# PLEASE REMOVE IT FIRST OR USE A DIFFERENT NAME!\n");
        exit(1);
    }
    printf("OK\n");

    // storage for delta and vs
    int bdim = wdata_get_blocklength(&md); // get block size

    double complex *delta;
    delta = (double complex*)malloc(sizeof(double complex) * bdim); // allocate memory
    if (delta == NULL)
    {
        printf("Cannot allocate memory for the order parameter!\n");
        exit(1);
    }

    double *vs;
    vs = (double *)malloc(sizeof(double) * bdim * 3); // allocate memory
    if (vs == NULL)
    {
        printf("Cannot allocate memory for the velocity field!\n");
        exit(1);
    }
    double *vx = vs + 0 * bdim; // pointer to x-component
    double *vy = vs + 1 * bdim; // pointer to y-component
    double *vz = vs + 2 * bdim; // pointer to z-component

    // memory for derivatives
    double complex *delta_dx, *delta_dy, *delta_dz;
    delta_dx = (double complex*)malloc(sizeof(double complex) * bdim); // allocate memory
    delta_dy = (double complex*)malloc(sizeof(double complex) * bdim); // allocate memory
    delta_dz = (double complex*)malloc(sizeof(double complex) * bdim); // allocate memory
    if (delta_dx == NULL || delta_dy == NULL || delta_dz == NULL)
    {
        printf("Cannot allocate memory for the derivatives!\n");
        exit(1);
    }

    // add variable to metadata file
    wdata_variable _vs={"vs", "vector8", "none", "wdat"};
    sprintf(_vs.name, vsname); // set name

    if(precision=='x') // use the same as for delta
    {
        wdata_variable _vdelta;
        wdata_get_variable(&md, deltaname, &_vdelta);
        precision=wdata_get_var_precision(&_vdelta);
    }

    if(precision=='d') // double
        sprintf(_vs.type, "vector8"); // double complex
    else if(precision=='f') // float
        sprintf(_vs.type, "vector4"); // double complex
    else
    {
        printf("Error: Unknown precision `%c`!\n", precision);
        exit(1);
    }

    printf("# CREATING VARIABLE `%s` OF TYPE `%s`\n", _vs.name, _vs.type);
    wdata_add_variable(&md, &_vs);  // add to metadata structure

    // variables for loop
    int icycle, ixyz;
    int idx=0; // cycle index
    double epsilon=1.0e-10; // small number to avoid division by zero

    if(md.datadim==1) // 1D case
    {
        // initialize wderiv lib
        wderiv_init_1d(md.nx, md.dx);

        printf("# COMPUTING SUPERFLUID VELOCITY FOR 1D DATA [NUMBER OF CYCLES: %d]...\n", md.cycles);
        for (icycle = 0; icycle < md.cycles; icycle++) // for each cycle
        {
            // load data
            ierr = wdata_read_cycle_d(&md, deltaname, icycle, (double *)delta);
            if (ierr != 0)
            {
                printf("ERROR: Cannot read %s [ierr=%d]!\n", deltaname, ierr);
                return 1;
            }

            wderiv_gradient_1d_c(delta, delta_dx);

            for (ixyz = 0; ixyz < bdim; ixyz++) // for each lattice point
            {
                // compute vs = (hbar/2m) Im( delta^* * grad(delta) ) / |delta|^2
                // 2m accounts for mass of Cooper pair
                vx[ixyz] = (0.5) * cimag( conj(delta[ixyz])*delta_dx[ixyz] ) / (pow(cabs(delta[ixyz]),2) + epsilon);
                vy[ixyz] = (0.5) * 0.0; // no y -component
                vz[ixyz] = (0.5) * 0.0; // no z -component
            }

            // write to binary file
            ierr = wdata_write_cycle_d(&md, vsname, vs);
            if (ierr != 0)
            {
                printf("ERROR: Cannot write %s [ierr=%d]!\n", vsname, ierr);
                return 1;
            }
        }
    }
    else if(md.datadim==2) // 2D case
    {
        // initialize wderiv lib
        wderiv_init_2d(md.nx, md.ny, md.dx, md.dy);

        printf("# COMPUTING SUPERFLUID VELOCITY FOR 2D DATA [NUMBER OF CYCLES: %d]...\n", md.cycles);
        for (icycle = 0; icycle < md.cycles; icycle++) // for each cycle
        {
            // load data
            ierr = wdata_read_cycle_d(&md, deltaname, icycle, (double *)delta);
            if (ierr != 0)
            {
                printf("ERROR: Cannot read %s [ierr=%d]!\n", deltaname, ierr);
                return 1;
            }

            wderiv_gradient_2d_c(delta, delta_dx, delta_dy);

            for (ixyz = 0; ixyz < bdim; ixyz++) // for each lattice point
            {
                // compute vs = (hbar/2m) Im( delta^* * grad(delta) ) / |delta|^2
                // 2m accounts for mass of Cooper pair
                vx[ixyz] = (0.5) * cimag( conj(delta[ixyz])*delta_dx[ixyz] ) / (pow(cabs(delta[ixyz]),2) + epsilon);
                vy[ixyz] = (0.5) * cimag( conj(delta[ixyz])*delta_dy[ixyz] ) / (pow(cabs(delta[ixyz]),2) + epsilon);
                vz[ixyz] = (0.5) * 0.0; // no z -component
            }

            // write to binary file
            ierr = wdata_write_cycle_d(&md, vsname, vs);
            if (ierr != 0)
            {
                printf("ERROR: Cannot write %s [ierr=%d]!\n", vsname, ierr);
                return 1;
            }
        }
    }
    else if(md.datadim==3) // 3D case
    {
        // initialize wderiv lib
        wderiv_init_3d(md.nx, md.ny, md.nz, md.dx, md.dy, md.dz);

        printf("# COMPUTING SUPERFLUID VELOCITY FOR 3D DATA [NUMBER OF CYCLES: %d]...\n", md.cycles);
        for (icycle = 0; icycle < md.cycles; icycle++) // for each cycle
        {
            // load data
            ierr = wdata_read_cycle_d(&md, deltaname, icycle, (double *)delta);
            if (ierr != 0)
            {
                printf("ERROR: Cannot read %s [ierr=%d]!\n", deltaname, ierr);
                return 1;
            }

            wderiv_gradient_2d_c(delta, delta_dx, delta_dy);

            for (ixyz = 0; ixyz < bdim; ixyz++) // for each lattice point
            {
                // compute vs = (hbar/2m) Im( delta^* * grad(delta) ) / |delta|^2
                // 2m accounts for mass of Cooper pair
                vx[ixyz] = (0.5) * cimag( conj(delta[ixyz])*delta_dx[ixyz] ) / (pow(cabs(delta[ixyz]),2) + epsilon);
                vy[ixyz] = (0.5) * cimag( conj(delta[ixyz])*delta_dy[ixyz] ) / (pow(cabs(delta[ixyz]),2) + epsilon);
                vz[ixyz] = (0.5) * cimag( conj(delta[ixyz])*delta_dz[ixyz] ) / (pow(cabs(delta[ixyz]),2) + epsilon);
            }

            // write to binary file
            ierr = wdata_write_cycle_d(&md, vsname, vs);
            if (ierr != 0)
            {
                printf("ERROR: Cannot write %s [ierr=%d]!\n", vsname, ierr);
                return 1;
            }
        }
    }

    // add variable to metadata file
    printf("# ADDING VARIABLE `%s` TO METADATA FILE `%s`\n", _vs.name, wtxt);
    wdata_add_var_to_metadata_file(wtxt, &_vs); // add entry to metadata file

    printf("# ALL DONE!\n");

    return 0;
}
