/**
 * W-SLDA Toolkit
 * 
 * This code converts existing wdata set into new one defined on lattice with different resolution. 
 * */

// Standard libraries
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <getopt.h>
#include <libgen.h>
#include <math.h>
#include <complex.h>
#include <time.h>

// W-tools libs
#include "wdata.h"
#include "wdata-tools-utils.h"
#include "winterp.h"

void print_help(char *pname)
{
    printf("WDATA TOOL FOR INTERPOLATING EXISTING DATASET TO NEW RESOLUTION\n");
    printf("Usage: %s -w wtxt -o outprefix ...\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]\n");
    printf("\t -o, --outprefix: prefix for output dataset and new metadata file outprefix.wtxt [REQUIRED]\n");
    printf("\t -x, --nx: number of points along x direction for new representation\n");
    printf("\t           dx will be ajusted such that to keep length along x fixed. Default=input.nx\n");
    printf("\t -y, --ny: number of points along y direction for new representation\n");
    printf("\t           dy will be ajusted such that to keep length along y fixed. Default=input.ny\n");
    printf("\t -z, --nz: number of points along z direction for new representation\n");
    printf("\t           dz will be ajusted such that to keep length along z fixed. Default=input.nz\n");
    printf("\t -h, --help: print help\n");
}

int main(int argc, char **argv)
{

    int ierr;
    int c, i;
    char wtxt[1024]; // input
    char wout[1024]; // output
    int set_wtxt=0, set_out=0;
    int NX = -1;
    int NY = -1;
    int NZ = -1;

    // Construct the execution command
    char cmd[1024] = "";
    for (i = 0; i < argc; i++)
    {
        strcat(cmd, argv[i]);
        if (i < argc - 1)
            strcat(cmd, " ");
    }

    // Get the current date & time
    char date_time[256];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    strftime(date_time, sizeof(date_time), "%Y-%m-%d %H:%M:%S", t);

    while (1)
    {
        static struct option long_options[] =
        {
            /* These options don’t set a flag.
                We distinguish them by their indices. */
            {"wtxt",     required_argument,      0, 'w'},
            {"outprefix",   required_argument,      0, 'o'},
            {"nx",   required_argument,      0, 'x'},
            {"ny",   required_argument,      0, 'y'},
            {"nz",   required_argument,      0, 'z'},
            {"help",     no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "w:o:x:y:z:h",
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
                set_wtxt=1;
                break;

            case 'o':
                sprintf(wout, "%s", optarg); // copy string
                set_out=1;
                break;

            case 'x':
                NX=atoi(optarg);
                break;

            case 'y':
                NY=atoi(optarg);
                break;

            case 'z':
                NZ=atoi(optarg);
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort ();
        }
    }

    // check provided parameters
    if(set_wtxt==0)
    {
        printf("Error: You must provide --wtxt filename\n");
        print_help(argv[0]);
        exit(1);
    }
    if(set_out==0)
    {
        printf("Error: You must provide --outprefix filename\n");
        print_help(argv[0]);
        exit(1);
    }

    // create metadata handler
    wdata_metadata wdmd;

    printf("# READING INPUT DATA: `%s`\n", wtxt);
    file_operationl(wdata_parse_metadata_file(wtxt, &wdmd));

    // Check
    int inNX = wdata_getnx(&wdmd), inNY = wdata_getny(&wdmd), inNZ = wdata_getnz(&wdmd);
    double inDX = wdata_getdx(&wdmd), inDY = wdata_getdy(&wdmd), inDZ = wdata_getdz(&wdmd);
    double inLX = inDX * inNX, inLY = inDY * inNY, inLZ = inDZ * inNZ;
    printf("# **********************  INPUT LATTICE **********************\n");
    printf("# LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
    printf("# SPACING: %f x %f x %f\n", inDX, inDY, inDZ);
    printf("# VOLUME : %f x %f x %f\n", inLX, inLY, inLZ);

    if (NX == -1)
        NX = wdmd.nx;
    if (wdmd.datadim >= 1 && NX <= 1)
    {
        printf("# ERROR: you must provide nx!\n");
        exit(1);
    }

    if (NY == -1)
        NY = wdmd.ny;
    if (wdmd.datadim >= 2 && NY <= 1)
    {
        printf("# ERROR: you must provide ny!\n");
        exit(1);
    }

    if (NZ == -1)
        NZ = wdmd.nz;
    if (wdmd.datadim >= 3 && NZ <= 1)
    {
        printf("# ERROR: you must provide nz!\n");
        exit(1);
    }
    double DX = inLX / NX, DY = inLY / NY, DZ = inLZ / NZ;
    double LX = inLX, LY = inLY, LZ = inLZ;
    printf("# ********************** OUTPUT LATTICE **********************\n");
    printf("# LATTICE: %d x %d x %d\n", NX, NY, NZ);
    printf("# SPACING: %f x %f x %f\n", DX, DY, DZ);
    printf("# VOLUME : %f x %f x %f\n", LX, LY, LZ);

    if (wdmd.datadim >= 1 && fabs(LX - inLX) > 0.001)
    {
        printf("# ERROR: LX FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n");
        exit(1);
    }
    if (wdmd.datadim >= 2 && fabs(LY - inLY) > 0.001)
    {
        printf("# ERROR: LY FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n");
        exit(1);
    }
    if (wdmd.datadim >= 3 && fabs(LZ - inLZ) > 0.001)
    {
        printf("# ERROR: LZ FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n");
        exit(1);
    }

    // storage for input variable
    int bdim = wdata_get_blocklength(&wdmd); // get block size
    double *indata;
    cppmallocl(indata, 3 * bdim, double);

    char ctmp1[512], ctmp2[512];
    char basedir[512];
    getcwd(basedir, 512);
    strcpy(ctmp1, wtxt);
    strcpy(ctmp2, wout);
    char *indir = dirname(ctmp1);
    char *outdir = dirname(ctmp2);
    printf("# WORKING DIR: `%s` --> `%s`\n", indir, outdir);

    // prepare set output set
    wdata_metadata wdmdo = wdmd;
    wdata_setnx(&wdmdo, NX);
    wdata_setny(&wdmdo, NY);
    wdata_setnz(&wdmdo, NZ);
    wdata_setdx(&wdmdo, DX);
    wdata_setdy(&wdmdo, DY);
    wdata_setdz(&wdmdo, DZ);
    wdata_setprefix(&wdmdo, wout);
    wdmdo.issetwrkdir = 0;
    wdata_set_working_dir(&wdmd, indir);
    wdata_set_working_dir(&wdmdo, outdir);
    wdata_clear_database(&wdmdo);
    wdmdo.cycles = wdmd.cycles; // there will be same number of cycles as in dataset-1.wtxt
    wdmdo.issetwrkdir = 0;

    int bdimo = wdata_get_blocklength(&wdmdo); // get block size
    double *outdata;
    cppmallocl(outdata, 3 * bdimo, double);
        
    // interpolate for each variable
    int ivar, icycle;

    for (ivar = 0; ivar < wdmdo.nvar; ivar++)
    {
        printf("# INTERPOLATING `%s`...\n", wdmdo.var[ivar].name);
        for (icycle = 0; icycle < wdmdo.cycles; icycle++) // for each cycle
        {
            file_operationl(wdata_read_cycle_d(&wdmd, wdmdo.var[ivar].name, icycle, indata));
            if (wdmdo.datadim == 3)
                winterp_interpolation_3d(wdmdo.var[ivar].type[0], inNX, inNY, inNZ, indata, NX, NY, NZ, outdata);
            else if (wdmdo.datadim == 2)
                winterp_interpolation_2d(wdmdo.var[ivar].type[0], inNX, inNY, indata, NX, NY, outdata);
            else
                winterp_interpolation_1d(wdmdo.var[ivar].type[0], inNX, indata, NX, outdata);
            file_operationl(wdata_write_cycle_d(&wdmdo, wdmdo.var[ivar].name, outdata));
        }
    }

    // do not copy txt files
    wdmdo.ntxt=0;

    // write metadata file (with default name)
    char file_name[256];
    sprintf(file_name, "%s.wtxt", wout);
    printf("# WRITING `%s`\n", file_name);
    //     wdata_write_metadata_to_file(&wdmdo, file_name);
    wdata_write_metadata_to_file(&wdmdo, "");
    wdata_add_comment_to_metadata_file(file_name, "Genereted by wdata tool:");
    wdata_add_comment_to_metadata_file(file_name, cmd);
    printf("# DONE.\n");

    /* Arrays will be cleared automatically */
    return (EXIT_SUCCESS);
}
