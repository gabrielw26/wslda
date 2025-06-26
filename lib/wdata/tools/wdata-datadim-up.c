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
#include <getopt.h>
#include <math.h>
#include <complex.h>
#include <time.h>

// W-DATA Format
// Must be before wslda_toolkit.h !
#include "wdata.h"

#include "wdata-tools-utils.h"

// W-SLDA Toolkit API
#include "wslda_resize.h"


void print_help(char *pname)
{
    printf("WDATA TOOL: DATADIM MODIFIER\n");
    printf("Usage: %s -w wtxt -o outprefix ...\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]\n");
    printf("\t -o, --outprefix: prefix for output dataset and new metadata file outprefix.wtxt [REQUIRED]\n");
    printf("\t -d, --dim: target datadim of output set, default=input datadim+1\n");
    printf("\t            uniformity along added dimension is assumed.\n");
    printf("\t -h, --help: print help\n");
}

int main( int argc , char ** argv ) 
{
    int ierr;
    int c, i;
    char wtxt[1024]; // input
    char wout[1024]; // output
    int set_wtxt=0, set_out=0;
    int newdim = -1;

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
            {"dim",   required_argument,      0, 'd'},
            {"help",     no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "w:o:d:h",
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

            case 'd':
                newdim=atoi(optarg);
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
    file_operationl( wdata_parse_metadata_file(wtxt, &wdmd) );
    
    // Check
    int inNX=wdata_getnx(&wdmd), inNY=wdata_getny(&wdmd), inNZ=wdata_getnz(&wdmd);
    double inDX=wdata_getdx(&wdmd), inDY=wdata_getdy(&wdmd), inDZ=wdata_getdz(&wdmd);
    double inLX=inDX*inNX, inLY=inDY*inNY, inLZ=inDZ*inNZ;
    if(newdim==-1) newdim = wdmd.datadim+1;
    if(newdim<=0 || newdim>3 || newdim<=wdmd.datadim)
    {
        printf("Error: Incorrect new dimension for the dataset.\n");
        if(newdim<=0)
            printf("Error: dim<=0 !\n");
        else if(newdim>3)
            printf("Error: dim>3 !\n");
        else if(newdim<=wdmd.datadim)
            printf("Error: requested dim=%d, while input.datadim=%d!\n", newdim, wdmd.datadim);
        exit(1);
    }

    printf("# ************************ LATTICE ***************************\n");
    printf("# LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
    printf("# SPACING: %f x %f x %f\n", inDX, inDY, inDZ);
    printf("# VOLUME : %f x %f x %f\n", inLX, inLY, inLZ);
    printf("# DIM-IN : %d\n", wdmd.datadim);
    printf("# DIM-OUT: %d\n", newdim);
    
    // storage for input variable
    int bdim = wdata_get_blocklength(&wdmd); // get block size
    double *indata;
    cppmallocl(indata, 3*bdim, double);
    
    // prepare set output set
    wdata_metadata wdmdo = wdmd;
    wdata_setprefix(&wdmdo,wout);
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
        printf("# RESIZING `%s`...\n", wdmdo.var[ivar].name);
        for(icycle=0; icycle<wdmdo.cycles; icycle++) // for each cycle
        {
            file_operationl( wdata_read_cycle_d(&wdmd, wdmdo.var[ivar].name, icycle, indata) );
            if     (wdmd.datadim==1 && wdmdo.datadim==3) wslda_resize_array_1d_to_3d(wdmdo.var[ivar].type[0], inNX, indata, inNY, inNZ, outdata);
            else if(wdmd.datadim==2 && wdmdo.datadim==3) wslda_resize_array_2d_to_3d(wdmdo.var[ivar].type[0], inNX, inNY, indata, inNZ, outdata);
            else                                         wslda_resize_array_1d_to_2d(wdmdo.var[ivar].type[0], inNX, indata, inNY, outdata);
            file_operationl( wdata_write_cycle_d(&wdmdo, wdmdo.var[ivar].name, outdata) );
        }
    }
    
    // write metadata file (with default name)
    char file_name[512];
    char file_name2[512];
    for(ivar=0; ivar<wdmd.ntxt; ivar++)
    {
        wdata_get_txt_fullname(&wdmd, wdmd.txt[ivar].filename, file_name);
        wdata_get_txt_fullname(&wdmdo, wdmd.txt[ivar].filename, file_name2);
        copy_txt(file_name, file_name2);
    }

    sprintf(file_name, "%s.wtxt", wout);
    printf("# WRITING `%s`\n", file_name);
//     wdata_write_metadata_to_file(&wdmdo, file_name);
    wdata_write_metadata_to_file(&wdmdo, "");
    wdata_add_comment_to_metadata_file(file_name, "Genereted by wdata tool:");
    wdata_add_comment_to_metadata_file(file_name,cmd);
    printf("# DONE.\n");
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
