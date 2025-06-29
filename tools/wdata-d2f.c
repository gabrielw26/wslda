/**
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2025
 *
 * A tool for downgrading the storage precision from double to float.
 * */

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "wdata-tools-utils.h"

void print_help(char *pname)
{
#ifdef MODE_FD
    printf("A TOOL FOR UPGRADING THE STORAGE PRECISION FROM FLOAT TO DOUBLE OF GIVEN WDATA SET.\n");
#else
    printf("A TOOL FOR DOWNGRADING THE STORAGE PRECISION FROM DOUBLE TO FLOAT OF GIVEN WDATA SET.\n");
#endif
    printf("Usage: %s -w wtxt -o outprefix ...\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]\n");
    printf("\t -o, --outprefix: prefix for output dataset [REQUIRED]\n");
    printf("\t -h, --help: print help\n");
}

int main(int argc, char **argv)
{
    int c, i;
    char wtxt[1024]; // input
    char wout[1024]; // output
    int set_wtxt=0, set_out=0;

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
            {"help",     no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "w:o:h",
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

#ifdef MODE_FD
    printf("# %s [FLOAT] --> %s.wtxt [DOUBLE]\n", wtxt, wout);
#else
    printf("# %s [DOUBLE] --> %s.wtxt [FLOAT]\n", wtxt, wout);
#endif

    char ctmp1[512], ctmp2[512];
    char basedir[512];
    getcwd(basedir, 512);
    strcpy(ctmp1, wtxt);
    strcpy(ctmp2, wout);
    char *indir = dirname(ctmp1);
    char *outdir = dirname(ctmp2);

    // create metadata handler
    wdata_metadata md;
    wdata_metadata mdout;

    // read metadata from file
    int ierr;
    printf("# READING FILE: `%s`\n", wtxt);
    ierr = wdata_parse_metadata_file(wtxt, &md);
    if (ierr != 0)
    {
        printf("# Cannot read metadata file!\n");
        exit(1);
    }

    mdout = md;
    mdout.cycles = 0;
    strcpy(mdout.prefix, basename(wout));

    // set working dirs
    wdata_set_working_dir(&md, indir);
    wdata_set_working_dir(&mdout, outdir);
    wdata_clear_database(&mdout);

    int icycle, ivar;

    // change type of variables
    for (ivar = 0; ivar < md.nvar; ivar++)
    {
        sprintf(mdout.var[ivar].format, "wdat"); // use wdat format
        char prec = wdata_get_var_precision(&mdout.var[ivar]);
#ifdef MODE_FD
        if(prec=='f') // upgrade
        {
            if(mdout.var[ivar].type[0]=='r')
            {
                sprintf(mdout.var[ivar].type, "real8"); // use wdat format
            }
            else if(mdout.var[ivar].type[0]=='c')
            {
                sprintf(mdout.var[ivar].type, "complex16"); // use wdat format
            }
            else if(mdout.var[ivar].type[0]=='v')
            {
                char _type[32];
                int _dim;
                wdata_vectorvar_decompose(&mdout.var[ivar], _type, &_dim);
                sprintf(mdout.var[ivar].type, "vector8(%d)", _dim); // use wdat format
            }
        }
#else
        if(prec=='d') // downgrade
        {
            if(mdout.var[ivar].type[0]=='r')
            {
                sprintf(mdout.var[ivar].type, "real4"); // use wdat format
            }
            else if(mdout.var[ivar].type[0]=='c')
            {
                sprintf(mdout.var[ivar].type, "complex8"); // use wdat format
            }
            else if(mdout.var[ivar].type[0]=='v')
            {
                char _type[32];
                int _dim;
                wdata_vectorvar_decompose(&mdout.var[ivar], _type, &_dim);
                sprintf(mdout.var[ivar].type, "vector4(%d)", _dim); // use wdat format
            }
        }
#endif
    }

    float *data;
    int bdim = wdata_get_blocklength(&md); // get block size
    cppmallocl(data, bdim * 3, float);

    for (ivar = 0; ivar < md.nvar; ivar++)
    {
        printf("#\t %20s [%10s] --> [%10s]\n", md.var[ivar].name, md.var[ivar].type, mdout.var[ivar].type);

        for (icycle = 0; icycle < md.cycles; icycle++)
        {
            ierr = wdata_read_cycle_f(&md, md.var[ivar].name, icycle, data);
            if (ierr != 0)
            {
                printf("# ERROR: Cannot read datablock!\n");
                exit(1);
            }

            ierr = wdata_write_cycle_f(&mdout, mdout.var[ivar].name, data);
            if (ierr != 0)
            {
                printf("# ERROR: Cannot write datablock!\n");
                exit(1);
            }

            if(ivar == 0) wdata_add_cycle(&mdout);
        }
    }



    // wtxt
    char file_name[512];
    char file_name2[512];

    for (ivar = 0; ivar < md.ntxt; ivar++)
    {
        wdata_get_txt_fullname(&md, md.txt[ivar].filename, file_name);
        wdata_get_txt_fullname(&mdout, md.txt[ivar].filename, file_name2);
        copy_txt(file_name, file_name2);
    }

    // Write new metadata file
    sprintf(file_name, "%s.wtxt", wout);
    printf("# WRITING NEW METADAFILE: `%s`\n", file_name);
    FILE *f = fopen(file_name, "w");
    fprintf(f, "# Generated by: %s\n\n", cmd);
    wdata_print_metadata(&mdout, f);
    fclose(f);

    printf("# DONE.\n");

    return 0;
}
