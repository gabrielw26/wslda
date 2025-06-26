/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Tool for extracting subset from wdata set
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
#include "wdata-subset-fun.c"

void print_help_copy(char *pname)
{
    printf("TOOL FOR COPYING WDATA SET\n");
    printf("Usage: %s -w wtxt -o outprefix ...\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]\n");
    printf("\t -o, --outprefix: prefix for output dataset [REQUIRED]\n");
    printf("\t -h, --help: print help\n");
}

void print_help_subset(char *pname)
{
    printf("TOOL FOR SUBSET EXTRACTING FROM WDATA SET\n");
    printf("Usage: %s -w wtxt -o outprefix ...\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file if input dataset [REQUIRED]\n");
    printf("\t -o, --outprefix: prefix for output dataset [REQUIRED]\n");
    printf("\t -b, --begin: first cycle of subset to extraction, default=0\n");
    printf("\t -e, --end: last cycle of subset to extraction, default=-1 (the last cycle)\n");
    printf("\t -s, --stride: every `stride` frame will be taken only, default=1\n");
    printf("\t -h, --help: print help\n");
}

void print_help(char *pname)
{
#ifdef WDATA_SUBSET_MODE
    print_help_subset(pname);
#else
    print_help_copy(pname);
#endif
}

int main(int argc, char **argv)
{
    int c, i;
    char wtxt[1024]; // input
    char wout[1024]; // output
    int set_wtxt=0, set_out=0;
    int begin=0;
    int end=-1;
    int stride=1;

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
#ifdef WDATA_SUBSET_MODE
            {"begin",   required_argument,      0, 'b'},
            {"end",   required_argument,      0, 'e'},
            {"stride",   required_argument,      0, 's'},
#endif
            {"help",     no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "w:o:b:e:s:h",
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

#ifdef WDATA_SUBSET_MODE
            case 'b':
                begin=atoi(optarg);
                break;

            case 'e':
                end=atoi(optarg);
                break;

            case 's':
                stride=atoi(optarg);
                break;
#endif

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

    wdata_generic_cut(cmd, wtxt, wout, begin, end, stride);

    return 0;
}
