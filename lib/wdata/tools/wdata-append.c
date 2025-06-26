/**
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 *
 * This tool appends a wdata set to the existing one.
 * */

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <libgen.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "wdata-tools-utils.h"


void print_help(char *pname)
{
    printf("TOOL FOR APPENDING WDATA SET TO EXISTING ONE\n");
    printf("Usage: %s -w wtxt -a wtxt -t\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file to which new data will be appended  [REQUIRED]\n");
    printf("\t -a, --append: wdata descriptor file which will be appended [REQUIRED]\n");
    printf("\t -t, --text: update text files, like wlog, stdout\n");
    printf("\t -e, --eps: tolerance between least entry times in wdata sets, default=1.0e-6\n");
    printf("\t -y, --yes: apply yes answer to prompt\n");
    printf("\t -h, --help: print help\n");
}

void append_files(const char *filename1, const char *filename2) {
    FILE *file1, *file2;
    char buffer[1024]; // Buffer to hold lines of text
    printf("# APPENDING: %s <-- %s\n", filename1, filename2);

    // Open filename1 in append mode
    file1 = fopen(filename1, "a");
    if (file1 == NULL) {
        printf("# COULD NOT OPEN FILE %s FOR APPENDING. SKIPPIG!\n", filename1);
        return;
    }

    // Open filename2 in read mode
    file2 = fopen(filename2, "r");
    if (file2 == NULL) {
        printf("# COULD NOT OPEN FILE %s FOR READING. SKIPPING!\n", filename2);
        fclose(file1); // Ensure file1 is closed before exiting
        return;
    }

    // Add note
    fprintf(file1, "# >>>> --------- APPENDED %s BY wdata-append TOOL --------- <<<<\n", filename2);

    // Read content from filename2 and append it to filename1
    while (fgets(buffer, sizeof(buffer), file2) != NULL) {
        fputs(buffer, file1);
    }

    // Close both files
    fclose(file1);
    fclose(file2);

    // printf("# CONTENT OF %s HAS BEEN APPENDED TO %s.\n", filename2, filename1);

    return;
}

void replace_suffix(const char *filename, const char *suffix, char *result, size_t result_size) {
    const char *target = ".wtxt"; // Substring to replace
    char *position;

    // Find the position of ".wtxt" in the filename
    position = strstr(filename, target);
    if (position == NULL) {
        // If ".wtxt" is not found, copy the original filename to the result
        snprintf(result, result_size, "%s", filename);
        return;
    }

    // Calculate the length of the part before ".wtxt"
    size_t prefix_length = position - filename;

    // Construct the new string with the suffix
    snprintf(result, result_size, "%.*s%s", (int)prefix_length, filename, suffix);
}

int main(int argc, char **argv)
{
    int c, ierr,i;
    int ivar, icycle;
    char wtxt[1024]; // input
    char wapp[1024]; // to be appended
    int utext=0, set_wtxt=0, set_append=0;
    double eps=1.0e-6;
    int prompt_yes=0;

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
            {"append",   required_argument,      0, 'a'},
            {"eps",   required_argument,      0, 'e'},
            {"text",     no_argument,       0, 't'},
            {"yes",     no_argument,       0, 'y'},
            {"help",     no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "w:a:e:tyh",
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

            case 'a':
                sprintf(wapp, "%s", optarg); // copy string
                set_append=1;
                break;

            case 't':
                utext=1;
                break;

            case 'y':
                prompt_yes=1;
                break;

            case 'e':
                eps=atof(optarg);
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
    if(set_append==0)
    {
        printf("Error: You must provide --append filename\n");
        print_help(argv[0]);
        exit(1);
    }

    printf("# WDATA SET THAT WILL BE UPDATED : %s\n", wtxt);
    printf("# WDATA SET THAT WILL BE APPENDED: %s\n", wapp);
    // create metadata handler
    wdata_metadata mdw; // wtxt
    wdata_metadata mda; // append

    // read metadata from file
    ierr = wdata_parse_metadata_file(wtxt, &mdw);
    if (ierr != 0)
    {
        printf("# Cannot read metadata file %s!\n", wtxt);
        exit(1);
    }

    ierr = wdata_parse_metadata_file(wapp, &mda);
    if (ierr != 0)
    {
        printf("# Cannot read metadata file %s!\n", wapp);
        exit(1);
    }

    // make checks before starting updates
    if(mdw.nx!=mda.nx || mdw.ny!=mda.ny || mdw.nz!=mda.nz || mdw.dx!=mda.dx || mdw.dy!=mda.dy || mdw.dz!=mda.dz || mdw.datadim!=mda.datadim)
    {
        printf("# ERROR: Inconsistent lattices!\n");
        exit(1);
    }
    else printf("# CHECK: LATTICES ARE CONSISTENT...\n");

    if(mdw.dt!=mda.dt)
    {
        printf("# ERROR: Inconsistent dt!\n");
        exit(1);
    }
    else printf("# CHECK: TIMESTEPS ARE CONSISTENT...\n");

    double timew0=mdw.t0; // first frame
    double timew1=mdw.t0+mdw.dt*(mdw.cycles-1); // last frame
    double timea0=mda.t0; // first frame
    double timea1=mda.t0+mda.dt*(mda.cycles-1); // last frame
    printf("# CHECK: %s: TIME IN [%f,%f]\n", wtxt, timew0, timew1);
    printf("# CHECK: %s: TIME IN  [%f,%f]\n", wapp, timea0, timea1);
    if(timea0-timew1>mda.dt+eps)
    {
        printf("# ERROR: Last time frame in %s differs from first time frame %s by more than time step dt=%f!\n", wtxt, wapp, mda.dt);
        exit(1);
    }
    else printf("# CHECK: TIME INTERVALS ARE CONSISTENT...\n");

    for (ivar = 0; ivar < mdw.nvar; ivar++)
    {
        int testv=wdata_has_variable(&mda, mdw.var[ivar].name);
        if(testv==0)
        {
            printf(" ERROR: %s does not contain variable %s!\n", wapp, mdw.var[ivar].name);
            exit(1);
        }
    }
    printf("# CHECK: VARIABLES ARE CONSISTENT...\n");
    printf("# CHECKS ARE DONE!\n");

    // make simulation of appending
    double ltime;
    int start_cycle=-1;
    for (icycle = 0; icycle < mda.cycles; icycle++)
    {
        wdata_get_time(&mda, icycle, &ltime);
        if(ltime>timew1+eps)
        {
            start_cycle=icycle;
            break;
        }
    }

    if(start_cycle==-1)
    {
        printf("# NOTHING TO BE APPENDED FROM %s!\n", wapp);
        exit(1);
    }

    printf("# TO ORIGIN       : CYCLES %8d, TIME IN [%12.8g,%12.8g]\n", mdw.cycles, timew0, timew1);
    printf("# WILL BE APPENDED: CYCLES %8d, TIME IN [%12.8g,%12.8g]  (start_cycle=%d)\n", mda.cycles-start_cycle, timea0+start_cycle*mda.dt, timea1,start_cycle);
    printf("# AFTER UPDATE    : CYCLES %8d, TIME IN [%12.8g,%12.8g]\n", mdw.cycles+mda.cycles-start_cycle, timew0, timea1);

    if(prompt_yes==0)
    {
        // Print the prompt
        printf("# DO YOU WANT TO CONTINUE: y/[n]? ");
        char answer[10]; // Buffer to store user input (large enough for safety)
        char default_answer = 'n'; // Default value is 'n'

        // Read the user's input
        if (fgets(answer, sizeof(answer), stdin) != NULL) {
            // Remove the newline character if present
            answer[strcspn(answer, "\n")] = '\0';

            // Check if the user just pressed Enter
            if (strlen(answer) == 0) {
                answer[0] = default_answer; // Set default to 'n'
                answer[1] = '\0';          // Null-terminate the string
            }
        }

        // Check the user's answer
        if (answer[0] != 'y' && answer[0] != 'Y')
        {
            printf("# TERMINATING!\n");
            exit(1);
        }
    }
    printf("# APPENDING TO WDAT FILES...\n");

    // allocate memory
    double *data;
    int bdim = wdata_get_blocklength(&mdw); // get block size
    cppmallocl(data, bdim * 3, double);

    for (icycle = start_cycle; icycle < mda.cycles; icycle++)
    {
        for (ivar = 0; ivar < mdw.nvar; ivar++)
        {
            ierr = wdata_read_cycle(&mda, mdw.var[ivar].name, icycle, data);
            if (ierr != 0)
            {
                printf("# ERROR: Cannot read datablock!\n");
                return 1;
            }

            ierr = wdata_write_cycle(&mdw, mdw.var[ivar].name, data);
            if (ierr != 0)
            {
                printf("# ERROR: Cannot write datablock!\n");
                return 1;
            }
        } // ivar
        wdata_add_cycle(&mdw);
    } // icycle

    // Write new metadata file
    printf("# UPDATING METADAFILE: `%s`\n", wtxt);
    FILE *f = fopen(wtxt, "w");
    wdata_print_metadata(&mdw, f);
    fprintf(f, "# THIS WDATA SET HAS BEEN UPDATED BY wdata-append TOOL.\n");
    fprintf(f, "# EXECUTION COMMAND    : %s\n", cmd);
    fprintf(f, "# EXECUTION DATE & TIME: %s\n", date_time);
    fclose(f);

    if(utext==1)
    {
        char resultw[1024]; // Buffer to store the result
        char resulta[1024]; // Buffer to store the result

        replace_suffix(wtxt, "_input.txt", resultw, sizeof(resultw));
        replace_suffix(wapp, "_input.txt", resulta, sizeof(resulta));
        append_files(resultw, resulta);

        replace_suffix(wtxt, ".wlog", resultw, sizeof(resultw));
        replace_suffix(wapp, ".wlog", resulta, sizeof(resulta));
        append_files(resultw, resulta);

        replace_suffix(wtxt, ".stdout", resultw, sizeof(resultw));
        replace_suffix(wapp, ".stdout", resulta, sizeof(resulta));
        append_files(resultw, resulta);

        replace_suffix(wtxt, "_check.stamp", resultw, sizeof(resultw));
        replace_suffix(wapp, "_check.stamp", resulta, sizeof(resulta));
        append_files(resultw, resulta);
    }

    printf("# DONE.\n");

    return 0;
}
