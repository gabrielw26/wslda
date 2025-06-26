/**
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2023
 *
 * Tool for generating section along line for selected variable
 *
 * */

// winterp lib
#include "winterp.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <getopt.h>

// W-DATA Format
// Must be before wslda_toolkit.h !
#include "wdata.h"

#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "ERROR: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        return -1;                                                          \
    }

#define file_operationl(cmd)                                             \
    {                                                                    \
        ierr = cmd;                                                      \
        if (ierr)                                                        \
        {                                                                \
            fprintf(stderr, "FILE ERROR:: cannot execute: %s\n", #cmd);  \
            fprintf(stderr, "file=`%s`, line=%d\n", __FILE__, __LINE__); \
            fprintf(stderr, "Error=%d\nExiting!\n", ierr);               \
            return (EXIT_FAILURE);                                       \
        }                                                                \
    }

void print_help(char *pname)
{
    printf("TOOL FOR GENERATING SECTION ALONG LINE (x1,y1,z1)--(x2,y2,z2) FOR SELECTED VARIABLE\n");
    printf("Usage: %s -w wtxt -o output -v variable --x1 0 --y1 0 --z1 0 --x2 10 --y2 10 --z2 10 -p 100 -c 0\n", pname);
    printf("\t -w, --wtxt: wdata descriptor file [REQUIRED]\n");
    printf("\t -v, --var: variable name [REQUIRED]\n");
    printf("\t -o, --output: name of output file, default: prefix.var.cycle.txt\n");
    printf("\t -x, --x1: x-coordinate for starting point, default=0\n");
    printf("\t -y, --y1: y-coordinate for starting point, default=0, Ignore for 1D data\n");
    printf("\t -z, --z1: z-coordinate for starting point, default=0, Ignore for 1D and 2D data\n");
    printf("\t -i, --x2: x-coordinate for final point, default=nx*dx\n");
    printf("\t -j, --y2: y-coordinate for final point, default=ny*dy, Ignore for 1D data\n");
    printf("\t -k, --z2: z-coordinate for final point, default=nz*dz, Ignore for 1D and 2D data\n");
    printf("\t -p, --points: number of sampling points along line (x1,y1,z1)-(x2,y2,z2), default=100\n");
    printf("\t -c, --cycle: cycle id, default=0. Negative means take form the end, -1 is the last one.\n");
    printf("\t -s, --screen: print on screen the cross-section values\n");
    printf("\t -h, --help: print help\n");
}

char output[256];
int silent_mode;
void wprintf( const char * format, ... )
{
    va_list args;
    if(silent_mode==0)
    {
        va_start (args, format);
        vprintf (format, args);
        va_end (args);
    }

    FILE * f = fopen(output, "a");
    if(f==NULL) { printf("PROBLEM!\n"); fflush(stdout);}
    va_start (args, format);
    vfprintf (f, format, args);
    va_end (args);
    fclose(f);

  fflush(stdout);
}

int main( int argc , char ** argv )
{

    // parse command line: https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html
    int c, ierr,i;
    char wtxt[256];   int set_wtxt_default_name=1;
    /*char output[256];*/ int set_out_default_name=1;
    char var[256];    int set_var_default_name=1;
    double x1=-1.0, y1=-1.0, z1=-1.0;
    double x2=-1.0, y2=-1.0, z2=-1.0;
    int points=100;
    int idx = 0;
    silent_mode=1;

    while (1)
    {
        static struct option long_options[] =
        {
            /* These options don’t set a flag.
                We distinguish them by their indices. */
            {"wtxt",     required_argument,      0, 'w'},
            {"output",   required_argument,      0, 'o'},
            {"cycle",    required_argument,      0, 'c'},
            {"x1",       required_argument,      0, 'x'},
            {"y1",       required_argument,      0, 'y'},
            {"z1",       required_argument,      0, 'z'},
            {"x2",       required_argument,      0, 'i'},
            {"y2",       required_argument,      0, 'j'},
            {"z2",       required_argument,      0, 'k'},
            {"points",   required_argument,      0, 'p'},
            {"var",      required_argument,      0, 'v'},
            {"help",     no_argument,       0, 'h'},
            {"screen",     no_argument,       0, 's'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "w:o:c:x:y:z:i:j:k:p:v:hs",
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

            case 's':
                silent_mode=0;
                break;

            case 'w':
                sprintf(wtxt, "%s", optarg); // copy string
                set_wtxt_default_name=0;
                break;

            case 'o':
                sprintf(output, "%s", optarg); // copy string
                set_out_default_name=0;
                break;

            case 'v':
                sprintf(var, "%s", optarg); // copy string
                set_var_default_name=0;
                break;

            case 'x':
                x1 = atof(optarg);
                break;

            case 'y':
                y1 = atof(optarg);
                break;

            case 'z':
                z1 = atof(optarg);
                break;

            case 'i':
                x2 = atof(optarg);
                break;

            case 'j':
                y2 = atof(optarg);
                break;

            case 'k':
                z2 = atof(optarg);
                break;

            case 'p':
                points = atoi(optarg);
                break;

            case 'c':
                idx = atoi(optarg);
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                abort ();
        }
    }

    // check provided parameters
    if(set_wtxt_default_name)
    {
        printf("Error: You must provide --wtxt filename\n");
        print_help(argv[0]);
        exit(1);
    }
    if(set_var_default_name)
    {
        printf("Error: You must provide --var variable\n");
        print_help(argv[0]);
        exit(1);
    }



    // create metadata handler
    wdata_metadata wdmd;
    file_operationl( wdata_parse_metadata_file(wtxt, &wdmd) );

    // Check
    int inNX=wdata_getnx(&wdmd), inNY=wdata_getny(&wdmd), inNZ=wdata_getnz(&wdmd);
    double inDX=wdata_getdx(&wdmd), inDY=wdata_getdy(&wdmd), inDZ=wdata_getdz(&wdmd);
    double inLX=inDX*inNX, inLY=inDY*inNY, inLZ=inDZ*inNZ;
    double eF = wdata_getconst_value(&wdmd, "eF");
    double kF = wdata_getconst_value(&wdmd, "kF");
    if(idx<0) idx=wdmd.cycles+idx; // take from the end

    if(wdata_has_variable(&wdmd, var)==0)
    {
        printf("Error: Variable `%s` not present in `%s`\n", var, wtxt);
        exit(1);
    }

    if(set_out_default_name)
    {
        sprintf(output,"%s.%s.%06d.txt", wdmd.prefix, var, idx);
    }
    // create empty file
    printf("# CREATING OUTPUT FILE `%s`\n", output);
    FILE *outf = fopen(output, "w");
    if(outf==NULL)
    {
        printf("Error: Cannot create output file `%s`\n", output);
        exit(1);
    }
    fclose(outf);
    wprintf("# INPUT DATA: `%s`\n", wtxt);
    wprintf("# ************************ LATTICE ***************************\n");
    wprintf("# LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
    wprintf("# SPACING: %f x %f x %f\n", inDX, inDY, inDZ);
    wprintf("# VOLUME : %f x %f x %f\n", inLX, inLY, inLZ);
    wprintf("# DIM    : %d\n", wdmd.datadim);
    wprintf("# CYCLES : %d\n", wdmd.cycles);
    wprintf("# ********************* CROSS-SECTION *************************\n");
    wprintf("# CYCLE: %d\n", idx);
    wprintf("# POINTS: %d\n", points);

    if(x1<0.0) x1=0.0;
    if(y1<0.0) y1=0.0;
    if(z1<0.0) z1=0.0;
    if(x2<0.0) x2=wdmd.nx*wdmd.dx;
    if(y2<0.0) y2=wdmd.ny*wdmd.dy;
    if(z2<0.0) z2=wdmd.nz*wdmd.dz;
    double dx = (x2-x1)/points;
    double dy = (y2-y1)/points;
    double dz = (z2-z1)/points;
    double x,y,z;

    if(wdmd.datadim==3) // 3D-data
    {
        wdata_variable myvar;
        wdata_get_variable(&wdmd, var, &myvar);
        int bdim=inNX*inNY*inNZ; // buffer size

        wprintf("# 3D SECTION: (%f,%f,%f) -- (%f,%f,%f)\n", x1,y1,z1,x2,y2,z2);
        wprintf("# VARIBALE: %s [%s]\n", var,myvar.type);
        wprintf("# COLUMNS: \n");
        wprintf("# 1: x\n");
        wprintf("# 2: y\n");
        wprintf("# 3: z\n");

        if(myvar.type[0]=='r') // real
        {
            wprintf("# 4: %s(x,y,z)\n", var);
            printf("# LOADING: %s\n", var);
            double *dreal;
            cppmallocl(dreal, bdim, double);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dreal) );
            winterp_interpolator idreal;
            winterp_create_interpolator_3d_r(inNX, inNY, inNZ, inDX, inDY, inDZ, dreal, &idreal);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                y=y1+dy*i;
                z=z1+dz*i;
                double vreal=0;
                winterp_getvalue_3d_r(&idreal, x,y,z, &vreal);
                wprintf("%16.8g %16.8g %16.8g %16.8g\n", x,y,z,vreal);
            }
        }

        if(myvar.type[0]=='c') // complex
        {
            wprintf("# 4: Re[%s(x,y,z)]\n", var);
            wprintf("# 5: Im[%s(x,y,z)]\n", var);
            wprintf("# 6: Abs[%s(x,y,z)]\n", var);
            wprintf("# 7: Arg[%s(x,y,z)]\n", var);
            printf("# LOADING: %s\n", var);
            double complex *dcomplex;
            cppmallocl(dcomplex, bdim, double complex);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dcomplex) );
            winterp_interpolator idcomplex;
            winterp_create_interpolator_3d_c(inNX, inNY, inNZ, inDX, inDY, inDZ, dcomplex, &idcomplex);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                y=y1+dy*i;
                z=z1+dz*i;
                double complex vcomplex;
                winterp_getvalue_3d_c(&idcomplex, x,y,z, &vcomplex);
                wprintf("%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", x,y,z,creal(vcomplex),cimag(vcomplex),cabs(vcomplex),carg(vcomplex));
            }
        }

        if(myvar.type[0]=='v') // vector type
        {
            wprintf("# 4: %s_x(x,y,z)\n", var);
            wprintf("# 5: %s_y(x,y,z)\n", var);
            wprintf("# 6: %s_z(x,y,z)\n", var);
            wprintf("# 7: magnitude[%s(x,y,z)]\n", var);
            printf("# LOADING: %s\n", var);
            double *dvector;
            cppmallocl(dvector, bdim*3, double);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dvector) );
            winterp_interpolator idvector;
            winterp_create_interpolator_3d_v(inNX, inNY, inNZ, inDX, inDY, inDZ, dvector, &idvector);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                y=y1+dy*i;
                z=z1+dz*i;
                double vvector_x=0;
                double vvector_y=0;
                double vvector_z=0;
                winterp_getvalue_3d_v(&idvector, x,y,z, &vvector_x,&vvector_y,&vvector_z);
                wprintf("%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", x,y,z,vvector_x,vvector_y,vvector_z,
                    sqrt(pow(vvector_x,2)+pow(vvector_y,2)+pow(vvector_z,2))
                );
            }
        }
    }

    if(wdmd.datadim==2) // 2D-data
    {
        wdata_variable myvar;
        wdata_get_variable(&wdmd, var, &myvar);
        int bdim=inNX*inNY; // buffer size

        wprintf("# 2D SECTION: (%f,%f) -- (%f,%f)\n", x1,y1,x2,y2);
        wprintf("# VARIBALE: %s [%s]\n", var,myvar.type);
        wprintf("# COLUMNS: \n");
        wprintf("# 1: x\n");
        wprintf("# 2: y\n");

        if(myvar.type[0]=='r') // real
        {
            wprintf("# 3: %s(x,y)\n", var);
            printf("# LOADING: %s\n", var);
            double *dreal;
            cppmallocl(dreal, bdim, double);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dreal) );
            winterp_interpolator idreal;
            winterp_create_interpolator_2d_r(inNX, inNY, inDX, inDY, dreal, &idreal);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                y=y1+dy*i;
                double vreal=0;
                winterp_getvalue_2d_r(&idreal, x,y, &vreal);
                wprintf("%16.8g %16.8g %16.8g\n", x,y,vreal);
            }
        }

        if(myvar.type[0]=='c') // complex
        {
            wprintf("# 3: Re[%s(x,y)]\n", var);
            wprintf("# 4: Im[%s(x,y)]\n", var);
            wprintf("# 5: Abs[%s(x,y)]\n", var);
            wprintf("# 6: Arg[%s(x,y)]\n", var);
            printf("# LOADING: %s\n", var);
            double complex *dcomplex;
            cppmallocl(dcomplex, bdim, double complex);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dcomplex) );
            winterp_interpolator idcomplex;
            winterp_create_interpolator_2d_c(inNX, inNY, inDX, inDY, dcomplex, &idcomplex);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                y=y1+dy*i;
                double complex vcomplex;
                winterp_getvalue_2d_c(&idcomplex, x,y, &vcomplex);
                wprintf("%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", x,y,creal(vcomplex),cimag(vcomplex),cabs(vcomplex),carg(vcomplex));
            }
        }

        if(myvar.type[0]=='v') // vector type
        {
            wprintf("# 3: %s_x(x,y)\n", var);
            wprintf("# 4: %s_y(x,y)\n", var);
            wprintf("# 5: %s_z(x,y)\n", var);
            wprintf("# 6: magnitude[%s(x,y)]\n", var);
            printf("# LOADING: %s\n", var);
            double *dvector;
            cppmallocl(dvector, bdim*3, double);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dvector) );
            winterp_interpolator idvector;
            winterp_create_interpolator_2d_v(inNX, inNY, inDX, inDY, dvector, &idvector);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                y=y1+dy*i;
                double vvector_x=0;
                double vvector_y=0;
                double vvector_z=0;
                winterp_getvalue_2d_v(&idvector, x,y, &vvector_x,&vvector_y,&vvector_z);
                wprintf("%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", x,y,vvector_x,vvector_y,vvector_z,
                    sqrt(pow(vvector_x,2)+pow(vvector_y,2)+pow(vvector_z,2))
                );
            }
        }
    }

    if(wdmd.datadim==1) // 1D-data
    {
        wdata_variable myvar;
        wdata_get_variable(&wdmd, var, &myvar);
        int bdim=inNX; // buffer size

        wprintf("# 1D SECTION: (%f) -- (%f)\n", x1,x2);
        wprintf("# VARIBALE: %s [%s]\n", var,myvar.type);
        wprintf("# COLUMNS: \n");
        wprintf("# 1: x\n");

        if(myvar.type[0]=='r') // real
        {
            wprintf("# 2: %s(x)\n", var);
            printf("# LOADING: %s\n", var);
            double *dreal;
            cppmallocl(dreal, bdim, double);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dreal) );
            winterp_interpolator idreal;
            winterp_create_interpolator_1d_r(inNX, inDX, dreal, &idreal);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                double vreal=0;
                winterp_getvalue_1d_r(&idreal, x, &vreal);
                wprintf("%16.8g %16.8g\n", x,vreal);
            }
        }

        if(myvar.type[0]=='c') // complex
        {
            wprintf("# 2: Re[%s(x)]\n", var);
            wprintf("# 3: Im[%s(x)]\n", var);
            wprintf("# 4: Abs[%s(x)]\n", var);
            wprintf("# 5: Arg[%s(x)]\n", var);
            printf("# LOADING: %s\n", var);
            double complex *dcomplex;
            cppmallocl(dcomplex, bdim, double complex);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dcomplex) );
            winterp_interpolator idcomplex;
            winterp_create_interpolator_1d_c(inNX, inDX, dcomplex, &idcomplex);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                double complex vcomplex;
                winterp_getvalue_1d_c(&idcomplex, x, &vcomplex);
                wprintf("%16.8g %16.8g %16.8g %16.8g %16.8g\n", x,creal(vcomplex),cimag(vcomplex),cabs(vcomplex),carg(vcomplex));
            }
        }

        if(myvar.type[0]=='v') // vector type
        {
            wprintf("# 2: %s_x(x)\n", var);
            wprintf("# 3: %s_y(x)\n", var);
            wprintf("# 4: %s_z(x)\n", var);
            wprintf("# 5: magnitude[%s(x)]\n", var);
            printf("# LOADING: %s\n", var);
            double *dvector;
            cppmallocl(dvector, bdim*3, double);
            file_operationl( wdata_read_cycle(&wdmd, var, idx, dvector) );
            winterp_interpolator idvector;
            winterp_create_interpolator_1d_v(inNX, inDX, dvector, &idvector);

            for(i=0; i<=points; i++) // iterate over points
            {
                x=x1+dx*i;
                double vvector_x=0;
                double vvector_y=0;
                double vvector_z=0;
                winterp_getvalue_1d_v(&idvector, x, &vvector_x,&vvector_y,&vvector_z);
                wprintf("%16.8g %16.8g %16.8g %16.8g %16.8g\n", x,vvector_x,vvector_y,vvector_z,
                    sqrt(pow(vvector_x,2)+pow(vvector_y,2)+pow(vvector_z,2))
                );
            }
        }
    }
    printf("# OUPUT WRITTEN TO `%s`.\n", output);
    printf("# DONE.\n");

    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
