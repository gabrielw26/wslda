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
#include <libgen.h>
#include <math.h>
#include <complex.h>

// W-tools libs
#include "wdata.h"
#include "winterp.h"

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

int main(int argc, char **argv)
{

    int ierr;
    printf("# WDATA SET INTERPOLATOR\n");

    if (argc != 6 && argc != 5 && argc != 4)
    {
        printf("Usage: %s input.wtxt outprefix new-nx new-ny new-nz\n", argv[0]);
        printf("\tinput.wtxt    -metadata file (INPUT)\n");
        printf("\toutprefix    - new metadata file will be written to outprefix.wtxt (OUTPUT)\n");
        printf("\tnew-nx       - number of points along x direction for new representation (nx*dx is fixed)\n");
        printf("\tnew-ny       - number of points along y direction for new representation (ny*dy is fixed)\n");
        printf("\t               required only if datadim>=2\n");
        printf("\tnew-nz       - number of points along z direction for new representation (nz*dz is fixed)\n");
        printf("\t               required only if datadim =3\n");
        return 0;
    }

    // create metadata handler
    wdata_metadata wdmd;

    printf("# READING INPUT DATA: `%s`\n", argv[1]);
    file_operationl(wdata_parse_metadata_file(argv[1], &wdmd));

    // Check
    int inNX = wdata_getnx(&wdmd), inNY = wdata_getny(&wdmd), inNZ = wdata_getnz(&wdmd);
    double inDX = wdata_getdx(&wdmd), inDY = wdata_getdy(&wdmd), inDZ = wdata_getdz(&wdmd);
    double inLX = inDX * inNX, inLY = inDY * inNY, inLZ = inDZ * inNZ;
    printf("# **********************  INPUT LATTICE **********************\n");
    printf("# LATTICE: %d x %d x %d\n", inNX, inNY, inNZ);
    printf("# SPACING: %f x %f x %f\n", inDX, inDY, inDZ);
    printf("# VOLUME : %f x %f x %f\n", inLX, inLY, inLZ);

    int NX = 1;
    if (argc >= 4)
        NX = atoi(argv[3]);
    if (wdmd.datadim >= 1 && NX == 1)
    {
        printf("# ERROR: you must provide new-nx!\n");
        return 0;
    }
    int NY = 1;
    if (argc >= 5)
        NY = atoi(argv[4]);
    if (wdmd.datadim >= 2 && NY == 1)
    {
        printf("# ERROR: you must provide new-ny!\n");
        return 0;
    }
    int NZ = 1;
    if (argc >= 6)
        NZ = atoi(argv[5]);
    if (wdmd.datadim >= 3 && NZ == 1)
    {
        printf("# ERROR: you must provide new-nz!\n");
        return 0;
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
        return (EXIT_FAILURE);
    }
    if (wdmd.datadim >= 2 && fabs(LY - inLY) > 0.001)
    {
        printf("# ERROR: LY FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n");
        return (EXIT_FAILURE);
    }
    if (wdmd.datadim >= 3 && fabs(LZ - inLZ) > 0.001)
    {
        printf("# ERROR: LZ FOR INPUT AND TARGET LATTICE INCOMPATIBLE!\n");
        return (EXIT_FAILURE);
    }

    // storage for input variable
    int bdim = wdata_get_blocklength(&wdmd); // get block size
    double *indata;
    cppmallocl(indata, 3 * bdim, double);

    char ctmp1[512], ctmp2[512];
    char basedir[512];
    getcwd(basedir, 512);
    strcpy(ctmp1, argv[1]);
    strcpy(ctmp2, argv[2]);
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
    wdata_setprefix(&wdmdo, argv[2]);
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
        if( strcmp(wdmdo.var[ivar].type,"real4") == 0 || strcmp(wdmdo.var[ivar].type,"complex8") == 0 || strcmp(wdmdo.var[ivar].type,"vector4") == 0)
        {
            printf("ERROR: Unsupported type %s for variable %s!\n", wdmdo.var[ivar].type, wdmdo.var[ivar].name);
            return (EXIT_FAILURE);
        }
    }

    for (ivar = 0; ivar < wdmdo.nvar; ivar++)
    {
        printf("# INTERPOLATING `%s`...\n", wdmdo.var[ivar].name);
        for (icycle = 0; icycle < wdmdo.cycles; icycle++) // for each cycle
        {
            file_operationl(wdata_read_cycle(&wdmd, wdmdo.var[ivar].name, icycle, indata));
            if (wdmdo.datadim == 3)
                winterp_interpolation_3d(wdmdo.var[ivar].type[0], inNX, inNY, inNZ, indata, NX, NY, NZ, outdata);
            else if (wdmdo.datadim == 2)
                winterp_interpolation_2d(wdmdo.var[ivar].type[0], inNX, inNY, indata, NX, NY, outdata);
            else
                winterp_interpolation_1d(wdmdo.var[ivar].type[0], inNX, indata, NX, outdata);
            file_operationl(wdata_write_cycle(&wdmdo, wdmdo.var[ivar].name, outdata));
        }
    }

    // write metadata file (with default name)
    char file_name[256];
    sprintf(file_name, "%s.wtxt", argv[2]);
    printf("# WRITING `%s`\n", file_name);
    //     wdata_write_metadata_to_file(&wdmdo, file_name);
    wdata_write_metadata_to_file(&wdmdo, "");
    wdata_add_comment_to_metadata_file(file_name, "Genereted by wdata tool:");
    char cmd[1024];
    if (argc == 4)
        sprintf(cmd, "\t%s %s %s %s", argv[0], argv[1], argv[2], argv[3]);
    else if (argc == 5)
        sprintf(cmd, "\t%s %s %s %s %s", argv[0], argv[1], argv[2], argv[3], argv[4]);
    else if (argc == 6)
        sprintf(cmd, "\t%s %s %s %s %s %s", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
    wdata_add_comment_to_metadata_file(file_name, cmd);
    printf("# DONE.\n");

    /* Arrays will be cleared automatically */
    return (EXIT_SUCCESS);
}
