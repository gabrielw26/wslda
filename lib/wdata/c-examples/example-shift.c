/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Code shits domain by vector [xs,ys,zs]
 * 
 * To compile:
 * g++ shift.c -o shift -I/home/gabrielw/MyProjects/wdata/c -L/home/gabrielw/MyProjects/wdata -lwdata -lm
 * */

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>

typedef std::complex<double> Complex;

#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "error: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "error: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        return -1;                                                          \
    }

#define XS 1
#define YS 1
#define ZS 1

int main()
{
    int ierr;

    // create metadata handler
    wdata_metadata md;

    // read metadata from file
    ierr = wdata_parse_metadata_file("test-40.wtxt", &md);
    if (ierr != 0)
    {
        printf("Cannot read metadata file!\n");
        return 1;
    }

    // print metadata file
    wdata_print_metadata(&md, stdout);

    wdata_metadata md2 = md; // copy it
    wdata_setprefix(&md2, "test-s1-40");
    wdata_clear_database(&md2);

    // allocate memory for variables
    int bdim = wdata_get_blocklength(&md); // get block size

    double *dataR, *dataR2;  // real data
    Complex *dataC, *dataC2; // complex data
    double *dataV, *dataV2;  // vector data
    cppmallocl(dataR, bdim, double);
    cppmallocl(dataC, bdim, Complex);
    cppmallocl(dataV, bdim * 3, double); // factor 3 accounts for three compoments of vector variable
    cppmallocl(dataR2, bdim, double);
    cppmallocl(dataC2, bdim, Complex);
    cppmallocl(dataV2, bdim * 3, double); // factor 3 accounts for three compoments of vector variable

    int icycle, ivar;
    int ix, iy, iz, ixyz;
    int ix2, iy2, iz2, ixyz2;
    for (icycle = 0; icycle < md.cycles; icycle++)
    {
        printf("Processing cycle %d\n", icycle);
        for (ivar = 0; ivar < md.nvar; ivar++)
        {
            if (md.var[ivar].type[0] == 'r')
            {
                printf("Shifting %s...\n", md.var[ivar].name);
                ierr = wdata_read_cycle(&md, md.var[ivar].name, icycle, dataR);
                if (ierr != 0)
                {
                    printf("ERROR: Cannot read!\n");
                    return 1;
                }

                ixyz = 0;
                for (ix = 0; ix < md.nx; ix++)
                    for (iy = 0; iy < md.ny; iy++)
                        for (iz = 0; iz < md.nz; iz++)
                        {
                            ix2 = (ix + XS) % md.nx;
                            iy2 = (iy + YS) % md.ny;
                            iz2 = (iz + ZS) % md.nz;
                            ixyz2 = iz2 + md.nz * iy2 + md.nz * md.ny * ix2;

                            dataR2[ixyz2] = dataR[ixyz];

                            ixyz++;
                        }

                // and add cycle to binary sets
                wdata_add_cycle(&md2);
                ierr = wdata_write_cycle(&md2, md.var[ivar].name, dataR2);
                if (ierr != 0)
                {
                    printf("ERROR: Cannot add!\n");
                    return 1;
                }
            }

            if (md.var[ivar].type[0] == 'c')
            {
                printf("Shifting %s...\n", md.var[ivar].name);
                ierr = wdata_read_cycle(&md, md.var[ivar].name, icycle, dataC);
                if (ierr != 0)
                {
                    printf("ERROR: Cannot read!\n");
                    return 1;
                }

                ixyz = 0;
                for (ix = 0; ix < md.nx; ix++)
                    for (iy = 0; iy < md.ny; iy++)
                        for (iz = 0; iz < md.nz; iz++)
                        {
                            ix2 = (ix + XS) % md.nx;
                            iy2 = (iy + YS) % md.ny;
                            iz2 = (iz + ZS) % md.nz;
                            ixyz2 = iz2 + md.nz * iy2 + md.nz * md.ny * ix2;

                            dataC2[ixyz2] = dataC[ixyz];

                            ixyz++;
                        }

                // and add cycle to binary sets
                wdata_add_cycle(&md2);
                ierr = wdata_write_cycle(&md2, md.var[ivar].name, dataC2);
                if (ierr != 0)
                {
                    printf("ERROR: Cannot add!\n");
                    return 1;
                }
            }

            if (md.var[ivar].type[0] == 'v')
            {
                printf("Shifting %s...\n", md.var[ivar].name);
                ierr = wdata_read_cycle(&md, md.var[ivar].name, icycle, dataV);
                if (ierr != 0)
                {
                    printf("ERROR: Cannot read!\n");
                    return 1;
                }

                ixyz = 0;
                for (ix = 0; ix < md.nx; ix++)
                    for (iy = 0; iy < md.ny; iy++)
                        for (iz = 0; iz < md.nz; iz++)
                        {
                            ix2 = (ix + XS) % md.nx;
                            iy2 = (iy + YS) % md.ny;
                            iz2 = (iz + ZS) % md.nz;
                            ixyz2 = iz2 + md.nz * iy2 + md.nz * md.ny * ix2;

                            dataV2[ixyz2 + 0 * bdim] = dataV[ixyz + 0 * bdim];
                            dataV2[ixyz2 + 1 * bdim] = dataV[ixyz + 1 * bdim];
                            dataV2[ixyz2 + 2 * bdim] = dataV[ixyz + 2 * bdim];

                            ixyz++;
                        }

                // and add cycle to binary sets
                wdata_add_cycle(&md2);
                ierr = wdata_write_cycle(&md2, md.var[ivar].name, dataV2);
                if (ierr != 0)
                {
                    printf("ERROR: Cannot add!\n");
                    return 1;
                }
            }
        }
    }

    // write metadata file (with default name)
    wdata_write_metadata_to_file(&md2, "");

    return 0;
}
