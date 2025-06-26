/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Code demonstrates how to write wdata files
 * */

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

typedef double complex Complex;

#define SQRT_2PI 2.506628274631
double function_x(double x, double sigma)
{
    return 1. / (sigma * SQRT_2PI) * exp(-0.5 * x * x / (sigma * sigma));
}

double function_xyz(double x, double y, double z, double time)
{
    double val = 0.0;
    double sigmax = 2.0 + 0.20 * time;
    double sigmay = 3.0 + 0.15 * time;
    double sigmaz = 4.0 + 0.10 * time;

    val = function_x(x, sigmax) * function_x(y, sigmay) * function_x(z, sigmaz);

    return val;
}

double function_time(double t0, const int icycle)
{
    return t0 + 0.1 * pow(icycle, 2);
}

// Generates a random double in the range (a, b)
double random_in_range(double a, double b)
{
    // Ensure a < b
    if (a >= b) {
        fprintf(stderr, "Invalid range: a must be less than b\n");
        exit(EXIT_FAILURE);
    }
    double r = ((double) rand() / (double) RAND_MAX); // [0, 1]
    // To strictly exclude endpoints, add a small epsilon
    double epsilon = 1e-9;
    r = r * (1.0 - 2.0 * epsilon) + epsilon; // (0, 1)
    return a + r * (b - a); // (a, b)
}

#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "error: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "error: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        return -1;                                                          \
    }

int main()
{
    // create metadata handler
    wdata_metadata md;

    // Lattice
    md.datadim = 3; // 3D data
    md.nx = 24;
    md.ny = 28;
    md.nz = 32;
    md.dx = -1.0; // varying x
    md.dy = -1.0; // varying y
    md.dz = -1.0; // varying z
    md.x0 = 0.0;
    md.y0 = 0.0;
    md.z0 = 0.0;

    strcpy(md.prefix, "test");
    md.t0 = 0.0;
    md.dt = -1.0; // varying t

    // add variables to data set
    // for each variable binary file of name `prefix_`varname`.wdat will be created
    wdata_variable vdensity_a = {"density_a", "real", "none", "wdat"};
    wdata_add_variable(&md, &vdensity_a);

    wdata_variable vdelta = {"delta", "complex", "none", "wdat"};
    wdata_add_variable(&md, &vdelta);

    wdata_variable vcurrent_a = {"current_a", "vector", "none", "wdat"};
    wdata_add_variable(&md, &vcurrent_a);

    // add links to data sets
    // links are alternative names of the same variable
    wdata_link ldensity_b = {"density_b", "density_a"};
    wdata_add_link(&md, &ldensity_b);

    wdata_link lcurrent_b = {"current_b", "current_a"};
    wdata_add_link(&md, &lcurrent_b);

    // add constants
    wdata_const lconst_eF = {"eF", 0.5, "MeV"};
    wdata_add_const(&md, &lconst_eF);

    wdata_const lconst_kF = {"kF", 1.0, "1/fm"};
    wdata_add_const(&md, &lconst_kF);

    // just in case - clear data sets
    // it removes binary files if they alredy exists
    wdata_clear_database(&md);

    // generate some artifical data
    int bdim = wdata_get_blocklength(&md); // get block size
    double time, x, y, z;

    double *dataR;  // real data
    Complex *dataC; // complex data
    double *dataV;  // vector data
    cppmallocl(dataR, bdim, double);
    cppmallocl(dataC, bdim, Complex);
    cppmallocl(dataV, bdim * 3, double); // factor 3 accounts for three compoments of vector variable

    int ncycles = 10; // number of cycles to be generated

    // varying time
    for (int icycle = 0; icycle < ncycles; icycle++)
    {
        time=function_time(md.t0,icycle);
        wdata_add_time(&md, icycle, &time);
    }

    // varying x
    for (int ix = 0; ix < md.nx; ix++)
    {
        x =  (ix - md.nx / 2) + md.dx*random_in_range(-0.49,0.49);
        wdata_add_x(&md, &x);
    }

    // varying y
    for (int iy = 0; iy < md.ny; iy++)
    {
        y = (iy - md.ny / 2) + md.dy*random_in_range(-0.49,0.49);
        wdata_add_y(&md, &y);
    }

    // varying z
    for (int iz = 0; iz < md.nz; iz++)
    {
        z = (iz - md.nz / 2) + md.dz*random_in_range(-0.49,0.49);
        wdata_add_z(&md, &z);
    }

    for (int icycle = 0; icycle < ncycles; icycle++)
    {
        // time decomposition
        wdata_get_time(&md, icycle, &time);
        printf("CURRENT TIME: %lf\n", time);

        int ixyz = 0;
        for (int ix = 0; ix < md.nx; ix++)
            for (int iy = 0; iy < md.ny; iy++)
                for (int iz = 0; iz < md.nz; iz++)
                {
                    // coordinate decomposition
                    wdata_get_x(&md, ix, &x);
                    wdata_get_y(&md, iy, &y);
                    wdata_get_z(&md, iz, &z);

                    // add data to array
                    dataR[ixyz] = function_xyz(x, y, z, time);
                    dataC[ixyz] = dataR[ixyz] + I * 0.0;

                    // vector - solid body rotation around z axis
                    dataV[ixyz + 0 * bdim] = -1.0 * y;
                    dataV[ixyz + 1 * bdim] = 1.0 * x;
                    dataV[ixyz + 2 * bdim] = 0.0;

                    ixyz++;
                }

        int ierr;

        // add info about new cycle to metadata
        wdata_add_cycle(&md);

        // and add cycle to binary sets
        ierr = wdata_write_cycle(&md, "density_a", dataR);
        // alternatively you can use:
        // ierr=wdata_add_datablock(&md, &vdensity_a, dataR);
        if (ierr != 0)
        {
            printf("ERROR: Cannot add density_a!\n");
            return 1;
        }

        ierr = wdata_write_cycle(&md, "delta", dataC);
        if (ierr != 0)
        {
            printf("ERROR: Cannot add delta!\n");
            return 1;
        }

        ierr = wdata_write_cycle(&md, "current_a", dataV);
        if (ierr != 0)
        {
            printf("ERROR: Cannot add delta!\n");
            return 1;
        }
    }

    // write metadata file (with default name)
    wdata_write_metadata_to_file(&md, "");
    return 0;
}
