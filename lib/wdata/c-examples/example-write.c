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
    wdata_reset_metadata(&md);

    // Lattice
    md.datadim = 3; // 3D data
    md.nx = 24;
    md.ny = 28;
    md.nz = 32;
    md.dx = 1.0;
    md.dy = 1.0;
    md.dz = 1.0;
    md.x0 = -md.nx/2;
    md.y0 = -md.ny/2;
    md.z0 = -md.nz/2;

    strcpy(md.prefix, "test");
    md.t0 = 0.0;
    md.dt = 1.0;

    // add variables to data set
    // for each variable binary file of name `prefix_`varname`.wdat will be created
    wdata_variable vdensity_a = {"density_a", "real", "none", "wdat"};
    wdata_add_variable(&md, &vdensity_a);

    wdata_variable vdensity_f = {"density_f", "real4", "none", "wdat"}; // density_a but in float precision
    wdata_add_variable(&md, &vdensity_f);

    wdata_variable vdelta = {"delta", "complex", "none", "wdat"};
    wdata_add_variable(&md, &vdelta);

    wdata_variable vdelta_f = {"delta_f", "complex8", "none", "wdat"}; // delta but in float precision
    wdata_add_variable(&md, &vdelta_f);

    wdata_variable vcurrent_a = {"current_a", "vector", "none", "wdat"};
    wdata_add_variable(&md, &vcurrent_a);

    wdata_variable vcurrent_f = {"current_f", "vector4", "none", "wdat"}; // current_a but in float precision
    wdata_add_variable(&md, &vcurrent_f);

    wdata_variable vcurrent_2 = {"current_2", "vector(2)", "none", "wdat"}; // current_a but in 2d, since third compoment is always zero
    wdata_add_variable(&md, &vcurrent_2);

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
    double time;

    double *dataR;  // real data
    Complex *dataC; // complex data
    double *dataV;  // vector data
    cppmallocl(dataR, bdim, double);
    cppmallocl(dataC, bdim, Complex);
    cppmallocl(dataV, bdim * 3, double); // factor 3 accounts for three compoments of vector variable

    int ncycles = 10; // number of cycles to be generated

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
                    double x = md.dx * (ix - md.nx / 2);
                    double y = md.dy * (iy - md.ny / 2);
                    double z = md.dz * (iz - md.nz / 2);

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

        // use native function, it assumes that the user provides binary data in a format consistent with the variable type
        wdata_write_cycle(&md, "density_a", dataR);
        wdata_write_cycle(&md, "delta", dataC);
        wdata_write_cycle(&md, "current_a", dataV);

        // input array is provided in double precision, it will be downgraded to a float precision variable if type requests this
        wdata_write_cycle_d(&md, "density_f", dataR);
        wdata_write_cycle_d(&md, "delta_f", (double *)dataC);
        wdata_write_cycle_d(&md, "current_f", dataV);

        // write vector in 2d
        wdata_write_cycle_d(&md, "current_2", dataV); // saves only x and y coordinates
    }

    // write metadata file (with default name)
    wdata_write_metadata_to_file(&md, "");
    return 0;
}
