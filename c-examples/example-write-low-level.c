/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Code demonstrates how to write wdata files using only standard I/O functions.
 *
 * Compile command
 *  gcc example-write-low-level.c -o example-write-low-level -lm
 * */

// standard libs
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

    // Settings
    int datadim = 3; // 3D data
    int nx = 24;
    int ny = 28;
    int nz = 32;
    double dx = 1.0;
    double dy = 1.0;
    double dz = 1.0;

    char prefix[256] = "test";

    double t0 = 0.0;
    double dt = 1.0;

    // Allocate memory for variables: real, complex and vector
    double *dataR;  // real data
    Complex *dataC; // complex data
    double *dataV;  // vector data
    cppmallocl(dataR, nx*ny*nz, double);
    cppmallocl(dataC, nx*ny*nz, Complex);
    cppmallocl(dataV, nx*ny*nz*3, double); // factor 3 accounts for three compoments of vector variable

    // Fill array with (artificial) data data write to file
    int ncycles = 10; // number of cycles to be generated
    for (int icycle = 0; icycle < ncycles; icycle++)
    {
        // time decomposition
        double time=t0+dt*icycle;
        printf("CURRENT TIME: %lf\n", time);

        // Fill array with (artificial) data
        int ixyz = 0;
        for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++)
                for (int iz = 0; iz < nz; iz++)
                {
                    // coordinate decomposition
                    double x = dx * (ix - nx / 2);
                    double y = dy * (iy - ny / 2);
                    double z = dz * (iz - nz / 2);

                    // add data to array
                    dataR[ixyz] = function_xyz(x, y, z, time);
                    dataC[ixyz] = dataR[ixyz] + I * 0.0;

                    // vector - solid body rotation around z axis
                    dataV[ixyz + 0 * nx*ny*nz] = -1.0 * y;
                    dataV[ixyz + 1 * nx*ny*nz] = 1.0 * x;
                    dataV[ixyz + 2 * nx*ny*nz] = 0.0;

                    ixyz++;
                }

        // write it to binary file
        // Note: to simplify notation, I skip checking for errors.
        FILE *pFile;
        char file_name[256];

        sprintf(file_name, "%s_density.wdat", prefix); // real variable
        pFile= fopen (file_name, "ab"); // open file in `append` and `binary` mode
        fwrite (dataR, nx*ny*nz*sizeof(double), 1, pFile); // add datablock
        fclose(pFile);

        sprintf(file_name, "%s_delta.wdat", prefix); // complex variable
        pFile= fopen (file_name, "ab"); // open file in `append` and `binary` mode
        fwrite (dataC, nx*ny*nz*sizeof(Complex), 1, pFile); // add datablock
        fclose(pFile);

        sprintf(file_name, "%s_current.wdat", prefix); // vector variable
        pFile= fopen (file_name, "ab"); // open file in `append` and `binary` mode
        fwrite (dataV, nx*ny*nz*sizeof(double)*3, 1, pFile); // add datablock
        fclose(pFile);
    }

    // Add descriptor file to binary files
    char wtxt_name[256];
    sprintf(wtxt_name, "%s.wtxt", prefix);
    printf("CREATING DESCRIPTOR %s.\n", wtxt_name);
    FILE *wtxt = fopen(wtxt_name, "w");
    fprintf(wtxt,"nx      %d  # lattice\n", nx);
    fprintf(wtxt,"ny      %d  # lattice\n", ny);
    fprintf(wtxt,"nz      %d  # lattice\n", nz);
    fprintf(wtxt,"dx      %lf # lattice\n", dx);
    fprintf(wtxt,"dy      %lf # lattice\n", dy);
    fprintf(wtxt,"dz      %lf # lattice\n", dz);
    fprintf(wtxt,"datadim %d  # 3D data\n", 3);
    fprintf(wtxt,"prefix  %s  # prefix for files belonging to this data set\n", prefix);
    fprintf(wtxt,"cycles  %d  # number of saved cycles\n", ncycles);
    fprintf(wtxt,"t0      %lf # for time extraction\n", t0);
    fprintf(wtxt,"dt      %lf # for time extraction\n", dt);
    fprintf(wtxt,"# variables\n");
    fprintf(wtxt,"# tag       name      type     unit   format\n");
    fprintf(wtxt,"var      density      real     none     wdat\n");
    fprintf(wtxt,"var        delta   complex     none     wdat\n");
    fprintf(wtxt,"var      current    vector     none     wdat\n");
    fclose(wtxt);

    return 0;
}
