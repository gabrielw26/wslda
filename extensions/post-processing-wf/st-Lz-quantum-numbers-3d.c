/**
 * W-SLDA Toolkit
 * 
 * This script extracts expectation value of angular momentum for each quasi-particle orbital.
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 st-Lz-quantum-numbers-3d.c -I. -I$WSLDA/hpc-engine -o st-Lz-quantum-numbers-3d -lm -lfftw3 -I$WSLDA/lib/wderiv/c -L$WSLDA/lib/wderiv -lwderiv
 * 
 * */ 

// Standard libraries
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

// W-SLDA Toolkit API
#include "wslda_toolkit.h"

#include "wderiv.h"
#include <fftw3.h>
int compute_derivative(int nx, int ny, double dx, double dy, double complex *f, double complex *df_dx, double complex *wrk, fftw_plan plan_f, fftw_plan plan_b, int axis);


int main( int argc , char ** argv ) 
{ 
    // command line
    if(argc!=3)
    {
        printf("Usage: %s prefix iogroups\n", argv[0]);
        return( EXIT_FAILURE ) ;        
    }
    
    char *prefix = argv[1];
    int iogrp= atoi(argv[2]);
    printf("# PREFIX FOR FILES: %s\n, IOGRP=%d", prefix, iogrp);
    
    // lattice settings
    int nx, ny, nz;
    double dx, dy, dz;
    
    // other variables
    int ix, iy, iz, ixyz; // lattice iterators
    int ierr;             // error flag
    char file_name[256];
    double kF, ec, beta;
    double mu[2];
    int nwf;
    double testNu=0.0, testNv=0.0;
    
    sprintf(file_name, "%s/s3dpca.info", prefix);
    printf("# READING INFO FILE: %s\n", file_name);
    file_operationl(read_checkpoint_info_pca(file_name, &nwf, &nx, &ny, &nz, &dx, &dy, &dz, &kF, mu, &ec, &beta));
    double eF = kF * kF / 2.0;

    printf("# LATTICE: %d x %d x %d\n", nx, ny, nz);
    printf("# SPACING: %.2f x %.2f x %.2f\n", dx, dy, dz);
    printf("# VOLUME : %.2f x %.2f x %.2f\n", dx * nx, dy * ny, dz * nz);
     
    // states information for each kz
    int igrp, iwf;
    double ei;
    double complex *u, *v;
    cppmallocl(u, nx * ny * nz, double complex);
    cppmallocl(v, nx * ny * nz, double complex);
    
    double complex *du_dx, *dv_dx, *du_dy, *dv_dy;
    cppmallocl(du_dx, nx * ny * nz, double complex);
    cppmallocl(du_dy, nx * ny * nz, double complex);
    cppmallocl(dv_dx, nx * ny * nz, double complex);
    cppmallocl(dv_dy, nx * ny * nz, double complex);
    
    // file pointers 
    FILE *pFile_en, *pFile_kkz, *pFile_wfu, *pFile_wfv;

    sprintf(file_name, "%s_Lz.txt", prefix);
    printf("# CREATING FILE `%s`\n", file_name);
    FILE *pFile = fopen(file_name, "w");
    fprintf(pFile, "# kF    %16.8g\n", kF);
    fprintf(pFile, "# eF    %16.8g\n", eF);
    fprintf(pFile, "# COLUMNS:\n");
    fprintf(pFile, "# 1: ei/eF\n");
    fprintf(pFile, "# 2: <u|u>\n");
    fprintf(pFile, "# 3: <v|v>\n");
    fprintf(pFile, "# 4: <u|Lz|u>\n");
    fprintf(pFile, "# 5: <v|Lz|v>\n");
    fprintf(pFile, "# 6: <u|Lz|u>/<u|u>\n");
    fprintf(pFile, "# 7: <v|Lz|v>/<v|v>\n");

    // wderiv
    wderiv_init_3d(nx, ny, nz, dx, dy, dz);
    
    for (igrp = 0; igrp < iogrp; igrp++)
    {
        // read file header
        sprintf(file_name, "%s/s3dpca.%04d.info", prefix, igrp);
        printf("# OPENING: %s\n", file_name);
        file_operationl(read_checkpoint_info_pca(file_name, &nwf, &nx, &ny, &nz, &dx, &dy, &dz, &kF, mu, &ec, &beta));

        // open files
        sprintf(file_name, "%s/s3dpca.%04d.en", prefix, igrp);
        pFile_en = fopen(file_name, "rb");
        if(pFile_en==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s3dpca.%04d.wfu", prefix, igrp);
        pFile_wfu = fopen(file_name, "rb");
        if(pFile_wfu==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s3dpca.%04d.wfv", prefix, igrp);
        pFile_wfv = fopen(file_name, "rb");
        if(pFile_wfv==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}

        // read wave-functions
        printf("# PROCESSING WF[%d] FOR IOGROUP=%d\n", nwf, igrp);
        for (iwf = 0; iwf < nwf; iwf++)
        {
            fread(&ei, sizeof(double), 1, pFile_en);  // eigen energy
            fread(u, sizeof(double complex) * nx * ny * nz, 1, pFile_wfu); // u-component
            fread(v, sizeof(double complex) * nx * ny * nz, 1, pFile_wfv); // v-compoment

            // normalize
            ixyz = 0;
            for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++) for (iz = 0; iz < nz; iz++)
            {
                u[ixyz] /= sqrt(dx * dy * dz);
                v[ixyz] /= sqrt(dx * dy * dz);
                ixyz++;
            }
            
            // integrate |u|^2 and |v|^2
            double u2 = 0.0, v2 = 0.0;
            ixyz = 0;
            for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++) for (iz = 0; iz < nz; iz++)
            {
                u2 += (pow(creal(u[ixyz]), 2) + pow(cimag(u[ixyz]), 2)) * dx * dy * dz;
                v2 += (pow(creal(v[ixyz]), 2) + pow(cimag(v[ixyz]), 2)) * dx * dy * dz;
                ixyz++;
            }

            testNu+=u2; // for testing purposes
            testNv+=v2; // for testing purposes

            // compute derivatives
            wderiv_derivative_3d_c(WDERIV_DX, 1, u, du_dx);
            wderiv_derivative_3d_c(WDERIV_DY, 1, u, du_dy);
            wderiv_derivative_3d_c(WDERIV_DX, 1, v, dv_dx);
            wderiv_derivative_3d_c(WDERIV_DY, 1, v, dv_dy);

            // integrate <u|Lz|u> and <v|Lz|v>
            // Lz = x*py - y*px = -i*x*d/dy + i*y*d/dx
            // and take with negative sign ("-=")
            double uLz = 0.0, vLz = 0.0;
            double _x, _y;
            ixyz = 0;
            for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++) for (iz = 0; iz < nz; iz++)
            {
                _x = dx * (ix - nx / 2);
                _y = dy * (iy - ny / 2);

                uLz -= creal(conj(u[ixyz]) * (I * _y * du_dx[ixyz] - I * _x * du_dy[ixyz])) * dx * dy * dz;
                vLz -= creal(conj(v[ixyz]) * (I * _y * dv_dx[ixyz] - I * _x * dv_dy[ixyz])) * dx * dy * dz;

                ixyz++;
            }

            // write to the file
            fprintf(pFile, "%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", ei / eF, u2, v2, uLz, vLz, uLz / u2, vLz / v2);
    
        }
        
        // close files
        fclose(pFile_en);
        fclose(pFile_wfu);
        fclose(pFile_wfv);
    }


    fclose(pFile);
    wderiv_clean_3d();

    printf("TEST: Nu=%f, Nv=%f\n", testNu, testNv);
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}

