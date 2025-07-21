/**
 * W-SLDA Toolkit
 * 
 * This script extracts expectation value of angular momentum for each quasi-particle orbital.
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 st-Lz-quantum-numbers-2d.c -I. -I$WSLDA/hpc-engine -o st-Lz-quantum-numbers-2d -lm -lfftw3
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

#include <fftw3.h>
int compute_derivative(int nx, int ny, double dx, double dy, double complex *f, double complex *df_dx, double complex *wrk, fftw_plan plan_f, fftw_plan plan_b, int axis);


int main( int argc , char ** argv ) 
{ 
    // command line
    if(argc!=2)
    {
        printf("Usage: %s prefix\n", argv[0]);
        return( EXIT_FAILURE ) ;        
    }
    
    char *prefix = argv[1];
    printf("# PREFIX FOR FILES: %s\n", prefix);
    
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
    
    sprintf(file_name, "%s/s2dpca.info", prefix);
    printf("# READING INFO FILE: %s\n", file_name);
    file_operationl(read_checkpoint_info_pca(file_name, &nwf, &nx, &ny, &nz, &dx, &dy, &dz, &kF, mu, &ec, &beta));

    printf("# LATTICE: %d x %d x %d\n", nx, ny, nz);
    printf("# SPACING: %.2f x %.2f x %.2f\n", dx, dy, dz);
    printf("# VOLUME : %.2f x %.2f x %.2f\n", dx * nx, dy * ny, dz * nz);
     
    // states information for each kz
    int ikz, iwf;
    double ei, kz;
    double complex *u, *v;
    cppmallocl(u, nx * ny, double complex);
    cppmallocl(v, nx * ny, double complex);
    
    double complex *du_dx, *dv_dx, *du_dy, *dv_dy;
    cppmallocl(du_dx, nx * ny, double complex);
    cppmallocl(du_dy, nx * ny, double complex);
    cppmallocl(dv_dx, nx * ny, double complex);
    cppmallocl(dv_dy, nx * ny, double complex);
    
    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    double complex *wrk;
    cppmallocl(wrk, nx * ny, double complex);
    fftw_plan plan_f = fftw_plan_dft_2d(nx, ny, wrk, wrk, FFTW_FORWARD, USE_FFTW_PLANNER);
    fftw_plan plan_b = fftw_plan_dft_2d(nx, ny, wrk, wrk, FFTW_BACKWARD, USE_FFTW_PLANNER);
    
    // file pointers 
    FILE *pFile_en, *pFile_kkz, *pFile_wfu, *pFile_wfv;
    
    for (ikz = 0; ikz < nz / 2; ikz++)
    {
        // read file header
        sprintf(file_name, "%s/s2dpca.%04d.info", prefix, ikz);
        printf("# OPENING: %s\n", file_name);
        file_operationl(read_checkpoint_info_pca(file_name, &nwf, &nx, &ny, &nz, &dx, &dy, &dz, &kF, mu, &ec, &beta));
        
        // open files
        sprintf(file_name, "%s/s2dpca.%04d.en", prefix, ikz);
        pFile_en = fopen(file_name, "rb");
        if(pFile_en==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s2dpca.%04d.kkz", prefix, ikz);
        pFile_kkz = fopen(file_name, "rb");
        if(pFile_kkz==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s2dpca.%04d.wfu", prefix, ikz);
        pFile_wfu = fopen(file_name, "rb");
        if(pFile_wfu==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s2dpca.%04d.wfv", prefix, ikz);
        pFile_wfv = fopen(file_name, "rb");
        if(pFile_wfv==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        
        sprintf(file_name, "%s_Lz.%04d.txt", prefix, ikz);
        printf("# CREATING FILE `%s`\n", file_name);      
        double eF = kF * kF / 2.0;
        FILE *pFile = fopen(file_name, "w");
        fprintf(pFile, "# kF    %16.8g\n", kF);
        fprintf(pFile, "# eF    %16.8g\n", eF);
        fprintf(pFile, "# COLUMNS:\n");
        fprintf(pFile, "# 1: kz/kF\n");
        fprintf(pFile, "# 2: ei/eF\n");
        fprintf(pFile, "# 3: <u|u>\n");
        fprintf(pFile, "# 4: <v|v>\n");
        fprintf(pFile, "# 5: <u|Lz|u>\n");
        fprintf(pFile, "# 6: <v|Lz|v>\n");
        fprintf(pFile, "# 7: <u|Lz|u>/<u|u>\n");
        fprintf(pFile, "# 8: <v|Lz|v>/<v|v>\n");
        
        // read wave-functions
        printf("# PROCESSING WF[%d] FOR ikz=%d\n", nwf, ikz);
        for (iwf = 0; iwf < nwf; iwf++)
        {
            fread(&ei, sizeof(double), 1, pFile_en);  // eigen energy
            fread(&kz, sizeof(double), 1, pFile_kkz); // kz value
            fread(u, sizeof(double complex) * nx * ny, 1, pFile_wfu); // u-component
            fread(v, sizeof(double complex) * nx * ny, 1, pFile_wfv); // v-compoment
            
            // normalize
            ixyz = 0;
            for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++)
            {
                u[ixyz] /= sqrt(dx * dy);
                v[ixyz] /= sqrt(dx * dy);
                ixyz++;
            }
            
            
            // integrate |u|^2 and |v|^2
            double u2 = 0.0, v2 = 0.0;
            ixyz = 0;
            for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++)
            {
                u2 += (pow(creal(u[ixyz]), 2) + pow(cimag(u[ixyz]), 2)) * dx * dy;
                v2 += (pow(creal(v[ixyz]), 2) + pow(cimag(v[ixyz]), 2)) * dx * dy;
                ixyz++;
            }
            
            // compute derivatives
            compute_derivative(nx, ny, dx, dy, u, du_dx, wrk, plan_f, plan_b, 0);
            compute_derivative(nx, ny, dx, dy, u, du_dy, wrk, plan_f, plan_b, 1);
            compute_derivative(nx, ny, dx, dy, v, dv_dx, wrk, plan_f, plan_b, 0);
            compute_derivative(nx, ny, dx, dy, v, dv_dy, wrk, plan_f, plan_b, 1);
            
            // integrate <u|Lz|u> and <v|Lz|v>
            // Lz = x*py - y*px = -i*x*d/dy + i*y*d/dx
            // and take with negative sign ("-=")
            double uLz = 0.0, vLz = 0.0;
            double _x, _y;
            ixyz = 0;
            for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++)
            {
                _x = dx * (ix - nx / 2);
                _y = dy * (iy - ny / 2);

                uLz -= creal(conj(u[ixyz]) * (I * _y * du_dx[ixyz] - I * _x * du_dy[ixyz])) * dx * dy;
                vLz -= creal(conj(v[ixyz]) * (I * _y * dv_dx[ixyz] - I * _x * dv_dy[ixyz])) * dx * dy;

                ixyz++;
            }
            
            // write to the file
            fprintf(pFile, "%12.6f %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n", kz / kF, ei / eF, u2, v2, uLz, vLz, uLz / u2, vLz / v2);
    
        }
        
        
            
        // close files
        fclose(pFile_en);
        fclose(pFile_kkz);
        fclose(pFile_wfu);
        fclose(pFile_wfv);
        
        fclose(pFile);
        
    }
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}

// derivative computation
int compute_derivative(int nx, int ny, double dx, double dy, double complex *f, double complex *df_dx, double complex *wrk, fftw_plan plan_f, fftw_plan plan_b, int axis)
{
    int ixyz;
    int ix, iy;
    double kx, ky;

    // copy data to working array
    for (ixyz = 0; ixyz < nx * ny; ixyz++)
        wrk[ixyz] = f[ixyz];

    fftw_execute(plan_f);

    // multiply by momentum
    ixyz = 0;
    for (ix = 0; ix < nx; ix++)
    {
        for (iy = 0; iy < ny; iy++)
        {
            // extract momentum
            if (ix < nx / 2)
                kx = 2. * M_PI / ((double)nx * dx) * (double)(ix);
            else
                kx = 2. * M_PI / ((double)nx * dx) * (double)(ix - nx);
            if (ix == nx / 2)
                kx = 0.0;

            if (iy < ny / 2)
                ky = 2. * M_PI / ((double)ny * dy) * (double)(iy);
            else
                ky = 2. * M_PI / ((double)ny * dy) * (double)(iy - ny);
            if (iy == ny / 2)
                ky = 0.0;

            if (axis == 0)
                wrk[ixyz] *= I * kx / (nx * ny); // note: normalization factor is included
            else
                wrk[ixyz] *= I * ky / (nx * ny); // note: normalization factor is included

            ixyz++;
        }
    }

    fftw_execute(plan_b);

    // copy data to result table
    for (ixyz = 0; ixyz < nx * ny; ixyz++)
        df_dx[ixyz] = wrk[ixyz];

    return 0;
}
