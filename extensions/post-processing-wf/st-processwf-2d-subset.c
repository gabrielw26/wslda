/**
 * W-SLDA Toolkit
 * 
 * This code computes density arising from subset wave-functions.
 * Assumptions:
 *   - system is spin symmetric.
 * 
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 st-processwf-2d-subset.c -I. -I$WSLDA/hpc-engine -o st-processwf-2d-subset -lm -I$WSLDA/lib/wdata/c -L$WSLDA/lib/wdata -lwdata
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

// wdata lib
#include "wdata.h"


int main( int argc , char ** argv ) 
{ 
    // command line
    if(argc!=5)
    {
        printf("Usage: %s prefix e_min/eF e_max/eF outprefix\n", argv[0]);
        return( EXIT_FAILURE ) ;        
    }
    
    char *prefix = argv[1];
    char *outprefix = argv[4];
    printf("# PREFIX FOR FILES: %s\n", prefix);
    double emin = atof(argv[2]);
    double emax = atof(argv[3]);
    printf("# COMPUTING DENSITY FROM STATES IN THE ENERGY RANGE E_n/eF in [%f,%f]\n", emin, emax);
    
    // lattice settings
    int nx, ny, nz;
    double dx, dy, dz;
    
    // other variables
    int ix, iy, iz, ixyz; // lattice iterators
    int ierr;             // error flag
    char file_name[256];
    double kF, ec, beta, eF;
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

    double *rho; // buffer for density
    cppmallocl(rho, nx * ny, double);
    for(ixyz=0; ixyz<nx*ny; ixyz++) rho[ixyz]=0.0; // reset buffer
    int subset_nwf=0, total_nwf = 0;
    
    // file pointers 
    FILE *pFile_en, *pFile_kkz, *pFile_wfu, *pFile_wfv;
    
    for (ikz = 0; ikz < nz / 2; ikz++)
    {
        // read file header
        sprintf(file_name, "%s/s2dpca.%04d.info", prefix, ikz);
        printf("# OPENING: %s\n", file_name);
        file_operationl(read_checkpoint_info_pca(file_name, &nwf, &nx, &ny, &nz, &dx, &dy, &dz, &kF, mu, &ec, &beta));
        eF = kF*kF / 2.0;

        double weight=2.0; // (+kz, -kz)
        if(ikz==0) weight=1.0; // special case
        
        // open files
        sprintf(file_name, "%s/s2dpca.%04d.en", prefix, ikz);
        pFile_en = fopen(file_name, "rb");
        if(pFile_en==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s2dpca.%04d.kkz", prefix, ikz);
        pFile_kkz = fopen(file_name, "rb");
        if(pFile_en==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s2dpca.%04d.wfu", prefix, ikz);
        pFile_wfu = fopen(file_name, "rb");
        if(pFile_en==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        sprintf(file_name, "%s/s2dpca.%04d.wfv", prefix, ikz);
        pFile_wfv = fopen(file_name, "rb");
        if(pFile_en==NULL) {printf("ERROR: Cannot open %s\n", file_name); return( EXIT_FAILURE ) ;}
        
        // read wave-functions
        printf("# PROCESSING WF[%d] FOR ikz=%d\n", nwf, ikz);
        for (iwf = 0; iwf < nwf; iwf++)
        {
            fread(&ei, sizeof(double), 1, pFile_en);  // eigen energy
            fread(&kz, sizeof(double), 1, pFile_kkz); // kz value
            fread(u, sizeof(double complex) * nx * ny, 1, pFile_wfu); // u-component
            fread(v, sizeof(double complex) * nx * ny, 1, pFile_wfv); // v-compoment

            if(ei<0.0)
            {
                printf("# ERROR: E_n=%f < 0 indictes data for spinimbalanced case. NOT SUPPORTED!\n", ei);
                return( EXIT_FAILURE ) ;
            }

            double fbEn=fbeta(ei, beta);
            double fbmEn = 1.0 - fbEn;
            
            // normalize
            ixyz = 0;
            for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++)
            {
                u[ixyz] /= sqrt(dx * dy);
                v[ixyz] /= sqrt(dx * dy);
                ixyz++;
            }

            ixyz = 0;
            // printf("AAAA: %f %f %f %d %d %d\n", ei, ei/eF,eF,ei>emin*eF,ei<emax*eF,ei>emin*eF && ei<emax*eF);
            if(ei>emin*eF && ei<emax*eF) // add contribution to the density if state in requested energy range
            {
                for (ix = 0; ix < nx; ix++) for (iy = 0; iy < ny; iy++)
                {
                    // add contribution to the density
                    #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
                    rho[ixyz]+=(cnorm(v[ixyz])*fbmEn + cnorm(u[ixyz])*fbEn) / (dz * nz) * weight;

                    ixyz++;
                }

                subset_nwf += (int)weight;
            }
            total_nwf += (int)weight;
        }
        
        
        // close files
        fclose(pFile_en);
        fclose(pFile_kkz);
        fclose(pFile_wfu);
        fclose(pFile_wfv);
        
    }

    printf("# NUMBER OF STATES CONTRIBUTIG TO THE SUBSET DENSITY=%d FROM TOTAL NUMBER=%d\n", subset_nwf, total_nwf);

    // create wdata set
    wdata_metadata mdset;
    // Lattice
    mdset.datadim = 2; // 2D data
    mdset.nx = nx;
    mdset.ny = ny;
    mdset.nz = nz;
    mdset.dx = dx;
    mdset.dy = dy;
    mdset.dz = dz;

    strcpy(mdset.prefix, outprefix);
    mdset.t0 = 0.0;
    mdset.dt = 1.0;

    wdata_variable rho_a = {"rho_a", "real", "none", "wdat"};
    wdata_add_variable(&mdset, &rho_a);
    wdata_link rho_b = {"rho_b", "rho_a"};
    wdata_add_link(&mdset, &rho_b);

    wdata_setconst(&mdset, "eF", eF);
    wdata_setconst_unit(&mdset, "emin", emin, "eF");
    wdata_setconst_unit(&mdset, "emax", emax, "eF");
    wdata_setconst(&mdset, "subset_nwf", subset_nwf);
    wdata_setconst(&mdset, "total_nwf", total_nwf);


    // just in case - clear data sets
    // it removes binary files if they exist
    wdata_clear_database(&mdset);

    wdata_add_cycle(&mdset);
    wdata_write_cycle(&mdset, "rho_a", rho);

    // write metadata file (with default name)
    wdata_write_metadata_to_file(&mdset, "");
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
