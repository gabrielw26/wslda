
/**
 * Toolkit for ferron analysis
 * 
 * COMPILATION:
 * GW laptop
 *      gcc tkpca_cloud_size.c -o tkpca_cloud_size -O3 -lm
 * */

#define file_operationl( cmd )                                                  \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "FILE ERROR:: cannot execute: %s\n" , #cmd);          \
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        return( EXIT_FAILURE ) ;                                                \
    } }

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include "pca_settings.h"
#include "pca_macro.h"
#include "pca_utils.h"

inline double fbeta(double E, double beta)
{
    double bE=beta*E;
    if(bE>50.) return 0.0; // to avoid numerical problems
    else if(bE<-50.) return 1.0; // to avoid numerical problems
    else return 1.0/(exp(bE)+1.0);
}
#include "pca_io.h"

int main( int argc , char ** argv ) 
{
    char inprefix[256] =  "/home/gabrielw/repository/higgs/trap_ho/g0.05.r1";
    char file_name[256];
    int nx, ny, nz, nxyz;
    int ix, iy, iz, ixyz;
    int ierr;
    int nom;
    double dx, dy, dz, dt, t0, eF;
    
    // read delta file
    sprintf(file_name, "%s_delta.dpca", inprefix);
    printf("# TKPCA: Reading file `%s`\n", file_name);
    file_operationl( read_measurement_file_header(file_name,  &nx, &ny, &nz, &dx, &dy, &dz, &eF, &t0, &dt, &nom) );
    printf("# TKPCA: Lattice for input data: %d x %d x %d with lattice spacing %.2f x %.2f x %.2f\n", nx, ny, nz, dx, dy, dz);
    
    nxyz = nx*ny*nz;
    
    // arrays
    
    // order parameter
    double complex *delta;
    cppmallocl(delta, nxyz, double complex);

    // density of spin-up particles
    double *rho_a; 
    cppmallocl(rho_a, nxyz, double);
    
    // density of spin-down particles
    double *rho_b; 
    cppmallocl(rho_b, nxyz, double);
    
    int inom;
    for(inom=0; inom<nom; inom++) // for each frame
    {
        double time = t0+dt*inom; // time
        
        sprintf(file_name, "%s_delta.dpca", inprefix);
        file_operationl( read_measurement_entry(file_name, inom, delta, sizeof(double complex)*nxyz) ); 
        sprintf(file_name, "%s_density_a.dpca", inprefix);
        file_operationl( read_measurement_entry(file_name, inom, rho_a, sizeof(double)*nxyz) ); 
        sprintf(file_name, "%s_density_b.dpca", inprefix);
        file_operationl( read_measurement_entry(file_name, inom, rho_b, sizeof(double)*nxyz) );  
        
        ixyz=0;
        double N=0.0; // particle number
        
        for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++)
        {
            double x = dx*(ix-nx/2);
            double y = dy*(iy-ny/2);
            double z = dz*(iz-nz/2);
            
            
            N += (rho_a[ixyz]+rho_b[ixyz])*dx*dy*dz;
            
            // add here X(t), Y(t), Z(t)
            
            ixyz++;
            
        }
        
        printf("%12.6f %12.6f\n", time*eF, N);

    }
    
    return 0;
}
