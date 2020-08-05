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

typedef std::complex<double> Complex;


#define SQRT_2PI 2.506628274631
double function_x(double x, double sigma)
{
    return 1./(sigma*SQRT_2PI) * exp(-0.5*x*x/(sigma*sigma));
}

double function_xyz(double x, double y, double z, double time)
{
    double val=0.0;
    double sigmax = 2.0 + 0.20*time;
    double sigmay = 3.0 + 0.15*time;
    double sigmaz = 4.0 + 0.10*time;
    
    val = function_x(x, sigmax)*function_x(y, sigmay)*function_x(z, sigmaz);
    
    return val;
}

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }
    
int main()
{    
    // create metadata handler
    wdata_metadata md;
    
    // Lattice
    md.datadim=3; // 3D data
    md.NX=24; md.NY=28; md.NZ=32;
    md.DX=1.0; md.DY=1.0; md.DZ=1.0;
    
    strcpy(md.prefix, "test");
    md.t0=0.0; md.dt=1.0;
    
    // add variables to data set
    // for each variable binary file of name `prefix_`varname`.wdat will be created
    wdata_variable vdensity_a = {"density_a", "real", "none"};
    wdata_add_variable(&md, &vdensity_a);
    
    wdata_variable vdelta = {"delta", "complex", "none"};
    wdata_add_variable(&md, &vdelta);
    
    wdata_variable vcurrent_a = {"current_a", "vector", "none"};
    wdata_add_variable(&md, &vcurrent_a);
    
    // add links to data sets
    // links are alternative names of the same variable
    wdata_link ldensity_b = {"density_b", "density_a"};
    wdata_add_link(&md, &ldensity_b);
    
    wdata_link lcurrent_b = {"current_b", "current_a"};
    wdata_add_link(&md, &lcurrent_b);
    
    // add constants 
    wdata_const lconst_eF = {"eF", 0.5};
    wdata_add_const(&md, &lconst_eF);
    
    wdata_const lconst_kF = {"kF", 1.0};
    wdata_add_const(&md, &lconst_kF);
    
    // just in case - clear data sets
    // it removes binary files if they alredy exists
    wdata_clear_database(&md);

    // generate some artifical data 
    int bdim = wdata_get_blocklength(&md); // get block size
    
    double *dataR;         // real data
    Complex *dataC; // complex data
    double *dataV;         // vector data 
    cppmallocl(dataR,bdim,double);
    cppmallocl(dataC,bdim,Complex);
    cppmallocl(dataV,bdim*3,double); // factor 3 accounts for three compoments of vector variable
    
    
    int ncycles=10; // number of cycles to be generated

    for(int i=0; i<ncycles; i++)
    {
        int ixyz=0;
        for(int ix=0; ix<md.NX; ix++) for(int iy=0; iy<md.NY; iy++) for(int iz=0; iz<md.NZ; iz++)
        {
            // coordinate decomposition
            double x = md.DX*(ix-md.NX/2);
            double y = md.DY*(iy-md.NY/2);
            double z = md.DZ*(iz-md.NZ/2);
            
            // time decomposition
            double time = md.t0 + md.dt*i;
            
            // add data to array
            
            dataR[ixyz] = function_xyz(x,y,z,time);
            dataC[ixyz] = Complex(dataR[ixyz], 0.0);
            
            // vector - solid body rotation around z axis 
            dataV[ixyz+0*bdim] = -1.0*y;
            dataV[ixyz+1*bdim] =  1.0*x;
            dataV[ixyz+2*bdim] =  0.0;
            
            ixyz++;
        }
        
        int ierr;
        
        // add info about new cycle to metadata
        wdata_add_cycle(&md);
        
        // and add cycle to binary sets
        ierr = wdata_write_cycle(&md, "density_a", dataR);
        // alternatively you can use:
        // ierr=wdata_add_datablock(&md, &vdensity_a, dataR);
        if(ierr!=0) { printf("ERROR: Cannot add density_a!\n"); return 1;}

        ierr = wdata_write_cycle(&md, "delta", dataC);
        if(ierr!=0) { printf("ERROR: Cannot add delta!\n"); return 1;}

        ierr = wdata_write_cycle(&md, "current_a", dataV);
        if(ierr!=0) { printf("ERROR: Cannot add delta!\n"); return 1;}
    }
    
    // write metadata file (with default name)
    wdata_write_metadata_to_file(&md, "");
    
    return 0;
}


