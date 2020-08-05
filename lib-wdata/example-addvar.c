/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Code demonstrates how to add variable to existing data set
 * 
 * New variable is: 
 *   w_a = current_a / sqrt(density_a) (vector type)
 * */ 

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 


int main()
{
    int ierr;
    
    // create metadata handler
    wdata_metadata md;
    
    // read metadata from file
    ierr = wdata_parse_metadata_file("test.wtxt", &md);
    if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;}
    
    // storage for real and vector data
    int bdim = wdata_get_blocklength(&md); // get block size
    
    double *density_a;                          
    density_a = (double *) malloc(sizeof(double)*bdim); // allocate memory
    if(density_a==NULL) {printf("Cannot allocate density_a!\n"); return 1;}
    
    double *current_a;                          
    current_a = (double *) malloc(sizeof(double)*bdim*3); // allocate memory
    if(current_a==NULL) {printf("Cannot allocate current_a!\n"); return 1;}
 
    // add variable to metadata file
    wdata_variable w_a = {"w_a", "vector", "none"};
    wdata_add_variable(&md, &w_a); // add to metadata structure
    wdata_add_var_to_metadata_file("test.wtxt", &w_a); // add entry to metadata file
 
    // create binary file for variable
    double *vx = current_a+0*bdim; // pointer to x-component
    double *vy = current_a+1*bdim; // pointer to y-component
    double *vz = current_a+2*bdim; // pointer to z-component
    
    int icycle, ixyz;
    for(icycle=0; icycle<md.cycles; icycle++) // for each cycle
    {
        // load data
        ierr = wdata_read_cycle(&md, "density_a", icycle, density_a);
        if(ierr!=0) { printf("ERROR: Cannot read density_a [ierr=%d]!\n", ierr); return 1;}
        
        ierr = wdata_read_cycle(&md, "current_a", icycle, current_a);
        if(ierr!=0) { printf("ERROR: Cannot read current_a [ierr=%d]!\n", ierr); return 1;} 
        
        for(ixyz=0; ixyz<bdim; ixyz++) // for each lattice point
        {
            // update in place
            vx[ixyz] = vx[ixyz]/sqrt(density_a[ixyz]+1.0e-16); // +1.0e-16 to avoid division by zero
            vy[ixyz] = vy[ixyz]/sqrt(density_a[ixyz]+1.0e-16); // +1.0e-16 to avoid division by zero
            vz[ixyz] = vz[ixyz]/sqrt(density_a[ixyz]+1.0e-16); // +1.0e-16 to avoid division by zero
        }
        
        // write to binary file
        ierr = wdata_write_cycle(&md, "w_a", current_a);
        if(ierr!=0) { printf("ERROR: Cannot write w_a [ierr=%d]!\n", ierr); return 1;} 
    }
    
    return 0;
}
