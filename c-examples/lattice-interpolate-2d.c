/**
 * W-interp Library 
 * @see https://gitlab.fizyka.pw.edu.pl/wtools/winterp
 * 
 * This exmaple shows how to interpolate data defined in the lattice to new resolution. 
 * 
 * It is assumes that the data is written in WDATA format. 
 * */ 

// winterp lib
#include "winterp.h"

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <complex.h>

int main()
{
    int ierr;
    
    // WDATA: prepare for reading
    
    // create metadata handler
    wdata_metadata md;
    
    // read metadata from file
    ierr = wdata_parse_metadata_file("sample.wtxt", &md);
    if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;}
    
    // allocate memory for variables
    int bdim = wdata_get_blocklength(&md); // get block size
    
    double *dataR = (double *)malloc(bdim * sizeof(double)); // real data
    Complex *dataC= (double complex *)malloc(bdim * sizeof(double complex)); // complex data
    double *dataV= (double *)malloc(3*bdim * sizeof(double)); // vector data, factor 3 accounts for three compoments of vector variable

    // WDATA: prepare output dataset
    wdata_metadata mdout;
    mdout=md; // take copy
    // ... and change lattice while keeping the volumne ...
    mdout.NX=160;
    mdout.NY=160;
    mdout.DX=(md.DX*md.NX)/mdout.NX;
    mdout.DY=(md.DY*md.NY)/mdout.NY;
    // ... and change prefix for file ...
    wdata_setprefix(&mdout, "sample-winterp");
    // ... and clear database if already exists ...
    wdata_clear_database(&mdout);
    
    printf("# **********************  INPUT LATTICE **********************\n");
    printf("# LATTICE: %d x %d x %d\n", md.NX, md.NY, md.NZ);
    printf("# SPACING: %f x %f x %f\n", md.DX, md.DY, md.DZ);
    printf("# VOLUME : %f x %f x %f\n", md.DX*md.NX, md.DY*md.NY, md.DZ*md.NZ);
    printf("# ********************** OUTPUT LATTICE **********************\n");
    printf("# LATTICE: %d x %d x %d\n", mdout.NX, mdout.NY, mdout.NZ);
    printf("# SPACING: %f x %f x %f\n", mdout.DX, mdout.DY, mdout.DZ);
    printf("# VOLUME : %f x %f x %f\n", mdout.DX*mdout.NX, mdout.DY*mdout.NY, mdout.DZ*mdout.NZ);
    printf("# ************************************************************\n");
    
    
    
    
//     int icycle;
//     for(icycle=0; icycle<md.cycles; icycle++)
//     {
//         printf("Processing cycle %d\n", icycle);
//         
//         ierr = wdata_read_cycle(&md, "density_a", icycle, dataR);
//         if(ierr!=0) { printf("ERROR: Cannot read density_a!\n"); return 1;}
//         
//         // you can also do like this (note density_b is link to density_a)
//         // ierr = wdata_read_cycle(&md, "density_b", icycle, dataR);
//         
//         ierr = wdata_read_cycle(&md, "delta", icycle, dataC);
//         if(ierr!=0) { printf("ERROR: Cannot read delta!\n"); return 1;}
//         
//         ierr = wdata_read_cycle(&md, "current_a", icycle, dataV);
//         if(ierr!=0) { printf("ERROR: Cannot read current_a!\n"); return 1;}        
//     }
// 
//     // extract constants
//     double eF = wdata_getconst_value(&md, "eF");
//     printf("Constant eF=%f\n", eF);
    
    return 0;
} 
