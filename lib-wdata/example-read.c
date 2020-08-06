/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Code demonstrates how to read wdata files
 * */

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <complex>

typedef std::complex<double> Complex; 

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }
    
int main()
{
    int ierr;
    
    // create metadata handler
    wdata_metadata md;
    
    // read metadata from file
    ierr = wdata_parse_metadata_file("test.wtxt", &md);
    if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;}
    
    // print metadata file
    wdata_print_metadata(&md, stdout);
    
    // allocate memory for variables
    int bdim = wdata_get_blocklength(&md); // get block size
    
    double *dataR;         // real data
    Complex *dataC;        // complex data
    double *dataV;         // vector data 
    cppmallocl(dataR,bdim,double);
    cppmallocl(dataC,bdim,Complex);
    cppmallocl(dataV,bdim*3,double); // factor 3 accounts for three compoments of vector variable
    
    int icycle;
    for(icycle=0; icycle<md.cycles; icycle++)
    {
        printf("Processing cycle %d\n", icycle);
        
        ierr = wdata_read_cycle(&md, "density_a", icycle, dataR);
        if(ierr!=0) { printf("ERROR: Cannot read density_a!\n"); return 1;}
        
        // you can also do like this (note density_b is link to density_a)
        // ierr = wdata_read_cycle(&md, "density_b", icycle, dataR);
        
        ierr = wdata_read_cycle(&md, "delta", icycle, dataC);
        if(ierr!=0) { printf("ERROR: Cannot read delta!\n"); return 1;}
        
        ierr = wdata_read_cycle(&md, "current_a", icycle, dataV);
        if(ierr!=0) { printf("ERROR: Cannot read current_a!\n"); return 1;}        
    }

    // extract constants
    double eF = wdata_getconst_value(&md, "eF");
    printf("Constant eF=%f\n", eF);
    
    return 0;
} 
