/**
 * W-SLDA Toolkit
 * 
 * This file provides simple code that can be used for testing 
 * of correctness of implementation of external potentials. 
 * 
 * The potentials are computed according problem-definition.h file
 * and saved to W-DATA format.
 * Next, they can be inspected by VisIt. 
 * 
 * Copy this file to your project folder and compile using (select CODEDIM from 1, 2 or 3):
 *    gcc -std=gnu99 st-test-external-potentials.c -I. -I$WSLDA/hpc-engine -I$WSLDA/lib/wdata/c -L$WSLDA/lib/wdata -lwdata -o st-test-external-potentials -lm -lfftw3 -DCODEDIM=2
 * 
 * NOTE: you need before generate wdata lib for C compiler:
 *    cd $WSLDA/lib/wdata
 *    make lib
 * */   

// Standard libraries
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

// W-DATA Format
// Must be before wslda_toolkit.h !
#include "wdata.h"

// W-SLDA Toolkit API
// Note: it loads predefines.h!
#include "wslda_toolkit.h"

// Load my problem-definition.h
#include "problem-definition.h"

int main( int argc , char ** argv ) 
{
    // Read of input parameters
    char execcmd[ 256 ];
    int i,j;
    wsldapid=0; // global variable
    
    strcpy( execcmd , argv[ 0 ] ) ;
    for( i = 1 ; i < argc ; i++ ) 
    {
        strcat( execcmd , " " ) ; 
        strcat( execcmd , argv[ i ] ) ;
    }
    
    // parse command line
    i = readcmd( argc , argv ) ;
    if( i == -1 )
    {
        printf( "TERMINATING! NO INPUT FILE.\n" ) ;
        return( EXIT_FAILURE ) ;
    }
        
    // Read input file
    // Info from file is loaded into metadata structure
    j = parse_input_file(argv[i]);
    if ( j == 0 )
    {
        printf("PROBLEM WITH INPUT FILE: `%s`.\n" , argv[ i ] ) ;
        return( EXIT_FAILURE ) ;      
    }
        
    // create metadata handler
    wdata_metadata wdmd;
    wdata_set_lattice(&wdmd, NX, NY, NZ, DX, DY, DZ);
    strcpy(wdmd.prefix, input->outprefix);
    wdmd.t0=0.0; wdmd.dt=1.0;
    // add variables to data set
    // for each variable binary file of name `prefix_`varname`.wdat will be created
    wdata_variable v1 = {"V_ext_a", "real", "none", "wdat"}; wdata_add_variable(&wdmd, &v1);
    wdata_variable v2 = {"V_ext_b", "real", "none", "wdat"}; wdata_add_variable(&wdmd, &v2);
    wdata_variable v3 = {"delta_ext", "complex", "none", "wdat"}; wdata_add_variable(&wdmd, &v3);
    wdata_variable v4 = {"velocity_ext_a", "vector", "none", "wdat"}; wdata_add_variable(&wdmd, &v4);
    wdata_variable v5 = {"velocity_ext_b", "vector", "none", "wdat"}; wdata_add_variable(&wdmd, &v5);
    // just in case - clear data sets
    // it removes binary files if they alredy exists
    wdata_clear_database(&wdmd);


    int lNX, lNY, lNZ;
#if CODEDIM==1
    printf("# TESTING OF POTENTIALS FOR CODE: ST-WSLDA-1D\n");
    lNX=NX; lNY=1; lNZ=1;
    wdmd.datadim=1; // 1D data
#elif CODEDIM==2
    printf("# TESTING OF POTENTIALS FOR CODE: ST-WSLDA-2D\n");
    lNX=NX; lNY=NY; lNZ=1;
    wdmd.datadim=2; // 2D data
#else
    printf("# TESTING OF POTENTIALS FOR CODE: ST-WSLDA-3D\n");
    lNX=NX; lNY=NY; lNZ=NZ;
    wdmd.datadim=3; // 3D data
#endif

    // set reference value
    double kF;
    double mu[2] = {1.0, 1.0};
    printf("# mu[SPINA]=%f, mu[SPINB]=%f", mu[SPINA], mu[SPINB]);
    if(input->referencekF>0.0) 
    {
        kF = input->referencekF;
        printf("# kF=%f (TAKEN FROM input)\n", kF);
    }
    else
    {
        kF = 1.0;
        printf("# kF=%f (DEFAULT VALUE!, YOU CAN SET IT VIA input)\n", kF);
    }
    
    // load extra data if needed
    void *extra_data = NULL;
    size_t extra_data_size;
    int ierr=0;
    extra_data_size = get_extra_data_size(input->params);
    if(extra_data_size>0)
    {
        printf("# EXTRA_DATA IS ACTIVE.\n");
        printf("# ALLOCATING EXTRA_DATA OF SIZE %ld B.\n", extra_data_size); fflush(stdout);
        if ( ( extra_data = (void *) malloc( extra_data_size ) ) == NULL  )
        {                                                             
            fprintf( stderr , "error: cannot malloc()! Exiting!\n") ; 
            fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; 
            /* Arrays will be cleared automatically */
            return( EXIT_FAILURE ) ; 
        }
        
        ierr = load_extra_data(extra_data_size, extra_data, input->params);
        if(ierr!=0) 
        {
            fprintf( stderr , "error[%d]: cannot load_extra_data! Exiting!\n", ierr) ; 
            /* Arrays will be cleared automatically */
            return( EXIT_FAILURE ) ;
        }
    }
    

    int it;
    int ix, iy, iz, ixyz;
    int lNXYZ = lNX*lNY*lNZ;
    double dc_params[MAX_USER_PARAMS];
    
    // allocate array for results
    double *V_ext_a, *V_ext_b;
    double *velocity_ext_a, *velocity_ext_b;
    double complex *delta_ext_;
    
    cppmallocl(V_ext_a, lNXYZ, double);
    cppmallocl(V_ext_b, lNXYZ, double);
    cppmallocl(velocity_ext_a, 3*lNXYZ, double);
    cppmallocl(velocity_ext_b, 3*lNXYZ, double);
    cppmallocl(delta_ext_, lNXYZ, double complex);
    
    for(it=0; it<input->maxiters; it++)
    {
        printf("# it=%d\n", it);
        for(i=0; i<MAX_USER_PARAMS; i++) dc_params[i]=input->params[i];
        process_params(dc_params, &kF, mu, extra_data_size, extra_data);
        
        // Fill with data 
        ixyz=0;
        for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++) 
        {
            V_ext_a[ixyz] = v_ext(ix, iy, iz, it, SPINA, dc_params, extra_data_size, extra_data);
            V_ext_b[ixyz] = v_ext(ix, iy, iz, it, SPINB, dc_params, extra_data_size, extra_data);
            delta_ext_[ixyz] = delta_ext(ix, iy, iz, it, 0.0+I*0.0, dc_params, extra_data_size, extra_data);
            velocity_ext_a[ixyz+0*lNXYZ] = velocity_ext(ix, iy, iz, it, SPINA, XAXIS, dc_params, extra_data_size, extra_data);
            velocity_ext_a[ixyz+1*lNXYZ] = velocity_ext(ix, iy, iz, it, SPINA, YAXIS, dc_params, extra_data_size, extra_data);
            velocity_ext_a[ixyz+2*lNXYZ] = velocity_ext(ix, iy, iz, it, SPINA, ZAXIS, dc_params, extra_data_size, extra_data);
            velocity_ext_b[ixyz+0*lNXYZ] = velocity_ext(ix, iy, iz, it, SPINB, XAXIS, dc_params, extra_data_size, extra_data);
            velocity_ext_b[ixyz+1*lNXYZ] = velocity_ext(ix, iy, iz, it, SPINB, YAXIS, dc_params, extra_data_size, extra_data);
            velocity_ext_b[ixyz+2*lNXYZ] = velocity_ext(ix, iy, iz, it, SPINB, ZAXIS, dc_params, extra_data_size, extra_data);            
            ixyz++;
        }
        
        // add info about new cycle to metadata
        wdata_add_cycle(&wdmd);
        file_operationl( wdata_write_cycle(&wdmd, "V_ext_a", V_ext_a) );
        file_operationl( wdata_write_cycle(&wdmd, "V_ext_b", V_ext_b) );
        file_operationl( wdata_write_cycle(&wdmd, "delta_ext", delta_ext_) );
        file_operationl( wdata_write_cycle(&wdmd, "velocity_ext_a", velocity_ext_a) );
        file_operationl( wdata_write_cycle(&wdmd, "velocity_ext_b", velocity_ext_b) );
    }
    
    // write metadata file (with default name)
    wdata_write_metadata_to_file(&wdmd, "");
    printf("# RESULTS WRITTEN TO: `%s.wtxt`\n", input->outprefix);
    
    /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
