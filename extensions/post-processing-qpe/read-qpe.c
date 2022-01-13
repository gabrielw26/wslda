/**
 * W-SLDA Toolkit
 * 
 * This script shows how to read data from qpe file.
 * To generate this file you neeed to activate option `STORE_QPE` in predefines
 * 
 * gcc read-qpe.c -o read-qpe -I$WSLDA/hpc-engine -lm 
 * */  
// Standard libraries
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>


int main( int argc , char ** argv ) 
{
    
    if(argc!=3)
    {
        printf("Usage: %s prefix nwf\n", argv[0]);
        printf("\t- prefix: of target file, code is reading from file `prefix`_qpe.dpca\n");
        printf("\t- nwf: number of wave-functions in the simulations, take it from `prefix`.wlog\n");
        
        return( EXIT_FAILURE ) ; 
    }    
    
    // open binary file
    int nwf=atoi(argv[2]);
    char file_name[512];
    sprintf(file_name,"%s_qpe.dpca", argv[1]);
    FILE * fqpe = fopen(file_name, "rb");
    if(fqpe==NULL)
    {
        printf("Error: Cannot open file %s!\n", file_name);
        return( EXIT_FAILURE ) ;
    }
    
    // read header
    int nx, ny, nz;
    double dx, dy, dz;
    double eF, t0, dt;
    int i,nom; // number of measurments
    
    fread(&nom,sizeof(int),   1,fqpe);
    fread(&nx ,sizeof(int),   1,fqpe);
    fread(&ny ,sizeof(int),   1,fqpe);
    fread(&nz, sizeof(int),   1,fqpe);
    fread(&dx ,sizeof(double),1,fqpe);
    fread(&dy ,sizeof(double),1,fqpe);
    fread(&dz, sizeof(double),1,fqpe);
    fread(&eF, sizeof(double),1,fqpe);
    fread(&t0, sizeof(double),1,fqpe);
    fread(&dt, sizeof(double),1,fqpe);
    fread(&nom,sizeof(int),   1,fqpe);
    fseek(fqpe,sizeof(double)*nwf, SEEK_CUR); // skip entry i==0, thera are only zeros there
    
    // allocate array for quasiparticle energies
    double *qpe;
    if ( ( qpe = (double *) malloc( (nwf) * sizeof( double ) ) ) == NULL )
    {
        printf("Error: Cannot allocate array qpe!\n");
        return( EXIT_FAILURE ) ;
    }
    
    for(i=1; i<nom; i++)
    {
        double time = t0+dt*i;
        fread(qpe,sizeof(double),nwf,fqpe);
        
        // here you can process your qpe
        // ...
        // example: count number of states with energy smaller than 0.5*eF
        int cnt=0;
        int iwf;
        for(iwf=0; iwf<nwf; iwf++) if(fabs(qpe[iwf])<0.5*eF) cnt++;
        printf("%12.6f %6d\n", time*eF, cnt);
        
        // ...
    }
    
    // close binary file
    fclose(fqpe);
    
     /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
