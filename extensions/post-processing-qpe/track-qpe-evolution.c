/**
 * W-SLDA Toolkit
 *
 * This script shows how to read data from qpe file.
 * To generate this file you neeed to activate option `STORE_QPE` in predefines
 *
 * gcc track-qpe-evolution.c -o track-qpe-evolution -lm
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

    if(argc!=4)
    {
        printf("Usage: %s prefix nwf idx\n", argv[0]);
        printf("\t- prefix: of target file, code is reading from file `prefix`_qpe.dpca\n");
        printf("\t- nwf: number of wave-functions in the simulations, take it from `prefix`.wlog\n");
        printf("\t- idx: state index to be traced\n");

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

    int idx = atoi(argv[3]);

    // read header
    int nx, ny, nz;
    double dx, dy, dz;
    double eF, t0, dt;
    int i,nom; // number of measurments
    double qpe, time;

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

    // header
    printf("# WORKING WITH FILE: `%s`\n", file_name);
    printf("# eF: %f\n", eF);
    printf("# nwf: %d\n", nwf);
    printf("# idx: %d\n", idx);
    printf("# %16s %16s %16s %16s\n", "cycle", "t*eF", "E_n", "E_n/eF");


    // go to correct record

    for(i=0; i<nom; i++)
    {
        fseek(fqpe,sizeof(int)*5+sizeof(double)*6+sizeof(double)*i*nwf+sizeof(double)*idx, SEEK_SET); // set pointer

        // read from file
        fread(&qpe,sizeof(double),1,fqpe);

        time=t0+dt*i;
        printf("%16d %16.10f %16.10f %16.10f\n", i, time*eF, qpe, qpe/eF);
    }

    // close binary file
    fclose(fqpe);

     /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
