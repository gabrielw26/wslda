/**
 * W-SLDA Toolkit
 *
 * This script shows how to read data from qpe file.
 * To generate this file you neeed to activate option `STORE_QPE` in predefines
 *
 * gcc read-qpe-below-treshold.c -o read-qpe-below-treshold -lm
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

    if(argc!=5)
    {
        printf("Usage: %s prefix nwf cycle treshold\n", argv[0]);
        printf("\t- prefix: of target file, code is reading from file `prefix`_qpe.dpca\n");
        printf("\t- nwf: number of wave-functions in the simulations, take it from `prefix`.wlog\n");
        printf("\t- cycle: energies for given cycle will be printed\n");
        printf("\t- treshold: only states with |E_n|<treshold will be printed\n");

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

    int cycle = atoi(argv[3]);
    double treshold = atof(argv[4]);

    // read header
    int nx, ny, nz;
    double dx, dy, dz;
    double eF, t0, dt;
    int i,nom; // number of measurments
    double qpe;

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

    if(cycle>=nom)
    {
        printf("ERROR: Data set has only %d cycles! You have requested cycle=%d\n", nom, cycle);
        return( EXIT_FAILURE ) ;
    }

    // go to correct record
    fseek(fqpe,sizeof(double)*nwf*cycle, SEEK_CUR); // skip entry i==0, thera are only zeros there

    // header
    printf("# WORKING WITH FILE: `%s`\n", file_name);
    printf("# eF: %f\n", eF);
    printf("# nwf: %d\n", nwf);
    printf("# cycle: %d\n", cycle);
    printf("# time=%f [time*eF=%f]\n", t0+cycle*dt, (t0+cycle*dt)*eF);
    printf("# PRINTING STATES WITH |En|<%f [|En/eF|<%f]\n", treshold, treshold/eF);
    printf("# %16s %16s %16s\n", "n", "E_n", "E_n/eF");
    for(i=0; i<nwf; i++)
    {
        // read from file
        fread(&qpe,sizeof(double),1,fqpe);

        if(fabs(qpe)<treshold) printf("%16d %16.10f %16.10f\n", i, qpe, qpe/eF);
    }

    // close binary file
    fclose(fqpe);

     /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
