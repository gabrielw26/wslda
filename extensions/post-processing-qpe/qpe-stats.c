/**
 * W-SLDA Toolkit
 * 
 * This script computes basic statistics for qpes.
 * It is used for checking stability of the ABM integrator.
 * 
 * gcc qpe-stats.c -o qpe-stats -I$WSLDA/hpc-engine -lm 
 * */  

#define TOLERANCE_DIFF 1.0e-6

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
    int cnt=0;
    int iwf;
    
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
    
    // allocate arrays
    double *qpe;
    if ( ( qpe = (double *) malloc( (nwf) * sizeof( double ) ) ) == NULL )
    {
        printf("Error: Cannot allocate array qpe!\n");
        return( EXIT_FAILURE ) ;
    }
    
    double *min;
    if ( ( min = (double *) malloc( (nwf) * sizeof( double ) ) ) == NULL )
    {
        printf("Error: Cannot allocate array min!\n");
        return( EXIT_FAILURE ) ;
    }
    
    double *max;
    if ( ( max = (double *) malloc( (nwf) * sizeof( double ) ) ) == NULL )
    {
        printf("Error: Cannot allocate array max!\n");
        return( EXIT_FAILURE ) ;
    }
    
    // reset buffers
    for(iwf=0; iwf<nwf; iwf++) min[iwf]= 9.9e12;
    for(iwf=0; iwf<nwf; iwf++) max[iwf]=-9.9e12;
    
    for(i=1; i<nom; i++)
    {
//         printf("# Processing cycle %d.\n", i);
        double time = t0+dt*i;
        fread(qpe,sizeof(double),nwf,fqpe);
        
        for(iwf=0; iwf<nwf; iwf++)
        {
            if(qpe[iwf]>max[iwf]) max[iwf]=qpe[iwf];
            if(qpe[iwf]<min[iwf]) min[iwf]=qpe[iwf];
        }
        
    }
    
    // read the firts entry again for reference
    fseek(fqpe,sizeof(int)*5 + sizeof(double)*6 + sizeof(double)*nwf, SEEK_SET); 
    fread(qpe,sizeof(double),nwf,fqpe);
    
    // close binary file
    fclose(fqpe);
    
    // make report
    printf("# ONLY STATES WHERE diff=|max[qpe(t)]-min[qpe(t)]|>%g ARE PRINTED\n", TOLERANCE_DIFF);
    printf("%6s %16s %16s %16s %16s\n", "iwf", "qpe(t=0)", "min[qpe(t)]", "max[qpe(t)]", "diff");
    double max_diff=0.0;
    for(iwf=0; iwf<nwf; iwf++)
    {
        double diff = fabs(max[iwf]-min[iwf]);
        if(diff>max_diff) max_diff=diff;
        if(diff>TOLERANCE_DIFF)
            printf("%6d %16.10g %16.10g %16.10g %16.10g\n", iwf, qpe[iwf], min[iwf], max[iwf], diff);
    }
    
    printf("# MAX DIFF=%16.10g\n", max_diff);
    
     /* Arrays will be cleared automatically */
    return( EXIT_SUCCESS ) ;
}
