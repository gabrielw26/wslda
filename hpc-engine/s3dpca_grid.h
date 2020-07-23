// Author: Gabriel Wlazlowski

// Structure that stores parameters about BALCS grid

#ifndef __S3DPCA_GRID__
#define __S3DPCA_GRID__

typedef struct
{
    int ip; // grid identifier
    int iq; // grid identifier
    int nip; // size of local matrix
    int niq; // size of local matrix
    int p; // grid size
    int q; // grid size
    int mb; // block size
    int nb; // block size
    
} metadata_s3dpca_grid;

#endif
