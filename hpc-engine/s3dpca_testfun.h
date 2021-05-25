// Author: Gabriel Wlazlowski

// This file provides set of test functions

int test_matrix_elements(double complex *A, int nip, int niq, int mb, int nb, int ip, int iq, int p, int q)
{
    int ci, ri; // column and row iterator - global index
    int li, lj, ij; // local indices
    int ZERO = 0, ti;   
    int ixyz1, ixyz2;
    
    for(lj=0; lj<niq; lj++) // column-major iteration: over local index (column)
    {
        for(li=0; li<nip; li++) // over local index (row)
        {
            ij=li + nip * lj; // local index of elemnt
            
            ixyz1 = li+1; ixyz2=lj+1; // conversion to fortran standard
            // find indices in global matrix: row and colummn
            ri = indxl2g_( &ixyz1, &mb, &ip, &ZERO, &p )-1; // back to C standard
            ci = indxl2g_( &ixyz2, &nb, &iq, &ZERO, &q )-1; // back to C standard            
            
            A[ij]=1.*ri+I*ci;
            
        } // for(li=0; li<nip; li++)
    } // for(lj=0; lj<niq; lj++)
    return 0; // no error
}

int test_matrix_print(double complex *A, int nip, int niq, int mb, int nb, int ip, int iq, int p, int q, int print)
{
    int ci, ri; // column and row iterator - global index
    int li, lj, ij; // local indices
    int ZERO = 0, ti;   
    int ixyz1, ixyz2;
    
    for(lj=0; lj<niq; lj++) // column-major iteration: over local index (column)
    {
        for(li=0; li<nip; li++) // over local index (row)
        {
            ij=li + nip * lj; // local index of elemnt
            
            ixyz1 = li+1; ixyz2=lj+1; // conversion to fortran standard
            // find indices in global matrix: row and colummn
            ri = indxl2g_( &ixyz1, &mb, &ip, &ZERO, &p )-1; // back to C standard
            ci = indxl2g_( &ixyz2, &nb, &iq, &ZERO, &q )-1; // back to C standard            
            
            if(print) wprintf("A[%4d,%4d]=(%8.2f,%8.2f)\n", ri, ci, creal(A[ij]),  cimag(A[ij]));
            
        } // for(li=0; li<nip; li++)
    } // for(lj=0; lj<niq; lj++)
    return 0; // no error
}
