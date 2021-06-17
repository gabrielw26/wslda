// Author: Gabriel Wlazlowski
// Date: 06-09-2016

// File contains functions for files managment

#ifndef __PCA_IO__
#define __PCA_IO__

#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#define PCA_PRECISION 8

// Deletes file with given name
// return 1 - file deleted, 0-file dosen't exist, -1 - problem with the file
int urm( const char * fn )
{
    int i ;
    i = unlink( ( const char * ) fn );

    if ( i == 0 )
        return 1;
    else if (i==-1 && errno== ENOENT)
        return 0;
    else
        return -1;

}

// Checks if file exists
int exists(const char *filename) 
{  
    return !access(filename, F_OK);  
}

int create_directory(const char *dirname)
{
    int err = mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    if(err==0 || errno==EEXIST) return WSLDA_OK;
    else return WSLDA_ERR_CANNOT_CREATE_DIR;
}

/**
 * return
 * 0 - OK, else - problem
 * */
int create_measurement_file_with_header(const char * file_name, 
                                        int nx, int ny, int nz, 
                                        double dx, double dy, double dz,
                                        double eF, double t0, double dt
                                       )
{
    // We do not overwrite!   
    if(exists(file_name)) 
    {
        if(md.overwrite==0) return -1;
        else urm(file_name);
    }
    
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_CREAT | O_RDWR  , fd_mode ) ) == -1 ) return -2;
    
    long int bytes_written ;
    int werr=0;
    int prec=PCA_PRECISION;
    int tmp=0;
    
#define mio_hentry_wrt(var,type,werrval) \
    if (werr==0 && ( bytes_written = write( fd , ( const void * ) &var , sizeof( type ) ) ) != sizeof( type ) ) werr=werrval;
    
    mio_hentry_wrt(prec,int,-3);
    mio_hentry_wrt(nx,int,-4);
    mio_hentry_wrt(ny,int,-5);
    mio_hentry_wrt(nz,int,-6);
    mio_hentry_wrt(dx,double,-7);
    mio_hentry_wrt(dy,double,-8);
    mio_hentry_wrt(dz,double,-9);
    mio_hentry_wrt(eF,double,-10);
    mio_hentry_wrt(t0,double,-11); // start time value 
    mio_hentry_wrt(dt,double,-12); // start time value 
    mio_hentry_wrt(tmp,int,-13); // This field shows current number of measurements in the file
 
    // Set shifts for measuremnts counter
#define MIO_CNT_INTS 4
#define MIO_CNT_double 6
    
    // Close file
    close( fd );
        
    return werr;    
}


/**
 * return
 * 0 - OK, else - problem
 * */
int read_measurement_file_header(const char * file_name, 
                                        int *nx, int *ny, int *nz, 
                                        double *dx, double *dy, double *dz,
                                        double *eF, double *t0, double *dt,
                                        int *number_of_entries
                                       )
{
    
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDONLY  , fd_mode ) ) == -1 ) return -1;
    
    int werr=0;
    long int bytes_read ;
    int prec;
    
#define mio_hentry_read(var,type,werrval) \
    if (werr==0 && ( bytes_read = read( fd , ( void * ) var , sizeof( type ) ) ) != sizeof( type ) ) werr=werrval;

    mio_hentry_read(&prec,int,-2);
    if(werr==0 && prec!=PCA_PRECISION) werr=-3;
    mio_hentry_read(nx,int,-4);
    mio_hentry_read(ny,int,-5);
    mio_hentry_read(nz,int,-6);
    mio_hentry_read(dx,double,-7);
    mio_hentry_read(dy,double,-8);
    mio_hentry_read(dz,double,-9);
    mio_hentry_read(eF,double,-10);
    mio_hentry_read(t0,double,-11);
    mio_hentry_read(dt,double,-12);
    mio_hentry_read(number_of_entries,int,-13);
    
    // Close file
    close(fd);
    
    return werr;    
}

/**
 * Adds single entry to the file
 * */
int add_measurement_entry(const char * file_name, void * array, size_t size)
{
    int fd ; /* file descriptor for handling the file or device */
    int werr=0;
    int nr_rec;
    long int bytes_rw ;
        
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDWR  , fd_mode ) ) == -1 ) return -1;
  
    // Update counter 
    int shift=MIO_CNT_INTS*sizeof(int) + MIO_CNT_double*sizeof(double);
    if (werr==0 && (lseek( fd, (long) shift, 0 ))!=shift) werr=-11;
    if (werr==0 && ( bytes_rw = read( fd , ( void * ) &nr_rec , sizeof( int ) ) ) != sizeof( int ) ) werr=-2;
    nr_rec=nr_rec+1;
    if (werr==0 && (lseek( fd, (long) shift, 0 ))!=shift) werr=-12;
    if (werr==0 && ( bytes_rw = write( fd , ( const void * ) &nr_rec , sizeof( int ) ) ) != sizeof( int ) ) werr=-3;
    
    // Write record
    lseek(fd, 0L, 2); // Move pointer to the end
    if (werr==0 && ( bytes_rw = write( fd , array , size ) ) != size ) werr=-5;
    
    // Close file
    close( fd );
    
    return werr;
}

/**
 * reads single entry from the file
 * */
int read_measurement_entry(const char * file_name, int entry_number, void * array, size_t size)
{
    int fd ; /* file descriptor for handling the file or device */
    int werr=0;
    int nr_rec;
    long int bytes_rw ;
        
    mode_t fd_mode = S_IRUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDONLY  , fd_mode ) ) == -1 ) return -1;
  
    // get counter 
    long int shift=MIO_CNT_INTS*sizeof(int) + MIO_CNT_double*sizeof(double);
    if (werr==0 && (lseek( fd, (long) shift, 0 ))!=shift) werr=-11;
    if (werr==0 && ( bytes_rw = read( fd , ( void * ) &nr_rec , sizeof( int ) ) ) != sizeof( int ) ) werr=-2;    
    
    if(entry_number>=nr_rec) werr=-2;
    shift=(MIO_CNT_INTS+1)*sizeof(int) + MIO_CNT_double*sizeof(double) + entry_number*size;
    
    if (werr==0 && (lseek( fd, (long) shift, 0 ))!=shift) werr=-12;
    if (werr==0 && ( bytes_rw = read( fd , array , size ) ) != size ) werr=-5;
    
    // Close file
    close( fd );
    
    return werr;    
}

/**
 * Creates header for checkpoint
 * */
int create_checkpoint_info(const char * file_name, 
                           int nwf, int nx, int ny, int nz, double dx, double dy, double dz,
                           double kF, double mu, double ec, double time
                          )
{ 
    if(exists(file_name)) 
    {
        if(md.overwrite==0) return -1; // We do not overwrite!  
        else urm(file_name);
    }
    
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_CREAT | O_RDWR  , fd_mode ) ) == -1 ) return -2; 
        
    long int bytes_written ;
    int werr=0;
    int prec=PCA_PRECISION;
    
#define cio_info_wrt(var,type,werrval) \
    if (werr==0 && ( bytes_written = write( fd , ( const void * ) &var , sizeof( type ) ) ) != sizeof( type ) ) werr=werrval;
    
    cio_info_wrt(prec,int,-3);
    cio_info_wrt(nwf,int,-13);
    cio_info_wrt(nx,int,-4);
    cio_info_wrt(ny,int,-5);
    cio_info_wrt(nz,int,-6);
    cio_info_wrt(dx,double,-7);
    cio_info_wrt(dy,double,-8);
    cio_info_wrt(dz,double,-9);
    cio_info_wrt(kF,double,-14);
    cio_info_wrt(mu,double,-10); 
    cio_info_wrt(ec,double,-11);
    cio_info_wrt(time,double,-12);
    
    // Close file
    close( fd );
        
    return werr;    
}

/**
 * Creates header for checkpoint
 * */
int create_checkpoint_info_pca(const char * file_name, 
                           int nwf, int nx, int ny, int nz, double dx, double dy, double dz,
                           double kF, double *mu, double ec, double beta
                          )
{ 
    if(exists(file_name)) 
    {
        if(md.overwrite==0) return -1; // We do not overwrite!  
        else urm(file_name);
    }
    
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_CREAT | O_RDWR  , fd_mode ) ) == -1 ) return -2; 
        
    long int bytes_written ;
    int werr=0;
    int prec=PCA_PRECISION;
    
#define cio_info_wrt(var,type,werrval) \
    if (werr==0 && ( bytes_written = write( fd , ( const void * ) &var , sizeof( type ) ) ) != sizeof( type ) ) werr=werrval;
    
    cio_info_wrt(prec,int,-3);
    cio_info_wrt(nwf,int,-13);
    cio_info_wrt(nx,int,-4);
    cio_info_wrt(ny,int,-5);
    cio_info_wrt(nz,int,-6);
    cio_info_wrt(dx,double,-7);
    cio_info_wrt(dy,double,-8);
    cio_info_wrt(dz,double,-9);
    cio_info_wrt(kF,double,-14);
    cio_info_wrt(mu[SPINA],double,-10); 
    cio_info_wrt(mu[SPINB],double,-10);
    cio_info_wrt(ec,double,-11);
    cio_info_wrt(beta,double,-12);
    
    // Close file
    close( fd );
        
    return werr;    
}

/**
 * Reads info from checkpoint file
 * */
int read_checkpoint_info(const char * file_name, 
                           int *nwf, int *nx, int *ny, int *nz, double *dx, double *dy, double *dz,
                           double *kF, double *mu, double *ec, double *time
                          )
{
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDONLY  , fd_mode ) ) == -1 ) return -1;
    
    int werr=0;
    long int bytes_read ;
    int prec;
    
#define cio_info_read(var,type,werrval) \
    if (werr==0 && ( bytes_read = read( fd , ( void * ) var , sizeof( type ) ) ) != sizeof( type ) ) werr=werrval;

    cio_info_read(&prec,int,-2);
    if(werr==0 && prec!=PCA_PRECISION) werr=-3;
    cio_info_read(nwf,int,-4);
    cio_info_read(nx,int,-5);
    cio_info_read(ny,int,-6);
    cio_info_read(nz,int,-7);
    cio_info_read(dx,double,-8);
    cio_info_read(dy,double,-9);
    cio_info_read(dz,double,-10);
    cio_info_read(kF,double,-14);
    cio_info_read(mu,double,-11);
    cio_info_read(ec,double,-12);
    cio_info_read(time,double,-13);
    
    // Close file
    close(fd);
    
    return werr;
}

/**
 * Reads info from checkpoint file
 * */
int read_checkpoint_info_pca(const char * file_name, 
                           int *nwf, int *nx, int *ny, int *nz, double *dx, double *dy, double *dz,
                           double *kF, double *mu, double *ec, double *beta
                          )
{
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDONLY  , fd_mode ) ) == -1 ) return WSLDA_ERR_CANNOT_OPEN_FILE;
    
    int werr=0;
    long int bytes_read ;
    int prec;
    
#define cio_info_read(var,type,werrval) \
    if (werr==0 && ( bytes_read = read( fd , ( void * ) var , sizeof( type ) ) ) != sizeof( type ) ) werr=werrval;

    cio_info_read(&prec,int,-2);
    if(werr==0 && prec!=PCA_PRECISION) werr=-3;
    cio_info_read(nwf,int,-4);
    cio_info_read(nx,int,-5);
    cio_info_read(ny,int,-6);
    cio_info_read(nz,int,-7);
    cio_info_read(dx,double,-8);
    cio_info_read(dy,double,-9);
    cio_info_read(dz,double,-10);
    cio_info_read(kF,double,-14);
    cio_info_read(&mu[SPINA],double,-11);
    cio_info_read(&mu[SPINB],double,-11);
    cio_info_read(ec,double,-12);
    cio_info_read(beta,double,-13);
    
    // Close file
    close(fd);
    
    return werr;
}

int touch_file(const char * file_name)
{
    if(exists(file_name)) 
    {
        if(md.overwrite==0) return WSLDA_ERR_CANNOT_OVERWRITE; // We do not overwrite!  
        else urm(file_name);
    }
    
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_CREAT | O_RDWR  , fd_mode ) ) == -1 ) return -2; 
    
    // Close file
    close( fd );
        
    return 0;
}

int checkpoint_save_u_and_delta(const char * file_name, int nxyz, double *u, double complex *delta)
{ 
    if(exists(file_name)) 
    {
        if(md.overwrite==0) return -1; // We do not overwrite!  
        else urm(file_name);
    }
    
    int fd ; /* file descriptor for handling the file or device */
    int werr=0;
    long int bytes_rw ;
    
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_CREAT | O_RDWR  , fd_mode ) ) == -1 ) return -2; 
    
    if (werr==0 && ( bytes_rw = write( fd , ( const void * ) u , nxyz*sizeof( double ) ) ) != nxyz*sizeof( double ) ) werr=-3;
    if (werr==0 && ( bytes_rw = write( fd , ( const void * ) delta , nxyz*sizeof( double complex ) ) ) != nxyz*sizeof( double complex ) ) werr=-4;

    // Close file
    close( fd );
        
    return werr;
}

int checkpoint_read_u_and_delta(const char * file_name, int nxyz, double *u, double complex *delta)
{     
    int fd ; /* file descriptor for handling the file or device */
    int werr=0;
    long int bytes_rw ;
    
    mode_t fd_mode = S_IRUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDONLY  , fd_mode ) ) == -1 ) return -2; 
    
    if (werr==0 && ( bytes_rw = read( fd , ( void * ) u , nxyz*sizeof( double ) ) ) != nxyz*sizeof( double ) ) werr=-3;
    if (werr==0 && ( bytes_rw = read( fd , ( void * ) delta , nxyz*sizeof( double complex ) ) ) != nxyz*sizeof( double complex ) ) werr=-4;

    // Close file
    close( fd );
        
    return werr;
}

int read_binary_file(const char * file_name, unsigned long int size, unsigned long int bshift, void * data)
{     
    
    FILE *pFile;
    
    pFile= fopen (file_name, "r");
    if (pFile==NULL)  return -1; // cannot open
    
    if(fseek ( pFile, bshift, SEEK_SET ) != 0 ) return -2; // cannot seek pointer
    
    size_t test_ele = fread (data , size, 1, pFile);
    if(test_ele!=1) return -3; // data read
    
    fclose(pFile);
        
    return 0;
}

int append_to_binary_file(const char * file_name, size_t size, void * data)
{     
    
    FILE *pFile;
    
    pFile= fopen (file_name, "ab");
    if (pFile==NULL)  return -1; // cannot open    
        
    size_t test_ele = fwrite (data , size, 1, pFile);
    if(test_ele!=1) return -3; // data not written 
    
    fclose(pFile);
    
    return 0;
}

/**
 * Function add entries to the file with wave-functions and kkz
 * */
double fbeta(double E, double beta);
#if CODEDIM==1
#define _BSHIFT NX 
#else
#define _BSHIFT NX*NY
#endif
/**
 * Function writes data from kzSLpca solver
 * Saves only kkz and en and counts states
 * */
int append_wf_from_kzpcaSL_part1(char * prefix, double *En, double complex *psi, double ecut, double beta, double ky, double kz, int ikz, int nwftwrt, int *nwf)
{
    char file_name_u[512];
    char file_name_v[512];
    char file_name_kkz[512];
    char file_name_fbeta[512];
    
#if CODEDIM==1
    sprintf(file_name_u, "%s/s1dpca.%04d.wfu", prefix, ikz);
    sprintf(file_name_v, "%s/s1dpca.%04d.wfv", prefix, ikz);
    sprintf(file_name_kkz, "%s/s1dpca.%04d.kkyz", prefix, ikz);
    sprintf(file_name_fbeta, "%s/s1dpca.%04d.en", prefix, ikz);
#else
    sprintf(file_name_u, "%s/s2dpca.%04d.wfu", prefix, ikz);
    sprintf(file_name_v, "%s/s2dpca.%04d.wfv", prefix, ikz);
    sprintf(file_name_kkz, "%s/s2dpca.%04d.kkz", prefix, ikz);
    sprintf(file_name_fbeta, "%s/s2dpca.%04d.en", prefix, ikz);
#endif
    
    // open files
    FILE *fkkz = fopen(file_name_kkz, "ab");
    FILE *ffbeta = fopen(file_name_fbeta, "ab");
    
    if (fkkz==NULL)  return -3; // cannot open  
    if (ffbeta==NULL)  return -4; // cannot open  
    
    int ien;
    size_t test_ele;
    double fbEn;
    
    for(ien=0; ien<nwftwrt; ien++) // for each eigen-energy 
    {
        if(fabs(En[ien])>ecut) continue; // above cut-off - skip!!!
        if(md.spinsymmetry==1 && En[ien]<0.0) continue; // spin symmetric mode - take only positive states
        
            
        (*nwf)++; // we have new state
#if CODEDIM==1
        test_ele = fwrite(&ky, sizeof(double), 1, fkkz);
        if(test_ele!=1) return -7; // data not written 
#endif
        test_ele = fwrite(&kz, sizeof(double), 1, fkkz);
        if(test_ele!=1) return -7; // data not written 
        
//         fbEn=fbeta(En[ien], beta);
        fbEn=En[ien];
        test_ele = fwrite(&fbEn, sizeof(double), 1, ffbeta);
        if(test_ele!=1) return -8; // data not written
    }
    
    // close files
    fclose(fkkz);
    fclose(ffbeta);
    
    return 0;
}

/**
 * Function writes data from kzSLpca solver
 * Saves only u components
 * */
int append_wf_from_kzpcaSL_part2(char * prefix, double *En, double complex *psi, double ecut, double beta, double ky, double kz, int ikz, int nwftwrt, int *nwf)
{
    char file_name_u[512];
    char file_name_v[512];
    char file_name_kkz[512];
    char file_name_fbeta[512];
    
#if CODEDIM==1
    sprintf(file_name_u, "%s/s1dpca.%04d.wfu", prefix, ikz);
    sprintf(file_name_v, "%s/s1dpca.%04d.wfv", prefix, ikz);
    sprintf(file_name_kkz, "%s/s1dpca.%04d.kkyz", prefix, ikz);
    sprintf(file_name_fbeta, "%s/s1dpca.%04d.en", prefix, ikz);
#else
    sprintf(file_name_u, "%s/s2dpca.%04d.wfu", prefix, ikz);
    sprintf(file_name_v, "%s/s2dpca.%04d.wfv", prefix, ikz);
    sprintf(file_name_kkz, "%s/s2dpca.%04d.kkz", prefix, ikz);
    sprintf(file_name_fbeta, "%s/s2dpca.%04d.en", prefix, ikz);
#endif
    // open files
    FILE *fu = fopen(file_name_u, "ab");
    
    if (fu==NULL)  return -1; // cannot open  
    
    int ien;
    double complex *u, *v; // psi=(u,v)
    size_t test_ele;
    
    for(ien=0; ien<nwftwrt; ien++) // for each eigen-energy 
    {
        
        if(fabs(En[ien])>ecut) continue; // above cut-off - skip!!!
        if(md.spinsymmetry==1 && En[ien]<0.0) continue; // spin symmetric mode - take only positive states
        
        // docompose state
        u = psi + ien*2*_BSHIFT; 
        
        test_ele = fwrite(u, sizeof(double complex)*_BSHIFT, 1, fu);
        if(test_ele!=1) return -5; // data not written 
        
    }
    
    // close files
    fclose(fu);
    
    return 0;
}

/**
 * Function writes data from kzSLpca solver
 * Saves only v components
 * */
int append_wf_from_kzpcaSL_part3(char * prefix, double *En, double complex *psi, double ecut, double beta, double ky, double kz, int ikz, int nwftwrt, int *nwf)
{
    char file_name_u[512];
    char file_name_v[512];
    char file_name_kkz[512];
    char file_name_fbeta[512];
    
#if CODEDIM==1
    sprintf(file_name_u, "%s/s1dpca.%04d.wfu", prefix, ikz);
    sprintf(file_name_v, "%s/s1dpca.%04d.wfv", prefix, ikz);
    sprintf(file_name_kkz, "%s/s1dpca.%04d.kkyz", prefix, ikz);
    sprintf(file_name_fbeta, "%s/s1dpca.%04d.en", prefix, ikz);
#else
    sprintf(file_name_u, "%s/s2dpca.%04d.wfu", prefix, ikz);
    sprintf(file_name_v, "%s/s2dpca.%04d.wfv", prefix, ikz);
    sprintf(file_name_kkz, "%s/s2dpca.%04d.kkz", prefix, ikz);
    sprintf(file_name_fbeta, "%s/s2dpca.%04d.en", prefix, ikz);
#endif
    
    // open files
    FILE *fv = fopen(file_name_v, "ab");
    
    if (fv==NULL)  return -2; // cannot open   
    
    int ien;
    double complex *u, *v; // psi=(u,v)
    size_t test_ele;
    
    for(ien=0; ien<nwftwrt; ien++) // for each eigen-energy 
    {
        if(fabs(En[ien])>ecut) continue; // above cut-off - skip!!!
        if(md.spinsymmetry==1 && En[ien]<0.0) continue; // spin symmetric mode - take only positive states
        
        
        // docompose state
        u = psi + ien*2*_BSHIFT;
        v = u + _BSHIFT;
        
        
        test_ele = fwrite(v, sizeof(double complex)*_BSHIFT, 1, fv);
        if(test_ele!=1) return -6; // data not written
        
    }
    
    // close files
    fclose(fv);
    
    return 0;
}

int checkpoint_save_u_and_delta_kzpca(const char * file_name, int nxyz, double *u, double complex *delta)
{ 
    if(exists(file_name)) 
    {
        if(md.overwrite==0) return -1; // We do not overwrite!  
        else urm(file_name);
    }
    
    int fd ; /* file descriptor for handling the file or device */
    int werr=0;
    long int bytes_rw ;
    
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_CREAT | O_RDWR  , fd_mode ) ) == -1 ) return -2; 
    
    if (werr==0 && ( bytes_rw = write( fd , ( const void * ) u , nxyz*sizeof( double )*2 ) ) != nxyz*sizeof( double )*2 ) werr=-3;
    if (werr==0 && ( bytes_rw = write( fd , ( const void * ) delta , nxyz*sizeof( double complex ) ) ) != nxyz*sizeof( double complex ) ) werr=-4;

    // Close file
    close( fd );
        
    return werr;
}

/**
 * Function determines number of wave-functions in each separate kz-file,
 * @param prefix for construction file names (INPUT)
 * @param nz size of kz points
 * @param nwf total number of wf, for checking correctness of input set (INPUT/OUTPUT)
 * @param nwf_per_kz number of wf for each kz, array of size NZ/2 (OUTPUT)
 * */
int scan_kzpca_info_files(const char * prefix, int nz, int *nwf, int *nwf_per_kz)
{
    char file_name[256];
    int ikz;
    FILE * pFile;
    int i, tnwf=0;
    
    int ikzadd=0;
#ifdef USE_CUBIC_CUTOFF
    ikzadd=1;
#endif
    
    for(ikz=0; ikz<nz/2+ikzadd; ikz++)
    {
        sprintf(file_name, "%s/s2dpca.%04d.info", prefix, ikz);
        
        pFile = fopen(file_name, "rb");
        if(pFile==NULL) return 1000+ikz;
        fread(&i          , sizeof(int)         , 1 , pFile); // percision
        fread(&i          , sizeof(int)         , 1 , pFile); // nwf
        fclose(pFile);
        
        nwf_per_kz[ikz]=i;
        if(ikz==0) tnwf+=i;
        else if(ikz==nz/2) tnwf+=i;
        else       tnwf+=2*i;
        
//         wprintf("# scan_kzpca_info_files: %4d %4d %4d\n", ikz, i, tnwf);
    }
    
    if(tnwf!=*nwf) return -1;
    
    tnwf=0;
    for(ikz=0; ikz<nz/2+ikzadd; ikz++)  tnwf+=nwf_per_kz[ikz];
    *nwf=tnwf;
    
    return 0;
}


/**
 * Function determines number of wave-functions in each kmode,
 * @param prefix for construction file names (INPUT)
 * @param codedim codes dimensonality that calls it
 * @param kvecs_to_consder number of k-modes
 * @param kvecs k-modes
 * @param nwf total number of wf, for checking correctness of input set (INPUT/OUTPUT)
 * @param nwf_per_kxy number of wf for each k-mode, of size kvecs_to_consder (OUTPUT)
 * */
int scan_stwslda1d_info_files(const char * prefix, int codedim, int kvecs_to_consder, wslda_kmode *kvecs, int *nwf, int *nwf_per_kyz)
{
    char file_name[256];
    int ikz;
    FILE * pFile;
    int i, tnwf=0;
    
    for(ikz=0; ikz<kvecs_to_consder; ikz++)
    {
        sprintf(file_name, "%s/s1dpca.%04d.info", prefix, ikz);
        
        pFile = fopen(file_name, "rb");
        if(pFile==NULL) return WSLDA_ERR_CANNOT_OPEN_FILE;
        fread(&i          , sizeof(int)         , 1 , pFile); // percision
        fread(&i          , sizeof(int)         , 1 , pFile); // nwf
        fclose(pFile);
        
        nwf_per_kyz[ikz]=i;
        tnwf+=i*kvecs[ikz].weight;        
//         wprintf("# scan_kzpca_info_files: %4d %4d %4d\n", ikz, i, tnwf);
    }
    
    if(tnwf!=*nwf) return WSLDA_ERR_BINARY_FILE_CORRUPTED;
    
    double *kytmp, *kztmp;
    cppmallocl(kytmp,NY*NZ,double);
    cppmallocl(kztmp,   NZ,double);
    
    tnwf=0;
    if(codedim==1)
    {
        for(ikz=0; ikz<kvecs_to_consder; ikz++)  tnwf+=nwf_per_kyz[ikz];
    }
    else if(codedim==2)
    {
        for(ikz=0; ikz<kvecs_to_consder; ikz++) 
        {
            i=wslda_kmodes_1d_getcnt2d(kvecs[ikz].ky,kvecs[ikz].kz, kytmp, kztmp);
            tnwf+=nwf_per_kyz[ikz]*i;
            
//             int ii;
//             for(ii=0; ii<i; ii++)
//                 wprintf("TTT: %6d %12.8f %12.8f %6d %6d %12.8f %12.8f\n", ikz, kvecs[ikz].ky,kvecs[ikz].kz, i, ii, kytmp[ii], kztmp[ii]); // TODO
        }
    }
    else if(codedim==3)
    {
        for(ikz=0; ikz<kvecs_to_consder; ikz++)  tnwf+=nwf_per_kyz[ikz]*kvecs[ikz].weight;
    }
    *nwf=tnwf;
    
    free(kytmp);
    free(kztmp);
    
    return WSLDA_OK;
}

/**
 * Function reads wf from kzSLpca standard
 * */
int read_kzSLpca_wf(const char * prefix, int nz, int *nwf_per_kz, int mylidx, int myuidx, 
                    double complex *h_wavefun, double *h_fbetaEn, double *h_kkz)
{
    char file_name[256];
    int ikz, iwf=0, ii;
    int nwfip = myuidx-mylidx;

    char file_name_u[512];
    char file_name_v[512];
    char file_name_kkz[512];
    char file_name_fbeta[512];
    
    FILE *fu;
    FILE *fv;
    FILE *fkkz;
    FILE *ffbeta;
    
    int ikzadd=0;
#ifdef USE_CUBIC_CUTOFF
    ikzadd=1;
#endif
    
//     wprintf("mylidx=%d, myuidx=%d\n", mylidx, myuidx);
    
    for(ikz=0; ikz<nz/2+ikzadd; ikz++)
    {
        // reset pointer to file
        fu=NULL;
        
        for(ii=0; ii<nwf_per_kz[ikz]; ii++)
        {
            if(iwf>=mylidx && iwf<myuidx)
            {
//                 wprintf("loading iwf=%d %d\n", iwf, ikz);
                
                if(fu==NULL) // open files
                {
//                     wprintf("OPENING iwf=%d, file=%d\n", iwf, ikz);
                    sprintf(file_name_u, "%s/s2dpca.%04d.wfu", prefix, ikz);
                    sprintf(file_name_v, "%s/s2dpca.%04d.wfv", prefix, ikz);
                    sprintf(file_name_kkz, "%s/s2dpca.%04d.kkz", prefix, ikz);
                    sprintf(file_name_fbeta, "%s/s2dpca.%04d.en", prefix, ikz);
                    
                    fu = fopen(file_name_u, "rb");
                    fv = fopen(file_name_v, "rb");
                    fkkz = fopen(file_name_kkz, "rb");
                    ffbeta = fopen(file_name_fbeta, "rb");
                    
                    if (fu==NULL)  return -1; // cannot open  
                    if (fv==NULL)  return -2; // cannot open  
                    if (fkkz==NULL)  return -3; // cannot open  
                    if (ffbeta==NULL)  return -4; // cannot open        
                    
                    // shift pointer to correct position
                    if(fseek ( fu, sizeof(double complex)*NXY*ii, SEEK_SET ) != 0 ) return -11; // cannot seek pointer
                    if(fseek ( fv, sizeof(double complex)*NXY*ii, SEEK_SET ) != 0 ) return -12; // cannot seek pointer
                    if(fseek ( fkkz, sizeof(double)*ii, SEEK_SET ) != 0 ) return -13; // cannot seek pointer
                    if(fseek ( ffbeta, sizeof(double)*ii, SEEK_SET ) != 0 ) return -14; // cannot seek pointer
                }
                
                if( fread(h_wavefun + (size_t)NXY*(iwf-mylidx)                    , sizeof(double complex)*NXY, 1 , fu) != 1) return -21;
                if( fread(h_wavefun + (size_t)NXY*(iwf-mylidx) + (size_t)NXY*nwfip, sizeof(double complex)*NXY, 1 , fv) != 1) return -22;
                if( fread(h_fbetaEn + (iwf-mylidx) , sizeof(double), 1 , ffbeta) != 1) return -23;
                if( fread(h_kkz     + (iwf-mylidx) , sizeof(double), 1 , fkkz)   != 1) return -23;
            }
            
            iwf++;
        }
        
        // close files if opened
        if(fu!=NULL)
        {
            fclose(fu);
            fclose(fv);
            fclose(fkkz);
            fclose(ffbeta);
//             wprintf("CLOSING iwf=%d, file=%d\n", iwf, ikz);
        }
    }
    
    return 0;
}

/**
 * Function reads wf from kzSLpca standard
 * */
int read_kzSLpca_wf_with_doubling(const char * prefix, int nz, int *nwf_per_kz, int mylidx, int myuidx, 
                    double complex *h_wavefun, double *h_fbetaEn, double *h_kkz)
{
    char file_name[256];
    int ikz, iwf=0, ii, dd, dcoeff;
    int nwfip = myuidx-mylidx;

    char file_name_u[512];
    char file_name_v[512];
    char file_name_kkz[512];
    char file_name_fbeta[512];
    
    FILE *fu;
    FILE *fv;
    FILE *fkkz;
    FILE *ffbeta;
    
    int ikzadd=0;
#ifdef USE_CUBIC_CUTOFF
    ikzadd=1;
#endif
    
//     wprintf("mylidx=%d, myuidx=%d\n", mylidx, myuidx);
    
    for(ikz=0; ikz<nz/2+ikzadd; ikz++)
    {
        // reset pointer to file
        fu=NULL;
        
        if(ikz==0) dcoeff=1;
        else if(ikz==nz/2) dcoeff=1;
        else dcoeff=2;
        
        for(dd=0; dd<dcoeff; dd++)
        {
            if(fu!=NULL && dd==1) // reset pointer to the beginning
            {
                // shift pointer to correct position
                if(fseek ( fu, 0, SEEK_SET ) != 0 ) return -111; // cannot seek pointer
                if(fseek ( fv, 0, SEEK_SET ) != 0 ) return -112; // cannot seek pointer
                if(fseek ( fkkz, 0, SEEK_SET ) != 0 ) return -113; // cannot seek pointer
                if(fseek ( ffbeta, 0, SEEK_SET ) != 0 ) return -114; // cannot seek pointer
            } 
            
            for(ii=0; ii<nwf_per_kz[ikz]; ii++)
            {
                if(iwf>=mylidx && iwf<myuidx)
                {
    //                 wprintf("loading iwf=%d %d\n", iwf, ikz);
                    
                    if(fu==NULL) // open files
                    {
    //                     wprintf("OPENING iwf=%d, file=%d\n", iwf, ikz);
                        sprintf(file_name_u, "%s/s2dpca.%04d.wfu", prefix, ikz);
                        sprintf(file_name_v, "%s/s2dpca.%04d.wfv", prefix, ikz);
                        sprintf(file_name_kkz, "%s/s2dpca.%04d.kkz", prefix, ikz);
                        sprintf(file_name_fbeta, "%s/s2dpca.%04d.en", prefix, ikz);
                        
                        fu = fopen(file_name_u, "rb");
                        fv = fopen(file_name_v, "rb");
                        fkkz = fopen(file_name_kkz, "rb");
                        ffbeta = fopen(file_name_fbeta, "rb");
                        
                        if (fu==NULL)  return -1; // cannot open  
                        if (fv==NULL)  return -2; // cannot open  
                        if (fkkz==NULL)  return -3; // cannot open  
                        if (ffbeta==NULL)  return -4; // cannot open        
                        
                        // shift pointer to correct position
                        if(fseek ( fu, sizeof(double complex)*NXY*ii, SEEK_SET ) != 0 ) return -11; // cannot seek pointer
                        if(fseek ( fv, sizeof(double complex)*NXY*ii, SEEK_SET ) != 0 ) return -12; // cannot seek pointer
                        if(fseek ( fkkz, sizeof(double)*ii, SEEK_SET ) != 0 ) return -13; // cannot seek pointer
                        if(fseek ( ffbeta, sizeof(double)*ii, SEEK_SET ) != 0 ) return -14; // cannot seek pointer
                    }
                    
                    if( fread(h_wavefun + (size_t)NXY*(iwf-mylidx)                    , sizeof(double complex)*NXY, 1 , fu) != 1) return -21;
                    if( fread(h_wavefun + (size_t)NXY*(iwf-mylidx) + (size_t)NXY*nwfip, sizeof(double complex)*NXY, 1 , fv) != 1) return -22;
                    if( fread(h_fbetaEn + (iwf-mylidx) , sizeof(double), 1 , ffbeta) != 1) return -23;
                    if( fread(h_kkz     + (iwf-mylidx) , sizeof(double), 1 , fkkz)   != 1) return -23;
                    
                    
                    if(dd==1) // revert sign of kz vector
                        h_kkz[iwf-mylidx]*=-1.0; 
                }
                
                iwf++;
            }
        }
        
        // close files if opened
        if(fu!=NULL)
        {
            fclose(fu);
            fclose(fv);
            fclose(fkkz);
            fclose(ffbeta);
//             wprintf("CLOSING iwf=%d, file=%d\n", iwf, ikz);
        }
    }
    
    return 0;
}

/**
 * Function reads wf from kzSLpca standard
 * */
int read_stwslda1d_wf(const char * prefix, int codedim, int kvecs_to_consder, wslda_kmode *kvecs, int *nwf_per_kyz, int mylidx, int myuidx, 
                    double complex *h_wavefun, double *h_fbetaEn, double *h_kkyz, int *h_cnt)
{
    char file_name[256];
    int ikz, iwf=0, ii;
    int nwfip = myuidx-mylidx;

    char file_name_u[512];
    char file_name_v[512];
    char file_name_kkz[512];
    char file_name_fbeta[512];
    
    FILE *fu;
    FILE *fv;
    FILE *fkkz;
    FILE *ffbeta;
    
    double complex *_u, *_v;
    double _ky, _kz, _en; 
    cppmallocl(_u,NX,double complex);
    cppmallocl(_v,NX,double complex);
    
    double *kytmp, *kztmp;
    cppmallocl(kytmp,NY*NZ,double);
    cppmallocl(kztmp,NY*NZ,double);
    
    int lNdim, dcoeff, dd;
    int ix, iy, iz, ixyz;
    
//     wprintf("mylidx=%d, myuidx=%d\n", mylidx, myuidx);
    
    for(ikz=0; ikz<kvecs_to_consder; ikz++)
    {
        // reset pointer to file
        fu=NULL;
        
        if(codedim==1) 
        {
            dcoeff=1;
            lNdim=NX;
        }
        else if(codedim==2)
        {
            dcoeff=wslda_kmodes_1d_getcnt2d(kvecs[ikz].ky, kvecs[ikz].kz, kytmp, kztmp);
            lNdim=NX*NY;
        }
        else if(codedim==3)
        {
            dcoeff=1;
            wslda_kmodes_1d_get_modes(kvecs[ikz].ky, kvecs[ikz].kz, &dcoeff, kytmp, kztmp);
            lNdim=NX*NY*NZ;
//             wprintf("CCC: %6d %6d %6d\n", ikz, dcoeff, kvecs[ikz].weight);
            if(dcoeff!=kvecs[ikz].weight) return -122;
        }

        
        for(dd=0; dd<dcoeff; dd++)
        {
            if(fu!=NULL && dd>=1) // reset pointer to the beginning
            {
                // shift pointer to correct position
                if(fseek ( fu, 0, SEEK_SET ) != 0 ) return -111; // cannot seek pointer
                if(fseek ( fv, 0, SEEK_SET ) != 0 ) return -112; // cannot seek pointer
                if(fseek ( fkkz, 0, SEEK_SET ) != 0 ) return -113; // cannot seek pointer
                if(fseek ( ffbeta, 0, SEEK_SET ) != 0 ) return -114; // cannot seek pointer
            } 

            
            for(ii=0; ii<nwf_per_kyz[ikz]; ii++)
            {
                if(iwf>=mylidx && iwf<myuidx)
                {
//                     wprintf("loading iwf=%d %d\n", iwf, ikz);
                    
                    if(fu==NULL) // open files
                    {
//                         wprintf("OPENING iwf=%d, file=%d\n", iwf, ikz);
                        sprintf(file_name_u, "%s/s1dpca.%04d.wfu", prefix, ikz);
                        sprintf(file_name_v, "%s/s1dpca.%04d.wfv", prefix, ikz);
                        sprintf(file_name_kkz, "%s/s1dpca.%04d.kkyz", prefix, ikz);
                        sprintf(file_name_fbeta, "%s/s1dpca.%04d.en", prefix, ikz);
                        
                        fu = fopen(file_name_u, "rb");
                        fv = fopen(file_name_v, "rb");
                        fkkz = fopen(file_name_kkz, "rb");
                        ffbeta = fopen(file_name_fbeta, "rb");
                        
                        if (fu==NULL)  return -1; // cannot open  
                        if (fv==NULL)  return -2; // cannot open  
                        if (fkkz==NULL)  return -3; // cannot open  
                        if (ffbeta==NULL)  return -4; // cannot open        
                        
                        // shift pointer to correct position
                        if(fseek ( fu, sizeof(double complex)*NX*ii, SEEK_SET ) != 0 ) return -11; // cannot seek pointer
                        if(fseek ( fv, sizeof(double complex)*NX*ii, SEEK_SET ) != 0 ) return -12; // cannot seek pointer
                        if(fseek ( fkkz, sizeof(double)*ii*2, SEEK_SET ) != 0 ) return -13; // cannot seek pointer
                        if(fseek ( ffbeta, sizeof(double)*ii, SEEK_SET ) != 0 ) return -14; // cannot seek pointer
                    }
                    
                    if( fread(_u, sizeof(double complex)*NX, 1 , fu) != 1) return -21;
                    if( fread(_v, sizeof(double complex)*NX, 1 , fv) != 1) return -22;
                    if( fread(&_en , sizeof(double), 1 , ffbeta) != 1) return -23;
                    if( fread(&_ky, sizeof(double), 1 , fkkz)   != 1) return -23;
                    if( fread(&_kz, sizeof(double), 1 , fkkz)   != 1) return -24;
                    
                    if(codedim==1)
                    {
                        ixyz=0;
                        for(ix=0; ix<NX; ix++) 
                        {
                            h_wavefun[(size_t)NX*(iwf-mylidx)                    + ixyz] = _u[ix];
                            h_wavefun[(size_t)NX*(iwf-mylidx) + (size_t)NX*nwfip + ixyz] = _v[ix];
                            ixyz++;
                        }
                        h_fbetaEn[iwf-mylidx]=_en;
                        h_kkyz[iwf-mylidx      ]=_ky;
                        h_kkyz[iwf-mylidx+nwfip]=_kz;
                        h_cnt [iwf-mylidx      ]= wslda_kmodes_1d_get_weight(_ky, _kz);
                    }
                    else if(codedim==2)
                    {
                        if(_ky!=kvecs[ikz].ky) return -44;
                        if(_kz!=kvecs[ikz].kz) return -45;
                        _ky=kytmp[dd];
                        ixyz=0;
                        for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) 
                        {
                            h_wavefun[(size_t)NXY*(iwf-mylidx)                     + ixyz] = _u[ix]*cexp(I*_ky*iy*DY)/sqrt(LY);
                            h_wavefun[(size_t)NXY*(iwf-mylidx) + (size_t)NXY*nwfip + ixyz] = _v[ix]*cexp(I*_ky*iy*DY)/sqrt(LY);
                            ixyz++;
                        }
                        h_fbetaEn[iwf-mylidx]=_en;
                        h_kkyz[iwf-mylidx]=kztmp[dd];
                    }
                    else if(codedim==3)
                    {
                        if(_ky!=kvecs[ikz].ky) return -54;
                        if(_kz!=kvecs[ikz].kz) return -55; 
                        // revert sign of ky vector
                        _ky=kytmp[dd];
                        _kz=kztmp[dd];
                        ixyz=0;
                        for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++) 
                        {
                            h_wavefun[(size_t)NXYZ*(iwf-mylidx)                      + ixyz] = _u[ix]*cexp(I*_ky*iy*DY)*cexp(I*_kz*iz*DZ)/sqrt(LY*LZ);
                            h_wavefun[(size_t)NXYZ*(iwf-mylidx) + (size_t)NXYZ*nwfip + ixyz] = _v[ix]*cexp(I*_ky*iy*DY)*cexp(I*_kz*iz*DZ)/sqrt(LY*LZ);
                            ixyz++;
                        }
                        h_fbetaEn[iwf-mylidx]=_en;
                    }
                    
                }
                
                iwf++;
            }
        }
        
        // close files if opened
        if(fu!=NULL)
        {
            fclose(fu);
            fclose(fv);
            fclose(fkkz);
            fclose(ffbeta);
//             wprintf("CLOSING iwf=%d, file=%d\n", iwf, ikz);
        }
    }
    
    free(_u); free(_v);
    free(kytmp); free(kztmp);
    
    return 0;
}
#undef _BSHIFT

// ---------------------------------------- s3dpca IO -------------------------------------------
/**
 * Function writes data from s3dpca solver
 * Saves only en and counts states
 * */
int append_wf_from_s3dpca_part1(char * prefix, double *En, double complex *psi, double ecut, double beta, int idgroup, int nwftwrt, int *nwf)
{
    char file_name_u[512];
    char file_name_v[512];
    char file_name_fbeta[512];
    
    sprintf(file_name_u, "%s/s3dpca.%04d.wfu", prefix, idgroup);
    sprintf(file_name_v, "%s/s3dpca.%04d.wfv", prefix, idgroup);
    sprintf(file_name_fbeta, "%s/s3dpca.%04d.en", prefix, idgroup);
    
    // open files
    FILE *ffbeta = fopen(file_name_fbeta, "ab");
    
    if (ffbeta==NULL)  return -4; // cannot open  
    
    int ien;
    size_t test_ele;
    double fbEn;
    
    for(ien=0; ien<nwftwrt; ien++) // for each eigen-energy 
    {
         
        if(fabs(En[ien])>ecut) continue; // above cut-off - skip!!!
        if(md.spinsymmetry==1 && En[ien]<0.0) continue; // spin symmetric mode - take only positive states
         
        (*nwf)++; // we have new state
        
//         fbEn=fbeta(En[ien], beta);
        fbEn=En[ien];
        test_ele = fwrite(&fbEn, sizeof(double), 1, ffbeta);
        if(test_ele!=1) return -8; // data not written
    }
    
    // close files
    fclose(ffbeta);
    
    return 0;
}

/**
 * Function writes data from s3dpca solver
 * Saves only u components
 * */
int append_wf_from_s3dpca_part2(char * prefix, double *En, double complex *psi, double ecut, double beta, int idgroup, int nwftwrt, int *nwf)
{
    char file_name_u[512];
    char file_name_v[512];
    char file_name_fbeta[512];
    
    sprintf(file_name_u, "%s/s3dpca.%04d.wfu", prefix, idgroup);
    sprintf(file_name_v, "%s/s3dpca.%04d.wfv", prefix, idgroup);
    sprintf(file_name_fbeta, "%s/s3dpca.%04d.en", prefix, idgroup);
    
    // open files
    FILE *fu = fopen(file_name_u, "ab");
    
    if (fu==NULL)  return -1; // cannot open  
    
    int ien;
    double complex *u, *v; // psi=(u,v)
    size_t test_ele;
    
    for(ien=0; ien<nwftwrt; ien++) // for each eigen-energy 
    {
        
        if(fabs(En[ien])>ecut) continue; // above cut-off - skip!!!
        if(md.spinsymmetry==1 && En[ien]<0.0) continue; // spin symmetric mode - take only positive states
        
        // docompose state
        u = psi + ien*2*NXYZ; 
        
        test_ele = fwrite(u, sizeof(double complex)*NXYZ, 1, fu);
        if(test_ele!=1) return -5; // data not written 
        
    }
    
    // close files
    fclose(fu);
    
    return 0;
}

/**
 * Function writes data from s3dpca solver
 * Saves only v components
 * */
int append_wf_from_s3dpca_part3(char * prefix, double *En, double complex *psi, double ecut, double beta, int idgroup, int nwftwrt, int *nwf)
{
    char file_name_u[512];
    char file_name_v[512];
    char file_name_fbeta[512];
    
    sprintf(file_name_u, "%s/s3dpca.%04d.wfu", prefix, idgroup);
    sprintf(file_name_v, "%s/s3dpca.%04d.wfv", prefix, idgroup);
    sprintf(file_name_fbeta, "%s/s3dpca.%04d.en", prefix, idgroup);
    
    // open files
    FILE *fv = fopen(file_name_v, "ab");
    
    if (fv==NULL)  return -2; // cannot open   
    
    int ien;
    double complex *u, *v; // psi=(u,v)
    size_t test_ele;
    
    for(ien=0; ien<nwftwrt; ien++) // for each eigen-energy 
    {
        if(fabs(En[ien])>ecut) continue; // above cut-off - skip!!!
        if(md.spinsymmetry==1 && En[ien]<0.0) continue; // spin symmetric mode - take only positive states
        
        
        // docompose state
        u = psi + ien*2*NXYZ;
        v = u + NXYZ;
        
        
        test_ele = fwrite(v, sizeof(double complex)*NXYZ, 1, fv);
        if(test_ele!=1) return -6; // data not written
        
    }
    
    // close files
    fclose(fv);
    
    return 0;
}

/**
 * Function determines number of wave-functions in each separate kz-file,
 * @param prefix for construction file names (INPUT)
 * @param number_of_files number of files to scan
 * @param nwf total number of wf, for checking correctness of input set (INPUT/OUTPUT)
 * @param nwf_per_file number of wf for each file, array of size number_of_files (OUTPUT)
 * */
int scan_s3dpca_info_files(const char * prefix, int number_of_files, int *nwf, int *nwf_per_file)
{
    char file_name[256];
    int ikz;
    FILE * pFile;
    int i, tnwf=0;
    
    for(ikz=0; ikz<number_of_files; ikz++)
    {
        sprintf(file_name, "%s/s3dpca.%04d.info", prefix, ikz);
        
        pFile = fopen(file_name, "rb");
        if(pFile==NULL) return WSLDA_ERR_S3DPCA_INFO_FILES_MISSING_FILE;
        fread(&i          , sizeof(int)         , 1 , pFile); // percision
        fread(&i          , sizeof(int)         , 1 , pFile); // nwf
        fclose(pFile);
        
        nwf_per_file[ikz]=i;
        tnwf+=i;
        
    }
    

    if(tnwf!=*nwf) return WSLDA_ERR_S3DPCA_INFO_FILES; // files are not consistent with excepted nwf
    
    return 0;
}

/**
 * Function reads wf from kzSLpca standard
 * */
int read_s3dpca_wf(const char * prefix, int number_of_files, int *nwf_per_file, int mylidx, int myuidx, 
                    double complex *h_wavefun, double *h_fbetaEn)
{
    char file_name[256];
    int ikz, iwf=0, ii;
    int nwfip = myuidx-mylidx;

    char file_name_u[512];
    char file_name_v[512];
    char file_name_fbeta[512];
    
    FILE *fu;
    FILE *fv;
    FILE *ffbeta;
    
//     wprintf("mylidx=%d, myuidx=%d\n", mylidx, myuidx);
    
    for(ikz=0; ikz<number_of_files; ikz++)
    {
        // reset pointer to file
        fu=NULL;
        
        for(ii=0; ii<nwf_per_file[ikz]; ii++)
        {
            if(iwf>=mylidx && iwf<myuidx)
            {
//                 wprintf("loading iwf=%d %d\n", iwf, ikz);
                
                if(fu==NULL) // open files
                {
//                     wprintf("OPENING iwf=%d, file=%d\n", iwf, ikz);
                    sprintf(file_name_u, "%s/s3dpca.%04d.wfu", prefix, ikz);
                    sprintf(file_name_v, "%s/s3dpca.%04d.wfv", prefix, ikz);
                    sprintf(file_name_fbeta, "%s/s3dpca.%04d.en", prefix, ikz);
                    
                    fu = fopen(file_name_u, "rb");
                    fv = fopen(file_name_v, "rb");
                    ffbeta = fopen(file_name_fbeta, "rb");
                    
                    if (fu==NULL)  return -1; // cannot open  
                    if (fv==NULL)  return -2; // cannot open  
                    if (ffbeta==NULL)  return -4; // cannot open        
                    
                    // shift pointer to correct position
                    if(fseek ( fu, sizeof(double complex)*NXYZ*ii, SEEK_SET ) != 0 ) return -11; // cannot seek pointer
                    if(fseek ( fv, sizeof(double complex)*NXYZ*ii, SEEK_SET ) != 0 ) return -12; // cannot seek pointer
                    if(fseek ( ffbeta, sizeof(double)*ii, SEEK_SET ) != 0 ) return -14; // cannot seek pointer
                }
                
                if( fread(h_wavefun + (size_t)NXYZ*(iwf-mylidx)                     , sizeof(double complex)*NXYZ, 1 , fu) != 1) return -21;
                if( fread(h_wavefun + (size_t)NXYZ*(iwf-mylidx) + (size_t)NXYZ*nwfip, sizeof(double complex)*NXYZ, 1 , fv) != 1) return -22;
                if( fread(h_fbetaEn + (iwf-mylidx) , sizeof(double), 1 , ffbeta) != 1) return -23;
            }
            
            iwf++;
        }
        
        // close files if opened
        if(fu!=NULL)
        {
            fclose(fu);
            fclose(fv);
            fclose(ffbeta);
//             wprintf("CLOSING iwf=%d, file=%d\n", iwf, ikz);
        }
    }
    
    return 0;
}

/**
 * Function add entries to check.stamp file
 * @param prefix is used for creation file name of form prefix_check.stamp
 * @param indens number of densities stored in densities array
 * @param ndens number of elements for each density
 * @param densities array with densities, total size is indens*ndens
 * @param ineregies number of entries in array energies
 * @param energies energies of the system
 * @return 0 - ok, otherwies error
 **/
int check_stamp_entry(const char *file_name, int idens, int ndens, double *densities, int ienergies, double *energies)
{
    // write
    FILE *check_stamp = fopen(file_name, "a");
    if(check_stamp==NULL) return 1;
    
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [20];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,20,"%x-%X",timeinfo);
    fprintf(check_stamp, "CHECK STAMP DATE: %s\n", buffer);
    
    double sum;
    int i, n;
    for(i=0; i<idens; i++)
    {
        sum=0.0;
        for(n=0; n<ndens; n++) sum+=densities[i*ndens + n];
        fprintf(check_stamp, "SUM(DENSITY[%2d])=%16.8g\n", i, sum);
    }
    
    for(i=0; i<ienergies; i++)
    {
        fprintf(check_stamp, "ENERGY[%2d])=%16.8f\n", i, energies[i]);
    }
    
    fclose(check_stamp);
                
    return 0;
}

/**
 * Function add entries to check.stamp file
 * @param prefix is used for creation file name of form prefix_check.stamp
 * @param indens number of densities stored in densities array
 * @param ndens number of elements for each density
 * @param densities array with densities, total size is indens*ndens
 * @param ineregies number of entries in array energies
 * @param energies energies of the system
 * @param dens_coeff sum of densities will be mutiplied by this coeff before writing the stamp
 * @return 0 - ok, otherwies error
 **/
int check_stamp_entry_coeff(const char *file_name, int idens, int ndens, double *densities, int ienergies, double *energies, double dens_coeff)
{
    // write
    FILE *check_stamp = fopen(file_name, "a");
    if(check_stamp==NULL) return 1;
    
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [20];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (buffer,20,"%x-%X",timeinfo);
    fprintf(check_stamp, "CHECK STAMP DATE: %s\n", buffer);
    
    double sum;
    int i, n;
    for(i=0; i<idens; i++)
    {
        sum=0.0;
        for(n=0; n<ndens; n++) sum+=densities[i*ndens + n];
        fprintf(check_stamp, "SUM(DENSITY[%2d])=%16.8g\n", i, sum*dens_coeff);
    }
    
    for(i=0; i<ienergies; i++)
    {
        fprintf(check_stamp, "ENERGY[%2d])=%16.8f\n", i, energies[i]);
    }
    
    fclose(check_stamp);
                
    return 0;
}

/**
 * Checks if output file exists
 * if yes, returns error
 * */
int check_if_can_overwrite_files()
{
    if(md.overwrite==1) return WSLDA_OK;
    
    char fname[1024];
    sprintf(fname, "%s_input.txt", md.outprefix); if(exists(fname)) return WSLDA_ERR_CANNOT_OVERWRITE;
    sprintf(fname, "%s.wlog", md.outprefix); if(exists(fname)) return WSLDA_ERR_CANNOT_OVERWRITE;
    sprintf(fname, "%s.wtxt", md.outprefix); if(exists(fname)) return WSLDA_ERR_CANNOT_OVERWRITE;
    sprintf(fname, "%s.stdout", md.outprefix); if(exists(fname)) return WSLDA_ERR_CANNOT_OVERWRITE;
    sprintf(fname, "%s_predefines.h", md.outprefix); if(exists(fname)) return WSLDA_ERR_CANNOT_OVERWRITE;
    sprintf(fname, "%s_problem-definition.h", md.outprefix); if(exists(fname)) return WSLDA_ERR_CANNOT_OVERWRITE;
    sprintf(fname, "%s_logger.h", md.outprefix); if(exists(fname)) return WSLDA_ERR_CANNOT_OVERWRITE;
    
    return WSLDA_OK;
}

#endif
