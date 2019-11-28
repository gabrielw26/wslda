// Author: Gabriel Wlazlowski
// Date: 06-04-2013

// This file implements I/O for measurements

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>

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
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
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
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
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
        
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
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
        
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDWR  , fd_mode ) ) == -1 ) return -1;
  
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
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
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
 * Reads info from checkpoint file
 * */
int read_checkpoint_info(const char * file_name, 
                           int *nwf, int *nx, int *ny, int *nz, double *dx, double *dy, double *dz,
                           double *kF, double *mu, double *ec, double *time
                          )
{
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
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

int touch_file(const char * file_name)
{
    if(exists(file_name)) 
    {
        if(md.overwrite==0) return -1; // We do not overwrite!  
        else urm(file_name);
    }
    
    int fd ; /* file descriptor for handling the file or device */
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
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
    
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
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
    
    mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
    if ( ( fd = open( file_name , O_RDONLY  , fd_mode ) ) == -1 ) return -2; 
    
    if (werr==0 && ( bytes_rw = read( fd , ( void * ) u , nxyz*sizeof( double ) ) ) != nxyz*sizeof( double ) ) werr=-3;
    if (werr==0 && ( bytes_rw = read( fd , ( void * ) delta , nxyz*sizeof( double complex ) ) ) != nxyz*sizeof( double complex ) ) werr=-4;

    // Close file
    close( fd );
        
    return werr;
}

int read_binary_file(const char * file_name, unsigned long int size, unsigned long int bshift, void * data)
{     
// //     NOTE - read function has problem if file is bigger than 20GB !!!
//     int fd ; /* file descriptor for handling the file or device */
//     long int bytes_rw ;
//     
//     mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
//     if ( ( fd = open( file_name , O_RDONLY  , fd_mode ) ) == -1 ) return -1; 
//     if( lseek( fd, bshift, 0 )!= bshift ) return -2;
//     if (( bytes_rw = read( fd , ( void * ) data , size ) ) != size ) return -3;
// 
//     // Close file
//     close( fd );
    
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
    // NOTE - write function has problem if file is bigger than 20GB !!!
//     int fd ; /* file descriptor for handling the file or device */
//     long int bytes_rw ;
//     
//     mode_t fd_mode = S_IRUSR | S_IWUSR | S_IRGRP; /* S_IRWXU ; S_IRGRP, S_IRWXG ; S_IROTH , S_IRWXO; etc. */
//     if ( ( fd = open( file_name , O_APPEND | O_WRONLY , fd_mode ) ) == -1 ) return -1; 
//     if (( bytes_rw = write( fd , ( const void * ) data , size ) ) != size ) return -3;
// 
//     // Close file
//     close( fd );
    
    // use std writing routines
//     printf("Adding to file `%s` data of size %ldB...\n", file_name, size);
    
    FILE *pFile;
    
    pFile= fopen (file_name, "ab");
    if (pFile==NULL)  return -1; // cannot open    
        
    size_t test_ele = fwrite (data , size, 1, pFile);
    if(test_ele!=1) return -3; // data not written 
    
    fclose(pFile);
    
    return 0;
}

