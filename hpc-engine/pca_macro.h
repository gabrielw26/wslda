#include "wslda_errors.h"

#ifdef VERBOSE

#define EXEC_START_INFO(cmd) \
    printf("VERBOSE[%3d]: LANCHING    %s [file=`%s`, line=%d]\n", ip, #cmd, __FILE__, __LINE__); fflush(stdout); 
    
#define EXEC_STOP_INFO(cmd) \
    printf("VERBOSE[%3d]: DONE %6d=%s [file=`%s`, line=%d]\n", ip, ierr, #cmd, __FILE__, __LINE__); fflush(stdout);
    
#else
    
#define EXEC_START_INFO(cmd)
#define EXEC_STOP_INFO(cmd)

#endif

// allocation of memory, not involving MPI
#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "ERROR: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    } 

#define ABORT \
    {                                                                           \
        MPI_Barrier( MPI_COMM_WORLD );                                          \
        fprintf( stderr , "ABORT!!! [ip=%d, file=`%s`, line=%d]\n", ip, __FILE__, __LINE__) ; \
        ierr = -1 ;                                                             \
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;                                    \
        return( EXIT_FAILURE ) ;                                                \
    }
    
#define ABORTip(ip) \
    {                                                                           \
        MPI_Barrier( MPI_COMM_WORLD );                                          \
        fprintf( stderr , "ABORT!!! [ip=%d, file=`%s`, line=%d]\n", ip, __FILE__, __LINE__) ; \
        ierr = -1 ;                                                             \
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;                                    \
        return( EXIT_FAILURE ) ;                                                \
    }
    
#define ABORT_NOBARRIER \
    {                                                                           \
        fprintf( stderr , "ABORT!!! [ip=%d, file=`%s`, line=%d]\n", ip, __FILE__, __LINE__) ; \
        ierr = -1 ;                                                             \
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;                                    \
        return( EXIT_FAILURE ) ;                                                \
    }
    
#define TESTLINE                                                                \
    { wprintf("# TESTLINE: PROCESS %4d REACHED LINE %d IN FILE %s\n", ip, __LINE__ , __FILE__); fflush(stdout); }
    
    
#define TESTLINE_BARRIER                                                                \
    { wprintf("# TESTLINE: PROCESS %4d REACHED LINE %d IN FILE %s\n", ip, __LINE__ , __FILE__); fflush(stdout); MPI_Barrier( MPI_COMM_WORLD ); }

    
// execution of function by GPU and CPU side.
// Check error and terminate if fail
#define gpu_exec( cmd )                                                         \
    { EXEC_START_INFO(cmd); ierr=cmd;  EXEC_STOP_INFO(cmd)                      \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "GPU ERROR: ip[%d]: cannot execute: %s\n" ,ip, #cmd) ;\
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        report_error(ierr, stderr);                                             \
        ierr=-1;                                                                \
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;                                    \
        return( EXIT_FAILURE ) ;                                                \
    } }
    
#define file_operation( cmd )                                                   \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "FILE ERROR: ip[%d]: cannot execute: %s\n" ,ip, #cmd);\
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        report_error(ierr, stderr);                                             \
        ierr=-1;                                                                \
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;                                    \
        return( EXIT_FAILURE ) ;                                                \
    } }
    
#define file_operationl(cmd)                                             \
    {                                                                    \
        ierr = cmd;                                                      \
        if (ierr)                                                        \
        {                                                                \
            fprintf(stderr, "FILE ERROR:: cannot execute: %s\n", #cmd);  \
            fprintf(stderr, "file=`%s`, line=%d\n", __FILE__, __LINE__); \
            fprintf(stderr, "Error=%d\nExiting!\n", ierr);               \
            return (EXIT_FAILURE);                                       \
        }                                                                \
    }
    
#define cpu_exec( cmd )                                                         \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "CPU ERROR: ip[%d]: cannot execute: %s\n" ,ip, #cmd) ;\
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        report_error(ierr, stderr);                                             \
        ierr=-1;                                                                \
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;                                    \
        return( EXIT_FAILURE ) ;                                                \
    } }
    
#define cpu_execl( cmd )                                                        \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "CPU ERROR: cannot execute: %s\n", #cmd) ;            \
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        report_error(ierr, stderr);                                             \
        ierr=-1;                                                                \
        return( EXIT_FAILURE ) ;                                                \
    } }
    
#define error_msg_mpi_abort(ip,msg)                                             \
    {                                                                           \
        fprintf( stderr , "ERROR [%d]: file=`%s`, line=%d: %s\n" ,ip, __FILE__,__LINE__, #msg) ;   \
        ierr=-1;                                                                \
        MPI_Abort( MPI_COMM_WORLD , ierr ) ;                                    \
        return( EXIT_FAILURE ) ;                                                \
    }
    
#define ixyz2ixiyiz(ixyz,_ix,_iy,_iz,i)     \
    i=ixyz;                                 \
    _ix=i/(NY*NZ);                          \
    i=i-_ix * NY * NZ;                      \
    _iy=i/NZ;                               \
    _iz=i-_iy * NZ;
    
#define ixyz2ix(ixyz,_ix,i)                 \
    i=ixyz;                                 \
    _ix=i/(NY*NZ);                          
    
#define ixyz2ixiy(ixyz,_ix,_iy,i)           \
    i=ixyz;                                 \
    _ix=i/(NY*NZ);                          \
    i=i-_ix * NY * NZ;                      \
    _iy=i/NZ;               
    
#define ixyz2ixiyizD2Z(ixyz,_ix,_iy,_iz,i)  \
    i=ixyz;                                 \
    _ix=i/(NY*(NZ/2+1));                    \
    i=i-_ix * NY * (NZ/2+1);                \
    _iy=i/(NZ/2+1);                         \
    _iz=i-_iy * (NZ/2+1);

// -------------------- 2D versions of decoding functions ------------------------------------
#define ixy2ixiy2d(ixy,_ix,_iy) \
    _ix=ixy/NY;                 \
    _iy=ixy-_ix * NY;   
    
#define ixy2ixiy2dD2Z(ixy,_ix,_iy) \
    _ix=ixy/(NY/2+1);              \
    _iy=ixy-_ix * (NY/2+1); 
    
