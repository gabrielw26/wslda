/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 15-10-2020
 * 
 * Declaration file for errors reporting
 * */ 


#include "wslda_errors.h"

void report_error(int errcode, FILE *stream)
{
    if(errcode==WSLDA_OK) return; // no reporting
    
    fprintf(stream, "==========================================================================\n");
    fprintf(stream, "WSLDA ERROR DESCRIPTION:\n");
    switch(errcode)
    {
        case WSLDA_ERR_CANNOT_CREATE_DIR:
            fprintf(stream, "\tCannot create directory for binary files.\n");
            fprintf(stream, "\tCheck your permissions.\n");
            break;
            
        case WSLDA_ERR_CANNOT_OPEN_FILE:
            fprintf(stream, "\tCannot open binary file.\n");
            fprintf(stream, "\tCheck if `inprefix` tag is correctly set in input file.\n");
            break;
            
        case WSLDA_ERR_CANNOT_OVERWRITE:
            fprintf(stream, "\tAttempt of overwriting existing file has taken.\n");
            fprintf(stream, "\tInput file tag overwrite=1 does not allow for this.\n");
            fprintf(stream, "\tChange overwrite tag or outprefix tag in input file.\n");
            break;
            
        case WSLDA_ERR_BINARY_FILE_CORRUPTED:
            fprintf(stream, "\tInput binary files do not satisfy expected sum-rules.\n");
            fprintf(stream, "\tThey seems to be damaged.\n");
            break;
            
        case WSLDA_ERR_INOCRRECT_PQ:
            fprintf(stream, "\tParameters p and q in input file are set incorrectly.\n");
            fprintf(stream, "\tCorrect input file.\n");
            break;
            
        case WSLDA_ERR_S3DPCA_INFO_FILES:
            fprintf(stream, "\tInput binary files do not satisfy expected sum-rules.\n");
            fprintf(stream, "\tThis situation may occur if value of iogroups used in static code (st)\n");
            fprintf(stream, "\t   is different from value of iogroups provided for time-dependent (td) code.\n");
            fprintf(stream, "\tCorrect iogropus in input file.\n");
            fprintf(stream, "\tIf iogroups is set correctly this error indicates that binary files\n");
            fprintf(stream, "\t   with wave-functions from static code may to be damaged.\n");
            break;
            
        case WSLDA_ERR_S3DPCA_INFO_FILES_MISSING_FILE:
            fprintf(stream, "\tInput binary files do not satisfy expected sum-rules.\n");
            fprintf(stream, "\tThis situation may occur if value of iogroups is set incorrectly.\n");
            fprintf(stream, "\tIn general it should have the same value as used in static code (st).\n");
            fprintf(stream, "\tIn some cases static code (st) cannot generate requested iogroups (problem is too small).\n");
            fprintf(stream, "\tCheck maximum value of number appearing in file names `inprefix`\\s3dpca.????.info\n");
            fprintf(stream, "\t  and in input file set iogroups as (maximum value+1).\n");
            fprintf(stream, "\tIf iogroups is set correctly this error indicates that binary files\n");
            fprintf(stream, "\t   with wave-functions from static code may to be damaged.\n");
            break;
            
        case WSLDA_ERR_ABDG_NOT_SET:
            fprintf(stream, "\tFor calculations with FUNCTIONAL==BDG\n");
            fprintf(stream, "\tit is required to set aBdG in input file!\n");
            break;
            
        case WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_DATA:
            fprintf(stream, "\tTD code cannot load data from checkpoint files!\n");
            fprintf(stream, "\tSuggestions that may help to solve the problem:\n");
            fprintf(stream, "\t\t- inprefix should point to folder with checkpoint files\n");
            fprintf(stream, "\t\t- you should have read permission to checkpoint files\n");
            fprintf(stream, "\t\t- make sure you use the same value of MPI_NP_PER_IO_GROUP as you used for writing\n");
            break;
            
        case WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_NWF:
            fprintf(stream, "\tTD code cannot read header of checkpoint files!\n");
            fprintf(stream, "\tSuggestions that may help to solve the problem:\n");
            fprintf(stream, "\t\t- inprefix should point to folder with checkpoint files\n");
            fprintf(stream, "\t\t- you should have read permission to checkpoint files\n");
            fprintf(stream, "\t\t- make sure you use the same value of MPI_NP_PER_IO_GROUP as you used for writing\n");
            break;
        
        default: 
            fprintf(stream, "\tThis error doesn not have description.\n");
    }
    fprintf(stream, "IF THIS INFORMATION IS NOT SUFFICIENT TO SOLVE YOUR PROBLEM\n");
    fprintf(stream, "\tCheck wiki pages: https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/home\n");
    fprintf(stream, "\tAsk for help WSLDA developers:\n");
    fprintf(stream, "\t\tusing Issues reporting system: https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/issues\n");
    fprintf(stream, "\t\tor by e-mail: wslda@fizyka.pw.edu.pl\n");
    fprintf(stream, "==========================================================================\n");
    
}
