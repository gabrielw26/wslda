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
            fprintf(stream, "\tInput file tag overwrite=0 does not allow for this.\n");
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
            
        case WSLDA_ERR_NAN_DETECTED:
            fprintf(stream, "\tNot a Number (NaN) has been detected!\n");
            fprintf(stream, "\tCheck settings of the code!\n");
            break;
            
        case WSLDA_ERR_INF_DETECTED:
            fprintf(stream, "\tInfinite (Inf) has been detected!\n");
            fprintf(stream, "\tCheck settings of the code!\n");
            break;
            
        case WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED:
            fprintf(stream, "\tUnsupported case of the interpolation!\n");
            fprintf(stream, "\tSupported cases for 3D -- ALL sizes bigger/smaller then the imput sizes.\n");
            break;
            
        case WSLDA_ERR_CANNOT_CREATE_CHECKPOINT_FILE:
            fprintf(stream, "\tCannot create the checkpoint file!\n");
            fprintf(stream, "\tCheck if you have write permission to the target location.\n");
            break;
    
        case WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE:
            fprintf(stream, "\tCannot add entry to the checkpoint file!\n");
            fprintf(stream, "\tCheck if you have disk space in the target location.\n");
            break;

        case WSLDA_ERR_CANNOT_OPEN_CHECKPOINT_FILE:
            fprintf(stream, "\tCannot open the checkpoint file!\n");
            fprintf(stream, "\tCheck if you have read permission to the target location.\n");
            break;    
            
        case WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE:
            fprintf(stream, "\tCannot read entry from the checkpoint file!\n");
            fprintf(stream, "\tCheck the correctness of the checkpoint file. It may be corrupted.\n");
            break;
            
        case WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE:
            fprintf(stream, "\tThe checkpoint file is incompatible with the code settings!\n");
            fprintf(stream, "\tYou cannot use provided checkpoint file to initialize the code!\n");
            break;
            
        case WSLDA_ERR_INTRISTIC_ERROR:
            fprintf(stream, "\tIt is intrinsic error of W-SLDA Toolkit.\n");
            fprintf(stream, "\tIt shouldn't have happened, but it did :-(\n");
            fprintf(stream, "\tPlease report this error to W-SLDA Teams and help us to improve the Toolkit.\n");
            break;
            
        default: 
            fprintf(stream, "\tThis error does not have description.\n");
    }
    fprintf(stream, "IF THIS INFORMATION IS NOT SUFFICIENT TO SOLVE YOUR PROBLEM\n");
    fprintf(stream, "\tCheck wiki pages: https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/home\n");
    fprintf(stream, "\tAsk for help WSLDA developers:\n");
    fprintf(stream, "\t\tusing Issues reporting system: https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/issues\n");
    fprintf(stream, "\t\tor by e-mail: wslda@fizyka.pw.edu.pl\n");
    fprintf(stream, "==========================================================================\n");
    
}

void warn_head(FILE *stream)
{
    fprintf(stream, "# !!! --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- !!!\n");
}
void warn_foot(FILE *stream)
{
    fprintf(stream, "# !!! --- ------- --- ------- --- ------- --- ------- --- ------- --- ------- --- !!!\n");
}

void report_warning(int errcode, FILE *stream)
{
    if(errcode==WSLDA_OK) return; // no reporting
    
    warn_head(stream);

    switch(errcode)
    {
        case WSLDA_WRN_SPINSYMMETRY0:
            fprintf(stream, "#\t Input file setting: `spinsymmetry 0` not comptible with predefines.h option SPINSYMMETRY_MODE!\n");
            fprintf(stream, "#\t Check if it is intended!\n");
            fprintf(stream, "#\t To avoid the code termination forcing: `spinsymmetry 1`!\n");
            break;
            
        case WSLDA_WRN_SPINSYMMETRY1:
            fprintf(stream, "#\t Input file setting: `spinsymmetry 1` not comptible with predefines.h option SPINSYMMETRY_MODE!\n");
            fprintf(stream, "#\t Check if it is intended!\n");
            fprintf(stream, "#\t To avoid the code termination forcing: `spinsymmetry 1`!\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_NOT_CONSITENT_BROYDEN:
            fprintf(stream, "#\t Binary data for Boyden algorithm not consistent with current settings!\n");
            fprintf(stream, "#\t The data will NOT be loaded!\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_DOINTERPOLATION:
            fprintf(stream, "#\t Resolution of the lattice has changed!\n");
            fprintf(stream, "#\t The code will interpolate given checkpoint data to the new resolution.\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_DORESIZE:
            fprintf(stream, "#\t Dimensonality of the lattice has changed!\n");
            fprintf(stream, "#\t The code will change dimensionality of given checkpoint data to the new lattice.\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_UNPREDICTED:
            fprintf(stream, "#\t Checkpoint files is incompatible with code settings!\n");
            fprintf(stream, "#\t The code will upload the data however the result of this operation may be unpredictable!\n");
            fprintf(stream, "#\t Make sure you understand what you are doing!\n");
            break;
            
        default: 
            fprintf(stream, "#\tThis warning does not have description.\n");
    }
    warn_foot(stream);
 
}

void print_warning(int errcode)
{
    report_warning(errcode, stdout); // to output 
    report_warning(errcode, stderr); // and to error process
}
