/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 15-10-2020
 * 
 * Declaration file for errors reporting
 * */ 


#include "wslda_errors.h"

void wprintf( const char * format, ... );
void wfprintf(FILE *stream,  const char * format, ... );

#include "jdb.h"
extern int wsldapid; // process id - global variable
void something_to_cheer_you_up(FILE *stream)
{
    // wfprintf(stream, "# ==========================================================================\n");
    wfprintf(stream, "# SOMETHING FUNNY FROM W-SLDA FOR A GOOD START:\n");
    jdb_message(stream);
    // wfprintf(stream, "# ==========================================================================\n");
}
void something_to_cheer_you_up_pid0(FILE *stream)
{
    if(wsldapid==0) something_to_cheer_you_up(stream);
}

void report_error(int errcode, FILE *stream)
{
    if(errcode==WSLDA_OK) return; // no reporting
    
    wfprintf(stream, "==========================================================================\n");
    wfprintf(stream, "WSLDA ERROR DESCRIPTION:\n");
    switch(errcode)
    {
        case WSLDA_ERR_CANNOT_CREATE_DIR:
            wfprintf(stream, "\tCannot create directory for binary files.\n");
            wfprintf(stream, "\tCheck your permissions.\n");
            break;
            
        case WSLDA_ERR_CANNOT_OPEN_FILE:
            wfprintf(stream, "\tCannot open binary file.\n");
            wfprintf(stream, "\tCheck if `inprefix` tag is correctly set in input file.\n");
            wfprintf(stream, "\tYou can do it by executing\n");
            wfprintf(stream, "\t\tls inprefix\n");
            wfprintf(stream, "\tand checking if there are binary files *.wfu and *wfv.\n");
            wfprintf(stream, "\tNote that these files are created only if `writewf` is set to 1\n");
            wfprintf(stream, "\tin the static calculations.\n");
            break;
            
        case WSLDA_ERR_CANNOT_OVERWRITE:
            wfprintf(stream, "\tAn attempt of overwriting existing file has taken.\n");
            wfprintf(stream, "\tInput file tag overwrite=0 does not allow for this.\n");
            wfprintf(stream, "\tChange overwrite tag or outprefix tag in the input file.\n");
            break;
            
        case WSLDA_ERR_BINARY_FILE_CORRUPTED:
            wfprintf(stream, "\tInput binary files do not satisfy expected sum-rules.\n");
            wfprintf(stream, "\tThey seems to be damaged.\n");
            break;
            
        case WSLDA_ERR_INOCRRECT_PQ:
            wfprintf(stream, "\tParameters p and q in input file are set incorrectly.\n");
            wfprintf(stream, "\tCorrect input file.\n");
            break;
            
        case WSLDA_ERR_S3DPCA_INFO_FILES:
            wfprintf(stream, "\tInput binary files do not satisfy expected sum-rules.\n");
            wfprintf(stream, "\tThis situation may occur if value of iogroups used in static code (st)\n");
            wfprintf(stream, "\t   is different from value of iogroups provided for time-dependent (td) code.\n");
            wfprintf(stream, "\t   Solution: Correct iogropus in input file.\n");
            wfprintf(stream, "\tThe problem can also appear if `writeecut` has been applied in static calculations.\n");
            wfprintf(stream, "\t   Then not all wave-functions needed for time evolution were written to disk.\n");
            wfprintf(stream, "\t   Solution: Regenerate wave-functions with deactivated `writeecut`.\n");
            wfprintf(stream, "\tIf cases listed above do not apply this error indicates that binary files\n");
            wfprintf(stream, "\t   with wave-functions from static code may be damaged.\n");
            break;
            
        case WSLDA_ERR_SCAN_INFO_FILES_SUM_FAILED:
            wfprintf(stream, "\tInput binary files do not satisfy expected sum-rules.\n");
            wfprintf(stream, "\tThis situation may occur if value of iogroups is set too large.\n");
            wfprintf(stream, "\t   Solution: Decrease value of iogroups.\n");
            wfprintf(stream, "\tThe problem can also appear if `writeecut` has been applied in static calculations.\n");
            wfprintf(stream, "\t   Then not all wave-functions needed for time evolution were written to disk.\n");
            wfprintf(stream, "\t   Solution: Regenerate wave-functions with deactivated `writeecut`.\n");
            wfprintf(stream, "\tIf cases listed above do not apply this error indicates that binary files\n");
            wfprintf(stream, "\t   with wave-functions from static code may be damaged.\n");
            break;
            
        case WSLDA_ERR_S3DPCA_INFO_FILES_MISSING_FILE:
            wfprintf(stream, "\tInput binary files do not satisfy expected sum-rules.\n");
            wfprintf(stream, "\tThis situation may occur if value of iogroups is set incorrectly.\n");
            wfprintf(stream, "\tIn general it should have the same value as used in static code (st).\n");
            wfprintf(stream, "\tIn some cases static code (st) cannot generate requested iogroups (problem is too small).\n");
            wfprintf(stream, "\tCheck maximum value of number appearing in file names `inprefix`\\s3dpca.????.info\n");
            wfprintf(stream, "\t  and in input file set iogroups as (maximum value+1).\n");
            wfprintf(stream, "\tIf iogroups is set correctly this error indicates that binary files\n");
            wfprintf(stream, "\t   with wave-functions from static code may to be damaged.\n");
            break;
            
        case WSLDA_ERR_ABDG_NOT_SET:
            wfprintf(stream, "\tFor calculations with FUNCTIONAL in [BDG, SLDAE] \n");
            wfprintf(stream, "\tit is required to set sclgth in the input file!\n");
            break;
            
        case WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_DATA:
            wfprintf(stream, "\tTD code cannot load data from checkpoint files!\n");
            wfprintf(stream, "\tSuggestions that may help to solve the problem:\n");
            wfprintf(stream, "\t\t- inprefix should point to folder with checkpoint files\n");
            wfprintf(stream, "\t\t- you should have read permission to checkpoint files\n");
            wfprintf(stream, "\t\t- make sure you use the same value of MPI_NP_PER_IO_GROUP as you used for writing\n");
            break;
            
        case WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_NWF:
            wfprintf(stream, "\tTD code cannot read header of checkpoint files!\n");
            wfprintf(stream, "\tSuggestions that may help to solve the problem:\n");
            wfprintf(stream, "\t\t- inprefix should point to folder with checkpoint files\n");
            wfprintf(stream, "\t\t- you should have read permission to checkpoint files\n");
            wfprintf(stream, "\t\t- make sure you use the same value of MPI_NP_PER_IO_GROUP as you used for writing\n");
            break;
            
        case WSLDA_ERR_NAN_DETECTED:
            wfprintf(stream, "\tNot a Number (NaN) has been detected!\n");
            wfprintf(stream, "\tCheck settings of the code!\n");
            break;
            
        case WSLDA_ERR_INF_DETECTED:
            wfprintf(stream, "\tInfinite (Inf) has been detected!\n");
            wfprintf(stream, "\tCheck settings of the code!\n");
            break;
            
        case WSLDA_ERR_INTERPOLATION_NOT_IMPLEMENTED:
            wfprintf(stream, "\tUnsupported case of the interpolation!\n");
            wfprintf(stream, "\tSupported cases for 3D -- ALL sizes bigger/smaller then the imput sizes.\n");
            break;
            
        case WSLDA_ERR_CANNOT_CREATE_CHECKPOINT_FILE:
            wfprintf(stream, "\tCannot create the checkpoint file!\n");
            wfprintf(stream, "\tCheck if you have write permission to the target location.\n");
            break;
    
        case WSLDA_ERR_CANNOT_WRITETO_CHECKPOINT_FILE:
            wfprintf(stream, "\tCannot add entry to the checkpoint file!\n");
            wfprintf(stream, "\tCheck if you have disk space in the target location.\n");
            break;

        case WSLDA_ERR_CANNOT_OPEN_CHECKPOINT_FILE:
            wfprintf(stream, "\tCannot open the checkpoint file!\n");
            wfprintf(stream, "\tCheck if you have read permission to the target location.\n");
            break;    
            
        case WSLDA_ERR_CANNOT_READFROM_CHECKPOINT_FILE:
            wfprintf(stream, "\tCannot read entry from the checkpoint file!\n");
            wfprintf(stream, "\tCheck the correctness of the checkpoint file. It may be corrupted.\n");
            break;
            
        case WSLDA_ERR_INCOMPATIBLE_CHECKPOINT_FILE:
            wfprintf(stream, "\tThe checkpoint file is incompatible with the code settings!\n");
            wfprintf(stream, "\tYou cannot use provided checkpoint file to initialize the code!\n");
            break;
            
        case WSLDA_ERR_INTRISTIC_ERROR:
            wfprintf(stream, "\tIt is intrinsic error of W-SLDA Toolkit.\n");
            wfprintf(stream, "\tIt shouldn't have happened, but it did :-(\n");
            wfprintf(stream, "\tPlease report this error to W-SLDA Teams and help us to improve the Toolkit.\n");
            break;
            
        case WSLDA_ERR_SPINSYMMETRY1:
            wfprintf(stream, "\t The code is compiled with SPINSYMMETRY_MODE option in predefines.h,\n");
            wfprintf(stream, "\t while the provided wave-functions contain SPINA and SPINB components separately.\n");
            wfprintf(stream, "\t Recreate wave-functions with selected `spinsymmetry 1` and rerun td code again,\n");
            wfprintf(stream, "\t\t or\n");
            wfprintf(stream, "\t Recompile td code with commented out SPINSYMMETRY_MODE, and rerun it again.\n");
            break;
            
        case WSLDA_ERR_SPINSYMMETRY0:
            wfprintf(stream, "\t The code is compiled with commented out SPINSYMMETRY_MODE option in predefines.h,\n");
            wfprintf(stream, "\t while the provided wave-functions contain only SPINB component.\n");
            wfprintf(stream, "\t Recreate wave-functions with selected `spinsymmetry 0` and rerun td code again,\n");
            wfprintf(stream, "\t\t or\n");
            wfprintf(stream, "\t Recompile td code with selected SPINSYMMETRY_MODE, and rerun it again.\n");
            break;

        case WSLDA_ERR_INITSTATE_FILE_NOT_FOUND:
            wfprintf(stream, "\t The specified initial state file cannot be found!\n");
            wfprintf(stream, "\t Check if the path to the file is correct.\n");
            wfprintf(stream, "\t Execute command: `ls inprefix` to check if the file exists.\n");
            break;

        case WSLDA_ERR_INITSTATE_FOR_1D:
            wfprintf(stream, "#\t The selected initial state is for 1D calculations!\n");
            wfprintf(stream, "#\t It is not compatible with your code!\n");
            break;
        
        case WSLDA_ERR_INITSTATE_FOR_2D:
            wfprintf(stream, "#\t The selected initial state is for 2D calculations!\n");
            wfprintf(stream, "#\t It is not compatible with your code!\n");
            break; 

        case WSLDA_ERR_INITSTATE_FOR_3D:
            wfprintf(stream, "#\t The selected initial state is for 3D calculations!\n");
            wfprintf(stream, "#\t It is not compatible with your code!\n");
            break; 

        case WSLDA_ERR_MISSING_FILE_FOR_INITSTATE:
            wfprintf(stream, "#\t Some files required for the selected initial state are missing!\n");
            wfprintf(stream, "#\t Make sure that all required files are present in the specified directory.\n");
            wfprintf(stream, "#\t Execute command: `ls inprefix` to check available files.\n");
            break;

        case WSLDA_ERR_MISSING_FILE_FOR_INITSTATE_KZADD:
            wfprintf(stream, "#\t Some files required for the selected initial state are missing!\n");
            wfprintf(stream, "#\t Make sure that all required files are present in the specified directory.\n");
            wfprintf(stream, "#\t Execute command: `ls inprefix` to check available files.\n");
            wfprintf(stream, "#\t This error can also appear if there is mismatch of regularization schemes between static and time-dependent codes.\n");
            break;
            
        default: 
            wfprintf(stream, "\tThis error does not have description.\n");
    }
    wfprintf(stream, "IF THIS INFORMATION IS NOT SUFFICIENT TO SOLVE YOUR PROBLEM\n");
    wfprintf(stream, "\tCheck wiki pages: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/home\n");
    wfprintf(stream, "\tAsk for help WSLDA developers:\n");
    wfprintf(stream, "\t\tusing Issues reporting system: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/issues\n");
    wfprintf(stream, "\t\tor by e-mail: wslda@fizyka.pw.edu.pl\n");
}

void warn_head(FILE *stream)
{
    wfprintf(stream, "# !!! --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- WARNING --- !!!\n");
}
void warn_foot(FILE *stream)
{
    wfprintf(stream, "# !!! --- ------- --- ------- --- ------- --- ------- --- ------- --- ------- --- !!!\n");
}

void report_warning(int errcode, FILE *stream)
{
    if(errcode==WSLDA_OK) return; // no reporting
    
    warn_head(stream);

    switch(errcode)
    {
        case WSLDA_WRN_SPINSYMMETRY0:
            wfprintf(stream, "#\t Input file setting: `spinsymmetry 0` not comptible with predefines.h option SPINSYMMETRY_MODE!\n");
            wfprintf(stream, "#\t Check if it is intended!\n");
            wfprintf(stream, "#\t To avoid the code termination forcing: `spinsymmetry 1`!\n");
            break;
            
        case WSLDA_WRN_SPINSYMMETRY1:
            wfprintf(stream, "#\t Input file setting: `spinsymmetry 1` not comptible with commented-out predefines.h option SPINSYMMETRY_MODE!\n");
            wfprintf(stream, "#\t Check if it is intended!\n");
            wfprintf(stream, "#\t To avoid the code termination forcing: `spinsymmetry 0`!\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_NOT_CONSITENT_BROYDEN:
            wfprintf(stream, "#\t Binary data for Boyden algorithm not consistent with current settings!\n");
            wfprintf(stream, "#\t The data will NOT be loaded!\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_DOINTERPOLATION:
            wfprintf(stream, "#\t Resolution of the lattice has changed!\n");
            wfprintf(stream, "#\t The code will interpolate given checkpoint data to the new resolution.\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_DORESIZE:
            wfprintf(stream, "#\t Dimensionality of the lattice has changed!\n");
            wfprintf(stream, "#\t The code will change the dimensionality of the given checkpoint data to the new lattice.\n");
            break;
            
        case WSLDA_WRN_CHECKPOINT_UNPREDICTED:
            wfprintf(stream, "#\t Checkpoint files is incompatible with code settings!\n");
            wfprintf(stream, "#\t The code will upload the data however the result of this operation may be unpredictable!\n");
            wfprintf(stream, "#\t Make sure you understand what you are doing!\n");
            break;

        case WSLDA_WRN_DIFFERENT_DXDYDZ:
            wfprintf(stream, "#\t You are using code with different lattice spacing (DX,DY,DZ).\n");
            wfprintf(stream, "#\t This type of usage is not recommended for non-expert users,\n");
            wfprintf(stream, "#\t as it requires deep knowledge of the regularization procedure.\n");
            break; 
            
        default: 
            wfprintf(stream, "#\tThis warning does not have description.\n");
    }
    warn_foot(stream);
 
}

void print_warning(int errcode)
{
    report_warning(errcode, stdout); // to output 
    report_warning(errcode, stderr); // and to error process
}
