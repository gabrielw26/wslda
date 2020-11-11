/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 15-10-2020
 * 
 * Declaration file for errors reporting
 * */ 
#ifndef __WSLDA_ERRORS__
#define __WSLDA_ERRORS__

#include <stdio.h>
void report_error(int errcode, FILE *stream);

// no error
#define WSLDA_OK  0

// errors
#define WSLDA_ERR_CANNOT_CREATE_DIR 10001
#define WSLDA_ERR_CANNOT_OPEN_FILE 10002
#define WSLDA_ERR_CANNOT_OVERWRITE 10003
#define WSLDA_ERR_BINARY_FILE_CORRUPTED 10004
#define WSLDA_ERR_INOCRRECT_PQ 10005
#define WSLDA_ERR_S3DPCA_INFO_FILES 10006
#define WSLDA_ERR_S3DPCA_INFO_FILES_MISSING_FILE 10007
#define WSLDA_ERR_ABDG_NOT_SET 10008
#define WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_DATA 10009 
#define WSLDA_ERR_TD_CANNOT_LOAD_CHECKPOINT_NWF 10010

#endif
