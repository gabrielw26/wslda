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
        
        default: 
            fprintf(stream, "\tThis error doesn not have description.\n");
    }
    fprintf(stream, "IF THIS INFORMATION IS NOT SUFFICIENT TO SOLVE YOUR PROBLEM\n");
    fprintf(stream, "\tCheck wikipages: http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/home\n");
    fprintf(stream, "\tAsk for help WSLDA developer:\n");
    fprintf(stream, "\t\tusing Issues reporting system: http://git2.if.pw.edu.pl/gabrielw/cold-atoms/issues\n");
    fprintf(stream, "\t\tsend request by e-mail: wslda@pw.edu.pl\n");
    fprintf(stream, "==========================================================================\n");
    
}
