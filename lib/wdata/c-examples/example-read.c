/** 
 * W-DATA format
 * @author Gabriel Wlazlowski, Warsaw University of Technology, 2020
 * 
 * Code demonstrates how to read wdata files
 * */

// wdata lib
#include "wdata.h"

// other libs
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

typedef double complex Complex;
typedef float  complex Complexf;

#include "example-utils.h"

int main()
{
    int ierr;

    // create metadata handler
    wdata_metadata md;
    wdata_reset_metadata(&md);

    // read metadata from file
    ierr = wdata_parse_metadata_file("test.wtxt", &md);
    if (ierr != 0)
    {
        printf("Cannot read metadata file! ERROR: #%d\n", ierr);
        return 1;
    }

    // print metadata file
    wdata_print_metadata(&md, stdout);

    // allocate memory for variables
    int bdim = wdata_get_blocklength(&md); // get block size

    // buffers for double variables
    double *dataR;  // real data
    Complex *dataC; // complex data
    double *dataV;  // vector data
    cppmallocl(dataR, bdim, double);
    cppmallocl(dataC, bdim, Complex);
    cppmallocl(dataV, bdim * 3, double); // factor 3 accounts for three compoments of vector variable

    // buffers for float variables
    float *dataRf;  // real data
    Complexf *dataCf; // complex data
    float *dataVf;  // vector data
    cppmallocl(dataRf, bdim, float);
    cppmallocl(dataCf, bdim, Complexf);
    cppmallocl(dataVf, bdim * 3, float); // factor 3 accounts for three compoments of vector variable

    int icycle;
    for (icycle = 0; icycle < md.cycles; icycle++)
    {
        printf("Processing cycle %d\n", icycle);

        double current_time;
        wdata_get_time(&md, icycle, &current_time);
        printf("# Current time is: %lf\n", current_time);

        // standard method of reading
        wdata_read_cycle(&md, "density_a", icycle, dataR);
        wdata_read_cycle(&md, "delta", icycle, dataC);
        wdata_read_cycle(&md, "current_a", icycle, dataV);

        // read variable in float
        wdata_read_cycle(&md, "density_f", icycle, dataRf);
        wdata_read_cycle(&md, "delta_f", icycle, dataCf);
        wdata_read_cycle(&md, "current_f", icycle, dataVf);

        // check correctness
        test_array_diff_df(bdim, dataR, dataRf);
        test_array_diff_df(bdim*2, (double*)dataC, (float *)dataCf);
        test_array_diff_df(bdim*3, dataV, dataVf);

        // read float variables and cast them into doubles
        wdata_read_cycle_d(&md, "density_f", icycle, dataR);
        wdata_read_cycle_d(&md, "delta_f", icycle, (double*)dataC);
        wdata_read_cycle_d(&md, "current_f", icycle, dataV);

        // check correctness
        test_array_diff_df(bdim, dataR, dataRf);
        test_array_diff_df(bdim*2, (double*)dataC, (float *)dataCf);
        test_array_diff_df(bdim*3, dataV, dataVf);

        // read double variables and cast them into floats
        wdata_read_cycle_f(&md, "density_a", icycle, dataRf);
        wdata_read_cycle_f(&md, "delta", icycle, (float*)dataCf);
        wdata_read_cycle_f(&md, "current_a", icycle, dataVf);

        // check correctness
        test_array_diff_df(bdim, dataR, dataRf);
        test_array_diff_df(bdim*2, (double*)dataC, (float *)dataCf);
        test_array_diff_df(bdim*3, dataV, dataVf);
    }

    // extract constants
    double eF = wdata_getconst_value(&md, "eF");
    printf("Constant eF=%f\n", eF);

    return 0;
}
