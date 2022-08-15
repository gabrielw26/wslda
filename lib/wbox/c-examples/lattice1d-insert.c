#include <complex.h>
#include <math.h>

#include "wbox.h"

#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "ERROR: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        return -1;                                                          \
    }

double gaussian_1d(const int ix, const int nx)
{
#define SQRT_2PI1 2.506628274631
    const double sigma = (double)sqrt(nx) * 0.2;
    const int shift = (double)(nx - 1) / 2.;
    return 4. / (sigma * SQRT_2PI1) * exp(-0.5 * pow((double)(ix - shift), 2.) / (sigma * sigma));
}

int main()
{
    wbox_md_t wbmd_small;    // WBox variable that store information of SMALL box
    wbox_md_t wbmd_big;      // WBox variable that store information of BIG box
    wbox_insert_t wb_insert; // WBox variable that stores information about the ANCHOR of SMALL to BIG box

    int nxs = 10, nys = -1, nzs = -1; // LATTICE parameters of SMALL box
    int nxb = 20, nyb = -1, nzb = -1; // LATTICE parameters of BIG box
    int ax = 0, ay = -1, az = -1;     // ANCHOR parameters in `units` of BIG box. Iteration is aq: 0-> nqb-nqs, here `q` is a cartesian coordinate X, Y or Z

    double *wdata_small; // Pointer to SMALL box variable
    double *wdata_big;   // Pointer to BIG box variable

    // MEMORY ALLOCATION
    cppmallocl(wdata_small, nxs, double);
    cppmallocl(wdata_big, nxb, double);

    // Fill data in SMALL box
    for (size_t ixs = 0; ixs < nxs; ixs++)
        wdata_small[ixs] = gaussian_1d(ixs, nxs);

    // WBox initialization. See `wbox_init()` description
    wbox_init(&wbmd_small, nxs, nys, nzs, 'R', wdata_small);
    wbox_init(&wbmd_big, nxb, nyb, nzb, 'r', wdata_big);

    // WBox insert initialization. See `wbox_insert_init()` description
    wbox_insert_init(&wb_insert, ax, ay, az);

    // PRINT WBox STRUCT
    wbox_print_box_params(wbmd_small, "WBOX SMALL");
    wbox_print_box_params(wbmd_big, "WBOX BIG");

    // PRINT WBox_Insert STRUCT
    wbox_print_insert_params(wb_insert);

    // Insert SMALL box into BIG box
    wbox_insert(&wbmd_small, &wb_insert, &wbmd_big);

    /* PRINTS for checks. */
    wbox_print(&wbmd_small, "WBOX SMALL");
    wbox_print(&wbmd_big, "WBOX BIG");

    return EXIT_SUCCESS;
}