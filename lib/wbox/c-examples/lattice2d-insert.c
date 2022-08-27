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

double gaussian_2d(const int ix, const int iy, const int nx, const int ny)
{
#define SQRT_2PI2 4.4428829
    const double sigma = (double)(sqrt(nx * ny)) * 0.1;
    const int shift_x = (double)(nx - 1) / 2.;
    const int shift_y = (double)(ny - 1) / 2.;
    const double exponent = pow((double)(ix - shift_x), 2.) * pow((double)(iy - shift_y), 2.);
    return 4. / (sigma * SQRT_2PI2) * exp(-0.5 * exponent / (sigma * sigma));
}

int main()
{
    wbox_md_t wbmd_small;    // WBox variable that store information of SMALL box
    wbox_md_t wbmd_big;      // WBox variable that store information of BIG box
    wbox_insert_t wb_insert; // WBox variable that stores information about the ANCHOR of SMALL to BIG box

    int ixy;
    int nxs = 10, nys = 20, nzs = -1; // LATTICE parameters of SMALL box
    int nxb = 20, nyb = 20, nzb = -1; // LATTICE parameters of BIG box
    int ax = 5, ay = 3, az = -1;      // ANCHOR parameters in `units` of BIG box. Iteration is aq: 0-> nqb-nqs, here `q` is a cartesian coordinate X, Y or Z

    double *wdata_small; // Pointer to SMALL box variable
    double *wdata_big;   // Pointer to BIG box variable

    // MEMORY ALLOCATION
    cppmallocl(wdata_small, nxs * nys, double);
    cppmallocl(wdata_big, nxb * nyb, double);

    // Fill data in SMALL box
    for (size_t ixs = 0; ixs < nxs; ixs++)
        for (size_t iys = 0; iys < nys; iys++)
        {
            ixy = ixs * nys + iys;
            wdata_small[ixy] = gaussian_2d(ixs, iys, nxs, nys);
        }
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