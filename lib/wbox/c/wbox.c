#include <ctype.h>

#include "wbox.h"

int wbox_init(wbox_md_t *wbmd, const int nx, const int ny, const int nz, const char datatype, const void *data)
{
    const char dt = tolower(datatype);

    // INIT LATTICE
    wbmd->nx = nx;
    wbmd->ny = ny;
    wbmd->nz = nz;

    // INIT DIMENSION
    if (wbmd->nx <= 0 && wbmd->ny <= 0 && wbmd->nz <= 0)
    {
        wbmd->dim = 0;
        return WBOX_ERR_INIT_DIM;
    }
    else if (wbmd->ny <= 0 && wbmd->nz <= 0)
        wbmd->dim = 1;
    else if (wbmd->nz <= 0)
        wbmd->dim = 2;
    else
        wbmd->dim = 3;

    // INIT DATA
    if (dt == 'i')
        wbmd->data = (int *)data;
    else if (dt == 'r')
        wbmd->data = (double *)data;
    else if (dt == 'c')
        wbmd->data = (double complex *)data;
    else
        return WBOX_ERR_INIT_DATATYPE_NOT_IMPLEMENTED;

    wbmd->datatype = dt;

    return WBOX_SUCCESS_INIT;
}

int wbox_insert_1d(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big)
{
    if (wbox_insert_check(wbmd_small, wbox_insert, wbmd_big))
        return WBOX_ERR_INSERT_1D;

    // ITERATORS & SHIFTS
    int ixb = 0;
    int i_shift;
    const int nxs = wbmd_small->nx;
    const int nx_big = wbmd_big->nx;

    // LOCAL ANCHOR
    const int ax = wbox_insert->ax;

    // DATATYPE
    const char dt = wbmd_big->datatype;

    // INSERT LOOP
    for (ixb = 0; ixb < nx_big; ixb++)
    {
        i_shift = ixb - ax;
        if (ixb >= ax && ixb < ax + nxs) // FILL EDGE
        {
          if (dt == 'i')
          *((int *)wbmd_big->data + ixb) += *((int *)wbmd_small->data + i_shift);
          else if (dt == 'r')
          *((double *)wbmd_big->data + ixb) += *((double *)wbmd_small->data + i_shift);
          else if (dt == 'c')
          *((double complex *)wbmd_big->data + ixb) += *((double complex *)wbmd_small->data + i_shift);
        }
    }

    return WBOX_SUCCESS_INSERT;
}

int wbox_insert_2d(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big)
{
    if (wbox_insert_check(wbmd_small, wbox_insert, wbmd_big))
        return WBOX_ERR_INSERT_2D;

    // ITERATORS & SHIFTS
    int ixb = 0, iyb = 0;
    int ixy_big, i_shift;
    const int nxs = wbmd_small->nx;
    const int nys = wbmd_small->ny;
    const int nx_big = wbmd_big->nx;
    const int ny_big = wbmd_big->ny;

    // LOCAL ANCHORS
    const int ax = wbox_insert->ax;
    const int ay = wbox_insert->ay;

    // DATATYPES
    const char dt = wbmd_big->datatype;

    for (ixb = 0; ixb < nx_big; ixb++)
    {
        for (iyb = 0; iyb < ny_big; iyb++)
        {
            ixy_big = (ixb * ny_big) + iyb;
            i_shift = (ixb - ax) * nys + iyb - ay;

            if ((ixb >= ax && ixb < ax + nxs) &&
                (iyb >= ay && iyb < ay + nys))
            {
              if (dt == 'i')
                  *((int *)wbmd_big->data + ixy_big) += *((int *)wbmd_small->data + i_shift);
              else if (dt == 'r')
                  *((double *)wbmd_big->data + ixy_big) += *((double *)wbmd_small->data + i_shift);
              else if (dt == 'c')
                  *((double complex *)wbmd_big->data + ixy_big) += *((double complex *)wbmd_small->data + i_shift);
            }
        }
    }
    return WBOX_SUCCESS_INSERT;
}

int wbox_insert_3d(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big)
{
    if (wbox_insert_check(wbmd_small, wbox_insert, wbmd_big))
        return WBOX_ERR_INSERT_2D;

    // ITERATORS & SHIFTS
    int ixb = 0, iyb = 0, izb = 0;
    int ixyz_big, i_shift;
    const int nxs = wbmd_small->nx;
    const int nys = wbmd_small->ny;
    const int nzs = wbmd_small->nz;
    const int nxyzs = nxs * nys * nzs;
    const int nx_big = wbmd_big->nx;
    const int ny_big = wbmd_big->ny;
    const int nz_big = wbmd_big->nz;

    // LOCAL ANCHORS
    const int ax = wbox_insert->ax;
    const int ay = wbox_insert->ay;
    const int az = wbox_insert->az;

    // DATATYPES
    const char dt = wbmd_big->datatype;

    for (ixb = 0; ixb < nx_big; ixb++)
    {
        for (iyb = 0; iyb < ny_big; iyb++)
        {
            for (izb = 0; izb < nz_big; izb++)
            {
                ixyz_big = (ixb * ny_big * nz_big) + (iyb * nz_big) + izb;
                i_shift = ((ixb - ax) * nys * nzs) + ((iyb - ay) * nzs) + (izb - az);
                if ((ixb >= ax && ixb < ax + nxs) &&
                    (iyb >= ay && iyb < ay + nys) &&
                    (izb >= az && izb < az + nzs))
                {
                  if (dt == 'i')
                      *((int *)wbmd_big->data + ixyz_big) += *((int *)wbmd_small->data + i_shift);
                  else if (dt == 'r')
                      *((double *)wbmd_big->data + ixyz_big) += *((double *)wbmd_small->data + i_shift);
                  else if (dt == 'c')
                      *((double complex *)wbmd_big->data + ixyz_big) += *((double complex *)wbmd_small->data + i_shift);
                }
            }
        }
    }
    return WBOX_SUCCESS_INSERT;
}

int wbox_insert(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big)
{
    int wbox_insert_return;
    if (wbmd_small->dim == 1)
        wbox_insert_return = wbox_insert_1d(&wbmd_small, &wbox_insert, &wbmd_big);
    else if (wbmd_small->dim == 2)
        wbox_insert_return = wbox_insert_2d(&wbmd_small, &wbox_insert, &wbmd_big);
    else if (wbmd_small->dim == 3)
        wbox_insert_return = wbox_insert_3d(&wbmd_small, &wbox_insert, &wbmd_big);

    return wbox_insert_return;
}

int wbox_insert_check(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, const wbox_md_t *wbmd_big)
{
    fprintf(stdout, "# WBOX INSERT CHECKLIST:\n");

    // DIMENSION
    if (wbmd_small->dim == wbmd_big->dim)
        fprintf(stdout, "\t# %-12s: OK!\n", "DIMENSION");
    else
    {
        fprintf(stdout, "\t# !!!ERROR!!! \n");
        fprintf(stdout, "\t     SMALL BOX DIM: %1d, BIG BOX DIM: %1d\n", wbmd_small->dim, wbmd_big->dim);
        fprintf(stdout, "\t     FIX: Match number of positive values of lattice size (NX, NY, NZ)\n");
        fprintf(stdout, "\t          in SMALL and BIG box.");

        return WBOX_ERR_INSERT_CHECK_DIM;
    }

    // todo: check anchor point, if anchor + small > big fail
    fprintf(stdout, "\t# %-12s: OK! (TODO)\n", "ANCHOR");

    // todo: check if data is of the same type
    fprintf(stdout, "\t# %-12s: OK! (TODO)\n", "DATATYPE");

    return WBOX_SUCCESS_INSERT;
}

int wbox_insert_init(wbox_insert_t *wbox_insert, const int ax, const int ay, const int az)
{
    if (ax < 0)
        return WBOX_ERR_INSERT_INIT_AX_MONE;
    else if (ay < 0)
        wbox_insert->dim = 1;
    else if (az < 0)
        wbox_insert->dim = 2;
    else
        wbox_insert->dim = 3;

    wbox_insert->ax = ax;
    wbox_insert->ay = ay;
    wbox_insert->az = az;

    return WBOX_SUCCESS_INSERT_INIT;
}

/* ***************************************************************************************** */
/* ***************************************************************************************** */
// int wbox_merge1d()
// {
// }

// int wbox_merge2d()
// {
// }

// int wbox_merge3d()
// {
// }

// int wbox_merge()
// {
// }

/* ***************************************************************************************** */
/* ***************************************************************************************** */
// int wbox_rotate2d()
// {
// }

// int wbox_rotate3d()
// {
// }

// int wbox_rotate()
// {
// }

/* ***************************************************************************************** */
/* ***************************************************************************************** */
int wbox_print(const wbox_md_t *wbmd, const char *title)
{
    if (wbmd->dim == 1)
        wbox_print_1d(wbmd, title);
    else if (wbmd->dim == 2)
        wbox_print_2d(wbmd, title);
    else if (wbmd->dim == 3)
        wbox_print_3d(wbmd, title);
    else
        return WBOX_ERR_PRINT_DIM;

    return WBOX_SUCCESS_PRINT;
}

int wbox_print_1d(const wbox_md_t *wbmd, const char *title)
{
    const int nx = wbmd->nx;

    const char dt = wbmd->datatype;

    printf("\n# PRINT: %s\n", title);

    for (size_t ix = 0; ix < nx; ix++)
    {
        if (dt == 'i')
            printf("[%-6d] ", *((int *)wbmd->data + ix));
        else if (dt == 'r')
            printf("[%-8.3lf] ", *((double *)wbmd->data + ix));
        else if (dt == 'c')
            printf("[%-4.2lf + i%4.2lf] ", creal(*((double complex *)wbmd->data + ix)), cimag(*((double complex *)wbmd->data + ix)));
    }
    printf("\n");
}

int wbox_print_2d(const wbox_md_t *wbmd, const char *title)
{
    const int nx = wbmd->nx;
    const int ny = wbmd->ny;
    int ixy;

    const char dt = wbmd->datatype;

    printf("\n# PRINT: %s\n", title);

    for (size_t iy = 0; iy < ny; iy++)
    {
        for (size_t ix = 0; ix < nx; ix++)
        {
            ixy = ix * ny + iy;

            if (dt == 'i')
                printf("[%-6d] ", *((int *)wbmd->data + ixy));
            else if (dt == 'r')
                printf("[%-6.3lf] ", *((double *)wbmd->data + ixy));
            else if (dt == 'c')
                printf("[%-4.2lf + i%4.2lf] ", creal(*((double complex *)wbmd->data + ixy)), cimag(*((double complex *)wbmd->data + ixy)));
        }
        printf("\n");
    }
}

int wbox_print_3d(const wbox_md_t *wbmd, const char *title)
{
    const int nx = wbmd->nx;
    const int ny = wbmd->ny;
    const int nz = wbmd->nz;
    int ixyz;

    const char dt = wbmd->datatype;

    printf("\n# PRINT: %s\n", title);

    for (size_t iz = 0; iz < nz; iz++)
    {
        printf("# IZ : %3ld\n", iz);
        for (size_t iy = 0; iy < ny; iy++)
        {
            for (size_t ix = 0; ix < nx; ix++)
            {
                ixyz = ix * ny * nz + iy * nz + iz;
                if (dt == 'i')
                    printf("[%-6d] ", *((int *)wbmd->data + ixyz));
                else if (dt == 'r')
                    printf("[%-8.3lf] ", *((double *)wbmd->data + ixyz));
                else if (dt == 'c')
                    printf("[%-4.2lf + i%4.2lf] ", creal(*((double complex *)wbmd->data + ixyz)), cimag(*((double complex *)wbmd->data + ixyz)));
            }
            printf("\n");
        }
        printf("\n");
    }
}

int wbox_print_box_params(const wbox_md_t wbmd, const char *title)
{
    printf("%-16s  (NX: %2d, NY: %2d, NZ: %2d), DIM: %d, &DATA: %p, DATATYPE: %c\n",
           title, wbmd.nx, wbmd.ny, wbmd.nz, wbmd.dim, wbmd.data, wbmd.datatype);

    return WBOX_SUCCESS_PRINT;
}
int wbox_print_insert_params(const wbox_insert_t wb_insert)
{
    printf("WBOX INSERT (AX: %2d, AY: %2d, AZ: %2d), DIM: %d\n",
           wb_insert.ax, wb_insert.ay, wb_insert.az, wb_insert.dim);

    return WBOX_SUCCESS_PRINT;
}



int wbox_boundary_avg_int(wbox_md_t *wbmd)
{
  const int nxs = wbmd->nx;
  const int nys = wbmd->ny;
  const int nzs = wbmd->nz;

  const int dim = wbmd->dim;
  printf("%d\n", dim);

  int mean = 0;
  int plane_mean = 0;
  int nof_points;
  if (dim == 3)
  {
    nof_points = 2 * (nxs * nys + nys * nzs + nxs * nzs) -4 * (nxs + nys + nzs) + 8;
    // One need to avoid summind doubly edges and vertices three times.
    // This is why sums do not always go from 0 to N
    // Summing up for all X,Y at Z=0 and Z=NZ
    for (int ixs = 0; ixs < nxs; ixs++)
    for (int iys = 0; iys < nys; iys++)
    {
      int ixys;
      ixys = (ixs * nys  + iys) * nzs;              // IZ = 0
      plane_mean += *((int *)wbmd->data + ixys); // BOTTOM XY plane
      ixys = (ixs * nys + iys) * nzs + (nzs - 1);  // IZ = NZ - 1
      plane_mean += *((int *)wbmd->data + ixys); // TOP XY plane
    }
    // Summing up for all X, Z:1...NZ-1 at Y=0 and Y=NY
    for (int ixs = 0; ixs < nxs; ixs++)
    for (int izs = 1; izs < nzs-1; izs++)
    {
      int ixzs;
      ixzs = ixs * nys * nzs + izs;                    // IY = 0
      plane_mean += *((int *)wbmd->data + ixzs); // LEFT YZ plane
      ixzs = (ixs * nys + (nys - 1)) * nzs + izs;  // IY = NY - 1
      plane_mean += *((int *)wbmd->data + ixzs); // RIGHT YZ plane
    }
    // Summing up  Y:1...NY-1, Z:1...NZ-1 at X=0 and X=NX
    for (int iys = 1; iys < nys-1; iys++)
    for (int izs = 1; izs < nzs-1; izs++)
    {
      int iyzs;
      iyzs = iys * nzs + izs;                          // IX = 0
      plane_mean += *((int *)wbmd->data + iyzs); // FRONT XZ plane
      iyzs = ((nxs - 1) * nys + iys) * nzs + izs;  // IX = NX - 1
      plane_mean += *((int *)wbmd->data + iyzs); // BACK XZ plane
    }
  }
  else if (dim == 2)
  {
    nof_points = 2 * (nxs + nys - 2);
    // One need to avoid summind doubly edges
    // This is why sums do not always go from 0 to N
    for (int ixs = 0; ixs < nxs; ixs++)
    {
        plane_mean += *((int *)wbmd->data + ixs);             // Bottom X edge
        plane_mean += *((int *)wbmd->data + ixs + nxs * (nys - 1)); // Top X edge
    }
    for (int iys = 1; iys < nys-1; iys++)
    {
        plane_mean += *((int *)wbmd->data + iys * nxs + 1);       // Left Y edge
        plane_mean += *((int *)wbmd->data + (iys - 2)* nxs); // Right Y edge
    }

  }
  else if (dim == 1)
  {
    return (*((int *)wbmd->data) + *((int *)wbmd->data + wbmd->nx - 1)) / 2; // Mean value of TWO edge points
  }


  mean = plane_mean / nof_points; // Mean value of SIX planes
  return mean;
}

double wbox_boundary_avg_double(wbox_md_t *wbmd)
{
  const int nxs = wbmd->nx;
  const int nys = wbmd->ny;
  const int nzs = wbmd->nz;

  const int dim = wbmd->dim;

  double mean = 0;
  double plane_mean = 0;
  int nof_points;
  if (dim == 3)
  {
    nof_points = 2 * (nxs * nys + nys * nzs + nxs * nzs) -4 * (nxs + nys + nzs) + 8;
    // One need to avoid summind doubly edges and vertices three times.
    // This is why sums do not always go from 0 to N
    // Summing up for all X,Y at Z=0 and Z=NZ
    for (int ixs = 0; ixs < nxs; ixs++)
    for (int iys = 0; iys < nys; iys++)
    {
      int ixys;
      ixys = (ixs * nys  + iys) * nzs;              // IZ = 0
      plane_mean += *((double *)wbmd->data + ixys); // BOTTOM XY plane
      ixys = (ixs * nys + iys) * nzs + (nzs - 1);  // IZ = NZ - 1
      plane_mean += *((double *)wbmd->data + ixys); // TOP XY plane
    }
    // Summing up for all X, Z:1...NZ-1 at Y=0 and Y=NY
    for (int ixs = 0; ixs < nxs; ixs++)
    for (int izs = 1; izs < nzs-1; izs++)
    {
      int ixzs;
      ixzs = ixs * nys * nzs + izs;                    // IY = 0
      plane_mean += *((double *)wbmd->data + ixzs); // LEFT YZ plane
      ixzs = (ixs * nys + (nys - 1)) * nzs + izs;  // IY = NY - 1
      plane_mean += *((double *)wbmd->data + ixzs); // RIGHT YZ plane
    }
    // Summing up  Y:1...NY-1, Z:1...NZ-1 at X=0 and X=NX
    for (int iys = 1; iys < nys-1; iys++)
    for (int izs = 1; izs < nzs-1; izs++)
    {
      int iyzs;
      iyzs = iys * nzs + izs;                          // IX = 0
      plane_mean += *((double *)wbmd->data + iyzs); // FRONT XZ plane
      iyzs = ((nxs - 1) * nys + iys) * nzs + izs;  // IX = NX - 1
      plane_mean += *((double *)wbmd->data + iyzs); // BACK XZ plane
    }
  }
  else if (dim == 2)
  {
    nof_points = 2 * (nxs + nys - 2);
    // One need to avoid summind doubly edges
    // This is why sums do not always go from 0 to N
    for (int ixs = 0; ixs < nxs; ixs++)
    {
        plane_mean += *((double *)wbmd->data + ixs);             // Bottom X edge
        plane_mean += *((double *)wbmd->data + ixs + nxs * (nys - 1)); // Top X edge
    }
    for (int iys = 1; iys < nys-1; iys++)
    {
        plane_mean += *((double *)wbmd->data + iys * nxs + 1);       // Left Y edge
        plane_mean += *((double *)wbmd->data + (iys - 2)* nxs); // Right Y edge
    }

  }
  else if (dim == 1)
  {
    return (*((double *)wbmd->data) + *((double *)wbmd->data + wbmd->nx - 1)) / 2; // Mean value of TWO edge points
  }


  mean = plane_mean / nof_points; // Mean value of SIX planes
  return mean;
}

double complex wbox_boundary_avg_complex(wbox_md_t *wbmd)
{
  const int nxs = wbmd->nx;
  const int nys = wbmd->ny;
  const int nzs = wbmd->nz;

  const int dim = wbmd->dim;

  double complex mean = 0;
  double complex plane_mean = 0;
  int nof_points;
  if (dim == 3)
  {
    nof_points = 2 * (nxs * nys + nys * nzs + nxs * nzs) -4 * (nxs + nys + nzs) + 8;
    // One need to avoid summind doubly edges and vertices three times.
    // This is why sums do not always go from 0 to N
    // Summing up for all X,Y at Z=0 and Z=NZ
    for (int ixs = 0; ixs < nxs; ixs++)
    for (int iys = 0; iys < nys; iys++)
    {
      int ixys;
      ixys = (ixs * nys  + iys) * nzs;              // IZ = 0
      plane_mean += *((double complex*)wbmd->data + ixys); // BOTTOM XY plane
      ixys = (ixs * nys + iys) * nzs + (nzs - 1);  // IZ = NZ - 1
      plane_mean += *((double complex*)wbmd->data + ixys); // TOP XY plane
    }
    // Summing up for all X, Z:1...NZ-1 at Y=0 and Y=NY
    for (int ixs = 0; ixs < nxs; ixs++)
    for (int izs = 1; izs < nzs-1; izs++)
    {
      int ixzs;
      ixzs = ixs * nys * nzs + izs;                    // IY = 0
      plane_mean += *((double complex*)wbmd->data + ixzs); // LEFT YZ plane
      ixzs = (ixs * nys + (nys - 1)) * nzs + izs;  // IY = NY - 1
      plane_mean += *((double complex*)wbmd->data + ixzs); // RIGHT YZ plane
    }
    // Summing up  Y:1...NY-1, Z:1...NZ-1 at X=0 and X=NX
    for (int iys = 1; iys < nys-1; iys++)
    for (int izs = 1; izs < nzs-1; izs++)
    {
      int iyzs;
      iyzs = iys * nzs + izs;                          // IX = 0
      plane_mean += *((double complex*)wbmd->data + iyzs); // FRONT XZ plane
      iyzs = ((nxs - 1) * nys + iys) * nzs + izs;  // IX = NX - 1
      plane_mean += *((double complex*)wbmd->data + iyzs); // BACK XZ plane
    }
  }
  else if (dim == 2)
  {
    nof_points = 2 * (nxs + nys - 2);
    // One need to avoid summind doubly edges
    // This is why sums do not always go from 0 to N
    for (int ixs = 0; ixs < nxs; ixs++)
    {
        plane_mean += *((double complex*)wbmd->data + ixs);             // Bottom X edge
        plane_mean += *((double complex*)wbmd->data + ixs + nxs * (nys - 1)); // Top X edge
    }
    for (int iys = 1; iys < nys-1; iys++)
    {
        plane_mean += *((double complex*)wbmd->data + iys * nxs + 1);       // Left Y edge
        plane_mean += *((double complex*)wbmd->data + (iys - 2)* nxs); // Right Y edge
    }

  }
  else if (dim == 1)
  {
    return (*((double complex*)wbmd->data) + *((double complex*)wbmd->data + wbmd->nx - 1)) / 2; // Mean value of TWO edge points
  }


  mean = plane_mean / nof_points; // Mean value of SIX planes
  return mean;
}
