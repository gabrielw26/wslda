#include <complex.h> // todo: create gcc g++ compatibility. See wdata, winterp, wderiv for template
#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    int ax;
    int ay;
    int az;

    int dim;
} wbox_insert_t;

typedef struct
{
    int nx;
    int ny;
    int nz;

    int dim;

    char datatype;

    void *data;
} wbox_md_t;

// SUCCESS
#define WBOX_SUCCESS_INIT EXIT_SUCCESS

#define WBOX_SUCCESS_INSERT EXIT_SUCCESS
#define WBOX_SUCCESS_INSERT_INIT EXIT_SUCCESS
#define WBOX_SUCCESS_INSERT_CHECK EXIT_SUCCESS

#define WBOX_SUCCESS_MERGE EXIT_SUCCESS

#define WBOX_SUCCESS_ROTATE EXIT_SUCCESS

#define WBOX_SUCCESS_PRINT EXIT_SUCCESS

// ERRORS
#define WBOX_ERR_INIT_DIM 0001
#define WBOX_ERR_INIT_DATATYPE_NOT_IMPLEMENTED 0002

#define WBOX_ERR_INSERT 1000
#define WBOX_ERR_INSERT_INIT_AX_MONE 1001 // Error that is produced by `wbox_insert_init()` function and is connected to the case when `nx` parameter is equal `-1`
#define WBOX_ERR_INSERT_CHECK_DIM 1002    // Error that is produed by `wbox_insert_check()` function. Dimension of small and big box is different
#define WBOX_ERR_INSERT_1D 1003
#define WBOX_ERR_INSERT_2D 1004
#define WBOX_ERR_INSERT_3D 1005

#define WBOX_ERR_MERGE 2000

#define WBOX_ERR_ROTATE 3000

#define WBOX_ERR_PRINT_DIM 4000 // Error that is produced by `wbox_print()` function and is connected to not dimensionality of the data other than 1, 2 or 3.

/**
 * @brief Function to initialize lattice. Function initializes the lattice.
 *
 * @param wbmd [out] - pointer to WBox metadata structure. Each WBox structure is responsible for only one variable, like rho_a.
 * @param nx [in] - number of lattice points in OX direction. If `<0` then return `WBOX_ERR_INIT_DIM`
 * @param ny [in] - number of lattice points in OY direction. If `<0` then dimension is `1` and only iterates over OX axis
 * @param nz [in] - number of lattice points in OZ direction. If `<0` then dimension is `2` and iterates over OX and OY axis
 * @param datatype [in] - character describing type of data. Can be 'i'nteger, 'r'eal (double), 'c'omplex (double complex).
 * @param data [in] - pointer to a data of `datatype` type.
 * @return `WBOX_SUCCESS_INIT` or `WBOX_ERR_INIT` with help description.
 */
int wbox_init(wbox_md_t *wbmd, const int nx, const int ny, const int nz, const char datatype, const void *data);

/* ***************************************************************************************** */
/* ***************************************************************************************** */
/**
 * @brief Function that initializes insert parameters.
 *
 * @param wbox_insert [out] - `wbox_insert_t` type that stores the infromation about the anchor.
 * @param ax [in] - anchor for NX lattice value. If `-1` then ignored and `WBOX_ERR_INSERT_INIT` is returned.
 * @param ay [in] - anchor for NY lattice value. If `-1` then ignored and `wbox_insert_1d` should be executed.
 * @param az [in] - anchor for NZ lattice value. If `-1` then ignored and `wbox_insert_2d` should be executed.
 * @return `WBOX_SUCCESS_INSERT_INIT` or `WBOX_ERR_INSERT_INIT` whith help description.
 */
int wbox_insert_init(wbox_insert_t *wbox_insert, const int ax, const int ay, const int az);

/**
 * @brief Master function that inserts data from smaller box (`wmbd`) into the bigger box(`wmbd_big`) that is anchored in `wbox_insert`.
 * Function checks the anchor placements. First checks the dimension of box and places calls adequate insert function.
 *
 * @param wbmd_small [in] - `wbox_md_t` type that stores the information about data from small box.
 * @param wbox_insert [in] - `wbox_insert_t` type that stores the information about the anchor point.
 * @param wbmd_big [out] - `wbox_md_t` type that stores the information about data from big box.
 * @return `WBOX_SUCCESS_INSERT` or `WBOX_ERR_INSERT` with help description.
 */
int wbox_insert(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big);


/**
 * @brief Function that calculates average over the boundaries of the box (`wmbd`).
 * It works for 1D, 2D, and 3D.
 *
 * @param wbmd [in] - `wbox_md_t` type that stores the information about data from box.
 * @return the average value of integer data
 */
int wbox_boundary_avg_int(wbox_md_t *wbmd);

/**
 * @brief Function that calculates average over the boundaries of the box (`wmbd`).
 * It works for 1D, 2D, and 3D.
 *
 * @param wbmd [in] - `wbox_md_t` type that stores the information about data from box.
 * @return the average value of double data
 */
double wbox_boundary_avg_double(wbox_md_t *wbmd);

/**
 * @brief Function that calculates average over the boundaries of the box (`wmbd`).
 * It works for 1D, 2D, and 3D.
 *
 * @param wbmd [in] - `wbox_md_t` type that stores the information about data from box.
 * @return the average value of complex data
 */
double complex wbox_boundary_avg_complex(wbox_md_t *wbmd);


/**
 * @brief Function that inserts data from smaller 1D box (`wmbd`) into the 1D bigger box(`wmbd_big`) that is anchored in `wbox_insert`.
 * Function checks the anchor placements.
 *
 * @param wbmd_small [in] - `wbox_md_t` type that stores the information about data from small box.
 * @param wbox_insert [in] - `wbox_insert_t` type that stores the information about the anchor point.
 * @param wbmd_big [out] - `wbox_md_t` type that stores the information about data from big box.
 * @return `WBOX_SUCCESS_INSERT` or `WBOX_ERR_INSERT` with help description.
 */
int wbox_insert_1d(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big);

/**
 * @brief Function that inserts data from smaller 2D box (`wmbd`) into the 2D bigger box(`wmbd_big`) that is anchored in `wbox_insert`.
 * Function checks the anchor placements.
 *
 * @param wbmd_small [in] - `wbox_md_t` type that stores the information about data from small box.
 * @param wbox_insert [in] - `wbox_insert_t` type that stores the information about the anchor point.
 * @param wbmd_big [out] - `wbox_md_t` type that stores the information about data from big box.
 * @return `WBOX_SUCCESS_INSERT` or `WBOX_ERR_INSERT` with help description.
 */
int wbox_insert_2d(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big);

/**
 * @brief Function that inserts data from smaller 3D box (`wmbd`) into the 3D bigger box(`wmbd_big`) that is anchored in `wbox_insert`.
 * Function checks the anchor placements.
 *
 * @param wbmd_small [in] - `wbox_md_t` type that stores the information about data from small box.
 * @param wbox_insert [in] - `wbox_insert_t` type that stores the information about the anchor point.
 * @param wbmd_big [out] - `wbox_md_t` type that stores the information about data from big box.
 * @return `WBOX_SUCCESS_INSERT` or `WBOX_ERR_INSERT` with help description.
 */
int wbox_insert_3d(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, wbox_md_t *wbmd_big);

/**
 * @brief Function that checks the possible errors for inserting a data from small box into big box.
 *
 * @param wbmd_small [in] - `wbox_md_t` type that stores the information about data from small box.
 * @param wbox_insert [in] - `wbox_insert_t` type that stores the information about the anchor point.
 * @param wbmd_big [out] - `wbox_md_t` type that stores the information about data from big box.
 * @return `WBOX_SUCCESS_INSERT_CHECK` or `WBOX_ERR_INSERT_CHECK` with help description.
 */
int wbox_insert_check(const wbox_md_t *wbmd_small, const wbox_insert_t *wbox_insert, const wbox_md_t *wbmd_big);

/* ***************************************************************************************** */
/* ***************************************************************************************** */
int wbox_merge(const wbox_md_t *wbmd_one, const wbox_md_t *wbmd_two, wbox_md_t *wbmd_result);
int wbox_merge_1d(const wbox_md_t *wbmd_one, const wbox_md_t *wbmd_two, wbox_md_t *wbmd_result);
int wbox_merge_2d(const wbox_md_t *wbmd_one, const wbox_md_t *wbmd_two, wbox_md_t *wbmd_result);
int wbox_merge_3d(const wbox_md_t *wbmd_one, const wbox_md_t *wbmd_two, wbox_md_t *wbmd_result);

/* ***************************************************************************************** */
/* ***************************************************************************************** */
int wbox_rotate_2d(wbox_md_t *wbmd);
int wbox_rotate_3d(wbox_md_t *wbmd);
int wbox_rotate(wbox_md_t *wbmd);

/* ***************************************************************************************** */
/* ***************************************************************************************** */
/**
 * @brief Master function that prints data in terminal
 *
 * @param wbmd [in] - pointer to data that one wants to be printed
 * @param title [in] - char pointer with a title, eg. "WBOX SMALL"
 * @return `WBOX_SUCCESS_PRINT` or `WBOX_ERR_PRINT` with a help description
 */
int wbox_print(const wbox_md_t *wbmd, const char *title);

/**
 * @brief Function that prints data in terminal for 1D data
 *
 * @param wbmd [in] - pointer to data that one wants to be printed
 * @param title [in] - char pointer with a title, eg. "WBOX SMALL"
 * @return `WBOX_SUCCESS_PRINT` or `WBOX_ERR_PRINT` with a help description
 */
int wbox_print_1d(const wbox_md_t *wbmd, const char *title);

/**
 * @brief Function that prints data in terminal for 2D data
 *
 * @param wbmd [in] - pointer to data that one wants to be printed
 * @param title [in] - char pointer with a title, eg. "WBOX SMALL"
 * @return `WBOX_SUCCESS_PRINT` or `WBOX_ERR_PRINT` with a help description
 */
int wbox_print_2d(const wbox_md_t *wbmd, const char *title);

/**
 * @brief Function that prints data in terminal for 3D data
 *
 * @param wbmd [in] - pointer to data that one wants to be printed
 * @param title [in] - char pointer with a title, eg. "WBOX SMALL"
 * @return `WBOX_SUCCESS_PRINT` or `WBOX_ERR_PRINT` with a help description
 */
int wbox_print_3d(const wbox_md_t *wbmd, const char *title);

/**
 * @brief Function that prints `wbox_md_t` struct in terminal
 *
 * @param wbmd [in] - pointer to `wbox_md_t` that one wants to be printed
 * @param title [in] - char pointer with a title, eg. "WBOX SMALL"
 * @return `WBOX_SUCCESS_PRINT` or `WBOX_ERR_PRINT` with a help description
 */
int wbox_print_box_params(const wbox_md_t wbmd, const char *title);

/**
 * @brief Function that prints `wbox_insert_t` struct in terminal
 *
 * @param wb_insert [in] - pointer to `wbox_insert_t` struct that one wants to be printed
 * @return `WBOX_SUCCESS_PRINT` or `WBOX_ERR_PRINT` with a help description
 */
int wbox_print_insert_params(const wbox_insert_t wb_insert);
