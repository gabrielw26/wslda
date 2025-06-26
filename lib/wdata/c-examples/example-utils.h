
#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "error: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "error: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        return -1;                                                          \
    }

/**
 * You can use this function to check diff between two arrays
 * */
void test_array_diff_df(int N, double *a, float *b)
{
    int ixyz=0;

    double d,d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ixyz=0; ixyz<N; ixyz++)
    {
        d = a[ixyz]-(double)b[ixyz];

        d2=d*d;
        sumd2+=d2;
        if(d2>maxd2) maxd2=d2;
    }

    if(sqrt(maxd2)>1.0e-9)
    {
        printf("#           |max[a-b]| : %16.8g\n", sqrt(maxd2));
        printf("#         SUM[(a-b)^2] : %16.8g\n", sumd2);
        printf("# SQRT(SUM[(a-b)^2])/N : %16.8g\n", sqrt(sumd2)/N);
    }
}

/**
 * You can use this function to check diff between two arrays
 * */
void test_array_diff_dd(int N, double *a, double *b)
{
    int ixyz=0;

    double d,d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ixyz=0; ixyz<N; ixyz++)
    {
        d = a[ixyz]-b[ixyz];

        d2=d*d;
        sumd2+=d2;
        if(d2>maxd2) maxd2=d2;
    }

    if(sqrt(maxd2)>1.0e-9)
    {
        printf("#           |max[a-b]| : %16.8g\n", sqrt(maxd2));
        printf("#         SUM[(a-b)^2] : %16.8g\n", sumd2);
        printf("# SQRT(SUM[(a-b)^2])/N : %16.8g\n", sqrt(sumd2)/N);
    }
}

/**
 * You can use this function to check diff between two arrays
 * */
void test_array_diff_ff(int N, float *a, float *b)
{
    int ixyz=0;

    double d,d2;
    double maxd2 = 0.0;
    double sumd2 = 0.0;
    for(ixyz=0; ixyz<N; ixyz++)
    {
        d = (double)a[ixyz]-(double)b[ixyz];

        d2=d*d;
        sumd2+=d2;
        if(d2>maxd2) maxd2=d2;
    }

    if(sqrt(maxd2)>1.0e-9)
    {
        printf("#           |max[a-b]| : %16.8g\n", sqrt(maxd2));
        printf("#         SUM[(a-b)^2] : %16.8g\n", sumd2);
        printf("# SQRT(SUM[(a-b)^2])/N : %16.8g\n", sqrt(sumd2)/N);
    }
}
