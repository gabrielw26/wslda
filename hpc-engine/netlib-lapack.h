 /**
 * WSLDA Toolkit
 * Author: Gabriel Wlazlowski
 * 
 * File created to get code compatible with some systems like Summit
 */

#ifndef __WSLDA_NETLIBSCALAPACK__
#define __WSLDA_NETLIBSCALAPACK__

#ifdef FOTRAN_NO_UNDERSCORE

#define numroc_ numroc
#define descinit_ descinit
#define pzheevr_ pzheevr
#define pdsyevr_ pdsyevr
#define pzheevd_ pzheevd
#define pdsyevd_ pdsyevd
#define pzheev_ pzheev
#define pdsyev_ pdsyev
#define indxl2g_ indxl2g
#define dgetri_ dgetri
#define dgetrf_ dgetrf

#endif

#endif

