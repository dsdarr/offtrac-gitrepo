/********+*********+*********+*********+*********+*********+*********+*
 *    This header file defines the metric arrays for HIM to be of     *
 *  appropriate sizes.                                                *
 ********+*********+*********+*********+*********+*********+*********+*/

#ifdef MAINFILE
#define DOUBLE double
#else
#define DOUBLE extern double
#endif

#ifndef CARTESIAN
# if defined(XMETRIC_J) && defined(XMETRIC_I)
#  define XMETRIC_I_J
# endif 
# if defined(YMETRIC_J) && defined(YMETRIC_I)
#  define YMETRIC_I_J
# endif 

# if defined(XMETRIC_I_J) || defined(YMETRIC_I_J) || \
    (defined(XMETRIC_J) && defined(YMETRIC_I)) || \
    (defined(XMETRIC_I) && defined(YMETRIC_J))
#  define XYMETRICS_I_J
#  else
#   if defined(XMETRIC_J) || defined(YMETRIC_J)
#    define XYMETRICS_J
#   else
#    if defined(XMETRIC_I) || defined(YMETRIC_I)
#     define XYMETRICS_I
#    else
#     define CARTESIAN
#   endif /* XMETRIC_I or YMETRIC_I */
#  endif /* XMETRIC_J YMETRIC_J */
# endif /* XMETRIC_I_J or YMETRIC_I_J */
#endif /* CARTESIAN */

#ifdef CARTESIAN
  DOUBLE dx, dy, Idx, Idy, dxdy, Idxdy;

# define DXh(i,j) dx
# define IDXh(i,j) Idx
# define DXq(i,j) dx
# define IDXq(i,j) Idx
# define DXu(i,j) dx
# define IDXu(i,j) Idx
# define DXv(i,j) dx
# define IDXv(i,j) Idx

# define DYh(i,j) dy
# define IDYh(i,j) Idy
# define DYq(i,j) dy
# define IDYq(i,j) Idy
# define DYu(i,j) dy
# define IDYu(i,j) Idy
# define DYv(i,j) dy
# define IDYv(i,j) Idy

# define DXDYh(i,j) dxdy
# define IDXDYh(i,j) Idxdy
# define DXDYq(i,j) dxdy
# define IDXDYq(i,j) Idxdy

#else

# ifdef XMETRIC_I_J
   DOUBLE dxh[NXMEM][NYMEM], Idxh[NXMEM][NYMEM]; /* dxp is the metric term in the */
   DOUBLE dxq[NXMEM][NYMEM], Idxq[NXMEM][NYMEM]; /* x-direction at p points, in m.*/
   DOUBLE dxu[NXMEM][NYMEM], Idxu[NXMEM][NYMEM]; /* Idxp is the inverse of dxp,   */
   DOUBLE dxv[NXMEM][NYMEM], Idxv[NXMEM][NYMEM]; /* in m-1.                       */
#  define DXh(i,j) dxh[i][j]
#  define IDXh(i,j) Idxh[i][j]
#  define DXq(i,j) dxq[i][j]
#  define IDXq(i,j) Idxq[i][j]
#  define DXu(i,j) dxu[i][j]
#  define IDXu(i,j) Idxu[i][j]
#  define DXv(i,j) dxv[i][j]
#  define IDXv(i,j) Idxv[i][j]
# else
#  ifdef XMETRIC_J
    DOUBLE dxh[NYMEM], Idxh[NYMEM]; /* dxp is the metric term in the */
    DOUBLE dxq[NYMEM], Idxq[NYMEM]; /* x-direction at p points, in m.*/
                                    /* Idxp is the inverse of dxp,   */
                                    /* in m-1.                       */
#    define DXh(i,j) dxh[j]
#    define IDXh(i,j) Idxh[j]
#    define DXq(i,j) dxq[j]
#    define IDXq(i,j) Idxq[j]
#    define DXu(i,j) dxh[j]
#    define IDXu(i,j) Idxh[j]
#    define DXv(i,j) dxq[j]
#    define IDXv(i,j) Idxq[j]
#  else
#   ifdef XMETRIC_I
     DOUBLE dxh[NXMEM], Idxh[NXMEM]; /* dxp is the metric term in the */
     DOUBLE dxq[NXMEM], Idxq[NXMEM]; /* x-direction at p points, in m.*/
                                     /* Idxp is the inverse of dxp,   */
i                                    /* in m-1.                       */
#    define DXh(i,j) dxh[i]
#    define IDXh(i,j) Idxh[i]
#    define DXq(i,j) dxq[i]
#    define IDXq(i,j) Idxq[i]
#    define DXu(i,j) dxq[i]
#    define IDXu(i,j) Idxq[i]
#    define DXv(i,j) dxh[i]
#    define IDXv(i,j) Idxh[i]
#   else
     DOUBLE dx, Idx;         /* dx is the metric term in the  */
                             /* x-direction, in m. Idx is the */
                             /* inverse of dx, in m-1.        */
#    define DXh(i,j) dx
#    define IDXh(i,j) Idx
#    define DXq(i,j) dx
#    define IDXq(i,j) Idx
#    define DXu(i,j) dx
#    define IDXu(i,j) Idx
#    define DXv(i,j) dx
#    define IDXv(i,j) Idx
#   endif /* XMETRIC_I */
#  endif /* XMETRIC_J */
# endif /* XMETRIC_I_J */

# ifdef YMETRIC_I_J
   DOUBLE dyh[NXMEM][NYMEM], Idyh[NXMEM][NYMEM]; /* dyp is the metric term in the */
   DOUBLE dyq[NXMEM][NYMEM], Idyq[NXMEM][NYMEM]; /* y-direction at p points, in m.*/
   DOUBLE dyu[NXMEM][NYMEM], Idyu[NXMEM][NYMEM]; /* Idyp is the inverse of dyp,   */
   DOUBLE dyv[NXMEM][NYMEM], Idyv[NXMEM][NYMEM]; /* in m-1.                       */
#  define DYh(i,j) dyh[i][j]
#  define IDYh(i,j) Idyh[i][j]
#  define DYq(i,j) dyq[i][j]
#  define IDYq(i,j) Idyq[i][j]
#  define DYu(i,j) dyu[i][j]
#  define IDYu(i,j) Idyu[i][j]
#  define DYv(i,j) dyv[i][j]
#  define IDYv(i,j) Idyv[i][j]
# else
#  ifdef YMETRIC_J
    DOUBLE dyh[NYMEM], Idyh[NYMEM]; /* dyp is the metric term in the */
    DOUBLE dyq[NYMEM], Idyq[NYMEM]; /* y-direction at p points, in m.*/
                                    /* Idyp is the inverse of dyp,   */
                                    /* in m-1.                       */
#   define DYh(i,j) dyh[j]
#   define IDYh(i,j) Idyh[j]
#   define DYq(i,j) dyq[j]
#   define IDYq(i,j) Idyq[j]
#   define DYu(i,j) dyh[j]
#   define IDYu(i,j) Idyh[j]
#   define DYv(i,j) dyq[j]
#   define IDYv(i,j) Idyq[j]
#  else
#   ifdef YMETRIC_I
     DOUBLE dyh[NXMEM], Idyh[NXMEM]; /* dyp is the metric term in the */
     DOUBLE dyq[NXMEM], Idyq[NXMEM]; /* y-direction at p points, in m.*/
                                   /* Idyp is the inverse of dyp,   */
                                   /* in m-1.                       */
#    define DYh(i,j) dyh[i]
#    define IDYh(i,j) Idyh[i]
#    define DYq(i,j) dyq[i]
#    define IDYq(i,j) Idyq[i]
#    define DYu(i,j) dyq[i]
#    define IDYu(i,j) Idyq[i]
#    define DYv(i,j) dyh[i]
#    define IDYv(i,j) Idyh[i]
#   else
     DOUBLE dy, Idy;              /* dy is the metric term in the  */
                                  /* y-direction, in m. Idy is the */
                                  /* inverse of dy, in m-1.        */
#    define DYh(i,j) dy
#    define IDYh(i,j) Idy
#    define DYq(i,j) dy
#    define IDYq(i,j) Idy
#    define DYu(i,j) dy
#    define IDYu(i,j) Idy
#    define DYv(i,j) dy
#    define IDYv(i,j) Idy
#   endif /* YMETRIC_I */
#  endif /* YMETRIC_J */
# endif /* YMETRIC_I_J */


# ifdef XYMETRICS_I_J
   DOUBLE dxdyh[NXMEM][NYMEM];     /* dxdyp is the area of a cell   */
   DOUBLE Idxdyh[NXMEM][NYMEM];    /* centered at a p point, in m2. */
   DOUBLE dxdyq[NXMEM][NYMEM];     /* Idxdyp is the inverse of      */
   DOUBLE Idxdyq[NXMEM][NYMEM];    /* dxdyp, in m-2.                */
#  define DXDYh(i,j) dxdyh[i][j]
#  define IDXDYh(i,j) Idxdyh[i][j]
#  define DXDYq(i,j) dxdyq[i][j]
#  define IDXDYq(i,j) Idxdyq[i][j]
# else
#  ifdef XYMETRICS_J
    DOUBLE dxdyh[NYMEM];           /* dxdyp is the area of a cell   */
    DOUBLE Idxdyh[NYMEM];          /* centered at a p point, in m2. */
    DOUBLE dxdyq[NYMEM];           /* Idxdyp is the inverse of      */
    DOUBLE Idxdyq[NYMEM];          /* dxdyp, in m-2.                */
#   define DXDYh(i,j) dxdyh[j]
#   define IDXDYh(i,j) Idxdyh[j]
#   define DXDYq(i,j) dxdyq[j]
#   define IDXDYq(i,j) Idxdyq[j]
#  else
#   ifdef XYMETRICS_I
     DOUBLE dxdyh[NXMEM];          /* dxdyp is the area of a cell   */
     DOUBLE Idxdyh[NXMEM];         /* centered at a p point, in m2. */
     DOUBLE dxdyq[NXMEM];          /* Idxdyp is the inverse of      */
     DOUBLE Idxdyq[NXMEM];         /* dxdyp, in m-2.                */
#    define DXDYh(i,j) dxdyh[i]
#    define IDXDYh(i,j) Idxdyh[i]
#    define DXDYq(i,j) dxdyq[i]
#    define IDXDYq(i,j) Idxdyq[i]
#   else
     DOUBLE dxdy, Idxdy;
#    define DXDYh(i,j) dxdy
#    define IDXDYh(i,j) Idxdy
#    define DXDYq(i,j) dxdy
#    define IDXDYq(i,j) Idxdy
#   endif /* XYMETRICS_I */
#  endif /* XYMETRICS_J  */
# endif /* XYMETRICS_I_J */
#endif /* CARTESIAN */

#undef DOUBLE
