/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 Rhys Ulerich
 * Copyright (C) 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 * bspline.c: higher-level logic build atop the GSL B-spline routines
 * $Id$
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline.h>
#include <suzerain/error.h>
#include <suzerain/pre_gsl.h>

int
suzerain_bspline_linear_combination(
    const size_t nderiv,
    const double * coeffs,
    const size_t npoints,
    const double * points,
    double * values,
    const size_t ldvalues,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Parameter sanity checks */
    if (ldvalues && ldvalues < npoints) {
        SUZERAIN_ERROR("ldvalues too small for npoints", SUZERAIN_EINVAL);
    }

    /* ldvalues == 0 signals that we only want derivative nderiv */
    const size_t jstart = (ldvalues == 0) ? nderiv : 0;

    size_t istart, iend;
    for (size_t i = 0; i < npoints; ++i) {

        gsl_bspline_deriv_eval_nonzero(points[i], nderiv,
                dB, &istart, &iend, w, dw);

        const double * coeff_start = coeffs + istart;

        for (size_t j = jstart; j <= nderiv; ++j) {
            double dot = 0.0;
            for (size_t k = 0; k < w->k; ++k) {
                dot += coeff_start[k] * (dB->data + j)[k * dB->tda];
            }
            values[i + j * ldvalues] = dot;
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_linear_combination_complex(
    const size_t nderiv,
    const complex_double *coeffs,
    const size_t npoints,
    const double * points,
    complex_double *values,
    const size_t ldvalues,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Parameter sanity checks */
    if (ldvalues && ldvalues < npoints) {
        SUZERAIN_ERROR("ldvalues too small for npoints", SUZERAIN_EINVAL);
    }

    /* ldvalues == 0 signals that we only want derivative nderiv */
    const size_t jstart = (ldvalues == 0) ? nderiv : 0;

    size_t istart, iend;
    for (size_t i = 0; i < npoints; ++i) {

        gsl_bspline_deriv_eval_nonzero(points[i], nderiv,
                dB, &istart, &iend, w, dw);

        const complex_double * coeff_start = coeffs + istart;

        for (size_t j = jstart; j <= nderiv; ++j) {
            double dotr = 0.0, doti = 0.0;
            for (size_t k = 0; k < w->k; ++k) {
                const double basis_value = (dB->data + j)[k * dB->tda];
                dotr += creal(coeff_start[k]) * basis_value;
                doti += cimag(coeff_start[k]) * basis_value;
            }
            values[i + j*ldvalues] = complex_double(dotr, doti);
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_integration_coefficients(
    const size_t nderiv,
    double * coeffs,
    size_t inc,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Obtain an appropriate order Gauss-Legendre integration rule */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc((w->k - nderiv + 1)/2);
    if (tbl == NULL) {
        SUZERAIN_ERROR("failed to obtain Gauss-Legendre rule from GSL",
                       SUZERAIN_ESANITY);
    }

    /* Zero integration coefficient values */
    for (size_t i = 0; i < w->n; ++i) {
        coeffs[i * inc] = 0.0;
    }

    /* Accumulate the breakpoint-by-breakpoint contributions to coeffs */
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval */
        const double a = gsl_bspline_breakpoint(i,   w);
        const double b = gsl_bspline_breakpoint(i+1, w);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_deriv_eval_nonzero(xj, nderiv,
                    dB, &kstart, &kend, w, dw);

            /* Accumulate weighted basis evaluations into coeffs */
            for (size_t k = kstart; k <= kend; ++k) {
                coeffs[k * inc] += wj * gsl_matrix_get(dB,
                                                       k - kstart, nderiv);
            }
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    return SUZERAIN_SUCCESS;
}

double
suzerain_bspline_spacing_greville_abscissae(
    size_t i,
    gsl_bspline_workspace *w)
{
    /* Find nearest abscissae indices silently folding back into range */
    size_t im = (i == 0       ) ? 1        : i - 1;
    size_t ip = (i == w->n - 1) ? w->n - 2 : i + 1;

    /* Compute and return the minimum distance */
    double x  = gsl_bspline_greville_abscissa(i,  w);
    double xm = gsl_bspline_greville_abscissa(im, w);
    double xp = gsl_bspline_greville_abscissa(ip, w);

    double dxm = fabs(x  - xm);
    double dxp = fabs(xp - x);
    return dxm < dxp ? dxm : dxp;
}

double
suzerain_bspline_spacing_breakpoints(
    size_t i,
    gsl_bspline_workspace *w)
{
    // Make indices {first,...,last} span the required portion of the knots
    size_t stride  = w->knots->stride;
    size_t k       = w->k;
    double * first = w->knots->data + (i  )*stride;
    double * last  = w->knots->data + (i+k)*stride;

    // Determine the minimum non-zero spacing between adjacent knots
    double retval = GSL_DBL_MAX;
    while (first != last) {

        double left  = *first;
        first       += stride;
        double right = *first;

        double dist = right - left;
        if (dist > 0 && dist < retval) {
            retval = dist;
        }
    }
    return retval;
}

int
suzerain_bspline_htstretch2_evdeltascale(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh)
{
    /* Parameter sanity checks */
    if (nderiv < 0) {
        SUZERAIN_ERROR("Requested derivative nderiv must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (k < 1) {
        SUZERAIN_ERROR("B-spline order k must be strictly positive",
                       SUZERAIN_EINVAL);
    } else if (htdelta < 0) {
        SUZERAIN_ERROR("Stretching parameter htdelta must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (N < k) {
        SUZERAIN_ERROR("Number of degrees of freedom N < B-spline order k",
                       SUZERAIN_EINVAL);
    }

    /* Coefficient lookup table for orders 4 through 17 from model document */
    /* The final two columns are the low and higher estimated error bounds */
#define NDERIVLOW  ( 1)
#define NDERIVHIGH ( 2)
#define KLOW       ( 4)
#define KHIGH      (17)
#define XSTR(s) STR(s)
#define STR(s)  #s
    static const double p[NDERIVHIGH - NDERIVLOW + 1][KHIGH - KLOW + 1][10] = {
    /* nderiv = 1 */ {
        /* 4*/  {  97891.     /  48235.   ,  -1901.  /  11697.  ,  -130995.    /  21709.  ,  -30229.   /  24081.  ,  5619.   /  10310.  ,  1620498.  /  27289.  ,  142.    /  5545.   ,  2601.    /  7318.   ,  -0.849  ,  0.918  },
        /* 5*/  {  -110189.   /  4436.    ,  34.     /  12241.  ,  -4172419.   /  12740.  ,  27053.    /  7741.   ,  1.      /  4663.   ,  547.      /  9693.   ,  1697.   /  6499.   ,  -1920.   /  7939.   ,  -1.91   ,  0.456  },
        /* 6*/  {  -821.      /  6150.    ,  -2903.  /  20696.  ,  101950.     /  5167.   ,  -17849.   /  17027.  ,  3099.   /  6728.   ,  646487.   /  9039.   ,  389.    /  10287.  ,  629.     /  1944.   ,  -1.42   ,  1.54   },
        /* 7*/  {  37007.     /  18425.   ,  -1064.  /  7817.   ,  -354668.    /  23235.  ,  -11655.   /  12098.  ,  1930.   /  4347.   ,  334349.   /  4917.   ,  589.    /  14681.  ,  1354.    /  4403.   ,  -1.71   ,  2.32   },
        /* 8*/  {  2473.      /  823.     ,  -5295.  /  35063.  ,  -247406.    /  5933.   ,  -4605.    /  5438.   ,  1122.   /  2545.   ,  464735.   /  9127.   ,  85.     /  2268.   ,  2385.    /  7247.   ,  -2.26   ,  3.26   },
        /* 9*/  {  9617.      /  9991.    ,  -631.   /  7734.   ,  31763.      /  6158.   ,  -7387.    /  8084.   ,  2022.   /  4423.   ,  891260.   /  16093.  ,  256.    /  5423.   ,  2593.    /  9197.   ,  -2.71   ,  5.18   },
        /*10*/  {  -33713.    /  7386.    ,  -761.   /  5137.   ,  1926089.    /  10690.  ,  -6921.    /  9599.   ,  2317.   /  5160.   ,  291351.   /  8860.   ,  265.    /  6544.   ,  4608.    /  13657.  ,  -5.56   ,  6.78   },
        /*11*/  {  3079.      /  6057.    ,  -1434.  /  8125.   ,  61849.      /  2773.   ,  -3043.    /  4610.   ,  1640.   /  3927.   ,  390201.   /  10112.  ,  719.    /  21422.  ,  3256.    /  9915.   ,  -3.04   ,  5.68   },
        /*12*/  {  -1133.     /  7526.    ,  41.     /  1685.   ,  324882.     /  6137.   ,  -5771.    /  6511.   ,  3772.   /  8097.   ,  2982597.  /  56293.  ,  301.    /  6138.   ,  1820.    /  7983.   ,  -8.17   ,  8.79   },
        /*13*/  {  76489.     /  140702.  ,  -1306.  /  12921.  ,  394799.     /  14952.  ,  -5967.    /  8981.   ,  10198.  /  23773.  ,  3000311.  /  82205.  ,  763.    /  20232.  ,  5967.    /  20437.  ,  -3.93   ,  7.11   },
        /*14*/  {  1889.      /  2990.    ,  -28.    /  2747.   ,  106241.     /  4297.   ,  -5824.    /  8067.   ,  8194.   /  18519.  ,  272113.   /  6793.   ,  517.    /  12610.  ,  926.     /  3631.   ,  -5.56   ,  8.43   },
        /*15*/  {  12728.     /  8271.    ,  325.    /  6967.   ,  -338210.    /  13231.  ,  -3283.    /  4378.   ,  739.    /  1661.   ,  401689.   /  8928.   ,  627.    /  15115.  ,  2647.    /  11721.  ,  -7.34   ,  7.78   },
        /*16*/  {  -871534.   /  4585.    ,  23089.  /  8570.   ,  70672227.   /  5777.   ,  -137793.  /  11641.  ,  24788.  /  4143.   ,  102179.   /  30235.  ,  -4403.  /  6191.   ,  109963.  /  9789.   ,  -56.6   ,  53.8   },
        /*17*/  {  12184154.  /  189207.  ,  15079.  /  5376.   ,  -42085883.  /  9480.   ,  -57193.   /  6761.   ,  58893.  /  13753.  ,  50031.    /  12802.  ,  -8767.  /  16502.  ,  45583.   /  5015.   ,  -62.7   ,  47.4   },
    },
    /* nderiv = 2 */ {
        /* 4*/  {  -91411.   /  5050.   ,  -1.     /  2802.   ,  -1172302.  /  8721.   ,  36245.   /  10134.  ,  1.     /  5977.   ,  80.       /  5163.   ,  2087.   /  6715.   ,  -1052.  /  4629.   ,  -2.62   ,  1.14   },
        /* 5*/  {  -75193.   /  13574.  ,  -5.     /  5733.   ,  -125327.   /  1527.   ,  50159.   /  19102.  ,  1.     /  4713.   ,  279.      /  4237.   ,  9311.   /  34379.  ,  -1694.  /  6683.   ,  -0.691  ,  0.183  },
        /* 6*/  {  2225.     /  2291.   ,  -414.   /  13573.  ,  -10958.    /  4867.   ,  -23049.  /  11897.  ,  7061.  /  12654.  ,  543407.   /  10936.  ,  250.    /  19167.  ,  2443.   /  4449.   ,  -1.10   ,  3.90   },
        /* 7*/  {  34701.    /  20338.  ,  -371.   /  23883.  ,  -42818.    /  2689.   ,  -76859.  /  46329.  ,  6122.  /  10999.  ,  639841.   /  19837.  ,  -503.   /  9629.   ,  8051.   /  8750.   ,  -3.76   ,  13.4   },
        /* 8*/  {  4825.     /  1948.   ,  98.     /  23433.  ,  -325375.   /  8918.   ,  -12865.  /  8157.   ,  3097.  /  5423.   ,  249574.   /  8269.   ,  -1328.  /  22645.  ,  2256.   /  2155.   ,  -6.05   ,  12.9   },
        /* 9*/  {  -25910.   /  10437.  ,  -31.    /  35091.  ,  4141690.   /  45639.  ,  -8064.   /  5755.   ,  7497.  /  14038.  ,  96422.    /  3093.   ,  -382.   /  6359.   ,  8728.   /  8061.   ,  -9.99   ,  10.9   },
        /*10*/  {  3894.     /  6067.   ,  -211.   /  6112.   ,  24043.     /  3050.   ,  -4266.   /  3403.   ,  1052.  /  2499.   ,  2425053.  /  43118.  ,  -248.   /  8379.   ,  6251.   /  7412.   ,  -8.61   ,  19.0   },
        /*11*/  {  4689.     /  5320.   ,  -329.   /  4061.   ,  2688.      /  4481.   ,  -3256.   /  2927.   ,  3703.  /  7393.   ,  230007.   /  7703.   ,  -687.   /  18191.  ,  2706.   /  2707.   ,  -19.7   ,  17.8   },
        /*12*/  {  -413593.  /  16678.  ,  -213.   /  8032.   ,  5544667.   /  5189.   ,  -15275.  /  14748.  ,  6463.  /  19007.  ,  1467571.  /  22503.  ,  -203.   /  7771.   ,  12836.  /  14889.  ,  -20.2   ,  22.6   },
        /*13*/  {  39391.    /  13017.  ,  -2737.  /  20572.  ,  -664032.   /  6695.   ,  -3657.   /  4114.   ,  2765.  /  6824.   ,  1757113.  /  43536.  ,  -172.   /  6179.   ,  10441.  /  11628.  ,  -19.0   ,  26.4   },
        /*14*/  {  3112.     /  1799.   ,  -311.   /  1922.   ,  -149759.   /  3508.   ,  -6836.   /  8427.   ,  3394.  /  9043.   ,  260436.   /  5803.   ,  -205.   /  9414.   ,  2476.   /  2931.   ,  -21.2   ,  31.1   },
        /*15*/  {  35927.    /  12573.  ,  -2323.  /  16156.  ,  -731236.   /  6493.   ,  -2766.   /  3751.   ,  1699.  /  5232.   ,  816098.   /  15589.  ,  -77.    /  3170.   ,  4910.   /  5627.   ,  -25.5   ,  30.4   },
        /*16*/  {  -2839.    /  7152.   ,  -6234.  /  22253.  ,  222397.    /  2640.   ,  -2855.   /  4641.   ,  4229.  /  12834.  ,  351357.   /  8068.   ,  -194.   /  7885.   ,  2265.   /  2578.   ,  -26.6   ,  29.6   },
        /*17*/  {  54963.    /  11992.  ,  -2959.  /  11710.  ,  -806251.   /  3144.   ,  -3173.   /  5096.   ,  9227.  /  24975.  ,  992267.   /  28762.  ,  -487.   /  16865.  ,  4971.   /  5087.   ,  -30.3   ,  29.1   },
    }
    };

    /* Where are we in the coefficient table? */
    const int ndx = nderiv - NDERIVLOW;
    if (ndx < 0 || (size_t) ndx >= sizeof(p)/sizeof(p[0])) {
        SUZERAIN_ERROR("Derivative nderiv outside of implemented range: ["
                XSTR(NDERIVLOW) "," XSTR(NDERIVHIGH) "]", SUZERAIN_EINVAL);
    }
    const int kdx = k - KLOW;
    if (kdx < 0 || (size_t) kdx >= sizeof(p[0])/sizeof(p[0][0])) {
        SUZERAIN_ERROR("B-spline order k outside of implemented range: ["
                XSTR(KLOW) "," XSTR(KHIGH) "]", SUZERAIN_EINVAL);
    }

#undef STR
#undef XSTR
#undef KHIGH
#undef KLOW
#undef NDERIVHIGH
#undef NDERIVLOW

    /* Now that we know where we are, assign symbolic names to bits of data */
    const double a    = p[ndx][kdx][0];
    const double b    = p[ndx][kdx][1];
    const double c    = p[ndx][kdx][2];
    const double d    = p[ndx][kdx][3];
    const double e    = p[ndx][kdx][4];
    const double f    = p[ndx][kdx][5];
    const double g    = p[ndx][kdx][6];
    const double h    = p[ndx][kdx][7];
    const double errl = p[ndx][kdx][8]; // Relative error measured in percent
    const double errh = p[ndx][kdx][9]; // Relative error measured in percent

    /* Compute the empirical fit being careful about integer k, N */
    const double deltahat = tanh(htdelta)*sqrt(1. + htdelta);
    *C  = pow((double)k, d + e * deltahat);
    *C *= 1.0 + pow((f * (double)k) / ((double)N - (double)k + 1.0),
                    g * (double)k + h * deltahat);
    *C += a * (double)k;
    *C += b * deltahat;
    *C += c / sqrt((double)k);

    // Using percent error bounds on the fit, estimate true range
    if (Clow)  *Clow  = *C / (1 - errl/100);
    if (Chigh) *Chigh = *C / (1 - errh/100);

    return SUZERAIN_SUCCESS;
}
