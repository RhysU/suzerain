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
    static const double p[NDERIVHIGH - NDERIVLOW + 1][KHIGH - KLOW + 1][8] = {
    /* nderiv = 1 */ {
        /*  4 */  {  28199. / 18559. , 3927.  / 2984.  , 279.   / 12247. , 1.    / 176331. , 4065883270. / 15549. , 14177. / 11138. , -2.67  , 1.41  },
        /*  5 */  {  7757.  / 3927.  , 20542. / 17101. , 230.   / 9363.  , 1.    / 265878. , 2128419709. / 6255.  , 5259.  / 3941.  , -2.72  , 1.23  },
        /*  6 */  {  10865. / 5099.  , 1091.  / 904.   , 164.   / 6807.  , 1.    / 453003. , 5432575189. / 10284. , 21811. / 15439. , -2.80  , 2.30  },
        /*  7 */  {  13884. / 6317.  , 14388. / 11713. , 608.   / 29097. , 1.    / 255011. , 3769915654. / 13467. , 7944.  / 5285.  , -2.81  , 5.94  },
        /*  8 */  {  12229. / 5291.  , 4835.  / 3924.  , 75.    / 4861.  , 1.    / 58399.  , 795198202.  / 12991. , 9161.  / 5777.  , -4.19  , 9.87  },
        /*  9 */  {  25027. / 10442. , 15387. / 12422. , 29.    / 6615.  , 1.    / 35517.  , 392916968.  / 11105. , 13529. / 7975.  , -6.81  , 15.9  },
        /* 10 */  {  71736. / 26099. , 658.   / 549.   , -132.  / 11383. , 1.    / 10349.  , 55467055.   / 5618.  , 15653. / 8658.  , -10.0  , 20.5  },
        /* 11 */  {  8894.  / 3301.  , 3852.  / 3143.  , -501.  / 15115. , 1.    / 3321.   , 5203926.    / 1723.  , 79529. / 41257. , -13.91 , 25.8  },
        /* 12 */  {  36035. / 11981. , 2123.  / 1773.  , -343.  / 5598.  , 3.    / 9793.   , 11432778.   / 4049.  , 20121. / 9800.  , -18.7  , 29.6  },
        /* 13 */  {  15247. / 6749.  , 9584.  / 7237.  , -239.  / 2591.  , 3.    / 8105.   , 12399465.   / 5524.  , 31808. / 14733. , -23.3  , 31.0  },
        /* 14 */  {  18863. / 8900.  , 7703.  / 5665.  , -3485. / 26908. , 5.    / 11174.  , 15662698.   / 9109.  , 21891. / 9520.  , -28.7  , 34.2  },
        /* 15 */  {  32932. / 12113. , 18604. / 14609. , -1466. / 10879. , 23.   / 37123.  , 17525538.   / 14141. , 19083. / 8363.  , -27.0  , 25.0  },
        /* 16 */  {  7881.  / 2593.  , 7488.  / 6037.  , -1943. / 11766. , 451.  / 4271.   , 52824.      / 7273.  , 77819. / 33448. , -30.1  , 25.6  },
        /* 17 */  {  10021. / 4643.  , 4088.  / 2991.  , -1824. / 9431.  , 6721. / 22696.  , 10302.      / 4157.  , 3285.  / 1367.  , -33.9  , 36.4  },
    },
    /* nderiv = 2 */ {
        /*  4 */  {  9977.  / 7732.  , 10631. / 11380. , 161.   / 11392. , 1.    / 266113. , 3672286074. / 14915. , 10959. / 8753.  , -1.53  , 0.891 },
        /*  5 */  {  6539.  / 4135.  , 14393. / 13401. , 1.     / 3263.  , 1.    / 78781.  , 94818740.   / 1793.  , 15697. / 9319.  , -2.97  , 6.64  },
        /*  6 */  {  25651. / 16404. , 3609.  / 3068.  , 76.    / 8979.  , 1.    / 297667. , 922006952.  / 3639.  , 7060.  / 4537.  , -2.07  , 4.93  },
        /*  7 */  {  17713. / 9399.  , 8075.  / 7143.  , 138.   / 11605. , 1.    / 136002. , 755158481.  / 6247.  , 4541.  / 2910.  , -3.02  , 7.53  },
        /*  8 */  {  18994. / 9237.  , 3977.  / 3539.  , 73.    / 7282.  , 1.    / 111904. , 577913389.  / 5810.  , 9851.  / 6166.  , -3.98  , 9.66  },
        /*  9 */  {  5092.  / 2353.  , 22628. / 20033. , -5.    / 6396.  , 1.    / 16535.  , 145714069.  / 10323. , 11201. / 6507.  , -6.81  , 15.8  },
        /* 10 */  {  23815. / 11282. , 5713.  / 4909.  , -57.   / 4691.  , 1.    / 17887.  , 199358385.  / 13336. , 24985. / 13879. , -8.94  , 19.0  },
        /* 11 */  {  54551. / 23412. , 8530.  / 7467.  , -921.  / 28645. , 1.    / 7601.   , 59931439.   / 9921.  , 16247. / 8406.  , -12.9  , 24.0  },
        /* 12 */  {  4520.  / 1879.  , 9869.  / 8601.  , -190.  / 3569.  , 2.    / 6037.   , 45877661.   / 19725. , 29355. / 14462. , -16.3  , 25.9  },
        /* 13 */  {  10841. / 5071.  , 3286.  / 2717.  , -166.  / 2157.  , 3.    / 10337.  , 42626923.   / 16678. , 7193.  / 3383.  , -19.6  , 27.1  },
        /* 14 */  {  16455. / 6572.  , 4577.  / 3939.  , -1623. / 16238. , 2.    / 7185.   , 37896608.   / 14629. , 22877. / 10405. , -22.4  , 25.6  },
        /* 15 */  {  70913. / 29232. , 15681. / 13229. , -7270. / 65299. , 2.    / 4523.   , 16300312.   / 10573. , 28627. / 12682. , -22.8  , 22.2  },
        /* 16 */  {  47761. / 23309. , 21087. / 16793. , -541.  / 4410.  , 8.    / 10739.  , 4825216.    / 5155.  , 24509. / 10946. , -22.6  , 22.4  },
        /* 17 */  {  54429. / 16724. , 5829.  / 5306.  , -1289. / 8489.  , 1091. / 9833.   , 41415.      / 6419.  , 16851. / 7441.  , -25.8  , 28.8  },
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
    const double errl = p[ndx][kdx][6]; // Measured in percent
    const double errh = p[ndx][kdx][7]; // Measured in percent

    /* Compute the empirical fit using */
    *C  = a * pow(k, b);
    *C *= 1 + c*(htdelta/log(k)) + d*(k/N)*(1 + e*pow(htdelta,f));

    // Using percent error bounds on the fit, estimate true range
    if (Clow)  *Clow  = *C / (1 - errl/100);
    if (Chigh) *Chigh = *C / (1 - errh/100);

    return SUZERAIN_SUCCESS;
}
