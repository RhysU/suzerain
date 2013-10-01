/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012, 2013 Rhys Ulerich
 * Copyright (C) 2012, 2013 The PECOS Development Team
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
 */

/** @file
 * @copydoc bspline.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/bspline.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>
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
    if (SUZERAIN_UNLIKELY(ldvalues && ldvalues < npoints)) {
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
    if (SUZERAIN_UNLIKELY(ldvalues && ldvalues < npoints)) {
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

// Parameters necessary for linear_combination_function.
typedef struct {
    size_t                       nderiv;
    const double                *coeffs;
    gsl_matrix                  *dB;
    gsl_bspline_workspace       *w;
    gsl_bspline_deriv_workspace *dw;
    double                       offset; // Additive constant
} linear_combination_params;

// A GSL-ready way to evaluate a B-spline function plus offset at some point x
static
double linear_combination_function(
        double x, void * params)
{
    double retval = GSL_NAN;
    linear_combination_params * p = (linear_combination_params *) params;
    suzerain_bspline_linear_combination(
            p->nderiv, p->coeffs, 1U, &x, &retval, 0U, p->dB, p->w, p->dw);
    return retval + p->offset;
}

int
suzerain_bspline_crossing(
    const size_t nderiv,
    const double * coeffs,
    const double value,
    double * lower,
    double * upper,
    const size_t maxiter,
    const double epsabs,
    const double epsrel,
    double * location,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    // Wrap the incoming parameters into an gsl_function for evaluation
    linear_combination_params params = { nderiv, coeffs, dB, w, dw, -value };
    gsl_function f                   = { linear_combination_function, &params};

    // Initialize fsolver to use Brent-Dekker on [lower, upper]
    // Bracketing, rather than fdfsolver, avoids exiting user-specified region
    gsl_root_fsolver * const s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    if (SUZERAIN_UNLIKELY(s == NULL)) {
        SUZERAIN_ERROR("Could not obtain gsl_root_fsolver", SUZERAIN_ENOMEM);
    }
    int status = gsl_root_fsolver_set(s, &f, *lower, *upper);
    status = (GSL_SUCCESS == status) ? GSL_CONTINUE : status;

    // Proceed until success, failure, or maximum iterations reached
    // On success, overwrite *location as described in API documentation.
    *location = GSL_NAN;
    for (size_t iter = 1; iter < maxiter && status == GSL_CONTINUE; ++iter) {
        if (GSL_SUCCESS != (status = gsl_root_fsolver_iterate(s))) break;
        *lower = gsl_root_fsolver_x_lower(s);
        *upper = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(*lower, *upper, epsabs, epsrel);
    }
    if (status == GSL_SUCCESS) {
        *location = gsl_root_fsolver_root(s);
    }
    gsl_root_fsolver_free(s);
    return status;
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
    if (SUZERAIN_UNLIKELY(tbl == NULL)) {
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
suzerain_bspline_distance(
    const gsl_bspline_workspace *a,
    const gsl_bspline_workspace *b)
{
    double retval = GSL_DBL_MAX;
    const size_t a_nknot = a->knots->size;
    const size_t b_nknot = b->knots->size;
    if (a->k == b->k && a->n == b->n && a_nknot == b_nknot) {
        retval = 0;
        for (size_t i = 0; i < a_nknot; ++i) {
            const double dist_i = fabs(  gsl_vector_get(a->knots, i)
                                       - gsl_vector_get(b->knots, i));
            retval = GSL_MAX(retval, dist_i);
        }
    }
    return retval;
}

const double suzerain_bspline_distance_distinct = 3 * GSL_DBL_EPSILON;

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
suzerain_bspline_htstretch2_evdeltascale_greville_abscissae(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(nderiv < 0)) {
        SUZERAIN_ERROR("Requested derivative nderiv must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(k < 1)) {
        SUZERAIN_ERROR("B-spline order k must be strictly positive",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(htdelta < 0)) {
        SUZERAIN_ERROR("Stretching parameter htdelta must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(N < k)) {
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
    static const double p[NDERIVHIGH - NDERIVLOW + 1][KHIGH - KLOW + 1][11] = {
    /* nderiv = 1 */ {
        /* 4*/ {    8753. /  6138.,    184. /   3513.,  -   3868. /   6327.,  - 22756. /  5215.,   19435. / 28258.,   152882005. /  9731.,    1605. / 12298.,     814. /  5783.,   7755. / 19621.,  -1.52,   0.46 },
        /* 5*/ { -182723. /  7970.,    121. /   6080.,  -7061095. /  23119.,   115483. / 33477.,       1. /  5560.,         197. /  3683.,     761. /  2741.,  - 2307. /  6251.,   1845. /  5548.,  -1.12,   0.53 },
        /* 6*/ {    8193. /  7183.,     74. /   7487.,     10787. /   6960.,  - 12529. /  6209.,    3196. /  4771.,     4441318. / 12819.,     932. / 12135.,    3915. / 17923.,    870. /  2123.,  -1.39,   1.14 },
        /* 7*/ {    9063. /  5683.,    291. /   6739.,  -  26786. /   3839.,  -  2757. /  1352.,    2959. /  4613.,     2666741. /  4919.,     507. /  6947.,     498. /  2627.,   3572. /  8791.,  -1.50,   1.12 },
        /* 8*/ {   20956. / 18111.,     73. /   3043.,      3316. /   4207.,  - 17945. / 10877.,    6659. / 11166.,     2564433. /  9169.,     658. / 10813.,    2725. / 12594.,   1181. /  2847.,  -1.48,   1.18 },
        /* 9*/ {    4351. /  3996.,    204. /   4093.,     13726. /   5957.,  - 19388. / 12047.,   10397. / 18003.,     1200048. /  3823.,     246. /  4289.,    1395. /  6839.,   4244. / 10183.,  -1.54,   1.09 },
        /*10*/ {    3174. /  3221.,    340. /   4047.,    435835. /  78898.,  - 24457. / 14830.,    3057. /  5416.,     5464428. / 12505.,     166. /  3013.,     679. /  3691.,   2669. /  6422.,  -1.59,   1.11 },
        /*11*/ {    6141. /  6356.,    513. /   5051.,     35557. /   5371.,  -  5359. /  3384.,    6824. / 12279.,     3628675. /  8504.,    1012. / 19681.,     363. /  2008.,   1873. /  4520.,  -1.56,   0.94 },
        /*12*/ {    9953. /  8604.,    521. /   4751.,  -   2623. /   3140.,  - 12011. /  8190.,    9067. / 16522.,     1947584. /  5723.,     141. /  2959.,    1231. /  6652.,   8479. / 20445.,  -1.58,   1.12 },
        /*13*/ {    6428. /  6247.,     28. /    239.,     28437. /   6178.,  - 13912. / 10071.,    2717. /  5125.,     3512083. / 11082.,     542. / 12275.,     701. /  3776.,   7291. / 17441.,  -1.58,   1.08 },
        /*14*/ {   18221. / 17050.,    573. /   3767.,     37351. /  13944.,  - 12384. /  8881.,    2238. /  4193.,     3701821. / 10597.,     231. /  5389.,    2953. / 16836.,   3950. /  9507.,  -1.59,   1.10 },
        /*15*/ {    9251. /  9369.,   1001. /   5400.,   1033331. / 141480.,  - 13367. /  9358.,    7583. / 14554.,     2662722. /  5819.,     392. /  9519.,     555. /  3394.,   2181. /  5234.,  -1.56,   1.01 },
        /*16*/ {    8003. /  8200.,   2503. /  12108.,    263686. /  31457.,  -  5713. /  4107.,    4803. /  9290.,     8189082. / 18419.,      93. /  2371.,    1663. / 10299.,   2945. /  7073.,  -1.56,   1.08 },
        /*17*/ {   10594. /  8835.,   2615. /  11361.,  -  38602. /   5641.,  -  3464. /  2537.,    3372. /  6575.,     4481119. / 10024.,     355. /  9469.,     811. /  5107.,   1129. /  2714.,  -1.53,   0.97 }
    },
    /* nderiv = 2 */ {
        /* 4*/ { -132683. /  8121.,  -   8. /   7617.,  -1156372. /   9759.,    78080. / 22333.,       2. / 10139.,          83. /  5144.,     911. /  2744.,  - 4791. / 13862.,   1479. /  4412.,  -1.86,   1.81 },
        /* 5*/ { - 96752. / 19419.,     27. /  10108.,  - 360604. /   4723.,    24551. /  9533.,       1. /  5309.,         277. /  4330.,   13424. / 46671.,  - 3888. / 10265.,   1289. /  3770.,  -0.50,   0.32 },
        /* 6*/ {    2635. /  8541.,  - 299. /  15112.,    174058. /  24237.,  - 17342. / 13001.,    2894. / 10471.,      291982. /  9057.,      68. / 17203.,    5579. / 13080.,   6911. /  9355.,  -1.51,   2.62 },
        /* 7*/ {    7358. / 11623.,  - 120. /   6169.,     82091. /  23515.,  - 10147. /  9621.,    2735. / 11778.,      212175. /  7154.,      10. /  6231.,   10354. / 24479.,   2416. /  3043.,  -1.70,   3.97 },
        /* 8*/ {    4751. /  5346.,   2143. /  18802.,  -  46259. /  19493.,  -  4506. /  9193.,       4. /  6951.,       49249. /  1482.,  -   24. /  8359.,     959. /  3662.,   9327. /  7114.,  -4.50,   6.86 },
        /* 9*/ {   17841. / 19519.,    553. /   8538.,  -  20719. /   8083.,  -  3783. /  7057.,     856. / 10529.,      121124. /  5101.,  -  115. / 30303.,    2467. /  7116.,   4707. /  4139.,  -3.68,   7.17 },
        /*10*/ {   20122. /  5145.,  -   5. /   6272.,  -1003532. /  10415.,  - 11287. / 15030.,    3429. / 17696.,      258125. /  9817.,       7. /  7694.,    4286. /  9151.,   8126. /  9217.,  -5.97,   8.01 },
        /*11*/ {   32624. / 11799.,     23. /   3905.,  - 466211. /   6855.,  - 14065. / 10992.,    4693. / 11010.,      335615. /  7641.,     187. /  7994.,    4512. /  8461.,   1469. /  2419.,  -6.09,   6.61 },
        /*12*/ {   21510. / 19531.,    407. /   3029.,  -  40991. /   5317.,  - 34333. / 13993.,    7253. /  9181.,     1291629. /  5929.,     217. /  3456.,    2659. /  8092.,   4103. /  9334.,  -6.88,   6.21 },
        /*13*/ {    3889. /  5285.,   2389. /  15117.,     52947. /   6029.,  - 17456. /  5339.,    4061. /  4625.,    11135709. /  8452.,     488. /  7087.,    1474. /  6209.,    383. /   987.,  -6.88,   5.89 },
        /*14*/ {   13010. / 10487.,   1375. /  10628.,  - 201367. /  12385.,  - 24232. /  6743.,   11512. / 10801.,     6297203. /  2487.,     428. /  6171.,     648. /  3823.,   4064. / 11831.,  -6.54,   5.79 },
        /*15*/ {   25912. / 21893.,   1831. /  16699.,  -  58570. /   4063.,  - 14377. /  4017.,   18889. / 18281.,    16561349. /  4390.,     272. /  4179.,     367. /  2451.,    919. /  2686.,  -6.57,   5.36 },
        /*16*/ {    3826. /  5191.,    800. /  13057.,    213790. /  16397.,  - 16873. /  4447.,   10559. /  9248.,    84802177. / 11958.,     607. /  9491.,    1313. / 14052.,   2667. /  8263.,  -6.11,   5.04 },
        /*17*/ {    6127. /  6960.,  - 163. /   9791.,     22545. /   4961.,  -115738. / 31923.,   10232. / 10879.,   167869433. / 10467.,     439. /  8114.,   10720. / 81537.,   3302. / 10275.,  -5.35,   4.96 }
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
    const double a    = p[ndx][kdx][ 0];
    const double b    = p[ndx][kdx][ 1];
    const double c    = p[ndx][kdx][ 2];
    const double d    = p[ndx][kdx][ 3];
    const double e    = p[ndx][kdx][ 4];
    const double f    = p[ndx][kdx][ 5];
    const double g    = p[ndx][kdx][ 6];
    const double h    = p[ndx][kdx][ 7];
    const double i    = p[ndx][kdx][ 8];
    const double errl = p[ndx][kdx][ 9]; // Relative error measured in percent
    const double errh = p[ndx][kdx][10]; // Relative error measured in percent

    /* Compute the empirical fit being careful about integer k, N */
#pragma warning(push,disable:981)
    const double deltahat = tanh(htdelta)*pow(1. + htdelta, i);
#pragma warning(pop)
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

// This routine is mostly copy-n-paste from above but it maintained separately
// for both readability and to simplify adding new coefficients in the future.
int
suzerain_bspline_htstretch1_evdeltascale_greville_abscissae(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(nderiv < 0)) {
        SUZERAIN_ERROR("Requested derivative nderiv must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(k < 1)) {
        SUZERAIN_ERROR("B-spline order k must be strictly positive",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(htdelta < 0)) {
        SUZERAIN_ERROR("Stretching parameter htdelta must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(N < k)) {
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
    static const double p[NDERIVHIGH - NDERIVLOW + 1][KHIGH - KLOW + 1][11] = {
    /* nderiv = 1 */ {
        /* 4*/ {   4129. /  3316.,    1052. / 13345.,      5769. /  7262.,  - 29140. /  6539.,    3051. /  5753.,    271951060. /   6969.,   3691. / 26472.,    458. /  3255.,    3944. / 13269., -2.65,  1.39 },
        /* 5*/ {   3747. /  2863.,     722. /  6223.,  -   3064. / 13069.,  - 43215. / 13954.,    4867. /  8717.,     54938146. /  12669.,    607. /  5225.,    863. /  5306.,    1663. /  5196., -2.73,  1.43 },
        /* 6*/ {    793. /   577.,    2234. / 16675.,  -  12137. /  6441.,  - 19450. /  5421.,    8944. / 17333.,    201430625. /   6661.,    799. /  7654.,   1290. /  9821.,     974. /  3283., -3.00,  1.48 },
        /* 7*/ {  10401. /  9794.,    1356. /  8845.,     30469. / 10756.,  - 11151. /  3473.,    3321. /  7345.,   1029376385. /  36246.,    377. /  4242.,   1838. / 13631.,    2446. /  8123., -2.88,  1.50 },
        /* 8*/ {  13867. / 11300.,    1847. / 10004.,  -    899. /  1084.,  - 24247. /  8311.,    2025. /  4691.,    208465759. /   9898.,   1195. / 15107.,    131. /   957.,    3041. /  9933., -2.90,  1.47 },
        /* 9*/ {  11292. / 11417.,    2667. / 12475.,     71783. / 14593.,  - 23753. /  8884.,   15279. / 38327.,    170388497. /  10015.,    913. / 12919.,   2803. / 19841.,   16957. / 54130., -2.87,  1.41 },
        /*10*/ {   4941. /  9610.,     830. /  3611.,    279834. / 13813.,  -  4892. /  2245.,    1823. /  4351.,     11323441. /   2063.,    366. /  5827.,   1918. / 12879.,    2981. /  9362., -2.83,  1.42 },
        /*11*/ {  12725. / 12639.,    1233. /  4685.,     29103. /  5882.,  - 17641. /  7391.,     812. /  2251.,     93817310. /   6013.,    575. /  9852.,   1238. /  8665.,    2801. /  8822., -2.79,  1.32 },
        /*12*/ {   5215. /  5156.,    1009. /  3470.,     77977. / 15613.,  - 29719. / 12774.,    3890. / 10879.,     75236231. /   4566.,    579. / 10655.,   3977. / 28471.,    2948. /  9301., -2.86,  1.36 },
        /*13*/ {  16753. / 16972.,     794. /  2485.,    100717. / 15948.,  - 34434. / 15637.,    2493. /  7441.,     61372987. /   4090.,    479. /  9560.,    451. /  3156.,    1286. /  3991., -2.82,  1.34 },
        /*14*/ {   9761. /  9871.,    1461. /  4217.,       281. /    43.,  - 33703. / 15770.,    6378. / 19271.,    170777218. /  11389.,    543. / 11569.,   1880. / 13293.,    4685. / 14541., -2.82,  1.33 },
        /*15*/ {   7329. /  7286.,    2281. /  6042.,     48733. /  8398.,  -  8421. /  4222.,    1839. /  6014.,     87107021. /   7039.,   1578. / 36371.,   1153. /  7805.,    1189. /  3613., -2.72,  1.29 },
        /*16*/ {   7733. /  8163.,    3533. /  8665.,    121822. / 12527.,  -  6932. /  3573.,    2407. /  8036.,     57224710. /   4651.,    286. /  6989.,    377. /  2558.,    4584. / 13895., -2.72,  1.29 },
        /*17*/ {  10349. /  9079.,    4121. /  9390.,  -  49628. / 15233.,  - 16811. /  8933.,    1143. /  3880.,    103494425. /   8753.,   1589. / 41049.,    817. /  5544.,    2196. /  6637., -2.66,  1.24 }
    },
    /* nderiv = 2 */ {
        /* 4*/ { -37181. /  3006.,  -   35. / 19004.,  -1709799. / 12403.,    11429. /  3303.,       1. /  4633.,          207. /   9529.,   8025. / 25267.,  -1792. /  5289.,    1579. /  6332., -2.73,  2.66 },
        /* 5*/ {   6843. /  7364.,  -  113. /  5126.,  -  19508. / 11693.,  - 15410. /  9853.,    2215. /  5596.,       199827. /   5060.,    117. /  5482.,   3886. /  7977.,    3618. /  7763., -0.63,  1.38 },
        /* 6*/ {   1309. /  1894.,     472. /  6839.,     10796. / 10079.,  - 13268. / 14331.,    1942. / 14501.,       169657. /   5987.,      2. /  3287.,   1786. /  3447.,    9664. / 13071., -3.14,  5.22 },
        /* 7*/ {  62991. / 64619.,   10271. / 61625.,  - 146725. / 39021.,  -  6349. / 10362.,      23. /  3440.,        26371. /    794.,  -  17. / 10006.,    798. /  1681.,   13080. / 14621., -4.53,  9.31 },
        /* 8*/ {   5701. /  6404.,     543. / 10115.,  -  15926. / 10845.,  -  9214. / 11799.,    1960. /  8807.,       759503. /  31646.,      9. /  4685.,   3636. /  5393.,    5996. /  8919., -5.10,  7.39 },
        /* 9*/ {  57103. / 12153.,     349. /  6313.,  - 292899. /  2840.,  - 15197. / 10802.,    3317. /  6840.,      6115361. / 140104.,    431. /  8216.,   7073. / 12325.,   10286. / 21381., -6.57,  7.12 },
        /*10*/ {   4753. /  6483.,     215. /  1719.,     32327. /  6251.,  -  8563. /  4452.,    1981. /  2831.,      1180287. /  14438.,    599. /  7682.,   1695. /  4003.,    3259. /  8267., -8.44,  6.56 },
        /*11*/ {   3058. /  4871.,    2994. / 14713.,     34703. /  3382.,  - 15557. /  4053.,    6541. /  5144.,      5024384. /   5503.,   2720. / 22829.,    961. / 20910.,    2347. /  7629., -6.87,  7.13 },
        /*12*/ {  12242. / 10351.,     785. /  5066.,  - 108847. /  9907.,  - 38407. /  8945.,    8189. /  5423.,     54739986. /  21557.,   1903. / 16894.,  - 381. /  7330.,    6661. / 24933., -7.21,  7.11 },
        /*13*/ {  24826. / 20307.,    4749. / 40912.,  - 170233. / 12219.,  -  4785. /  1084.,   16494. / 11557.,    720241352. / 105099.,   2930. / 29579.,  - 807. / 19174.,    2725. / 10439., -7.16,  6.65 },
        /*14*/ {  13245. / 12827.,     163. /  5022.,  -  25753. /  4855.,  -176237. / 37876.,   21379. / 17208.,    315490949. /   7639.,   4100. / 49413.,      2. /  4983.,    1942. /  7671., -6.65,  6.31 },
        /*15*/ {   5358. /  5287.,  -  142. / 11949.,  -  16365. /  3631.,  -173358. / 39355.,   10908. / 10097.,    529204799. /   8657.,    119. /  1628.,    216. /  7747.,    1637. /  6397., -6.47,  5.75 },
        /*16*/ {   8995. /  9151.,  -  367. /  3056.,  -  13906. /  5041.,  - 40441. / 10087.,    4335. /  4703.,    366898466. /   4633.,    324. /  5153.,    440. /  8139.,    1633. /  6441., -5.99,  5.43 },
        /*17*/ {  19825. / 20211.,  - 3148. / 12971.,  -   7791. /  2918.,  - 27739. /  8000.,    4651. /  6107.,    516527143. /   8065.,    321. /  5962.,   1049. / 12527.,    1391. /  5436., -5.04,  5.04 }
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
    const double a    = p[ndx][kdx][ 0];
    const double b    = p[ndx][kdx][ 1];
    const double c    = p[ndx][kdx][ 2];
    const double d    = p[ndx][kdx][ 3];
    const double e    = p[ndx][kdx][ 4];
    const double f    = p[ndx][kdx][ 5];
    const double g    = p[ndx][kdx][ 6];
    const double h    = p[ndx][kdx][ 7];
    const double i    = p[ndx][kdx][ 8];
    const double errl = p[ndx][kdx][ 9]; // Relative error measured in percent
    const double errh = p[ndx][kdx][10]; // Relative error measured in percent

    /* Compute the empirical fit being careful about integer k, N */
#pragma warning(push,disable:981)
    const double deltahat = tanh(htdelta)*pow(1. + htdelta, i);
#pragma warning(pop)
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

int
suzerain_bspline_htstretch2_evdeltascale_breakpoints(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(nderiv < 0)) {
        SUZERAIN_ERROR("Requested derivative nderiv must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(k < 1)) {
        SUZERAIN_ERROR("B-spline order k must be strictly positive",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(htdelta < 0)) {
        SUZERAIN_ERROR("Stretching parameter htdelta must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(N < k)) {
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
    static const double p[NDERIVHIGH - NDERIVLOW + 1][KHIGH - KLOW + 1][11] = {
    /* nderiv = 1 */ {
        /*04*/  {  -994345./51933,  -139./11711,   -2110264./12191,  +28112./7633,    +2./10677,      +27./1829,           +2065./6124,    -3838./10157,   +1453./4742,    -2.68,  1.44  },
        /*05*/  {  -93077./15015,   +489./98290,   -345686./4859,    +29437./11373,   +1./5560,       +1108./20717,        +319./1149,     -3303./8950,    +4988./14999,   -1.13,  0.53  },
        /*06*/  {  +6279./10364,    +35./10096,    -47259./9020,     -25624./8609,    +4234./6261,    +2309833./5944,      +1157./14714,   +6707./31646,   +5623./13768,   -1.41,  1.15  },
        /*07*/  {  +5631./9478,     +46./7871,     -56332./7775,     -11902./4115,    +5599./8803,    +2255441./4813,      +935./13093,    +1442./7351,    +10872./26665,  -1.49,  1.12  },
        /*08*/  {  +2131./3722,     +33./9395,     -18268./2007,     -23766./9181,    +2393./4010,    +1290613./4587,      +853./13997,    +2345./10858,   +1333./3214,    -1.48,  1.19  },
        /*09*/  {  +8346./14407,    +21./2666,     -181851./15575,   -92355./35161,   +4731./8081,    +4026203./10967,     +885./14998,    +2650./13573,   +1162./2801,    -1.56,  1.10  },
        /*10*/  {  +12556./22955,   +47./4758,     -22051./1668,     -52441./19912,   +2852./5017,    +5531294./11845,     +189./3392,     +1473./8165,    +2714./6545,    -1.60,  1.11  },
        /*11*/  {  +6505./12074,    +53./5065,     -229994./14869,   -30809./12026,   +8299./14849,   +4589365./10364,     +817./15772,    +941./5273,     +4955./11977,   -1.57,  0.94  },
        /*12*/  {  +1027./2032,     +84./7387,     -97849./5855,     -18003./7204,    +4677./8414,    +2289877./5699,      +142./2909,     +2992./16927,   +5895./14276,   -1.59,  1.12  },
        /*13*/  {  +2257./4145,     +113./9434,    -290747./13769,   -24527./9972,    +3037./5584,    +9161227./22172,     +1063./23092,   +2065./12046,   +5725./13808,   -1.60,  1.08  },
        /*14*/  {  +1543./2888,     +13./1111,     -132939./5663,    -16089./6799,    +4452./8341,    +4863208./13921,     +575./13414,    +3896./22213,   +11870./28569,  -1.59,  1.10  },
        /*15*/  {  +22653./43537,   +35./2641,     -140002./5467,    -7540./3137,     +5961./11440,   +771535./1683,       +628./15247,    +1162./7109,    +2028./4867,    -1.57,  1.01  },
        /*16*/  {  +4214./7709,     +359./26018,   -285669./9440,    -24301./10259,   +1709./3305,    +6955255./15601,     +115./2931,     +998./6185,     +1727./4148,    -1.57,  1.06  },
        /*17*/  {  +2123./4028,     +153./10622,   -150472./4685,    -12858./5483,    +1010./1969,    +5960264./13293,     +432./11519,    +1059./6674,    +4727./11364,   -1.53,  0.97  }
    },
    /* nderiv = 2 */ {
        /*04*/  {  -27907./6471,    -1./2862,      -859927./17700,   +12097./4474,    +1./5071,       +132./8183,          +1402./4223,    -2741./7931,    +2897./8642,    -1.87,  1.82  },
        /*05*/  {  -47575./44741,   +8./12023,     -426749./20210,   +14057./8202,    +1./5306,       +199./3110,          +1454./5055,    -486./1283,     +3397./9936,    -0.50,  0.33  },
        /*06*/  {  -53385./10876,   -248./11847,   -918950./10271,   +16411./6994,    +1./3766,       +437./7105,          +10911./32239,  -2164./4157,    +4163./9164,    -3.82,  6.55  },
        /*07*/  {  -33501./7315,    -1201./48968,  -2598654./23501,  +11071./4992,    +3./11623,      +2140./26499,        +417./1345,     -3449./6311,    +4172./8755,    -4.37,  8.28  },
        /*08*/  {  +7690./13711,    +428./26279,   -133719./13165,   -11069./7762,    +5./7682,       +307111./9246,       -38./13071,     +502./1915,     +10037./7659,   -4.50,  6.86  },
        /*09*/  {  +3371./6391,     +35./4329,     -146754./12787,   -8916./6013,     +513./6302,     +53603./2257,        -21./5543,      +1495./4311,    +26657./23447,  -3.69,  7.18  },
        /*10*/  {  -25328./6131,    -244./4271,    -2506057./14800,  +11537./5822,    +3./8420,       +1237./5679,         +4033./13306,   -6484./8521,    +10813./19626,  -7.67,  6.64  },
        /*11*/  {  +18808./33329,   +1./2285,      -172691./9981,    -19452./8719,    +10135./23928,  +906281./20810,      +109./4736,     +14088./26381,  +9289./15254,   -6.11,  6.62  },
        /*12*/  {  +2267./4067,     +128./10473,   -301902./15317,   -26599./7786,    +14477./18417,  +2247545./10262,     +689./10975,    +2034./6199,    +7794./17687,   -6.89,  6.23  },
        /*13*/  {  +561./1193,      +69./5231,     -64062./3475,     -81411./19190,   +9189./10114,   +9619873./8086,      +1133./16048,   +4387./19624,   +2665./6837,    -6.88,  5.90  },
        /*14*/  {  +1416./2729,     +821./75258,   -196892./8403,    -136799./26769,  +3669./2908,    +11330951./2625,     +997./12735,    +323./3484,     +2035./6063,    -6.48,  5.74  },
        /*15*/  {  +4603./8165,     +125./14022,   -309439./10719,   -34303./6979,    +23489./18275,  +42428569./12008,    +221./2932,     +466./8579,     +6257./18589,   -6.48,  5.36  },
        /*16*/  {  +4366./12887,    +3./1381,      -176039./9963,    -38866./6841,    +1814./7671,    +17659379991./9881,  +1002./20813,   +2484./9679,    +4279./13457,   -6.21,  5.03  },
        /*17*/  {  +7109./9432,     -8./8887,      -242221./4975,    -77725./14824,   +4415./4977,    +616106394./5905,    +779./14545,    +563./3963,     +2044./6689,    -5.30,  4.95  }
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
    const double a    = p[ndx][kdx][ 0];
    const double b    = p[ndx][kdx][ 1];
    const double c    = p[ndx][kdx][ 2];
    const double d    = p[ndx][kdx][ 3];
    const double e    = p[ndx][kdx][ 4];
    const double f    = p[ndx][kdx][ 5];
    const double g    = p[ndx][kdx][ 6];
    const double h    = p[ndx][kdx][ 7];
    const double i    = p[ndx][kdx][ 8];
    const double errl = p[ndx][kdx][ 9]; // Relative error measured in percent
    const double errh = p[ndx][kdx][10]; // Relative error measured in percent

    /* Compute the empirical fit being careful about integer k, N */
#pragma warning(push,disable:981)
    const double deltahat = tanh(htdelta)*pow(1. + htdelta, i);
#pragma warning(pop)
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

// This routine is mostly copy-n-paste from above but it maintained separately
// for both readability and to simplify adding new coefficients in the future.
int
suzerain_bspline_htstretch1_evdeltascale_breakpoints(
    const int nderiv,
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(nderiv < 0)) {
        SUZERAIN_ERROR("Requested derivative nderiv must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(k < 1)) {
        SUZERAIN_ERROR("B-spline order k must be strictly positive",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(htdelta < 0)) {
        SUZERAIN_ERROR("Stretching parameter htdelta must be nonnegative",
                       SUZERAIN_EINVAL);
    } else if (SUZERAIN_UNLIKELY(N < k)) {
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
    static const double p[NDERIVHIGH - NDERIVLOW + 1][KHIGH - KLOW + 1][11] = {
    /* nderiv = 1 */ {
        /*04*/  {  -163797./7976,  -129./9529,    -1758235./13096,  +21164./5847,    +1./5226,       +106./6285,          +2814./9161,   -2449./7163,   +4175./18531,   -4.31,  1.57  },
        /*05*/  {  +5663./7242,    +593./19700,   -17774./3457,     -39874./9513,    +1649./3182,    +40507958./5253,     +1295./11043,  +2446./15393,  +3661./11453,   -2.74,  1.42  },
        /*06*/  {  +10265./15782,  +297./11086,   -23273./3947,     -42986./9589,    +3201./6197,    +484246172./16157,   +1281./12272,  +1548./11783,  +1851./6239,    -3.00,  1.49  },
        /*07*/  {  +57288./85633,  +322./12613,   -71132./8233,     -18617./4525,    +3575./7856,    +444110233./16461,   +724./8149,    +1225./9076,   +2373./7880,    -2.88,  1.50  },
        /*08*/  {  +8588./13239,   +95./3701,     -117987./10895,   -34861./9677,    +2442./4861,    +54119060./6243,     +259./3247,    +1924./14379,  +5498./18009,   -2.93,  1.44  },
        /*09*/  {  +1017./1747,    +688./25731,   -137439./11683,   -25133./6861,    +7192./17697,   +259647606./14077,   +213./2990,    +2885./20888,  +1409./4507,    -2.89,  1.39  },
        /*10*/  {  +20923./32505,  +167./6451,    -37233./2285,     -55810./16973,   +4232./9843,    +71707295./9002,     +302./4689,    +599./4306,    +454./1435,     -2.88,  1.36  },
        /*11*/  {  +22161./29455,  +320./11943,   -653677./28085,   -6501./1549,     +2139./7126,    +2128745627./9929,   +668./10895,   +3111./25325,  +2176./6937,    -2.88,  1.20  },
        /*12*/  {  +4646./7141,    +341./12725,   -156069./6856,    -23731./6285,    +12817./36191,  +1122612424./16775,  +401./7061,    +456./3751,    +931./2970,     -2.94,  1.25  },
        /*13*/  {  +2715./5072,    +88./3391,     -130763./6314,    -8653./2978,     +2033./5338,    +50503691./9068,     +417./8381,    +776./5319,    +14263./44203,  -2.82,  1.36  },
        /*14*/  {  +2926./5241,    +177./6796,    -360413./14557,   -22816./7957,    +1289./3401,    +35465564./6147,     +664./14173,   +4354./30573,  +6399./19859,   -2.83,  1.33  },
        /*15*/  {  +3969./7292,    +1003./37197,  -884245./32711,   -30830./10389,   +852./2777,     +169791343./13854,   +790./18201,   +1310./8877,   +1462./4443,    -2.72,  1.29  },
        /*16*/  {  +11651./19500,  +691./25426,   -265573./7916,    -21956./7537,    +1642./5449,    +104203394./8647,    +341./8328,    +787./5348,    +947./2871,     -2.72,  1.29  },
        /*17*/  {  +4658./8361,    +32./1195,     -307805./8984,    -86257./32806,   +4487./13460,   +84391415./18136,    +98./2553,     +2689./17830,  +417./1258,     -2.66,  1.26  }
    },
    /* nderiv = 2 */ {
        /*04*/  {  -25./2193,      -29./7386,     +21705./12154,    -16879./6734,    +2466./6163,    +519161./11194,      +173./8699,    +2408./5643,   +3539./7983,    -0.55,  1.40  },
        /*05*/  {  +406./4635,     -86./15761,    +31754./26355,    -10904./4445,    +12377./30774,  +235442./5763,       +256./10679,   +1933./4013,   +2401./5193,    -0.69,  1.41  },
        /*06*/  {  -255./3931,     +4./289,       +35710./11163,    -8181./4486,     +617./4620,     +207329./7312,       +3./4972,      +4834./9333,   +4369./5907,    -3.15,  5.22  },
        /*07*/  {  +145./3468,     +11./392,      +12625./7842,     -14327./9316,    +5./1154,       +153037./4538,       -19./10832,    +4709./9858,   +4091./4595,    -4.55,  9.30  },
        /*08*/  {  +353./12750,    +211./27483,   +15362./7525,     -11009./6412,    +1777./7985,    +1069317./44539,     +99./51035,    +1736./2575,   +201./299,      -5.10,  7.40  },
        /*09*/  {  +837./17816,    +101./14409,   +12500./7361,     -13912./5901,    +2173./4474,    +388075./8851,       +160./3031,    +9689./16929,  +5212./10837,   -6.57,  7.12  },
        /*10*/  {  +567./11569,    +169./8885,    +16629./10378,    -20250./6227,    +5989./8148,    +1014609./7487,      +1057./11463,  +1182./3541,   +1651./4190,    -8.27,  6.55  },
        /*11*/  {  +339./6184,     +131./6414,    +29995./22786,    -22258./4577,    +3982./3159,    +5343089./5085,      +2084./17565,  +271./5476,    +2004./6509,    -6.88,  7.13  },
        /*12*/  {  +709./17736,    +154./10635,   +38802./21445,    -68189./10448,   +30331./18247,  +474061690./23789,   +1247./10740,  -1041./12766,  +2209./8344,    -7.21,  7.09  },
        /*13*/  {  +777./16762,    +166./17161,   +14992./10401,    -61289./11316,   +3197./2249,    +14382851./1930,     +1004./10157,  -564./14039,   +1671./6391,    -7.17,  6.65  },
        /*14*/  {  +675./15224,    +35./12347,    +14493./10126,    -10967./2119,    +5541./4079,    +199291147./19442,   +506./5809,    -473./11713,   +3999./15815,   -6.61,  6.30  },
        /*15*/  {  +948./21929,    -1./28761,     +3081./2242,      -60199./12333,   +11503./8832,   +44740778./5103,     +915./11414,   -289./6281,    +544./2137,     -6.41,  5.74  },
        /*16*/  {  +355./8042,     -19./3033,     +3677./3091,      -49327./10489,   +10851./8434,   +344320415./33871,   +619./8433,    -1813./28867,  +4928./19725,   -5.89,  5.41  },
        /*17*/  {  +2325./54836,   -221./18262,   +5647./4841,      -217207./47365,  +5780./4629,    +75172547./5426,     +1650./24589,  -1219./17223,  +3585./14351,   -4.92,  5.02  }
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
    const double a    = p[ndx][kdx][ 0];
    const double b    = p[ndx][kdx][ 1];
    const double c    = p[ndx][kdx][ 2];
    const double d    = p[ndx][kdx][ 3];
    const double e    = p[ndx][kdx][ 4];
    const double f    = p[ndx][kdx][ 5];
    const double g    = p[ndx][kdx][ 6];
    const double h    = p[ndx][kdx][ 7];
    const double i    = p[ndx][kdx][ 8];
    const double errl = p[ndx][kdx][ 9]; // Relative error measured in percent
    const double errh = p[ndx][kdx][10]; // Relative error measured in percent

    /* Compute the empirical fit being careful about integer k, N */
#pragma warning(push,disable:981)
    const double deltahat = tanh(htdelta)*pow(1. + htdelta, i);
#pragma warning(pop)
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
