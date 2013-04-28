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
suzerain_bspline_htstretch1_evdeltascale(
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
