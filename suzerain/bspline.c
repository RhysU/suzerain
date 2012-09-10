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

static int
suzerain_bspline_htstretch2_evdeltascale1(
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh)
{

    /* Coefficient lookup table for orders 4 through 17 from model document */
    /* The final two columns are the low and higher estimated error bounds */
#define KLOW  ( 4)
#define KHIGH (17)
    static const double p[KHIGH - KLOW + 1][7] = {
        /* 4*/ {1037./641., 473./20764., 39./224212., 792520./93., 3155./2478., -2.66 , 1.41},
        /* 5*/ {739./455., 167./6799., 60./641219., 369130./27., 610./457., -2.72 , 1.26},
        /* 6*/ {1849./1135., 600./24907., 51./804298., 422796./23., 3585./2537., -2.80 , 2.30},
        /* 7*/ {4829./2957., 22./1053., 39./374987., 358703./34., 3643./2423., -2.81 , 5.94},
        /* 8*/ {3275./2003., 583./37716., 57./47897., 5279./6., 779./491., -4.19 , 9.87},
        /* 9*/ {9125./5576., 282./62887., 155./47511., 94118./309., 2235./1316., -6.82 , 15.9},
        /*10*/ {1369./836., -39./3434., 501./74569., 27691./196., 5078./2803., -10.0 , 20.5},
        /*11*/ {4046./2469., -252./7727., 693./53629., 279./4., 9532./4927., -13.9 , 25.8},
        /*12*/ {1186./723., -424./6971., 911./94585., 73428./823., 1443./701., -18.7 ,  29.6},
        /*13*/ {4945./3012., -653./7141., 625./41087., 45539./842., 2321./1071., -23.3 ,  31.0},
        /*14*/ {2523./1535., -531./4168., 601./16577., 17133./827., 508./219., -28.7 ,  34.2},
        /*15*/ {1876./1143., -937./7195., 338./4835., 5994./571., 5600./2413., -27.1 ,  25.0},
        /*16*/ {3177./1936., -617./3756., 217./1839., 5409./839., 4429./1898., -30.5 ,  25.6},
        /*17*/ {9839./6006., -172./891., 1199./3983., 5723./2355., 2442./1015., -33.9 ,  36.3}
    };

    /* Where are we in the coefficient table? */
    const int ndx = k - KLOW;
    if (ndx < 0 || (size_t) ndx >= sizeof(p)/sizeof(p[0])) {
#define XSTR(s) STR(s)
#define STR(s)  #s
        SUZERAIN_ERROR("B-spline order k outside of implemented range: ["
                       XSTR(KLOW) "," XSTR(KHIGH) "]",SUZERAIN_EINVAL);
#undef STR
#undef XSTR
    }
#undef KLOW
#undef KHIGH
    const double a    = p[ndx][0];
    const double b    = p[ndx][1];
    const double c    = p[ndx][2];
    const double d    = p[ndx][3];
    const double e    = p[ndx][4];
    const double errl = p[ndx][5]; // Measured in percent
    const double errh = p[ndx][6]; // Measured in percent

    // C^{(1)} &\approx k^a \left(
    //     1 + b \frac{\delta}{\ln{} k} + c \frac{k}{N_y}
    //     \left(1 + d \delta^e\right)
    // \right)
    *C  = pow(k, a);
    *C *= 1 + b*(htdelta/log(k)) + c*(k/N)*(1 + d*pow(htdelta,e));

    // Using percent error bounds on the fit, estimate true range
    if (Clow)  *Clow  = *C / (1 - errl/100);
    if (Chigh) *Chigh = *C / (1 - errh/100);

    return SUZERAIN_SUCCESS;
}

static int
suzerain_bspline_htstretch2_evdeltascale2(
    const int k,
    const double htdelta,
    const int N,
    double * const C,
    double * const Clow,
    double * const Chigh)
{
    /* Coefficient lookup table for orders 4 through 17 from model document */
    /* The final two columns are the low and higher estimated error bounds */
#define KLOW  ( 4)
#define KHIGH (17)
    static const double p[KHIGH - KLOW + 1][11] = {
        /* 4*/{ 570685./9., -590851./4., -4514312., -49033./14., -104./557., 4196./6889., -603./2624., -493./1455., -5736./793., -172 ,  84.0 },
        /* 5*/{ 4759./187., -1556./851., 37795./66., -897./211., -483./626., 6320./4273., -4206./937., 9881./2059., 3425./4039., -0.789 ,  3.14 },
        /* 6*/{ 8471./1718., -939./1994., 10683./5., -490./151., -15./49., 3728./2559., -486./109., 929./1033., 5515./2987., -0.387 ,  1.37 },
        /* 7*/{ 5437./950., 2886./553., 169467./77., -1522./431., -5537./906., 9485./6452., -1103./243., 469835./4., -75727., -0.248 ,  0.838 },
        /* 8*/{ 550042./11., -1599875./4., -2147483648./23., -423998./9., -571./4932., 8527./5767., -1873./406., 212./27., 79./92., -0.187 ,  0.622 },
        /* 9*/{ 2056./379., 2385./623., 128515./49., -1387./401., -412./1085., 4763./3209., -2399./512., 24009./478., -637./1051., -0.153 ,  0.508 },
        /*10*/{ 187./35., 1389./410., 11333./4., -905./266., -757./1996., 4344./2917., -2161./454., 4763./6571., 6062./4759., -0.130 ,  0.432 },
        /*11*/{ 5413./1024., 4073./1377., 36127./12., -1049./314., -2971./507., 611./409., -1592./329., 7858./1401., -1306./33., -0.113 ,  0.376 },
        /*12*/{ 8787./1690., 8389./3368., 25183./8., -1765./542., -3149./496., 3589./2395., -2245./456., 50386./281., -13591./28., -0.0998 ,  0.329 },
        /*13*/{ 4815./941., 2781./1688., 229996./71., -723./229., -965./2096., 2587./1721., -2552./509., 9183./7750., 4387./4524., -0.0873 ,  0.288 },
        /*14*/{ 5653./1135., -863./855., 153347./45., -2020./671., -780./1633., 9431./6255., -46./9., 3115./944., 4083./6946., -0.0893 ,  0.249 },
        /*15*/{ 2440./521., -620./137., 46790./13., -1724./607., -66./137., 9363./6191., -4987./956., 1542./703., 463./598., -0.125 ,  0.212 },
        /*16*/{ 929./179., 1147./349., 251877./50., -595./181., -3176./299., 7861./5182., -2879./540., 231./94., -45727./78., -0.164 ,  0.182 },
        /*17*/{ 3497./674., 2222./731., 325651./61., -75./23., -2724./215., 3836./2521., -2363./433., 8405./906., -4220./169., -0.207 ,  0.193 }
    };

    /* Where are we in the coefficient table? */
    const int ndx = k - KLOW;
    if (ndx < 0 || (size_t) ndx >= sizeof(p)/sizeof(p[0])) {
#define XSTR(s) STR(s)
#define STR(s)  #s
        SUZERAIN_ERROR("B-spline order k outside of implemented range: ["
                       XSTR(KLOW) "," XSTR(KHIGH) "]",SUZERAIN_EINVAL);
#undef STR
#undef XSTR
    }
#undef KLOW
#undef KHIGH
    const double a    = p[ndx][ 0];
    const double b    = p[ndx][ 1];
    const double c    = p[ndx][ 2];
    const double d    = p[ndx][ 3];
    const double e    = p[ndx][ 4];
    const double f    = p[ndx][ 5];
    const double g    = p[ndx][ 6];
    const double h    = p[ndx][ 7];
    const double i    = p[ndx][ 8];
    const double errl = p[ndx][ 9]; // Measured in percent
    const double errh = p[ndx][10]; // Measured in percent

    // C^{(2)} &\approx \left(a + \frac{b}{k} + c k^d\right)
    //     \left(
    //       1 + e \delta^f \left(\ln N_y\right)^g
    //       \left(1 + h k^i\right)
    //     \right)
    *C  = a + b/k + c*pow(k,d);
    *C *= 1 + e*pow(htdelta,f)*pow(log(N),g)*(1 + h*pow(k,i));

    // Using percent error bounds on the fit, estimate true range
    if (Clow)  *Clow  = *C / (1 - errl/100);
    if (Chigh) *Chigh = *C / (1 - errh/100);

    return SUZERAIN_SUCCESS;
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

    /* Dispatch to derivative-specific fits */
    switch (nderiv) {
        default:
            SUZERAIN_ERROR("Derivative nderiv unimplemented",
                           SUZERAIN_FAILURE);
        case 1:
            return suzerain_bspline_htstretch2_evdeltascale1(
                    k, htdelta, N, C, Clow, Chigh);
        case 2:
            return suzerain_bspline_htstretch2_evdeltascale2(
                    k, htdelta, N, C, Clow, Chigh);
    }
}
