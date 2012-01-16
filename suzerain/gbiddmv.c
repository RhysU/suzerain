/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * gbiddmv.c: provides custom, BLAS-like general band matrix-vector operations
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

// File iteration is used to generate bandwidth-specific routines.  See
// http://www.boost.org/libs/preprocessor/doc/topics/file_iteration.html for an
// overview.  We iterate first to generate internal, static routines and then
// continue on to the non-iterated logic.

#if !defined(BOOST_PP_IS_ITERATING) || !BOOST_PP_IS_ITERATING

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/gbiddmv.h>

#pragma warning(disable:1418 1572 2259)

// ---------------------------------------------------------
// Iterate on this file to generate fixed bandwidth routines
// ---------------------------------------------------------
#define BOOST_PP_ITERATION_LIMITS (0, 15)
#define BOOST_PP_FILENAME_1 <suzerain/gbiddmv.c>
#include BOOST_PP_ITERATE()

// -----------------------------------
// Generate general bandwidth routines
// -----------------------------------

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  suzerain_gbiddmv_internal_s
#define GBIDDMV_COMPONENT float
#define GBIDDMV_SCALAR    float
#define GBIDDMV_KL        const int kl,
#define GBIDDMV_KU        const int ku,
#include "gbiddmv.def"

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  suzerain_gbiddmv_internal_d
#define GBIDDMV_COMPONENT double
#define GBIDDMV_SCALAR    double
#define GBIDDMV_KL        const int kl,
#define GBIDDMV_KU        const int ku,
#include "gbiddmv.def"

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  suzerain_gbiddmv_internal_sc
#define GBIDDMV_COMPONENT float
#define GBIDDMV_SCALAR    float _Complex
#define GBIDDMV_KL        const int kl,
#define GBIDDMV_KU        const int ku,
#include "gbiddmv.def"

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  suzerain_gbiddmv_internal_dz
#define GBIDDMV_COMPONENT double
#define GBIDDMV_SCALAR    double _Complex
#define GBIDDMV_KL        const int kl,
#define GBIDDMV_KU        const int ku,
#include "gbiddmv.def"

// ------------------------------------------------------------------
// Provide externally callable logic dispatching to internal routines
// ------------------------------------------------------------------

int
suzerain_gbiddmv_s(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0,
        const float alpha1,
        const float *d1,
        const float alpha2,
        const float *d2,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
#define ARGS trans, n, alpha0, alpha1, d1, alpha2, d2, \
            a, lda, x, incx, beta, y, incy
            case  0: return suzerain_gbiddmv_internal_s0(ARGS);
            case  1: return suzerain_gbiddmv_internal_s1(ARGS);
            case  2: return suzerain_gbiddmv_internal_s2(ARGS);
            case  3: return suzerain_gbiddmv_internal_s3(ARGS);
            case  4: return suzerain_gbiddmv_internal_s4(ARGS);
            case  5: return suzerain_gbiddmv_internal_s5(ARGS);
            case  6: return suzerain_gbiddmv_internal_s6(ARGS);
            case  7: return suzerain_gbiddmv_internal_s7(ARGS);
            case  8: return suzerain_gbiddmv_internal_s8(ARGS);
            case  9: return suzerain_gbiddmv_internal_s9(ARGS);
            case 10: return suzerain_gbiddmv_internal_s10(ARGS);
            case 11: return suzerain_gbiddmv_internal_s11(ARGS);
            case 12: return suzerain_gbiddmv_internal_s12(ARGS);
            case 13: return suzerain_gbiddmv_internal_s13(ARGS);
            case 14: return suzerain_gbiddmv_internal_s14(ARGS);
            case 15: return suzerain_gbiddmv_internal_s15(ARGS);
#undef ARGS
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbiddmv_internal_s(trans, n, kl, ku,
                                       alpha0, alpha1, d1, alpha2, d2,
                                       a, lda, x, incx,
                                       beta,   y, incy);
}

int
suzerain_gbiddmv_d(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0,
        const double alpha1,
        const double *d1,
        const double alpha2,
        const double *d2,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
#define ARGS trans, n, alpha0, alpha1, d1, alpha2, d2, \
            a, lda, x, incx, beta, y, incy
            case  0: return suzerain_gbiddmv_internal_d0(ARGS);
            case  1: return suzerain_gbiddmv_internal_d1(ARGS);
            case  2: return suzerain_gbiddmv_internal_d2(ARGS);
            case  3: return suzerain_gbiddmv_internal_d3(ARGS);
            case  4: return suzerain_gbiddmv_internal_d4(ARGS);
            case  5: return suzerain_gbiddmv_internal_d5(ARGS);
            case  6: return suzerain_gbiddmv_internal_d6(ARGS);
            case  7: return suzerain_gbiddmv_internal_d7(ARGS);
            case  8: return suzerain_gbiddmv_internal_d8(ARGS);
            case  9: return suzerain_gbiddmv_internal_d9(ARGS);
            case 10: return suzerain_gbiddmv_internal_d10(ARGS);
            case 11: return suzerain_gbiddmv_internal_d11(ARGS);
            case 12: return suzerain_gbiddmv_internal_d12(ARGS);
            case 13: return suzerain_gbiddmv_internal_d13(ARGS);
            case 14: return suzerain_gbiddmv_internal_d14(ARGS);
            case 15: return suzerain_gbiddmv_internal_d15(ARGS);
#undef ARGS
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbiddmv_internal_d(trans, n, kl, ku,
                                       alpha0, alpha1, d1, alpha2, d2,
                                       a, lda, x, incx,
                                       beta,   y, incy);
}

int
suzerain_gbiddmv_sc(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha0[2],
        const float alpha1[2],
        const float *d1,
        const float alpha2[2],
        const float *d2,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy)
{
    float _Complex alpha0_c, alpha1_c, alpha2_c, beta_c;
    memcpy(&alpha0_c, alpha0, sizeof(float _Complex));
    memcpy(&alpha1_c, alpha1, sizeof(float _Complex));
    memcpy(&alpha2_c, alpha2, sizeof(float _Complex));
    memcpy(&beta_c,   beta,   sizeof(float _Complex));

    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
#define ARGS trans, n, alpha0_c, alpha1_c, (void *) d1,                 \
            alpha2_c, (void *) d2,                                      \
            (void *) a, lda, (void *) x, incx, beta_c, (void *) y, incy
            case  0: return suzerain_gbiddmv_internal_sc0(ARGS);
            case  1: return suzerain_gbiddmv_internal_sc1(ARGS);
            case  2: return suzerain_gbiddmv_internal_sc2(ARGS);
            case  3: return suzerain_gbiddmv_internal_sc3(ARGS);
            case  4: return suzerain_gbiddmv_internal_sc4(ARGS);
            case  5: return suzerain_gbiddmv_internal_sc5(ARGS);
            case  6: return suzerain_gbiddmv_internal_sc6(ARGS);
            case  7: return suzerain_gbiddmv_internal_sc7(ARGS);
            case  8: return suzerain_gbiddmv_internal_sc8(ARGS);
            case  9: return suzerain_gbiddmv_internal_sc9(ARGS);
            case 10: return suzerain_gbiddmv_internal_sc10(ARGS);
            case 11: return suzerain_gbiddmv_internal_sc11(ARGS);
            case 12: return suzerain_gbiddmv_internal_sc12(ARGS);
            case 13: return suzerain_gbiddmv_internal_sc13(ARGS);
            case 14: return suzerain_gbiddmv_internal_sc14(ARGS);
            case 15: return suzerain_gbiddmv_internal_sc15(ARGS);
#undef ARGS
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbiddmv_internal_sc(
            trans, n, kl, ku,
            alpha0_c, alpha1_c, (void *) d1, alpha2_c, (void *) d2,
            (void *) a, lda, (void *) x, incx,
            beta_c,          (void *) y, incy);
}

int
suzerain_gbiddmv_dz(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha0[2],
        const double alpha1[2],
        const double *d1,
        const double alpha2[2],
        const double *d2,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy)
{
    double _Complex alpha0_c, alpha1_c, alpha2_c, beta_c;
    memcpy(&alpha0_c, alpha0, sizeof(double _Complex));
    memcpy(&alpha1_c, alpha1, sizeof(double _Complex));
    memcpy(&alpha2_c, alpha2, sizeof(double _Complex));
    memcpy(&beta_c,   beta,   sizeof(double _Complex));

    // Dispatch to fixed bandwidth specialization for small bandwidth...
    if (kl == ku) {
        switch (kl) {
#define ARGS trans, n, alpha0_c, alpha1_c, (void *) d1,                 \
            alpha2_c, (void *) d2,                                      \
            (void *) a, lda, (void *) x, incx, beta_c, (void *) y, incy
            case  0: return suzerain_gbiddmv_internal_dz0(ARGS);
            case  1: return suzerain_gbiddmv_internal_dz1(ARGS);
            case  2: return suzerain_gbiddmv_internal_dz2(ARGS);
            case  3: return suzerain_gbiddmv_internal_dz3(ARGS);
            case  4: return suzerain_gbiddmv_internal_dz4(ARGS);
            case  5: return suzerain_gbiddmv_internal_dz5(ARGS);
            case  6: return suzerain_gbiddmv_internal_dz6(ARGS);
            case  7: return suzerain_gbiddmv_internal_dz7(ARGS);
            case  8: return suzerain_gbiddmv_internal_dz8(ARGS);
            case  9: return suzerain_gbiddmv_internal_dz9(ARGS);
            case 10: return suzerain_gbiddmv_internal_dz10(ARGS);
            case 11: return suzerain_gbiddmv_internal_dz11(ARGS);
            case 12: return suzerain_gbiddmv_internal_dz12(ARGS);
            case 13: return suzerain_gbiddmv_internal_dz13(ARGS);
            case 14: return suzerain_gbiddmv_internal_dz14(ARGS);
            case 15: return suzerain_gbiddmv_internal_dz15(ARGS);
#undef ARGS
        }
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_gbiddmv_internal_dz(
            trans, n, kl, ku,
            alpha0_c, alpha1_c, (void *) d1, alpha2_c, (void *) d2,
            (void *) a, lda, (void *) x, incx,
            beta_c,          (void *) y, incy);
}

#else

// ------------------------------------------------------------
// Generate fixed bandwidth routines using BOOST_PP_ITERATION()
// ------------------------------------------------------------
#define kl BOOST_PP_ITERATION()
#define ku BOOST_PP_ITERATION()

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbiddmv_internal_s, \
                                       BOOST_PP_ITERATION())
#define GBIDDMV_COMPONENT float
#define GBIDDMV_SCALAR    float
#define GBIDDMV_KL
#define GBIDDMV_KU
#include "gbiddmv.def"

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbiddmv_internal_d, \
                                       BOOST_PP_ITERATION())
#define GBIDDMV_COMPONENT double
#define GBIDDMV_SCALAR    double
#define GBIDDMV_KL
#define GBIDDMV_KU
#include "gbiddmv.def"

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbiddmv_internal_sc, \
                                       BOOST_PP_ITERATION())
#define GBIDDMV_COMPONENT float
#define GBIDDMV_SCALAR    float _Complex
#define GBIDDMV_KL
#define GBIDDMV_KU
#include "gbiddmv.def"

#define GBIDDMV_STATIC    static
#define GBIDDMV_FUNCTION  BOOST_PP_CAT(suzerain_gbiddmv_internal_dz, \
                                       BOOST_PP_ITERATION())
#define GBIDDMV_COMPONENT double
#define GBIDDMV_SCALAR    double _Complex
#define GBIDDMV_KL
#define GBIDDMV_KU
#include "gbiddmv.def"

#undef kl
#undef ku

#endif
