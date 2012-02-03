/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * sbmv.c: provides custom, BLAS-like general band matrix-vector operations
 * $Id$
 */

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
#include <suzerain/sbmv.h>

#pragma warning(disable:1418 1572 2259)

// ---------------------------------------------------------
// Iterate on this file to generate fixed bandwidth routines
// ---------------------------------------------------------
#define BOOST_PP_ITERATION_LIMITS (0, 15)
#define BOOST_PP_FILENAME_1 <suzerain/sbmv.c>
#include BOOST_PP_ITERATE()

// -----------------------------------
// Generate general bandwidth routines
// -----------------------------------

#define SBMV_STATIC    static
#define SBMV_FUNCTION  suzerain_sbmv_internal_s
#define SBMV_COMPONENT float
#define SBMV_SCALAR    float
#define SBMV_K         const int k,
#include "sbmv.def"

#define SBMV_STATIC    static
#define SBMV_FUNCTION  suzerain_sbmv_internal_d
#define SBMV_COMPONENT double
#define SBMV_SCALAR    double
#define SBMV_K         const int k,
#include "sbmv.def"

#define SBMV_STATIC    static
#define SBMV_FUNCTION  suzerain_sbmv_internal_sc
#define SBMV_COMPONENT float
#define SBMV_SCALAR    complex_float
#define SBMV_K         const int k,
#include "sbmv.def"

#define SBMV_STATIC    static
#define SBMV_FUNCTION  suzerain_sbmv_internal_dz
#define SBMV_COMPONENT double
#define SBMV_SCALAR    complex_double
#define SBMV_K         const int k,
#include "sbmv.def"

// ------------------------------------------------------------------
// Provide externally callable logic dispatching to internal routines
// ------------------------------------------------------------------

int
suzerain_sbmv_s(
        const char uplo,
        const int n,
        const int k,
        const float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    switch (k) {
#define ARGS uplo, n, alpha, a, lda, x, incx, beta, y, incy
        case  0: return suzerain_sbmv_internal_s0(ARGS);
        case  1: return suzerain_sbmv_internal_s1(ARGS);
        case  2: return suzerain_sbmv_internal_s2(ARGS);
        case  3: return suzerain_sbmv_internal_s3(ARGS);
        case  4: return suzerain_sbmv_internal_s4(ARGS);
        case  5: return suzerain_sbmv_internal_s5(ARGS);
        case  6: return suzerain_sbmv_internal_s6(ARGS);
        case  7: return suzerain_sbmv_internal_s7(ARGS);
        case  8: return suzerain_sbmv_internal_s8(ARGS);
        case  9: return suzerain_sbmv_internal_s9(ARGS);
        case 10: return suzerain_sbmv_internal_s10(ARGS);
        case 11: return suzerain_sbmv_internal_s11(ARGS);
        case 12: return suzerain_sbmv_internal_s12(ARGS);
        case 13: return suzerain_sbmv_internal_s13(ARGS);
        case 14: return suzerain_sbmv_internal_s14(ARGS);
        case 15: return suzerain_sbmv_internal_s15(ARGS);
#undef ARGS
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_s(uplo, n, k,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
suzerain_sbmv_d(
        const char uplo,
        const int n,
        const int k,
        const double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    switch (k) {
#define ARGS uplo, n, alpha, a, lda, x, incx, beta, y, incy
        case  0: return suzerain_sbmv_internal_d0(ARGS);
        case  1: return suzerain_sbmv_internal_d1(ARGS);
        case  2: return suzerain_sbmv_internal_d2(ARGS);
        case  3: return suzerain_sbmv_internal_d3(ARGS);
        case  4: return suzerain_sbmv_internal_d4(ARGS);
        case  5: return suzerain_sbmv_internal_d5(ARGS);
        case  6: return suzerain_sbmv_internal_d6(ARGS);
        case  7: return suzerain_sbmv_internal_d7(ARGS);
        case  8: return suzerain_sbmv_internal_d8(ARGS);
        case  9: return suzerain_sbmv_internal_d9(ARGS);
        case 10: return suzerain_sbmv_internal_d10(ARGS);
        case 11: return suzerain_sbmv_internal_d11(ARGS);
        case 12: return suzerain_sbmv_internal_d12(ARGS);
        case 13: return suzerain_sbmv_internal_d13(ARGS);
        case 14: return suzerain_sbmv_internal_d14(ARGS);
        case 15: return suzerain_sbmv_internal_d15(ARGS);
#undef ARGS
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_d(uplo, n, k,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
suzerain_sbmv_sc(
        const char uplo,
        const int n,
        const int k,
        const complex_float alpha,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    switch (k) {
#define ARGS uplo, n, alpha, a, lda, x, incx, beta, y, incy
        case  0: return suzerain_sbmv_internal_sc0(ARGS);
        case  1: return suzerain_sbmv_internal_sc1(ARGS);
        case  2: return suzerain_sbmv_internal_sc2(ARGS);
        case  3: return suzerain_sbmv_internal_sc3(ARGS);
        case  4: return suzerain_sbmv_internal_sc4(ARGS);
        case  5: return suzerain_sbmv_internal_sc5(ARGS);
        case  6: return suzerain_sbmv_internal_sc6(ARGS);
        case  7: return suzerain_sbmv_internal_sc7(ARGS);
        case  8: return suzerain_sbmv_internal_sc8(ARGS);
        case  9: return suzerain_sbmv_internal_sc9(ARGS);
        case 10: return suzerain_sbmv_internal_sc10(ARGS);
        case 11: return suzerain_sbmv_internal_sc11(ARGS);
        case 12: return suzerain_sbmv_internal_sc12(ARGS);
        case 13: return suzerain_sbmv_internal_sc13(ARGS);
        case 14: return suzerain_sbmv_internal_sc14(ARGS);
        case 15: return suzerain_sbmv_internal_sc15(ARGS);
#undef ARGS
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_sc(uplo, n, k,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
}

int
suzerain_sbmv_dz(
        const char uplo,
        const int n,
        const int k,
        const complex_double alpha,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    switch (k) {
#define ARGS uplo, n, alpha, a, lda, x, incx, beta, y, incy
        case  0: return suzerain_sbmv_internal_dz0(ARGS);
        case  1: return suzerain_sbmv_internal_dz1(ARGS);
        case  2: return suzerain_sbmv_internal_dz2(ARGS);
        case  3: return suzerain_sbmv_internal_dz3(ARGS);
        case  4: return suzerain_sbmv_internal_dz4(ARGS);
        case  5: return suzerain_sbmv_internal_dz5(ARGS);
        case  6: return suzerain_sbmv_internal_dz6(ARGS);
        case  7: return suzerain_sbmv_internal_dz7(ARGS);
        case  8: return suzerain_sbmv_internal_dz8(ARGS);
        case  9: return suzerain_sbmv_internal_dz9(ARGS);
        case 10: return suzerain_sbmv_internal_dz10(ARGS);
        case 11: return suzerain_sbmv_internal_dz11(ARGS);
        case 12: return suzerain_sbmv_internal_dz12(ARGS);
        case 13: return suzerain_sbmv_internal_dz13(ARGS);
        case 14: return suzerain_sbmv_internal_dz14(ARGS);
        case 15: return suzerain_sbmv_internal_dz15(ARGS);
#undef ARGS
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_dz(uplo, n, k,
                                     alpha, a, lda, x, incx,
                                     beta,          y, incy);
}

#else

// ------------------------------------------------------------
// Generate fixed bandwidth routines using BOOST_PP_ITERATION()
// ------------------------------------------------------------
#define k BOOST_PP_ITERATION()

#define SBMV_STATIC    static
#define SBMV_FUNCTION  BOOST_PP_CAT(suzerain_sbmv_internal_s, \
                                    BOOST_PP_ITERATION())
#define SBMV_COMPONENT float
#define SBMV_SCALAR    float
#define SBMV_K
#include "sbmv.def"

#define SBMV_STATIC    static
#define SBMV_FUNCTION  BOOST_PP_CAT(suzerain_sbmv_internal_d, \
                                    BOOST_PP_ITERATION())
#define SBMV_COMPONENT double
#define SBMV_SCALAR    double
#define SBMV_K
#include "sbmv.def"

#define SBMV_STATIC    static
#define SBMV_FUNCTION  BOOST_PP_CAT(suzerain_sbmv_internal_sc, \
                                    BOOST_PP_ITERATION())
#define SBMV_COMPONENT float
#define SBMV_SCALAR    complex_float
#define SBMV_K
#include "sbmv.def"

#define SBMV_STATIC    static
#define SBMV_FUNCTION  BOOST_PP_CAT(suzerain_sbmv_internal_dz, \
                                    BOOST_PP_ITERATION())
#define SBMV_COMPONENT double
#define SBMV_SCALAR    complex_double
#define SBMV_K
#include "sbmv.def"

#undef k

#endif
