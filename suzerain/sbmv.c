/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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
 * @copydoc sbmv.h
 */

// File iteration is used to generate bandwidth-specific routines.  See
// http://www.boost.org/libs/preprocessor/doc/topics/file_iteration.html for an
// overview.  We iterate first to generate internal, static routines and then
// continue on to the non-iterated logic.

#if !defined(BOOST_PP_IS_ITERATING) || !BOOST_PP_IS_ITERATING

#include <suzerain/sbmv.h>

#include <suzerain/common.h>

#pragma warning(disable:1418 1572 2259)

// ---------------------------------------------------------
// Iterate on this file to generate fixed bandwidth routines
// ---------------------------------------------------------
#define FIXEDBW_LOWER 0
#define FIXEDBW_UPPER 15
#define BOOST_PP_ITERATION_LIMITS (FIXEDBW_LOWER, FIXEDBW_UPPER)
#define BOOST_PP_FILENAME_1 <suzerain/sbmv.c>
#include BOOST_PP_ITERATE()

// -----------------------------------
// Generate general bandwidth routines
// -----------------------------------

#define STATIC   static
#define FUNCTION suzerain_sbmv_internal_s
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   float
#define K        const int k,
#include "sbmv.def"

#define STATIC   static
#define FUNCTION suzerain_sbmv_internal_ssc
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   complex_float
#define K        const int k,
#include "sbmv.def"

#define STATIC   static
#define FUNCTION suzerain_sbmv_internal_scc
#define TYPE_A   float
#define TYPE_X   complex_float
#define TYPE_Y   complex_float
#define K        const int k,
#include "sbmv.def"

#define STATIC   static
#define FUNCTION suzerain_sbmv_internal_d
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   double
#define K        const int k,
#include "sbmv.def"

#define STATIC   static
#define FUNCTION suzerain_sbmv_internal_ddz
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   complex_double
#define K        const int k,
#include "sbmv.def"

#define STATIC   static
#define FUNCTION suzerain_sbmv_internal_dzz
#define TYPE_A   double
#define TYPE_X   complex_double
#define TYPE_Y   complex_double
#define K        const int k,
#include "sbmv.def"

// ------------------------------------------------------------------
// Provide externally callable logic dispatching to internal routines
// ------------------------------------------------------------------

#define FIXEDBW_CASE(z,num,prefix)                       \
    case num: return BOOST_PP_CAT(prefix,num)(           \
        uplo, n, alpha, a, lda, x, incx, beta, y, incy);

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
        BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                FIXEDBW_CASE, suzerain_sbmv_internal_s)
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_s(uplo, n, k,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
suzerain_sbmv_ssc(
        const char uplo,
        const int n,
        const int k,
        const complex_float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    switch (k) {
        BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                FIXEDBW_CASE, suzerain_sbmv_internal_ssc)
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_ssc(uplo, n, k,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
}

int
suzerain_sbmv_scc(
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
        BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                FIXEDBW_CASE, suzerain_sbmv_internal_scc)
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_scc(uplo, n, k,
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
        BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                FIXEDBW_CASE, suzerain_sbmv_internal_d)
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_d(uplo, n, k,
                                    alpha, a, lda, x, incx,
                                    beta,          y, incy);
}

int
suzerain_sbmv_ddz(
        const char uplo,
        const int n,
        const int k,
        const complex_double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy)
{
    // Dispatch to fixed bandwidth specialization for small bandwidth...
    switch (k) {
        BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                FIXEDBW_CASE, suzerain_sbmv_internal_ddz)
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_ddz(uplo, n, k,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
}

int
suzerain_sbmv_dzz(
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
        BOOST_PP_REPEAT_FROM_TO(FIXEDBW_LOWER, BOOST_PP_INC(FIXEDBW_UPPER),
                                FIXEDBW_CASE, suzerain_sbmv_internal_dzz)
    }

    // ...otherwise employ a general bandwidth implementation
    return suzerain_sbmv_internal_dzz(uplo, n, k,
                                      alpha, a, lda, x, incx,
                                      beta,          y, incy);
}

#else

// ------------------------------------------------------------
// Generate fixed bandwidth routines using BOOST_PP_ITERATION()
// ------------------------------------------------------------
#define k BOOST_PP_ITERATION()

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_sbmv_internal_s, BOOST_PP_ITERATION())
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   float
#define K
#include "sbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_sbmv_internal_ssc, BOOST_PP_ITERATION())
#define TYPE_A   float
#define TYPE_X   float
#define TYPE_Y   complex_float
#define K
#include "sbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_sbmv_internal_scc, BOOST_PP_ITERATION())
#define TYPE_A   float
#define TYPE_X   complex_float
#define TYPE_Y   complex_float
#define K
#include "sbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_sbmv_internal_d, BOOST_PP_ITERATION())
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   double
#define K
#include "sbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_sbmv_internal_ddz, BOOST_PP_ITERATION())
#define TYPE_A   double
#define TYPE_X   double
#define TYPE_Y   complex_double
#define K
#include "sbmv.def"

#define STATIC   static
#define FUNCTION BOOST_PP_CAT(suzerain_sbmv_internal_dzz, BOOST_PP_ITERATION())
#define TYPE_A   double
#define TYPE_X   complex_double
#define TYPE_Y   complex_double
#define K
#include "sbmv.def"

#undef k

#endif
