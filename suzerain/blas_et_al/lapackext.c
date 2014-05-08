/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
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
 * @copydoc lapackext.h
 */

#include <suzerain/blas_et_al/lapackext.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al/blasext.h>
#include <suzerain/blas_et_al/blas.h>
#include <suzerain/blas_et_al/lapack.h>

// Shorthand
#define UNLIKELY(expr) SUZERAIN_UNLIKELY(expr)

// Many of the short methods have "inline" though their declarations do not.
// to allow inlining them later within this particular translation unit.

// Thank you captain obvious...
#pragma warning(disable:981)

inline int
suzerain_lapackext_sgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        int *ipiv,
        float *b,
        const int ldb)
{
    int info = 0;
    if (toupper(*fact) == 'N') {
        info  = suzerain_lapack_sgbtrf(n, n, kl, ku, ab, ldab, ipiv);
        *fact = 'F';
    }
    if (!info) {
        info = suzerain_lapack_sgbtrs(
                trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
    }
    return info;
}

inline int
suzerain_lapackext_dgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        int *ipiv,
        double *b,
        const int ldb)
{
    int info = 0;
    if (toupper(*fact) == 'N') {
        info = suzerain_lapack_dgbtrf(n, n, kl, ku, ab, ldab, ipiv);
        *fact = 'F';
    }
    if (!info) {
        info = suzerain_lapack_dgbtrs(
                trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
    }
    return info;
}

inline int
suzerain_lapackext_cgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        int *ipiv,
        complex_float *b,
        const int ldb)
{
    int info = 0;
    if (toupper(*fact) == 'N') {
        info = suzerain_lapack_cgbtrf(n, n, kl, ku, ab, ldab, ipiv);
        *fact = 'F';
    }
    if (!info) {
        info = suzerain_lapack_cgbtrs(
                trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
    }
    return info;
}

inline int
suzerain_lapackext_zgbsv(
        char * const fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        int *ipiv,
        complex_double *b,
        const int ldb)
{
    int info = 0;
    if (toupper(*fact) == 'N') {
        info = suzerain_lapack_zgbtrf(n, n, kl, ku, ab, ldab, ipiv);
        *fact = 'F';
    }
    if (!info) {
        info = suzerain_lapack_zgbtrs(
                trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb);
    }
    return info;
}

// Generate suzerain_lapackext_dsgbsvx(...) from template dsgbsvx.def
#define SINGLE    float                /* Type of single precision   */
#define SCHAR     s                    /* Abbreviation for single    */
#define DOUBLE    double               /* Type of double precision   */
#define DCHAR     d                    /* Abbreviation for double    */
#define NORM      double               /* Type of norm-based results */
#define DNRM2     suzerain_blas_dnrm2  /* Norm mapping to NORM       */
#include "dsgbsvx.def"

// Generate suzerain_lapackext_zcgbsvx(...) from template dsgbsvx.def
#define SINGLE    complex_float        /* Type of single precision   */
#define SCHAR     c                    /* Abbreviation for single    */
#define DOUBLE    complex_double       /* Type of double precision   */
#define DCHAR     z                    /* Abbreviation for double    */
#define NORM      double               /* Type of norm-based results */
#define DNRM2     suzerain_blas_dznrm2 /* Norm mapping to NORM       */
#include "dsgbsvx.def"
