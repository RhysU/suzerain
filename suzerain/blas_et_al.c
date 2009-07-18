/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * blas_et_al.c: wraps external implementations of BLAS, LAPACK, et al.
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BLAS_ET_AL_H__
#define __SUZERAIN_BLAS_ET_AL_H__

#include "config.h"

#ifdef HAVE_MKL
#include <mkl_types.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

#include <suzerain/blas_et_al.h>

void
suzerain_blas_dcopy(
        const int n,
        double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef HAVE_MKL
    MKL_INT _n    = n;
    MKL_INT _incx = incx;
    MKL_INT _incy = incy;
#else
#error "Sanity failure"
#endif

    dcopy(&_n, x, &_incx, y, &_incy);
}

double
suzerain_blas_ddot(
        const int n,
        const double *x,
        const int incx,
        const double *y,
        const int incy)
{
#ifdef HAVE_MKL
    MKL_INT _n    = n;
    MKL_INT _incx = incx;
    MKL_INT _incy = incy;
#else
#error "Sanity failure"
#endif

    return ddot(&_n, x, &_incx, y, &_incy);
}

double
suzerain_blas_dasum(
        const int n,
        const double *x,
        const int incx)
{
#ifdef HAVE_MKL
    MKL_INT _n    = n;
    MKL_INT _incx = incx;
#else
#error "Sanity failure"
#endif

    return dasum(&_n, x, &_incx);
}

void
suzerain_blas_dgbmv(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
    char _trans   = trans;
#ifdef HAVE_MKL
    MKL_INT _m    = m;
    MKL_INT _n    = n;
    MKL_INT _kl   = kl;
    MKL_INT _ku   = ku;
    MKL_INT _lda  = lda;
    MKL_INT _incx = incx;
    MKL_INT _incy = incy;
#else
#error "Sanity failure"
#endif
    double  _alpha = alpha;
    double  _beta  = beta;

    dgbmv(&_trans, &_m, &_n, &_kl, &_ku, &_alpha, a, &_lda,
          x, &_incx, &_beta, y, &_incy);
}

int
suzerain_lapack_dgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double *ab,
        const int ldab,
        int *ipiv)
{
#ifdef HAVE_MKL
    MKL_INT _m    = m;
    MKL_INT _n    = n;
    MKL_INT _kl   = kl;
    MKL_INT _ku   = ku;
    MKL_INT _ldab = ldab;
    MKL_INT _info = 0;
#else
#error "Sanity failure"
#endif

    dgbtrf(&_m, &_n, &_kl, &_ku, ab, &_ldab, ipiv, &_info);

    return _info;
}

int
suzerain_lapack_dgbtrs(
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
    char _trans   = trans;
#ifdef HAVE_MKL
    MKL_INT _n    = n;
    MKL_INT _kl   = kl;
    MKL_INT _ku   = ku;
    MKL_INT _nrhs = nrhs;
    MKL_INT _ldab = ldab;
    MKL_INT _ldb  = ldb;
    MKL_INT _info = 0;
#else
#error "Sanity failure"
#endif

    dgbtrs(&_trans, &_n, &_kl, &_ku, &_nrhs,
           ab, &_ldab, ipiv, b, &_ldb, &_info);

    return _info;
}

#endif /* __SUZERAIN_BLAS_ET_AL_H__ */
