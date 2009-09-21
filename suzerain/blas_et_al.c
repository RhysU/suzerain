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

#include <config.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MKL
#include <mkl_types.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

#include <suzerain/blas_et_al.h>
#include <suzerain/macros.h>

void *
suzerain_blas_malloc(size_t size)
{
#ifdef HAVE_MKL
    /* We do not use MKL_malloc to avoid later needing MKL_free calls. */
    /* Align at 16-byte boundaries per MKL user guide section 8. */
    const size_t alignment = 16;
#else
#error "Sanity failure"
#endif

    void * p = NULL;

    const int status = posix_memalign(&p, alignment, size);

    switch (status) {
        case 0:
            return p;
        case ENOMEM:
            return NULL;
        case EINVAL:
        default:
            fprintf(stderr, "Fatal posix_memalign error at %s:%d-- %s\n",
                    __FILE__, __LINE__, strerror(status));
            abort();
    }
}

void *
suzerain_blas_calloc(size_t nmemb, size_t size)
{
    const size_t total_bytes = nmemb * size;
    void * p = suzerain_blas_malloc(total_bytes);
    if (p != NULL) {
        memset(p, 0, total_bytes);
    }
    return p;
}

void
suzerain_blas_dcopy(
        const int n,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef HAVE_MKL
    const MKL_INT _n    = n;
    const MKL_INT _incx = incx;
    const MKL_INT _incy = incy;

    dcopy(&_n, x, &_incx, y, &_incy);
#else
#error "Sanity failure"
#endif
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

    return ddot(&_n, x, &_incx, y, &_incy);
#else
#error "Sanity failure"
#endif
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

    return dasum(&_n, x, &_incx);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_daxpy(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef HAVE_MKL
    const MKL_INT _n    = n;
    const MKL_INT _incx = incx;
    const MKL_INT _incy = incy;

    daxpy(&_n, &alpha, x, &_incx, y, &_incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_daxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy)
{
#ifdef HAVE_MKL
    /* Simulate daxpby since MKL lacks the routine. */
    if (SUZERAIN_UNLIKELY((alpha == 0.0 && beta == 1.0) || n <= 0)) return;

    const MKL_INT _n    = n;
    const MKL_INT _incx = incx;
    const MKL_INT _incy = incy;

    dscal(&_n, &beta, y, &_incy);
    daxpy(&_n, &alpha, x, &_incx, y, &_incy);
#else
#error "Sanity failure"
#endif
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
    const MKL_INT _m    = m;
    const MKL_INT _n    = n;
    const MKL_INT _kl   = kl;
    const MKL_INT _ku   = ku;
    const MKL_INT _lda  = lda;
    const MKL_INT _incx = incx;
    const MKL_INT _incy = incy;

    dgbmv(&_trans, &_m, &_n, &_kl, &_ku, &alpha, a, &_lda,
          x, &_incx, &beta, y, &_incy);
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *a,
        const int lda,
        const double beta,
        double *b,
        const int ldb)
{
#ifdef HAVE_MKL
    /* Simulate dgb_acc since MKL lacks the routine. */
    if (SUZERAIN_UNLIKELY((alpha == 0.0 && beta == 1.0) || m <= 0)) return;

    const int veclength = ku + 1 + kl;
    const double * const bj_end = b + n *ldb;
    const double *aj;
    double       *bj;

    for (aj = a, bj = b; bj < bj_end; aj += lda, bj += ldb) {
        suzerain_blas_daxpby(veclength, alpha, aj, 1, beta, bj, 1);
    }
#else
#error "Sanity failure"
#endif
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

    dgbtrf(&_m, &_n, &_kl, &_ku, ab, &_ldab, ipiv, &_info);

    return _info;
#else
#error "Sanity failure"
#endif
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

    dgbtrs(&_trans, &_n, &_kl, &_ku, &_nrhs,
           ab, &_ldab, ipiv, b, &_ldb, &_info);

    return _info;
#else
#error "Sanity failure"
#endif
}
