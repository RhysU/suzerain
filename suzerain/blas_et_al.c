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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop

#ifdef SUZERAIN_HAVE_MKL
#include <mkl_types.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#else
#error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

#include <suzerain/blas_et_al.h>

void *
suzerain_blas_malloc(size_t size)
{
#ifdef SUZERAIN_HAVE_MKL
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
suzerain_blas_sswap(
        const int n,
        float *x,
        const int incx,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        return sswap(&n, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        return sswap(&_n, x, &_incx, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dswap(
        const int n,
        double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        return dswap(&n, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        return dswap(&_n, x, &_incx, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_scopy(
        const int n,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        scopy(&n, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        scopy(&_n, x, &_incx, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dcopy(
        const int n,
        const double *x,
        const int incx,
        double *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        dcopy(&n, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        dcopy(&_n, x, &_incx, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}


float
suzerain_blas_sdot(
        const int n,
        const float *x,
        const int incx,
        const float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        return sdot(&n, x, &incx, y, &incy);
    } else {
        MKL_INT _n    = n;
        MKL_INT _incx = incx;
        MKL_INT _incy = incy;

        return sdot(&_n, x, &_incx, y, &_incy);
    }
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
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        return ddot(&n, x, &incx, y, &incy);
    } else {
        MKL_INT _n    = n;
        MKL_INT _incx = incx;
        MKL_INT _incy = incy;

        return ddot(&_n, x, &_incx, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}

float
suzerain_blas_sasum(
        const int n,
        const float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        return sasum(&n, x, &incx);
    } else {
        MKL_INT _n    = n;
        MKL_INT _incx = incx;

        return sasum(&_n, x, &_incx);
    }
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
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        return dasum(&n, x, &incx);
    } else {
        MKL_INT _n    = n;
        MKL_INT _incx = incx;

        return dasum(&_n, x, &_incx);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_saxpy(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        saxpy(&n, &alpha, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        saxpy(&_n, &alpha, x, &_incx, y, &_incy);
    }
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
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        daxpy(&n, &alpha, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        daxpy(&_n, &alpha, x, &_incx, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_saxpby(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate saxpby since MKL lacks the routine. */
    if (SUZERAIN_UNLIKELY((alpha == 0.0f && beta == 1.0f) || n <= 0)) return;

    if (sizeof(MKL_INT) == sizeof(int)) {
        if (beta != 1.0f)
            sscal(&n, &beta, y, &incy);
        saxpy(&n, &alpha, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        if (beta != 1.0f)
            sscal(&_n, &beta, y, &_incy);
        saxpy(&_n, &alpha, x, &_incx, y, &_incy);
    }
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
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate daxpby since MKL lacks the routine. */
    if (SUZERAIN_UNLIKELY((alpha == 0.0 && beta == 1.0) || n <= 0)) return;

    if (sizeof(MKL_INT) == sizeof(int)) {
        if (beta != 1.0)
            dscal(&n, &beta, y, &incy);
        daxpy(&n, &alpha, x, &incx, y, &incy);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        if (beta != 1.0)
            dscal(&_n, &beta, y, &_incy);
        daxpy(&_n, &alpha, x, &_incx, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_swaxpby(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        const float *y,
        const int incy,
        float *w,
        const int incw)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate swaxpby since MKL lacks the routine. */

    if (sizeof(MKL_INT) == sizeof(int)) {
        scopy(&n, y, &incy, w, &incw);
        if (beta != 1.0f)
            sscal(&n, &beta, w, &incw);
        saxpy(&n, &alpha, x, &incx, w, &incw);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;
        const MKL_INT _incw = incw;

        scopy(&_n, y, &_incy, w, &_incw);
        if (beta != 1.0f)
            sscal(&_n, &beta, w, &_incw);
        saxpy(&_n, &alpha, x, &_incx, w, &_incw);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dwaxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        const double *y,
        const int incy,
        double *w,
        const int incw)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate dwaxpby since MKL lacks the routine. */

    if (sizeof(MKL_INT) == sizeof(int)) {
        dcopy(&n, y, &incy, w, &incw);
        if (beta != 1.0)
            dscal(&n, &beta, w, &incw);
        daxpy(&n, &alpha, x, &incx, w, &incw);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;
        const MKL_INT _incw = incw;

        dcopy(&_n, y, &_incy, w, &_incw);
        if (beta != 1.0)
            dscal(&_n, &beta, w, &_incw);
        daxpy(&_n, &alpha, x, &_incx, w, &_incw);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_sscal(
        const int n,
        const float alpha,
        float *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        sscal(&n, &alpha, x, &incx);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;

        sscal(&_n, &alpha, x, &_incx);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_dscal(
        const int n,
        const double alpha,
        double *x,
        const int incx)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        dscal(&n, &alpha, x, &incx);
    } else {
        const MKL_INT _n    = n;
        const MKL_INT _incx = incx;

        dscal(&_n, &alpha, x, &_incx);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_sgbmv(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy)
{
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        sgbmv(&trans, &m, &n, &kl, &ku, &alpha, a, &lda,
              x, &incx, &beta, y, &incy);
    } else {
        const MKL_INT _m    = m;
        const MKL_INT _n    = n;
        const MKL_INT _kl   = kl;
        const MKL_INT _ku   = ku;
        const MKL_INT _lda  = lda;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        sgbmv(&trans, &_m, &_n, &_kl, &_ku, &alpha, a, &_lda,
              x, &_incx, &beta, y, &_incy);
    }
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
#ifdef SUZERAIN_HAVE_MKL
    if (sizeof(MKL_INT) == sizeof(int)) {
        dgbmv(&trans, &m, &n, &kl, &ku, &alpha, a, &lda,
              x, &incx, &beta, y, &incy);
    } else {
        const MKL_INT _m    = m;
        const MKL_INT _n    = n;
        const MKL_INT _kl   = kl;
        const MKL_INT _ku   = ku;
        const MKL_INT _lda  = lda;
        const MKL_INT _incx = incx;
        const MKL_INT _incy = incy;

        dgbmv(&trans, &_m, &_n, &_kl, &_ku, &alpha, a, &_lda,
              x, &_incx, &beta, y, &_incy);
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blas_sgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *a,
        const int lda,
        const float beta,
        float *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    /* Simulate sgb_acc since MKL lacks the routine. */
    if (SUZERAIN_UNLIKELY((alpha == 0.0f && beta == 1.0f) || m <= 0)) return;

    const int veclength = ku + 1 + kl;
    const float * const bj_end = b + n *ldb;
    const float *aj;
    float       *bj;

    for (aj = a, bj = b; bj < bj_end; aj += lda, bj += ldb) {
        suzerain_blas_saxpby(veclength, alpha, aj, 1, beta, bj, 1);
    }
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
#ifdef SUZERAIN_HAVE_MKL
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
suzerain_lapack_sgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float *ab,
        const int ldab,
        int *ipiv)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    if (sizeof(MKL_INT) == sizeof(int)) {
        int info = 0;
        sgbtrf((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
               ab, (int *)&ldab, ipiv, &info);
        return info;
    } else {
        MKL_INT _m    = m;
        MKL_INT _n    = n;
        MKL_INT _kl   = kl;
        MKL_INT _ku   = ku;
        MKL_INT _ldab = ldab;

        MKL_INT _info = 0;
        // FIXME: ipiv's contents are incompatible with LAPACK here
        sgbtrf(&_m, &_n, &_kl, &_ku, ab, &_ldab, ipiv, &_info);
        return _info;
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
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    if (sizeof(MKL_INT) == sizeof(int)) {
        int info = 0;
        dgbtrf((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
               ab, (int*)&ldab, ipiv, &info);
        return info;
    } else {
        MKL_INT _m    = m;
        MKL_INT _n    = n;
        MKL_INT _kl   = kl;
        MKL_INT _ku   = ku;
        MKL_INT _ldab = ldab;

        MKL_INT _info = 0;
        // FIXME: ipiv's contents are incompatible with LAPACK here
        dgbtrf(&_m, &_n, &_kl, &_ku, ab, &_ldab, ipiv, &_info);
        return _info;
    }
#else
#error "Sanity failure"
#endif
}

int
suzerain_lapack_sgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const float *ab,
        const int ldab,
        const int *ipiv,
        float *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    if (sizeof(MKL_INT) == sizeof(int)) {
        int info = 0;
        sgbtrs((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
            (float *)ab, (int*)&ldab, (int *)ipiv, b, (int*)&ldb, &info);
        return info;
    } else {
        MKL_INT _n    = n;
        MKL_INT _kl   = kl;
        MKL_INT _ku   = ku;
        MKL_INT _nrhs = nrhs;
        MKL_INT _ldab = ldab;
        MKL_INT _ldb  = ldb;

        MKL_INT _info = 0;
        // FIXME: ipiv's contents are incompatible with LAPACK here
        sgbtrs((char*)&trans, &_n, &_kl, &_ku, &_nrhs,
               (float *)ab, &_ldab, (int *)ipiv, b, &_ldb, &_info);
        return _info;
    }
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
        const double *ab,
        const int ldab,
        const int *ipiv,
        double *b,
        const int ldb)
{
#ifdef SUZERAIN_HAVE_MKL
    // Casts away const; MKL LAPACK API does not enforce its logical const-ness
    // software.intel.com/en-us/forums/intel-math-kernel-library/topic/70025/
    if (sizeof(MKL_INT) == sizeof(int)) {
        int info = 0;
        dgbtrs((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
              (double *)ab, (int*)&ldab, (int *)ipiv, b, (int*)&ldb, &info);
        return info;
    } else {
        MKL_INT _n    = n;
        MKL_INT _kl   = kl;
        MKL_INT _ku   = ku;
        MKL_INT _nrhs = nrhs;
        MKL_INT _ldab = ldab;
        MKL_INT _ldb  = ldb;

        MKL_INT _info = 0;
        // FIXME: ipiv's contents are incompatible with LAPACK here
        dgbtrs((char*)&trans, &_n, &_kl, &_ku, &_nrhs,
               (double *)ab, &_ldab, (int *)ipiv, b, &_ldb, &_info);
        return _info;
    }
#else
#error "Sanity failure"
#endif
}

void
suzerain_blasext_i2s_zaxpby2(
        const int m,
        const int n,
        const double * const alpha,
        const double * const x,
        const int incx,
        const int ldx,
        const double * const beta,
        double * const y_re,
        const int incy_re,
        const int ldy_re,
        double * const y_im,
        const int incy_im,
        const int ldy_im)
{
    /* alpha == NULL defaults to 1 */
    const double alpha_re = (alpha == NULL) ? 1.0 : alpha[0];
    const double alpha_im = (alpha == NULL) ? 0.0 : alpha[1];
    /* beta == NULL defaults to 0 */
    const double beta_re  = (beta == NULL)  ? 0.0 : beta[0];
    const double beta_im  = (beta == NULL)  ? 0.0 : beta[1];

    if (alpha_im == 0.0 && beta_im == 0.0) {

        // Avoid complex arithmetic for real-valued scaling
        for (int j = 0; j < m; ++j) {

            const double * p_x_re = x + 2*j*ldx;
            const double * p_x_im = x + 2*j*ldx + 1;
            double * p_y_re = y_re + j*ldy_re;
            double * p_y_im = y_im + j*ldy_im;

            for (int i = 0; i < n; ++i) {

                /* Compute complex-valued alpha*x */
                const double alpha_x_re = (*p_x_re)*alpha_re;
                const double alpha_x_im = (*p_x_im)*alpha_re;

                /* Compute complex-valued beta*y */
                const double beta_y_re =  (*p_y_re)*beta_re;
                const double beta_y_im =  (*p_y_im)*beta_re;

                /* Store result alpha*x + beta*y in y */
                *p_y_re = alpha_x_re + beta_y_re;
                *p_y_im = alpha_x_im + beta_y_im;

                /* Increment pointers for next iteration */
                p_x_re += 2*incx;
                p_x_im += 2*incx;
                p_y_re += incy_re;
                p_y_im += incy_im;
            }
        }

    } else {

        // Complex arithmetic required
        for (int j = 0; j < m; ++j) {

            const double * p_x_re = x + 2*j*ldx;
            const double * p_x_im = x + 2*j*ldx + 1;
            double * p_y_re = y_re + j*ldy_re;
            double * p_y_im = y_im + j*ldy_im;

            for (int i = 0; i < n; ++i) {

                /* Compute complex-valued alpha*x */
                const double alpha_x_re
                    = (*p_x_re)*alpha_re - (*p_x_im)*alpha_im;
                const double alpha_x_im
                    = (*p_x_re)*alpha_im + (*p_x_im)*alpha_re;

                /* Compute complex-valued beta*y */
                const double beta_y_re
                    =  (*p_y_re)*beta_re - (*p_y_im)*beta_im;
                const double beta_y_im
                    =  (*p_y_re)*beta_im + (*p_y_im)*beta_re;

                /* Store result alpha*x + beta*y in y */
                *p_y_re = alpha_x_re + beta_y_re;
                *p_y_im = alpha_x_im + beta_y_im;

                /* Increment pointers for next iteration */
                p_x_re += 2*incx;
                p_x_im += 2*incx;
                p_y_re += incy_re;
                p_y_im += incy_im;
            }
        }

    }
}
