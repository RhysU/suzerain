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
 * @copydoc lapack.h
 */

#include <suzerain/blas_et_al/lapack.h>

#include <suzerain/common.h>

// Odd looking, unused, anonymous enumerations are globally-scoped versions
// of "Compile Time Assertions" by Ralf Holly (http://drdobbs.com/184401873)
enum { // Required for strict aliasing workarounds
    assert_floatp  = 1/(sizeof(float*)  == sizeof(complex_float*) ),
    assert_doublep = 1/(sizeof(double*) == sizeof(complex_double*))
};

#ifdef SUZERAIN_HAVE_MKL
# include <mkl_types.h>
# include <mkl_lapack.h>
# define LAPACK_FUNC(name,NAME) name             /* Think AC_F77_WRAPPERS */
enum { // Required for binary interoperability with MKL
    assert_mkl_i = 1 /(sizeof(MKL_INT)       == sizeof(int)           ),
    assert_mkl_c = 1 /(sizeof(MKL_Complex8)  == sizeof(complex_float) ),
    assert_mkl_z = 1 /(sizeof(MKL_Complex16) == sizeof(complex_double))
};
#else
# error "No suitable BLAS and/or LAPACK library found during configuration"
#endif

// Shorthand
#define CAST_VOIDP(l)       SUZERAIN_CAST_VOIDP(l)
#define CONST_CAST_VOIDP(l) SUZERAIN_CONST_CAST_VOIDP(l)
#define UNLIKELY(expr)      SUZERAIN_UNLIKELY(expr)

// Many of the short methods have "inline" though their declarations do not.
// to allow inlining them later within this particular translation unit.

// Thank you captain obvious...
#pragma warning(disable:981)

inline void
suzerain_lapack_slacpy(
                char uplo,
                const int m,
                const int n,
                const float * a,
                const int lda,
                float *b,
                const int ldb)
{
    LAPACK_FUNC(slacpy,SLACPY)((char*)&uplo, (int*)&m, (int*)&n,
                               (float*)a, (int*)&lda,
                               (float*)b, (int*)&ldb);
}

inline void
suzerain_lapack_dlacpy(
                char uplo,
                const int m,
                const int n,
                const double * a,
                const int lda,
                double *b,
                const int ldb)
{
    LAPACK_FUNC(dlacpy,DLACPY)((char*)&uplo, (int*)&m, (int*)&n,
                               (double*)a, (int*)&lda,
                               (double*)b, (int*)&ldb);
}

inline void
suzerain_lapack_clacpy(
                char uplo,
                const int m,
                const int n,
                const complex_float * a,
                const int lda,
                complex_float *b,
                const int ldb)
{
    LAPACK_FUNC(clacpy,CLACPY)((char*)&uplo, (int*)&m, (int*)&n,
                               CAST_VOIDP(a), (int*)&lda,
                               CAST_VOIDP(b), (int*)&ldb);
}

inline void
suzerain_lapack_zlacpy(
                char uplo,
                const int m,
                const int n,
                const complex_double * a,
                const int lda,
                complex_double *b,
                const int ldb)
{
    LAPACK_FUNC(zlacpy,ZLACPY)((char*)&uplo, (int*)&m, (int*)&n,
                               CAST_VOIDP(a), (int*)&lda,
                               CAST_VOIDP(b), (int*)&ldb);
}

inline float
suzerain_lapack_slamch(char cmach)
{
    return LAPACK_FUNC(slamch,SLAMCH)((char*)&cmach);
}

inline double
suzerain_lapack_dlamch(char cmach)
{
    return LAPACK_FUNC(dlamch,DLAMCH)((char*)&cmach);
}

inline int
suzerain_lapack_sgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(sgbtrf,SGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               ab, (int*)&ldab, ipiv, &info);
    return info;
}

inline int
suzerain_lapack_dgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(dgbtrf,DGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               ab, (int*)&ldab, ipiv, &info);
    return info;
}

inline int
suzerain_lapack_cgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_float *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(cgbtrf,CGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               CAST_VOIDP(ab), (int *)&ldab, ipiv, &info);
    return info;
}

inline int
suzerain_lapack_zgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_double *ab,
        const int ldab,
        int *ipiv)
{
    int info = 0;
    LAPACK_FUNC(zgbtrf,ZGBTRF)((int*)&m, (int*)&n, (int*)&kl, (int*)&ku,
                               CAST_VOIDP(ab), (int *)&ldab, ipiv, &info);
    return info;
}

inline int
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
    int info = 0;
    LAPACK_FUNC(sgbtrs,SGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (float *)ab, (int*)&ldab,
                               (int *)ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
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
    int info = 0;
    LAPACK_FUNC(dgbtrs,DGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, (double *)ab, (int*)&ldab,
                               (int *)ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_cgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const complex_float *ab,
        const int ldab,
        const int *ipiv,
        complex_float *b,
        const int ldb)
{
    int info = 0;
    LAPACK_FUNC(cgbtrs,CGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, CAST_VOIDP(ab), (int*)&ldab,
                               (int *)ipiv, CAST_VOIDP(b ), (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_zgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const complex_double *ab,
        const int ldab,
        const int *ipiv,
        complex_double *b,
        const int ldb)
{
    int info = 0;
    LAPACK_FUNC(zgbtrs,ZGBTRS)((char*)&trans, (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, CAST_VOIDP(ab), (int*)&ldab,
                               (int *)ipiv, CAST_VOIDP(b ), (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_sgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float *ab,
        const int ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        float *work,
        int *iwork)
{
    int info = 0;
    LAPACK_FUNC(sgbcon,SGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               (float*)ab, (int*)&ldab, (int*) ipiv,
                               (float*)&anorm, rcond, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_dgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double *ab,
        const int ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        double *work,
        int *iwork)
{
    int info = 0;
    LAPACK_FUNC(dgbcon,DGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               (double*)ab, (int*)&ldab, (int*) ipiv,
                               (double*)&anorm, rcond, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_cgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_float *ab,
        const int ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        complex_float *work,
        float  *rwork)
{
    int info = 0;
    LAPACK_FUNC(cgbcon,CGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               CAST_VOIDP(ab), (int*)&ldab, (int*)ipiv,
                               (float*)&anorm, rcond, CAST_VOIDP(work), rwork,
                               &info);
    return info;
}

inline int
suzerain_lapack_zgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_double *ab,
        const int ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        complex_double *work,
        double  *rwork)
{
    int info = 0;
    LAPACK_FUNC(zgbcon,ZGBCON)((char*)&norm, (int*)&n, (int*)&kl, (int*)&ku,
                               CAST_VOIDP(ab), (int*)&ldab, (int*) ipiv,
                               (double*)&anorm, rcond, CAST_VOIDP(work), rwork,
                               &info);
    return info;
}

inline int
suzerain_lapack_sgbsv(
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
    int info;
    LAPACK_FUNC(sgbsv,SGBSV)((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
                             ab, (int*)&ldab, ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_dgbsv(
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
    int info;
    LAPACK_FUNC(dgbsv,DGBSV)((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
                             ab, (int*)&ldab, ipiv, b, (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_cgbsv(
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
    int info;
    LAPACK_FUNC(cgbsv,CGBSV)((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
                             CAST_VOIDP(ab), (int*)&ldab, ipiv,
                             CAST_VOIDP(b ), (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_zgbsv(
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
    int info;
    zgbsv((int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs,
          CAST_VOIDP(ab), (int*)&ldab, ipiv,
          CAST_VOIDP(b ), (int*)&ldb, &info);
    return info;
}

inline int
suzerain_lapack_sgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        float *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        float *r,
        float *c,
        float *b,
        const int ldb,
        float *x,
        const int ldx,
        float *rcond,
        float *ferr,
        float *berr,
        float *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(sgbsvx,SGBSVX)((char*)&fact, (char*)&trans,
                               (int*)&n, (int*)&kl, (int*)&ku, (int*)&nrhs, ab,
                               (int*)&ldab, afb, (int*)&ldafb, ipiv, equed, r,
                               c, b, (int*)&ldb, x, (int*)&ldx, rcond, ferr,
                               berr, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_dgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        double *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        double *r,
        double *c,
        double *b,
        const int ldb,
        double *x,
        const int ldx,
        double *rcond,
        double *ferr,
        double *berr,
        double *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(dgbsvx,DGBSVX)((char*)&fact, (char*)&trans, (int*)&n,
                               (int*)&kl, (int*)&ku, (int*)&nrhs, ab,
                               (int*)&ldab, afb, (int*)&ldafb, ipiv, equed, r,
                               c, b, (int*)&ldb, x, (int*)&ldx, rcond, ferr,
                               berr, work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_cgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        complex_float *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        float *r,
        float *c,
        complex_float *b,
        const int ldb,
        complex_float *x,
        const int ldx,
        float *rcond,
        float *ferr,
        float *berr,
        complex_float *work,
        float *rwork)
{
    int info;
    LAPACK_FUNC(cgbsvx,CGBSVX)(
            (char*)&fact, (char*)&trans, (int*)&n,
            (int*)&kl, (int*)&ku, (int*)&nrhs, CAST_VOIDP(ab),
            (int*)&ldab, CAST_VOIDP(afb), (int*)&ldafb, ipiv,
            equed, r, c, CAST_VOIDP(b  ), (int*)&ldb, CAST_VOIDP(x),
            (int*)&ldx, rcond, ferr, berr, CAST_VOIDP(work),
            rwork, &info);
    return info;
}

inline int
suzerain_lapack_zgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        complex_double *afb,
        const int ldafb,
        int *ipiv,
        char *equed,
        double *r,
        double *c,
        complex_double *b,
        const int ldb,
        complex_double *x,
        const int ldx,
        double *rcond,
        double *ferr,
        double *berr,
        complex_double *work,
        double *rwork)
{
    int info;
    LAPACK_FUNC(zgbsvx,ZGBSVX)(
            (char*)&fact, (char*)&trans,
            (int*)&n, (int*)&kl, (int*)&ku,
            (int*)&nrhs, CAST_VOIDP(ab), (int*)&ldab,
            CAST_VOIDP(afb), (int*)&ldafb, ipiv, equed,
            r, c, CAST_VOIDP(b), (int*)&ldb, CAST_VOIDP(x),
            (int*)&ldx, rcond, ferr, berr, CAST_VOIDP(work),
            rwork, &info);
    return info;
}

inline int
suzerain_lapack_sgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        float *afb,
        const int ldafb,
        int *ipiv,
        float *b,
        const int ldb,
        float *x,
        const int ldx,
        float *ferr,
        float *berr,
        float *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(sgbrfs,SGBRFS)((char*)&trans,
                               (int*)&n, (int*)&kl, (int*)&ku,
                               (int*)&nrhs, ab, (int*)&ldab, afb,
                               (int*)&ldafb, ipiv, b, (int*)&ldb, x,
                               (int*)&ldx, ferr, berr, work, iwork,
                               &info);
    return info;
}

inline int
suzerain_lapack_dgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        double *afb,
        const int ldafb,
        int *ipiv,
        double *b,
        const int ldb,
        double *x,
        const int ldx,
        double *ferr,
        double *berr,
        double *work,
        int *iwork)
{
    int info;
    LAPACK_FUNC(dgbrfs,DGBRFS)((char*)&trans, (int*)&n,
                               (int*)&kl, (int*)&ku, (int*)&nrhs, ab,
                               (int*)&ldab, afb, (int*)&ldafb, ipiv, b,
                               (int*)&ldb, x, (int*)&ldx, ferr, berr,
                               work, iwork, &info);
    return info;
}

inline int
suzerain_lapack_cgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        complex_float *afb,
        const int ldafb,
        int *ipiv,
        complex_float *b,
        const int ldb,
        complex_float *x,
        const int ldx,
        float *ferr,
        float *berr,
        complex_float *work,
        float *rwork)
{
    int info;
    LAPACK_FUNC(cgbrfs,CGBRFS)(
            (char*)&trans, (int*)&n,
            (int*)&kl, (int*)&ku, (int*)&nrhs,
            CAST_VOIDP(ab), (int*)&ldab, CAST_VOIDP(afb),
            (int*)&ldafb, ipiv, CAST_VOIDP(b), (int*)&ldb,
            CAST_VOIDP(x), (int*)&ldx, ferr, berr,
            CAST_VOIDP(work), rwork, &info);
    return info;
}

inline int
suzerain_lapack_zgbrfs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        complex_double *afb,
        const int ldafb,
        int *ipiv,
        complex_double *b,
        const int ldb,
        complex_double *x,
        const int ldx,
        double *ferr,
        double *berr,
        complex_double *work,
        double *rwork)
{
    int info;
    LAPACK_FUNC(zgbrfs,ZGBRFS)(
            (char*)&trans,
            (int*)&n, (int*)&kl, (int*)&ku,
            (int*)&nrhs, CAST_VOIDP(ab), (int*)&ldab,
            CAST_VOIDP(afb), (int*)&ldafb, ipiv, CAST_VOIDP(b),
            (int*)&ldb, CAST_VOIDP(x), (int*)&ldx, ferr,
            berr, CAST_VOIDP(work), rwork, &info);
    return info;
}

inline float
suzerain_lapack_slangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float *ab,
        const int ldab,
        float *work)
{
    return LAPACK_FUNC(slangb,SLANGB)((char*)&norm, (int*)&n, (int*)&kl,
                                      (int*)&ku, (float*)ab, (int*)&ldab,
                                      work);
}

inline double
suzerain_lapack_dlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double *ab,
        const int ldab,
        double *work)
{
    return LAPACK_FUNC(dlangb,DLANGB)((char*)&norm, (int*)&n, (int*)&kl,
                                      (int*)&ku, (double*)ab, (int*)&ldab,
                                      work);
}

inline float
suzerain_lapack_clangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_float *ab,
        const int ldab,
        float *work)
{
    return LAPACK_FUNC(clangb,CLANGB)(
            (char*)&norm, (int*)&n, (int*)&kl,
            (int*)&ku, CAST_VOIDP(ab), (int*)&ldab, work);
}

inline double
suzerain_lapack_zlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_double *ab,
        const int ldab,
        double *work)
{
    return LAPACK_FUNC(zlangb,ZLANGB)(
            (char*)&norm, (int*)&n, (int*)&kl,
            (int*)&ku, CAST_VOIDP(ab), (int*)&ldab, work);
}
