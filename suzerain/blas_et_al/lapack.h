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

#ifndef SUZERAIN_BLAS_ET_AL_LAPACK_H
#define SUZERAIN_BLAS_ET_AL_LAPACK_H

/*!\file
 * Wraps external LAPACK routines necessary for Suzerain.
 * Provided to insulate the library from potential variations in
 * type signatures as well as to consolidate all Fortran-from-C
 * parameter differences.
 *
 * \see A LAPACK reference for more details.
 */

#include <suzerain/common.h>
#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Copies all or part of a two-dimensional matrix A to another matrix B using
 * LAPACK's LACPY.
 *
 * @param uplo Specifies the part of the matrix A to be copied to matrix B.
 *             If 'U' or 'u' the upper triangular part is copied.
 *             If 'L' or 'l' the lower triangular part is copied.
 *             Otherwise, the entire matrix is copied.
 * @param m    Number of rows of the matrix A.
 * @param n    Number of columns of the matrix A.
 * @param a    The \c m by \c n matrix A.
 * @param lda  The leading dimension of \c a.
 * @param b    The \c m by \c n matrix B.
 * @param ldb  The leading dimension of \c b.
 */
void suzerain_lapack_slacpy(
                char uplo,
                const int m,
                const int n,
                const float * a,
                const int lda,
                float *b,
                const int ldb);

/*! \copydoc suzerain_lapack_slacpy */
void suzerain_lapack_dlacpy(
                char uplo,
                const int m,
                const int n,
                const double * a,
                const int lda,
                double *b,
                const int ldb);

/*! \copydoc suzerain_lapack_slacpy */
void suzerain_lapack_clacpy(
                char uplo,
                const int m,
                const int n,
                const complex_float * a,
                const int lda,
                complex_float *b,
                const int ldb);

/*! \copydoc suzerain_lapack_slacpy */
void suzerain_lapack_zlacpy(
                char uplo,
                const int m,
                const int n,
                const complex_double * a,
                const int lda,
                complex_double *b,
                const int ldb);

/**
 * Determines machine parameters for floating-point arithmetic.
 *
 * The input parameter \c cmach should specify which of the
 * following values is desired:
 * <dl>
 * <dt><tt>E</tt>, <tt>e</tt></dt>
 * <dd><tt>eps</tt>, relative machine precision</dd>
 * <dt><tt>S</tt>, <tt>s</tt></dt>
 * <dd><tt>sfmin</tt>, safe minimum,
 *     such that <tt>1/sfmin</tt> does not overflow</dd>
 * <dt><tt>B</tt>, <tt>b</tt></dt>
 * <dd><tt>base</tt>, base of the machine</dd>
 * <dt><tt>P</tt>, <tt>p</tt></dt>
 * <dd><tt>eps*base</tt></dd>
 * <dt><tt>N</tt>, <tt>n</tt></dt>
 * <dd><tt>t</tt>, number of (base) digits in the mantissa</dd>
 * <dt><tt>R</tt>, <tt>r</tt></dt>
 * <dd><tt>rnd</tt>, 1.0 when rounding occurs in addition, 0.0 otherwise</dd>
 * <dt><tt>M</tt>, <tt>m</tt></dt>
 * <dd><tt>emin</tt>, minimum exponent before (gradual) underflow</dd>
 * <dt><tt>U</tt>, <tt>u</tt></dt>
 * <dd><tt>rmin</tt>, underflow threshold - <tt>base**(emin-1)</tt></dd>
 * <dt><tt>L</tt>, <tt>l</tt></dt>
 * <dd><tt>emax</tt>, largest exponent before overflow</dd>
 * <dt><tt>O</tt>, <tt>o</tt></dt>
 * <dd><tt>rmax</tt>, overflow threshold  - <tt>(base**emax)*(1-eps)</tt></dd>
 * </dl>
 *
 * @param cmach The desired machine-specific quantity.
 *
 * @return the value requested per parameter \c cmach.
 */
float suzerain_lapack_slamch(char cmach);

/*! \copydoc suzerain_lapack_slamch */
double suzerain_lapack_dlamch(char cmach);

/*!
 * \brief Compute the LUP decomposition of a general banded matrix using
 * LAPACK's gbtrf.
 *
 * Stores the results back into the same matrix.  Note that the banded matrix
 * storage must have \c kl extra superdiagonals available to handle the
 * factorization fill in.
 *
 * \param m Number of rows in matrix \c ab.
 * \param n Number of columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in band storage of \c ab.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix computed in the decomposition.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrs for how to solve a linear system
 *      once you have decomposed the matrix.
 * \see A LAPACK reference for the <tt>lda >= 2*kl + ku + 1</tt> storage
 *      requirement and the resulting factored storage format.
 */
int
suzerain_lapack_sgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float *ab,
        const int ldab,
        int *ipiv);

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_dgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double *ab,
        const int ldab,
        int *ipiv);

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_cgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_float *ab,
        const int ldab,
        int *ipiv);

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_zgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        complex_double *ab,
        const int ldab,
        int *ipiv);

/*!
 * \brief Solve \f$ AX = B \f$ using the previously LUP decomposed general band
 * matrix \f$ A \f$ and LAPACK's gbtrs.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Number of rows and columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in nonfactored matrix \c ab.
 *      Note this is \e not the number of superdiagonals in the storage
 *      format of \c ab, but rather the number of superdiagonals required
 *      to store the non-factored matrix \c ab.  This is odd.
 * \param nrhs Number of right hand sides, or columns, in \c b.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix already computed in the decomposition.
 * \param b Matrix \f$ B \f$ containing right hand sides on invocation and
 *      solutions on return.
 * \param ldb Leading dimension of matrix \c b.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrf for how to decompose the matrix \f$ A \f$.
 */
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
        const int ldb);

/*! \copydoc suzerain_lapack_sgbtrs */
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
        const int ldb);

/*! \copydoc suzerain_lapack_sgbtrs */
int
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
        const int ldb);

/*! \copydoc suzerain_lapack_sgbtrs */
int
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
        const int ldb);

/*!
 * \brief Estimate the reciprocal of the condition number of a
 * real-valued general band matrix \f$ A \f$ using LAPACK's gbcon.
 *
 * \param norm One of '0', '1', or 'I' for the 1-norm, the 1-norm,
 *      or the infinity-norm, respectively.
 * \param n Number of rows and columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in nonfactored matrix \c ab.
 *      Note this is \e not the number of superdiagonals in the storage
 *      format of \c ab, but rather the number of superdiagonals required
 *      to store the non-factored matrix \c ab.  This is odd.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix already computed in the decomposition.
 * \param anorm The norm of the matrix according to the norm choice
 *      made in \c norm.
 * \param rcond The reciprocal of the condition number of the matrix
 *      \f$ A \f$ computed as <tt>1/(norm(A)*norm(inv(A)))</tt>.
 * \param work Work array of dimension <tt>3*n</tt>.
 * \param iwork Work array of dimension <tt>n</tt>.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrf for how to decompose the matrix \f$ A \f$.
 */
int
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
        int *iwork);

/*! \copydoc suzerain_lapack_sgbcon */
int
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
        int *iwork);

/*!
 * \brief Estimate the reciprocal of the condition number of a
 * complex-valued general band matrix \f$ A \f$ using LAPACK's gbcon.
 *
 * \param norm One of '0', '1', or 'I' for the 1-norm, the 1-norm,
 *      or the infinity-norm, respectively.
 * \param n Number of rows and columns in matrix \c ab.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in nonfactored matrix \c ab.
 *      Note this is \e not the number of superdiagonals in the storage
 *      format of \c ab, but rather the number of superdiagonals required
 *      to store the non-factored matrix \c ab.  This is odd.
 * \param ab General band storage of the matrix to factor.
 * \param ldab Leading dimension of \c ab.
 * \param ipiv Pivot matrix already computed in the decomposition.
 * \param anorm The norm of the matrix according to the norm choice
 *      made in \c norm.
 * \param rcond The reciprocal of the condition number of the matrix
 *      \f$ A \f$ computed as <tt>1/(norm(A)*norm(inv(A)))</tt>.
 * \param work Work array of dimension <tt>2*n</tt>.
 * \param rwork Work array of dimension <tt>n</tt>.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see gbtrf for how to decompose the matrix \f$ A \f$.
 */
int
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
        float  *rwork);

/*! \copydoc suzerain_lapack_cgbcon */
int
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
        double *rwork);

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an in-place LU factorization.
 *
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
suzerain_lapack_sgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float *ab,
        const int ldab,
        int *ipiv,
        float *b,
        const int ldb);

/*! \copydoc suzerain_lapack_sgbsv */
int
suzerain_lapack_dgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        double *ab,
        const int ldab,
        int *ipiv,
        double *b,
        const int ldb);

/*! \copydoc suzerain_lapack_sgbsv */
int
suzerain_lapack_cgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_float *ab,
        const int ldab,
        int *ipiv,
        complex_float *b,
        const int ldb);

/*! \copydoc suzerain_lapack_sgbsv */
int
suzerain_lapack_zgbsv(
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        complex_double *ab,
        const int ldab,
        int *ipiv,
        complex_double *b,
        const int ldb);

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an out-of-place LU factorization.  Error bounds on the solution
 * and a condition estimate are also provided.
 *
 * \param[in]     fact
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] equed
 * \param[in,out] r     Dimension \c n when fact == 'R' or 'B'
 * \param[in,out] c     Dimension \c n when fact == 'C' or 'B'
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    rcond
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>3*n</tt>
 * \param[out]    iwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
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
        int *iwork);

/*! \copydoc suzerain_lapack_sgbsvx */
int
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
        int *iwork);

/*!
 * \brief Compute the solution to a banded system of linear equations \f$ AX=B
 * \f$ using an out-of-place LU factorization.  Error bounds on the solution
 * and a condition estimate are also provided.
 *
 * \param[in]     fact
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] equed
 * \param[in,out] r     Dimension \c n when fact == 'R' or 'B'
 * \param[in,out] c     Dimension \c n when fact == 'C' or 'B'
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    rcond
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>2*n</tt>
 * \param[out]    rwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
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
        float *rwork);

/*! \copydoc suzerain_lapack_cgbsvx */
int
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
        double *rwork);

/*!
 * \brief Iteratively refine the solution to a banded system of linear
 * equations \f$AX=B\f$ and obtain error estimates.
 *
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>3*n</tt>
 * \param[out]    iwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
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
        int *iwork);

/*! \copydoc suzerain_lapack_sgbrfs */
int
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
        int *iwork);

/*!
 * \brief Iteratively refine the solution to a banded system of linear
 * equations \f$AX=B\f$ and obtain error estimates.
 *
 * \param[in]     trans
 * \param[in]     n
 * \param[in]     kl
 * \param[in]     ku
 * \param[in]     nrhs
 * \param[in,out] ab    Dimension (<tt>ldab</tt>,<tt>n</tt>)
 * \param[in]     ldab  Minimum <tt>kl + ku + 1</tt>
 * \param[in,out] afb   Dimension (<tt>ldafb</tt>,<tt>n</tt>)
 * \param[in]     ldafb Minimum <tt>2*kl + ku + 1</tt>
 * \param[in,out] ipiv  Dimension \c n
 * \param[in,out] b     Dimension (<tt>ldb</tt>,<tt>nrhs</tt>)
 * \param[in]     ldb   Minimum \c n
 * \param[in,out] x     Dimension (<tt>ldx</tt>,<tt>nrhs</tt>)
 * \param[in]     ldx   Minimum \c n
 * \param[out]    ferr  Dimension \c nrhs
 * \param[out]    berr  Dimension \c nrhs
 * \param[out]    work  Dimension <tt>2*n</tt>
 * \param[out]    rwork Dimension \c n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
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
        float *rwork);

/*! \copydoc suzerain_lapack_cgbrfs */
int
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
        double *rwork);

/*!
 * \brief Compute the one norm, or the Frobenius norm, or the infinity
 * norm, or the element of largest absolute value of square band matrix.
 */
float
suzerain_lapack_slangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float *ab,
        const int ldab,
        float *work);

/*! \copydoc suzerain_lapack_slangb */
double
suzerain_lapack_dlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double *ab,
        const int ldab,
        double *work);

/*! \copydoc suzerain_lapack_slangb */
float
suzerain_lapack_clangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_float *ab,
        const int ldab,
        float *work);

/*! \copydoc suzerain_lapack_slangb */
double
suzerain_lapack_zlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const complex_double *ab,
        const int ldab,
        double *work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BLAS_ET_AL_LAPACK_H */
