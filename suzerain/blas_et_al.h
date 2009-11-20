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
 * blas_et_al.h: wraps external implementations of BLAS, LAPACK, et al.
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BLAS_ET_AL_H__
#define __SUZERAIN_BLAS_ET_AL_H__

#include <suzerain/common.h>

/** \file
 * Wraps external BLAS and LAPACK routines necessary for Suzerain.
 * Provided to insulate the library from potential variations in
 * type signatures as well as to consolidate all Fortran-from-C
 * parameter differences.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*! \name Memory allocation
 * @{
 */

/*!
 * \brief Allocates memory aligned according to the underlying BLAS'
 * recommendations for performance and numerical stability.
 *
 * The memory is not cleared.  Any returned memory must later be freed using \c
 * free.
 *
 * \param size Number of bytes to allocate.
 *
 * \return On success, return a pointer to the allocated memory.  On failure,
 * return a NULL pointer.
 */
void *
suzerain_blas_malloc(size_t size);

/*!
 * \brief Allocate and clear memory aligned according to the underlying BLAS'
 * recommendations for performance and numerical stability.
 *
 * The memory is set to zero.  Any returned memory must later be freed using \c
 * free.
 *
 * \param nmemb Number of elements to allocate
 * \param size Size of each member in bytes.
 *
 * \return On success, return a pointer to the allocated memory.  On failure,
 * return a NULL pointer.
 */
void *
suzerain_blas_calloc(size_t nmemb, size_t size);

/*! @} */

/*! \name BLAS level 1 operations
 * @{
 */

/*!
 * \brief Perform \f$ y \leftarrow{} x \f$ using BLAS's copy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param x Source vector.
 * \param incx Source vector stride.
 * \param y Target vector.
 * \param incy Target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_scopy(
        const int n,
        const float *x,
        const int incx,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_scopy */
void
suzerain_blas_dcopy(
        const int n,
        const double *x,
        const int incx,
        double *y,
        const int incy);

/*!
 * \brief Compute \f$ x \cdot{} y \f$ using BLAS's dot.
 *
 * \param n Number of elements in \c x and \c y.
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param y Second source vector.
 * \param incy Second source vector stride.
 *
 * \return \f$ x \cdot{} y \f$.
 *
 * \see A BLAS reference for more details.
 */
float
suzerain_blas_sdot(
        const int n,
        const float *x,
        const int incx,
        const float *y,
        const int incy);

/*! \copydoc suzerain_blas_sdot */
double
suzerain_blas_ddot(
        const int n,
        const double *x,
        const int incx,
        const double *y,
        const int incy);

/*!
 * \brief Compute \f$ \left|\left| x \right|\right|_{1} \f$ using BLAS's asum.
 *
 * \param n Number of elements in \c x.
 * \param x Source vector.
 * \param incx Source vector stride.
 *
 * \return \f$ \left|\left| x \right|\right|_{1} \f$
 *
 * \see A BLAS reference for more details.
 */
float
suzerain_blas_sasum(
        const int n,
        const float *x,
        const int incx);

/*! \copydoc suzerain_blas_sasum */
double
suzerain_blas_dasum(
        const int n,
        const double *x,
        const int incx);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{}x + y \f$ using BLAS's axpy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param y Second source vector and target vector.
 * \param incy Second source vector and target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_saxpy(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpy */
void
suzerain_blas_daxpy(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        double *y,
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{}x + \beta{}y \f$
 * using BLAS's axpby.
 *
 * If axpby is not available, simulate it using scal and axpy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param y Second source vector and target vector.
 * \param incy Second source vector and target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_saxpby(
        const int n,
        const float alpha,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_saxpby */
void
suzerain_blas_daxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy);

/*!
 * \brief Compute \f$ w \leftarrow{} \alpha{}x + \beta{}y \f$
 * using BLAS's waxpby.
 *
 * If waxpby is not available, simulate it using copy, scal, and axpy.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param y Second source vector.
 * \param incy Second source vector.
 * \param w Target vector.
 * \param incw Target vector stride.
 *
 * \see A BLAS reference for more details.
 */
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
        const int incw);

/*! \copydoc suzerain_blas_swaxpby */
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
        const int incw);

/*! @} */

/*! \name BLAS level 2 operations
 * @{
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$
 * using BLAS's gbmv.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param m Number of rows in matrix \c a.
 * \param n Number of columns in matrix \c a.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y.
 *
 * \see A BLAS reference for more details, especially for general
 *      band storage matrix requirements.
 */
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
        const int incy);

/*! \copydoc suzerain_blas_sgbmv */
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
        const int incy);

/*!
 * \brief Compute \f$ B \leftarrow{} \alpha{}A + \beta{}B \f$ using
 * BLAS' gb_acc.
 *
 * Matrices \f$ A \f$ and \f$ B \f$ both have
 * band storage and must have the same shape and same number
 * of super- and subdiagonals.
 *
 * \param m Number of rows in matrices \f$ A \f$ and \f$ B \f$.
 * \param n Number of columns in matrices \f$ A \f$ and \f$ B \f$.
 * \param kl Number of subdiagonals in band storage of \c ab.
 * \param ku Number of superdiagonals in band storage of \c ab.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param a General band storage of the matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param ldb Leading dimension of \c b.
 *
 * \see A BLAS reference for more details, especially for general
 *      band storage matrix requirements.
 */
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
        const int ldb);

/*! \copydoc suzerain_blas_sgb_acc */
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
        const int ldb);

/*! @} */

/*! \name BLAS level 3 operations
 * @{
 */

/*! @} */

/*! \name LAPACK operations
 * @{
 */

/*!
 * \brief Compute the LUP decomposition of a general banded matrix using
 * LAPACK's gbtrf.
 *
 * Stores the results back into the same matrix.  Note that the matrix must
 * have extra superdiagonals available to handle the factorization fill in.
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
 * \see A LAPACK reference for more details, especially for the
 *      \c ku storage requirements and the resulting factored
 *      storage format.
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

/*!
 * \brief Solve \f$ AX = B \f$ using the previously LUP decomposed general band matrix
 * \f$ A \f$ and LAPACK's gbtrs.
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
 * \see A LAPACK reference for more details.
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

/*! @} */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif /* __SUZERAIN_BLAS_ET_AL_H__ */
