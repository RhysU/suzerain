/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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

#ifdef __cplusplus
extern "C" {
#endif

/*! \name Aligned memory allocation and deallocation
 * @{
 */

/*!
 * \brief Allocates memory aligned according to the underlying BLAS'
 * recommendations for performance and numerical stability.
 *
 * The memory is not cleared.  Any returned memory must later be freed using \c
 * suzerain_blas_free().
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
 * suzerain_blas_free.
 *
 * \param nmemb Number of elements to allocate
 * \param size Size of each member in bytes.
 *
 * \return On success, return a pointer to the allocated memory.  On failure,
 * return a NULL pointer.
 */
void *
suzerain_blas_calloc(size_t nmemb, size_t size);

/*!
 * \brief Free memory previously allocated by suzerain_blas_malloc() or
 * suzerain_blas_calloc().
 *
 * \param ptr Pointer to the previously allocated memory, which may be NULL.
 */
void
suzerain_blas_free(void *ptr);

/*! @} */

/*! \name BLAS level 1 operations
 * @{
 */

/*!
 * \brief Perform \f$ x \leftrightarrow{} y \f$ using BLAS's swap.
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
suzerain_blas_sswap(
        const int n,
        float *x,
        const int incx,
        float *y,
        const int incy);

/*! \copydoc suzerain_blas_sswap */
void
suzerain_blas_dswap(
        const int n,
        double *x,
        const int incx,
        double *y,
        const int incy);

/*! \copydoc suzerain_blas_sswap */
void
suzerain_blas_cswap(
        const int n,
        float (*x)[2],
        const int incx,
        float (*y)[2],
        const int incy);

/*! \copydoc suzerain_blas_sswap */
void
suzerain_blas_zswap(
        const int n,
        double (*x)[2],
        const int incx,
        double (*y)[2],
        const int incy);

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

/*! \copydoc suzerain_blas_scopy */
void
suzerain_blas_ccopy(
        const int n,
        const float (*x)[2],
        const int incx,
        float (*y)[2],
        const int incy);

/*! \copydoc suzerain_blas_scopy */
void
suzerain_blas_zcopy(
        const int n,
        const double (*x)[2],
        const int incx,
        double (*y)[2],
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
 * \brief Compute \f$ x \cdot{} y \f$ using BLAS's dotc.
 *
 * \param[in]  n Number of elements in \c x and \c y.
 * \param[in]  x First source vector.
 * \param[in]  incx First source vector stride.
 * \param[in]  y Second source vector.
 * \param[in]  incy Second source vector stride.
 * \param[out] dotc The complex-valued result \f$ x \cdot{} y \f$.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_cdotc(
        const int n,
        const float (*x)[2],
        const int incx,
        const float (*y)[2],
        const int incy,
        float dotc[2]);

/*! \copydoc suzerain_blas_cdotc */
void
suzerain_blas_zdotc(
        const int n,
        const double (*x)[2],
        const int incx,
        const double (*y)[2],
        const int incy,
        double dotc[2]);

/*!
 * \brief Compute \f$ \left|\left| x \right|\right|_{2} \f$ using BLAS's nrm2.
 *
 * \param n Number of elements in \c x.
 * \param x Source vector.
 * \param incx Source vector stride.
 *
 * \return \f$ \left|\left| x \right|\right|_{2} \f$
 *
 * \see A BLAS reference for more details.
 */
float
suzerain_blas_snrm2(
        const int n,
        const float *x,
        const int incx);

/*! \copydoc suzerain_blas_snrm2 */
double
suzerain_blas_dnrm2(
        const int n,
        const double *x,
        const int incx);

/*! \copydoc suzerain_blas_snrm2 */
float
suzerain_blas_scnrm2(
        const int n,
        const float (*x)[2],
        const int incx);

/*! \copydoc suzerain_blas_snrm2 */
double
suzerain_blas_dznrm2(
        const int n,
        const double (*x)[2],
        const int incx);

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

/*! \copydoc suzerain_blas_saxpy */
void
suzerain_blas_caxpy(
        const int n,
        const float alpha[2],
        const float (*x)[2],
        const int incx,
        float (*y)[2],
        const int incy);

/*! \copydoc suzerain_blas_saxpy */
void
suzerain_blas_zaxpy(
        const int n,
        const double alpha[2],
        const double (*x)[2],
        const int incx,
        double (*y)[2],
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

/*! \copydoc suzerain_blas_saxpby */
void
suzerain_blas_caxpby(
        const int n,
        const float alpha[2],
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy);

/*! \copydoc suzerain_blas_saxpby */
void
suzerain_blas_zaxpby(
        const int n,
        const double alpha[2],
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
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

/*!
 * \brief Compute \f$ x \leftarrow{} \alpha{}x \f$
 * using BLAS's scal.
 *
 * \param n Number of elements in \c x.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x Source and target vector.
 * \param incx Source and target vector stride.
 *
 * \see A BLAS reference for more details.
 */
void
suzerain_blas_sscal(
        const int n,
        const float alpha,
        float *x,
        const int incx);

/*! \copydoc suzerain_blas_sscal */
void
suzerain_blas_dscal(
        const int n,
        const double alpha,
        double *x,
        const int incx);

/*! \copydoc suzerain_blas_sscal */
void
suzerain_blas_cscal(
        const int n,
        const float alpha[2],
        float (*x)[2],
        const int incx);

/*! \copydoc suzerain_blas_sscal */
void
suzerain_blas_zscal(
        const int n,
        const double alpha[2],
        double (*x)[2],
        const int incx);

/*! @} */

/*! \name BLAS level 2 operations
 * @{
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$
 * using an external BLAS's gbmv.
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
suzerain_blas_sgbmv_external(
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

/*! \copydoc suzerain_blas_sgbmv_external */
void
suzerain_blas_dgbmv_external(
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
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$.
 * The computation may use internal, optimized kernels for some fixed bandwidth
 * cases while deferring to the BLAS for general use.
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
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$
 * for symmetric \f$ A \f$ using an external BLAS's sbmv.
 *
 * \param uplo One of 'U', 'u', 'L', or 'l' specifying if either the
 *             upper- or lower-triangular part of A is being supplied.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a Symmetric band storage for matrix \f$ A \f$.
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
suzerain_blas_ssbmv_external(
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv_external */
void
suzerain_blas_dsbmv_external(
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for symmetric
 * \f$ A \f$.  The computation may use internal, optimized kernels for some
 * fixed bandwidth cases while deferring to the BLAS for general use.
 *
 * \param uplo One of 'U', 'u', 'L', or 'l' specifying if either the
 *             upper- or lower-triangular part of A is being supplied.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a Symmetric band storage for matrix \f$ A \f$.
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
suzerain_blas_ssbmv(
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
        const int incy);

/*! \copydoc suzerain_blas_ssbmv */
void
suzerain_blas_dsbmv(
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

/*! \copydoc suzerain_blas_sgb_acc */
void
suzerain_blas_cgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float (*a)[2],
        const int lda,
        const float beta[2],
        float (*b)[2],
        const int ldb);

/*! \copydoc suzerain_blas_sgb_acc */
void
suzerain_blas_zgb_acc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double (*a)[2],
        const int lda,
        const double beta[2],
        double (*b)[2],
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

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_cgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        float (*ab)[2],
        const int ldab,
        int *ipiv);

/*! \copydoc suzerain_lapack_sgbtrf */
int
suzerain_lapack_zgbtrf(
        const int m,
        const int n,
        const int kl,
        const int ku,
        double (*ab)[2],
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

/*! \copydoc suzerain_lapack_sgbtrs */
int
suzerain_lapack_cgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const float (*ab)[2],
        const int ldab,
        const int *ipiv,
        float (*b)[2],
        const int ldb);

/*! \copydoc suzerain_lapack_sgbtrs */
int
suzerain_lapack_zgbtrs(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        const double (*ab)[2],
        const int ldab,
        const int *ipiv,
        double (*b)[2],
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
 * \see A LAPACK reference for more details.
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
 * \see A LAPACK reference for more details.
 */
int
suzerain_lapack_cgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const float (*ab)[2],
        const int ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        float (*work)[2],
        float  *rwork);

/*! \copydoc suzerain_lapack_cgbcon */
int
suzerain_lapack_zgbcon(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double (*ab)[2],
        const int ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        double (*work)[2],
        double  *rwork);

/*!
 * \brief Use the LU factorization to compute the solution to a real system
 * of linear equations \f$ AX=B \f$.
 *
 * @param[in]     fact
 * @param[in]     trans
 * @param[in]     n
 * @param[in]     kl
 * @param[in]     ku
 * @param[in]     nrhs
 * @param[in,out] ab    Dimension (ldab,n)
 * @param[in]     ldab  Minimum kl + ku + 1
 * @param[in,out] afb   Dimension (ldfab,n)
 * @param[in]     ldafb Minimum 2*kl + ku + 1
 * @param[in,out] ipiv  Dimension n
 * @param[in,out] equed
 * @param[in,out] r     Dimension n when fact == 'R' or 'B'
 * @param[in,out] c     Dimension n when fact == 'C' or 'B'
 * @param[in,out] b     Dimension (ldb,nrhs)
 * @param[in]     ldb   Minimum n
 * @param[in,out] x     Dimension (ldx,nrhs)
 * @param[in]     ldx   Minimum n
 * @param[out]    rcond
 * @param[out]    ferr  Dimension nrhs
 * @param[out]    berr  Dimension nrhs
 * @param[out]    work  Dimension 3*n
 * @param[out]    iwork Dimension n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see A LAPACK reference for the details.
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
 * \brief Use the LU factorization to compute the solution to a real system
 * of linear equations \f$ AX=B \f$.
 *
 * @param[in]     fact
 * @param[in]     trans
 * @param[in]     n
 * @param[in]     kl
 * @param[in]     ku
 * @param[in]     nrhs
 * @param[in,out] ab    Dimension (ldab,n)
 * @param[in]     ldab  Minimum kl + ku + 1
 * @param[in,out] afb   Dimension (ldfab,n)
 * @param[in]     ldafb Minimum 2*kl + ku + 1
 * @param[in,out] ipiv  Dimension n
 * @param[in,out] equed
 * @param[in,out] r     Dimension n when fact == 'R' or 'B'
 * @param[in,out] c     Dimension n when fact == 'C' or 'B'
 * @param[in,out] b     Dimension (ldb,nrhs)
 * @param[in]     ldb   Minimum n
 * @param[in,out] x     Dimension (ldx,nrhs)
 * @param[in]     ldx   Minimum n
 * @param[out]    rcond
 * @param[out]    ferr  Dimension nrhs
 * @param[out]    berr  Dimension nrhs
 * @param[out]    work  Dimension 2*n
 * @param[out]    rwork Dimension n
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 *
 * \see A LAPACK reference for the details.
 */
int
suzerain_lapack_cgbsvx(
        const char fact,
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const int nrhs,
        float (*ab)[2],
        const int ldab,
        float (*afb)[2],
        const int ldafb,
        int *ipiv,
        char *equed,
        float *r,
        float *c,
        float (*b)[2],
        const int ldb,
        float (*x)[2],
        const int ldx,
        float *rcond,
        float *ferr,
        float *berr,
        float (*work)[2],
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
        double (*ab)[2],
        const int ldab,
        double (*afb)[2],
        const int ldafb,
        int *ipiv,
        char *equed,
        double *r,
        double *c,
        double (*b)[2],
        const int ldb,
        double (*x)[2],
        const int ldx,
        double *rcond,
        double *ferr,
        double *berr,
        double (*work)[2],
        double *rwork);

/*!
 * \brief Compute the the one norm, or the Frobenius norm, or the infinity
 * norm, or the element of largest absolute value of square band matrix.
 *
 * \see A LAPACK reference for more details.
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
        const float (*ab)[2],
        const int ldab,
        float *work);

/*! \copydoc suzerain_lapack_slangb */
double
suzerain_lapack_zlangb(
        const char norm,
        const int n,
        const int kl,
        const int ku,
        const double (*ab)[2],
        const int ldab,
        double *work);

/*! @} */

/*! \name BLAS-like extensions
 * These are extensions built atop the BLAS, on vendor-specific BLAS-like routines,
 * and/or on custom coded loops.
 *
 * Some of these extensions refer to interleaved versus split complex storage.
 * Interleaved complex storage stores the imaginary component immediately after
 * the real component in memory, for example \c fftw_complex or
 * <tt>std::complex<T></tt>.  Split complex storage stores the real and
 * imaginary components in wholly separate locations.
 *
 * \see The FFTW manual for further discussion on
 * <a href="http://www.fftw.org/fftw3_doc/Interleaved-and-split-arrays.html">
 * interleaved versus split complex storage</a>.
 * @{
 */

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{}x + y \f$ where \f$x\f$ is
 * real-valued and \f$\alpha\f$ and \f$y\f$ are complex-valued.  Real-valued
 * strides are in units of <tt>double</tt> while complex-valued strides are in
 * units of <tt>double[2]</tt>.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride in units of <tt>double</tt>.
 * \param y Second source vector and target vector.
 * \param incy Second source vector and target vector stride in
 *             units of <tt>double[2]</tt>.
 */
void
suzerain_blasext_daxpzy(
        const int n,
        const double alpha[2],
        const double *x,
        const int incx,
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{}x + \beta{}y \f$ where \f$x\f$ is
 * real-valued and \f$\alpha\f$, \f$\beta\f$, and \f$y\f$ are complex-valued.
 * Real-valued strides are in units of <tt>double</tt> while complex-valued
 * strides are in units of <tt>double[2]</tt>.
 *
 * \param n Number of elements in \c x and \c y.
 * \param alpha Multiplicative scalar \f$ \alpha \f$
 * \param x First source vector.
 * \param incx First source vector stride in units of <tt>double</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param y Second source vector and target vector.
 * \param incy Second source vector and target vector stride in
 *             units of <tt>double[2]</tt>.
 */
void
suzerain_blasext_daxpzby(
        const int n,
        const double alpha[2],
        const double *x,
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$A\f$
 * using an external BLAS' gbmv.  Real-valued strides are in units of
 * <tt>float</tt> while complex-valued strides are in units of
 * <tt>float[2]</tt>.
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
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_sgbmzv_external(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$A\f$
 * using an external BLAS' gbmv.  Real-valued strides are in units of
 * <tt>double</tt> while complex-valued strides are in units of
 * <tt>double[2]</tt>.
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
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>double[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>double[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_dgbmzv_external(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$A\f$.
 * Real-valued strides are in units of <tt>float</tt> while complex-valued
 * strides are in units of <tt>float[2]</tt>.
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
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_sgbmzv(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$A\f$.
 * Real-valued strides are in units of <tt>double</tt> while complex-valued
 * strides are in units of <tt>double[2]</tt>.
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
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>double[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>double[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_dgbmzv(
        const char trans,
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$
 * using an external BLAS.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param d Contiguous storage for diagonal matrix \f$ D \f$.
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
suzerain_blasext_sgbdmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *d,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy);

/*! \copydoc suzerain_blasext_dgbdmv_external */
void
suzerain_blasext_dgbdmv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *d,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D\f$
 * and \f$A\f$ using an external BLAS.  Real-valued strides are in units of
 * <tt>float</tt> while complex-valued strides are in units of
 * <tt>float[2]</tt>.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param d Contiguous storage for diagonal matrix \f$ D \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_sgbdmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float *d,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D\f$
 * and \f$A\f$ using an external BLAS.  Real-valued strides are in units of
 * <tt>double</tt> while complex-valued strides are in units of
 * <tt>double[2]</tt>.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param d Contiguous storage for diagonal matrix \f$ D \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>double[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>double[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_dgbdmzv_external(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *d,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param d Contiguous storage for diagonal matrix \f$ D \f$.
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
suzerain_blasext_sgbdmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *d,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy);

/*! \copydoc suzerain_blasext_dgbdmv */
void
suzerain_blasext_dgbdmv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *d,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D\f$
 * and \f$A\f$.  Real-valued strides are in units of <tt>float</tt> while
 * complex-valued strides are in units of <tt>float[2]</tt>.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param d Contiguous storage for diagonal matrix \f$ D \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_sgbdmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha[2],
        const float *d,
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$D\f$
 * and \f$A\f$.  Real-valued strides are in units of <tt>double</tt> while
 * complex-valued strides are in units of <tt>double[2]</tt>.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param d Contiguous storage for diagonal matrix \f$ D \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>double[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>double[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_dgbdmzv(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *d,
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but symmetric, real-valued
 * \f$A\f$ using an external BLAS' sbmv.  Real-valued strides are in units of
 * <tt>float</tt> while complex-valued strides are in units of
 * <tt>float[2]</tt>.
 *
 * \param uplo One of 'U', 'u', 'L', or 'l' specifying if either the
 *             upper- or lower-triangular part of A is being supplied.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a Symmetric band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_ssbmzv_external(
        const char uplo,
        const int n,
        const int k,
        const float alpha[2],
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but symmetric, real-valued
 * \f$A\f$ using an external BLAS' sbmv.  Real-valued strides are in units of
 * <tt>double</tt> while complex-valued strides are in units of
 * <tt>double[2]</tt>.
 *
 * \param uplo One of 'U', 'u', 'L', or 'l' specifying if either the
 *             upper- or lower-triangular part of A is being supplied.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a Symmetric band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>double[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>double[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_dsbmzv_external(
        const char uplo,
        const int n,
        const int k,
        const double alpha[2],
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but symmetric, real-valued
 * \f$A\f$.  Real-valued strides are in units of <tt>float</tt> while
 * complex-valued strides are in units of <tt>float[2]</tt>.
 *
 * \param uplo One of 'U', 'u', 'L', or 'l' specifying if either the
 *             upper- or lower-triangular part of A is being supplied.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a Symmetric band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_ssbmzv(
        const char uplo,
        const int n,
        const int k,
        const float alpha[2],
        const float *a,
        const int lda,
        const float (*x)[2],
        const int incx,
        const float beta[2],
        float (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but symmetric, real-valued
 * \f$A\f$.  Real-valued strides are in units of <tt>double</tt> while
 * complex-valued strides are in units of <tt>double[2]</tt>.
 *
 * \param uplo One of 'U', 'u', 'L', or 'l' specifying if either the
 *             upper- or lower-triangular part of A is being supplied.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a Symmetric band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>double[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>double[2]</tt>.
 *
 * \see A BLAS reference for for general band storage matrix requirements.
 */
void
suzerain_blasext_dsbmzv(
        const char uplo,
        const int n,
        const int k,
        const double alpha[2],
        const double *a,
        const int lda,
        const double (*x)[2],
        const int incx,
        const double beta[2],
        double (*y)[2],
        const int incy);

/*!
 * \brief Compute \f$ B \leftarrow{} \alpha{}A + \beta{}B \f$ using where
 * \f$B\f$, \f$\alpha\f$, and \f$\beta\f$ are complex-valued and \f$A\f$ is
 * real-valued.  Real-valued strides are in units of <tt>double</tt> while
 * complex-valued strides are in units of <tt>double[2]</tt>.
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
 * \param lda Leading dimension of \c a in units of <tt>double</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$
 * \param b General band storage of the matrix \f$ B \f$.
 * \param ldb Leading dimension of \c b in units of <tt>double[2]</tt>.
 *
 * \see A BLAS reference for general band storage matrix requirements.
 */
void
suzerain_blasext_zgb_dacc(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double alpha[2],
        const double *a,
        const int lda,
        const double beta[2],
        double (*b)[2],
        const int ldb);

/*!
 * \brief Compute two-strided \f$ Y \leftarrow{} \alpha{}X+\beta{}Y \f$
 *        where \f$ X \f$ and \f$ Y \f$ have split and interleaved
 *        complex storage, respectively.
 *
 * Accumulated interleaved complex column-major matrices into
 * two split complex column-major matrices, one each for the real and
 * imaginary parts.  The input and output arguments must not coincide.
 * Note that \c incx and \c ldx are specified in terms of complex strides,
 * and not real-valued strides.
 *
 * \param m Number of rows in column-major matrices \f$ X \f$ and \f$ Y \f$.
 * \param n Number of columns in column-major matrices \f$ X \f$ and \f$ Y \f$.
 * \param alpha Scale factor \f$ \alpha \f$ to apply to matrix \f$ X \f$.
 * \param x Beginning of interleaved complex storage for matrix \f$ X \f$.
 *          Real and imaginary components must be in adjacent locations.
 * \param incx Increment between values in a column of \f$ X \f$,
 *             specified in terms of complex elements like <tt>double[2]</tt>.
 * \param ldx  Leading distance between values in a row of \f$ X \f$,
 *             specified in terms of complex elements like <tt>double[2]</tt>.
 * \param beta Scale factor \f$ \beta \f$ to apply to matrix \f$ Y \f$
 * \param y_re Beginning of the storage for \f$\operatorname{real} Y \f$.
 * \param incy_re Increment between values in a column of
 *                \f$ \operatorname{real} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 * \param ldy_re Leading distance between values in a row of
 *                \f$ \operatorname{real} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 * \param y_im Beginning of the storage for \f$\operatorname{imag} Y \f$.
 * \param incy_im Increment between values in a column of
 *                \f$ \operatorname{imag} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 * \param ldy_im Leading distance between values in a row of
 *                \f$ \operatorname{imag} Y \f$, specified in terms of
 *                <tt>sizeof(double)</tt>.
 */
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
        const int ldy_im);

/*!
 * \brief Compute \f$ \left|\left|A\right|\right| \f$ for a general
 *        banded matrix.
 *
 * \param m Number of rows in matrix \c a.
 * \param n Number of columns in matrix \c a.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param norm1 The one norm of the matrix.
 *
 * \see A BLAS reference for more details, especially for general
 *      band storage matrix requirements.
 *
 * \return Zero on successful execution.  Nonzero otherwise.
 */
int
suzerain_blasext_sgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float *a,
        const int lda,
        float *norm1);

/*! \copydoc suzerain_blasext_sgbnorm1 */
int
suzerain_blasext_dgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double *a,
        const int lda,
        double *norm1);

/*! \copydoc suzerain_blasext_sgbnorm1 */
int
suzerain_blasext_cgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const float (*a)[2],
        const int lda,
        float *norm1);

/*! \copydoc suzerain_blasext_sgbnorm1 */
int
suzerain_blasext_zgbnorm1(
        const int m,
        const int n,
        const int kl,
        const int ku,
        const double (*a)[2],
        const int lda,
        double *norm1);

/*! @} */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_BLAS_ET_AL_H__ */
