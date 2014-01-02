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

#ifndef SUZERAIN_SBMV_H
#define SUZERAIN_SBMV_H

/** @file
 * Provides custom, BLAS-like symmetric band matrix-vector operations.
 * Includes mixed real/complex operations and fixed bandwidth kernels.
 */

#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$
 * for symmetric, banded \f$A\f$.
 *
 * \param uplo Either 'U'/'u' or 'L'/'l' if the upper or lower triangular
 *      part of \f$A\f$ is supplied in \c a, respectively.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of super-diagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y.
 *
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for more details, especially for general
 *      band storage matrix requirements.
 */
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
        const int incy);

/*! \copydoc suzerain_sbmv_s */
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$\beta\f$, and \f$y\f$ but real-valued \f$x\f$ and
 * symmetric, real-valued \f$A\f$.
 *
 * \param uplo Either 'U'/'u' or 'L'/'l' if the upper or lower triangular
 *      part of \f$A\f$ is supplied in \c a, respectively.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of super-diagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for for general band storage matrix requirements.
 */
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
        const int incy);

/*! \copydoc suzerain_sbmv_ssc */
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
        const int incy);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} A x + \beta{} y \f$ for complex
 * \f$\alpha{}\f$, \f$x\f$, \f$\beta\f$, and \f$y\f$ but symmetric, real-valued
 * \f$A\f$.
 *
 * \param uplo Either 'U'/'u' or 'L'/'l' if the upper or lower triangular
 *      part of \f$A\f$ is supplied in \c a, respectively.
 * \param n Number of rows and columns in matrix \c a.
 * \param k Number of super-diagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a in units of <tt>float</tt>.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x in units of <tt>float[2]</tt>.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y in units of <tt>float[2]</tt>.
 *
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for for general band storage matrix requirements.
 */
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
        const int incy);

/*! \copydoc suzerain_sbmv_scc */
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
        const int incy);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_SBMV_H */
