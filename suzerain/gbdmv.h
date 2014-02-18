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

#ifndef SUZERAIN_GBDMV_H
#define SUZERAIN_GBDMV_H

/** @file
 * Provides custom, BLAS-like diagonal matrix times general band matrix-vector
 * operations.  Includes mixed real/complex operations and fixed bandwidth
 * kernels.
 */

#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D \begin{bmatrix} I & 0 \\ 0 & 1 +
 * \gamma \end{bmatrix} A x + \beta{} y \f$.
 *
 * Transposes of \f$ A \f$ can be taken using the \c trans parameter.
 *
 * \param trans One of 'N', 'T', or 'C' for no transpose, a transpose,
 *      or a conjugate transpose, respectively.
 * \param n Dimension of all matrices and vectors.
 * \param kl Number of subdiagonals in band storage of \c a.
 * \param ku Number of superdiagonals in band storage of \c a.
 * \param alpha Multiplicative scalar \f$ \alpha \f$.
 * \param d Storage for diagonal matrix \f$ D \f$.
 * \param ldd Leading dimension of \c d.
 * \param a General band storage for matrix \f$ A \f$.
 * \param lda Leading dimension of \c a.
 * \param x Vector to be multiplied.
 * \param incx Stride of vector \c x.
 * \param beta Multiplicative scalar \f$ \beta \f$.
 * \param y Vector to be added to product and to contain result.
 * \param incy Stride of vector \c y.
 * \param gamma Controls factor \f$ 1+\gamma \f$ scaling only the final row.
 *
 * \return Zero on success and a BLAS-like error code otherwise.
 * \see A BLAS reference for more details, especially for general
 *      band storage matrix requirements.
 */
int
suzerain_gbdmv_s(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const float beta,
        float *y,
        const int incy,
        const float gamma);

/*! \copydoc suzerain_gbdmv_s */
int
suzerain_gbdmv_d(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const double beta,
        double *y,
        const int incy,
        const double gamma);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D \begin{bmatrix} I & 0 \\ 0 & 1 +
 * \gamma \end{bmatrix} A x + \beta{} y \f$ for complex \f$\alpha{}\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$ \gamma \f$, \f$D\f$, \f$A\f$,
 * and \f$x\f$.
 *
 * \copydetails suzerain_gbdmv_s
 */
int
suzerain_gbdmv_ssc(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy,
        const float gamma);

/*! \copydoc suzerain_gbdmv_ssc */
int
suzerain_gbdmv_ddz(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy,
        const double gamma);

/*!
 * \brief Compute \f$ y \leftarrow{} \alpha{} D \begin{bmatrix} I & 0 \\ 0 & 1 +
 * \gamma \end{bmatrix} A x + \beta{} y \f$ for complex \f$\alpha{}\f$, \f$x\f$,
 * \f$\beta\f$, and \f$y\f$ but real-valued \f$\gamma\f$, \f$D\f$ and \f$A\f$.
 *
 * \copydetails suzerain_gbdmv_s
 */
int
suzerain_gbdmv_scc(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_float alpha,
        const float *d,
        const int ldd,
        const float *a,
        const int lda,
        const complex_float *x,
        const int incx,
        const complex_float beta,
        complex_float *y,
        const int incy,
        const float gamma);

/*! \copydoc suzerain_gbdmv_scc */
int
suzerain_gbdmv_dzz(
        const char trans,
        const int n,
        const int kl,
        const int ku,
        const complex_double alpha,
        const double *d,
        const int ldd,
        const double *a,
        const int lda,
        const complex_double *x,
        const int incx,
        const complex_double beta,
        complex_double *y,
        const int incy,
        const double gamma);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_GBDMV_H */
