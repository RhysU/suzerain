/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011, 2012 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *--------------------------------------------------------------------------
 * sbmv.h: provides custom, BLAS-like general band matrix-vector operations
 * $Id$
 */

#ifndef __SUZERAIN_SBMV_H
#define __SUZERAIN_SBMV_H

#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/** \file
 * Provides custom, BLAS-like general band matrix-vector operations.
 * Includes mixed real/complex operations and fixed bandwidth kernels.
 */

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

#endif /* __SUZERAIN_SBMV_H */
