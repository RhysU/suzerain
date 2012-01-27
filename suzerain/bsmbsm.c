/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
 * bsmbsm.c: routines for blocked square matrices with banded submatrices
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/bsmbsm.h>
#include <suzerain/blas_et_al.h>

// Shorthand
#define UNLIKELY(expr) SUZERAIN_UNLIKELY(expr)

void
suzerain_bsmbsm_saPxpby(
    char trans,
    int S,
    int n,
    const float alpha,
    const float * restrict x,
    int incx,
    const float beta,
    float * restrict y,
    int incy)
{
    if (UNLIKELY(S <  0)) return suzerain_blas_xerbla(__func__,  2);
    if (UNLIKELY(n <  0)) return suzerain_blas_xerbla(__func__,  3);
    if (UNLIKELY(x == y)) return suzerain_blas_xerbla(__func__, 58);

    // Adjust for P versus P^T operation
    switch (toupper(trans)) {
        case 'N': break;
        case 'T': S ^= n; n ^= S; S ^= n;  // (S <=> n) ==> (q <=> q^-1)
                  break;
        default:  return suzerain_blas_xerbla(__func__, 1);
    }

#pragma warning(push,disable:1572)
    const _Bool alpha_is_one  = (alpha == 1.0f);
    const _Bool beta_is_zero  = (beta  == 0.0f);
    const _Bool beta_is_one   = (beta  == 1.0f);
#pragma warning(pop)

    // Compute vector length
    const int N = S*n;

    // Dispatch based on stride characteristics
    if (incx == 1 && incy == 1) {

        // Dispatch to alpha- and beta-specific loops
        if (alpha_is_one) {
            if        (beta_is_zero) {                  // y := P x
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = x[ix];
                }
            } else if (beta_is_one) {                   // y := P x + y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] += x[ix];
                }
            } else {                                    // y := P x + beta y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = x[ix] + beta*y[i];
                }
            }
        } else {
            if        (beta_is_zero) {                  // y := alpha P x
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = alpha*x[ix];
                }
            } else if (beta_is_one) {                   // y := alpha P x + y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] += alpha*x[ix];
                }
            } else {                                    // y := alpha P x + beta y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = alpha*x[ix] + beta*y[i];
                }
            }
        }

    } else if (incx == 1) {

        // Adjust for possibly negative incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;

        // Dispatch to alpha- and beta-specific loops
        if (alpha_is_one) {
            if        (beta_is_zero) {                  // y := P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix];
                }
            } else if (beta_is_one) {                   // y := P x + y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] += x[ix];
                }
            } else {                                    // y := P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix] + beta*y[iy];
                }
            }
        } else {
            if        (beta_is_zero) {                  // y := alpha P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix];
                }
            } else if (beta_is_one) {                   // y := alpha P x + y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] += alpha*x[ix];
                }
            } else {                                    // y := alpha P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix] + beta*y[iy];
                }
            }
        }

    } else { // general strides

        // Adjust for possibly negative incx and incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;
        if (incx < 0) x += (1 - N)*incx;

        // Dispatch to alpha- and beta-specific loops
        if (alpha_is_one) {
            if        (beta_is_zero) {                  // y :=       P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix];
                }
            } else if (beta_is_one) {                   // y :=       P x +      y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] += x[ix];
                }
            } else {                                    // y :=       P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix] + beta*y[iy];
                }
            }
        } else {
            if        (beta_is_zero) {                  // y := alpha P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix];
                }
            } else if (beta_is_one) {                   // y := alpha P x +      y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] += alpha*x[ix];
                }
            } else {                                    // y := alpha P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix] + beta*y[iy];
                }
            }
        }

    }
}

void
suzerain_bsmbsm_daPxpby(
    char trans,
    int S,
    int n,
    const double alpha,
    const double * restrict x,
    int incx,
    const double beta,
    double * restrict y,
    int incy)
{
    if (UNLIKELY(S <  0)) return suzerain_blas_xerbla(__func__,  2);
    if (UNLIKELY(n <  0)) return suzerain_blas_xerbla(__func__,  3);
    if (UNLIKELY(x == y)) return suzerain_blas_xerbla(__func__, 58);

    // Adjust for P versus P^T operation
    switch (toupper(trans)) {
        case 'N': break;
        case 'T': S ^= n; n ^= S; S ^= n;  // (S <=> n) ==> (q <=> q^-1)
                  break;
        default:  return suzerain_blas_xerbla(__func__, 1);
    }

#pragma warning(push,disable:1572)
    const _Bool alpha_is_one  = (alpha == 1.0);
    const _Bool beta_is_zero  = (beta  == 0.0);
    const _Bool beta_is_one   = (beta  == 1.0);
#pragma warning(pop)

    // Compute vector length
    const int N = S*n;

    // Dispatch based on stride characteristics
    if (incx == 1 && incy == 1) {

        // Dispatch to alpha- and beta-specific loops
        if (alpha_is_one) {
            if        (beta_is_zero) {                  // y := P x
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = x[ix];
                }
            } else if (beta_is_one) {                   // y := P x + y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] += x[ix];
                }
            } else {                                    // y := P x + beta y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = x[ix] + beta*y[i];
                }
            }
        } else {
            if        (beta_is_zero) {                  // y := alpha P x
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = alpha*x[ix];
                }
            } else if (beta_is_one) {                   // y := alpha P x + y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] += alpha*x[ix];
                }
            } else {                                    // y := alpha P x + beta y
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i] = alpha*x[ix] + beta*y[i];
                }
            }
        }

    } else if (incx == 1) {

        // Adjust for possibly negative incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;

        // Dispatch to alpha- and beta-specific loops
        if (alpha_is_one) {
            if        (beta_is_zero) {                  // y := P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix];
                }
            } else if (beta_is_one) {                   // y := P x + y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] += x[ix];
                }
            } else {                                    // y := P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix] + beta*y[iy];
                }
            }
        } else {
            if        (beta_is_zero) {                  // y := alpha P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix];
                }
            } else if (beta_is_one) {                   // y := alpha P x + y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] += alpha*x[ix];
                }
            } else {                                    // y := alpha P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix] + beta*y[iy];
                }
            }
        }

    } else { // general strides

        // Adjust for possibly negative incx and incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;
        if (incx < 0) x += (1 - N)*incx;

        // Dispatch to alpha- and beta-specific loops
        if (alpha_is_one) {
            if        (beta_is_zero) {                  // y :=       P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix];
                }
            } else if (beta_is_one) {                   // y :=       P x +      y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] += x[ix];
                }
            } else {                                    // y :=       P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = x[ix] + beta*y[iy];
                }
            }
        } else {
            if        (beta_is_zero) {                  // y := alpha P x
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix];
                }
            } else if (beta_is_one) {                   // y := alpha P x +      y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] += alpha*x[ix];
                }
            } else {                                    // y := alpha P x + beta y
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy] = alpha*x[ix] + beta*y[iy];
                }
            }
        }

    }

}

void
suzerain_bsmbsm_caPxpby(
    char trans,
    int S,
    int n,
    const float alpha[2],
    const float (* restrict x)[2],
    int incx,
    const float beta[2],
    float (* restrict y)[2],
    int incy)
{
    if (UNLIKELY(S <  0)) return suzerain_blas_xerbla(__func__,  2);
    if (UNLIKELY(n <  0)) return suzerain_blas_xerbla(__func__,  3);
    if (UNLIKELY((void*)x == (void*)y))
                          return suzerain_blas_xerbla(__func__, 58);

    // Adjust for P versus P^T operation
    switch (toupper(trans)) {
        case 'N': break;
        case 'T': S ^= n; n ^= S; S ^= n;  // (S <=> n) ==> (q <=> q^-1)
                  break;
        default:  return suzerain_blas_xerbla(__func__, 1);
    }

#pragma warning(push,disable:1572)
    const _Bool alpha_is_real = (alpha[1] == 0.0f);
    const _Bool beta_is_real  = (beta[1]  == 0.0f);
    const _Bool alpha_is_one  = (alpha_is_real && alpha[0] == 1.0f);
    const _Bool beta_is_one   = (beta_is_real  && beta[0]  == 1.0f);
    const _Bool beta_is_zero  = (beta_is_real  && beta[0]  == 0.0f);
#pragma warning(pop)

    // Compute vector length
    const int N = S*n;

    // Dispatch based on stride characteristics
    if (incx == 1 && incy == 1) {

        // Dispatch to alpha- and beta-specific loops
        if        (alpha_is_one) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] = x[ix][0];                               // P x
                    y[i][1] = x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] += x[ix][0];                              // P x
                    y[i][1] += x[ix][1];                              // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] *= beta[0];                               // beta y
                    y[i][1] *= beta[0];
                    y[i][0] += x[ix][0];                              // P x
                    y[i][1] += x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]  = beta [0]*y[i][0] - beta [1]*y[i][1];    // beta y
                    tmp[1]  = beta [0]*y[i][1] + beta [1]*y[i][0];
                    y[i][0] = tmp[0] + x[ix][0];                      // P x
                    y[i][1] = tmp[1] + x[ix][1];
                }
            }
        } else if (alpha_is_real) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] = alpha[0]*x[ix][0];                      // alpha P x
                    y[i][1] = alpha[0]*x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[i][1] += alpha[0]*x[ix][1];                     // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] *= beta[0];                               // beta y
                    y[i][1] *= beta[0];
                    y[i][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[i][1] += alpha[0]*x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]  = beta [0]*y[i][0] - beta [1]*y[i][1];    // beta y
                    tmp[1]  = beta [0]*y[i][1] + beta [1]*y[i][0];
                    y[i][0] = tmp[0] + alpha[0]*x[ix][0];             // alpha P x
                    y[i][1] = tmp[1] + alpha[0]*x[ix][1];
                }
            }
        } else {// alpha_is_complex
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[i][1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[i][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0]; // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] *= beta[0];                               // beta y
                    y[i][1] *= beta[0];
                    y[i][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[i][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    tmp[1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                    tmp[0]  += beta [0]*y[i][0] - beta [1]*y[i][1];   // beta y
                    tmp[1]  += beta [0]*y[i][1] + beta [1]*y[i][0];
                    y[i][0] = tmp[0];
                    y[i][1] = tmp[1];
                }
            }
        }

    } else if (incx == 1) {

        // Adjust for possibly negative incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;

        // Dispatch to alpha- and beta-specific loops
        if        (alpha_is_one) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] = x[ix][0];                               // P x
                    y[iy][1] = x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];                              // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + x[ix][0];                      // P x
                    y[iy][1] = tmp[1] + x[ix][1];
                }
            }
        } else if (alpha_is_real) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];                     // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + alpha[0]*x[ix][0];             // alpha P x
                    y[iy][1] = tmp[1] + alpha[0]*x[ix][1];
                }
            }
        } else {// alpha_is_complex
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];  // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]   = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    tmp[1]   = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                    tmp[0]  += beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]  += beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0];
                    y[iy][1] = tmp[1];
                }
            }
        }

    } else { // general strides

        // Adjust for possibly negative incx and incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;
        if (incx < 0) x += (1 - N)*incx;

        // Dispatch to alpha- and beta-specific loops
        if        (alpha_is_one) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] = x[ix][0];                               // P x
                    y[iy][1] = x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];                              // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + x[ix][0];                      // P x
                    y[iy][1] = tmp[1] + x[ix][1];
                }
            }
        } else if (alpha_is_real) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];                     // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + alpha[0]*x[ix][0];             // alpha P x
                    y[iy][1] = tmp[1] + alpha[0]*x[ix][1];
                }
            }
        } else {// alpha_is_complex
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];  // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    float tmp[2];
                    tmp[0]   = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    tmp[1]   = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                    tmp[0]  += beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]  += beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0];
                    y[iy][1] = tmp[1];
                }
            }
        }

    }
}

void
suzerain_bsmbsm_zaPxpby(
    char trans,
    int S,
    int n,
    const double alpha[2],
    const double (* restrict x)[2],
    int incx,
    const double beta[2],
    double (* restrict y)[2],
    int incy)
{
    if (UNLIKELY(S <  0)) return suzerain_blas_xerbla(__func__,  2);
    if (UNLIKELY(n <  0)) return suzerain_blas_xerbla(__func__,  3);
    if (UNLIKELY((void*)x == (void*)y))
                          return suzerain_blas_xerbla(__func__, 58);

    // Adjust for P versus P^T operation
    switch (toupper(trans)) {
        case 'N': break;
        case 'T': S ^= n; n ^= S; S ^= n;  // (S <=> n) ==> (q <=> q^-1)
                  break;
        default:  return suzerain_blas_xerbla(__func__, 1);
    }

#pragma warning(push,disable:1572)
    const _Bool alpha_is_real = (alpha[1] == 0.0);
    const _Bool beta_is_real  = (beta[1]  == 0.0);
    const _Bool alpha_is_one  = (alpha_is_real && alpha[0] == 1.0);
    const _Bool beta_is_one   = (beta_is_real  && beta[0]  == 1.0);
    const _Bool beta_is_zero  = (beta_is_real  && beta[0]  == 0.0);
#pragma warning(pop)

    // Compute vector length
    const int N = S*n;

    // Dispatch based on stride characteristics
    if (incx == 1 && incy == 1) {

        // Dispatch to alpha- and beta-specific loops
        if        (alpha_is_one) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] = x[ix][0];                               // P x
                    y[i][1] = x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] += x[ix][0];                              // P x
                    y[i][1] += x[ix][1];                              // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] *= beta[0];                               // beta y
                    y[i][1] *= beta[0];
                    y[i][0] += x[ix][0];                              // P x
                    y[i][1] += x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]  = beta [0]*y[i][0] - beta [1]*y[i][1];    // beta y
                    tmp[1]  = beta [0]*y[i][1] + beta [1]*y[i][0];
                    y[i][0] = tmp[0] + x[ix][0];                      // P x
                    y[i][1] = tmp[1] + x[ix][1];
                }
            }
        } else if (alpha_is_real) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] = alpha[0]*x[ix][0];                      // alpha P x
                    y[i][1] = alpha[0]*x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[i][1] += alpha[0]*x[ix][1];                     // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] *= beta[0];                               // beta y
                    y[i][1] *= beta[0];
                    y[i][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[i][1] += alpha[0]*x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]  = beta [0]*y[i][0] - beta [1]*y[i][1];    // beta y
                    tmp[1]  = beta [0]*y[i][1] + beta [1]*y[i][0];
                    y[i][0] = tmp[0] + alpha[0]*x[ix][0];             // alpha P x
                    y[i][1] = tmp[1] + alpha[0]*x[ix][1];
                }
            }
        } else {// alpha_is_complex
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[i][1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[i][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0]; // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[i][0] *= beta[0];                               // beta y
                    y[i][1] *= beta[0];
                    y[i][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[i][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    tmp[1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                    tmp[0]  += beta [0]*y[i][0] - beta [1]*y[i][1];   // beta y
                    tmp[1]  += beta [0]*y[i][1] + beta [1]*y[i][0];
                    y[i][0] = tmp[0];
                    y[i][1] = tmp[1];
                }
            }
        }

    } else if (incx == 1) {

        // Adjust for possibly negative incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;

        // Dispatch to alpha- and beta-specific loops
        if        (alpha_is_one) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] = x[ix][0];                               // P x
                    y[iy][1] = x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];                              // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + x[ix][0];                      // P x
                    y[iy][1] = tmp[1] + x[ix][1];
                }
            }
        } else if (alpha_is_real) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];                     // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + alpha[0]*x[ix][0];             // alpha P x
                    y[iy][1] = tmp[1] + alpha[0]*x[ix][1];
                }
            }
        } else {// alpha_is_complex
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];  // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]   = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    tmp[1]   = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                    tmp[0]  += beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]  += beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0];
                    y[iy][1] = tmp[1];
                }
            }
        }

    } else { // general strides

        // Adjust for possibly negative incx and incy
        int iy = (incy < 0) ? (1 - N)*incy : 0;
        if (incx < 0) x += (1 - N)*incx;

        // Dispatch to alpha- and beta-specific loops
        if        (alpha_is_one) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] = x[ix][0];                               // P x
                    y[iy][1] = x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];                              // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += x[ix][0];                              // P x
                    y[iy][1] += x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + x[ix][0];                      // P x
                    y[iy][1] = tmp[1] + x[ix][1];
                }
            }
        } else if (alpha_is_real) {
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];                     // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0];                     // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]   = beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]   = beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0] + alpha[0]*x[ix][0];             // alpha P x
                    y[iy][1] = tmp[1] + alpha[0]*x[ix][1];
                }
            }
        } else {// alpha_is_complex
            if        (beta_is_zero) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0]  = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1]  = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else if (beta_is_one)  {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];  // y
                }
            } else if (beta_is_real) {
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    y[iy][0] *= beta[0];                               // beta y
                    y[iy][1] *= beta[0];
                    y[iy][0] += alpha[0]*x[ix][0] - alpha[1]*x[ix][1]; // alpha P x
                    y[iy][1] += alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                }
            } else {// beta_is_complex
                for (int i = 0; i < N; ++i, iy += incy) {
                    const int ix = incx*suzerain_bsmbsm_q(S, n, i);
                    double tmp[2];
                    tmp[0]   = alpha[0]*x[ix][0] - alpha[1]*x[ix][1];  // alpha P x
                    tmp[1]   = alpha[0]*x[ix][1] + alpha[1]*x[ix][0];
                    tmp[0]  += beta [0]*y[iy][0] - beta [1]*y[iy][1];  // beta y
                    tmp[1]  += beta [0]*y[iy][1] + beta [1]*y[iy][0];
                    y[iy][0] = tmp[0];
                    y[iy][1] = tmp[1];
                }
            }
        }

    }
}

gsl_permutation *
suzerain_bsmbsm_permutation(int S, int n)
{
    assert(S > 0);
    assert(n > 0);

    const ptrdiff_t N = S*n;
    gsl_permutation * const p = gsl_permutation_alloc(N);
    if (p) {
        size_t * const q = p->data;
        for (ptrdiff_t i = 0; i < N; ++i) {
            q[i] = suzerain_bsmbsm_q(S, n, i);
        }
    }
    return p;
}
