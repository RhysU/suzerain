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
    const float *x,
    int incx,
    const float beta,
    float *y,
    int incy)
{
    assert(0); // FIXME Implement
}

void
suzerain_bsmbsm_daPxpby(
    char trans,
    int S,
    int n,
    const double alpha,
    const double *x,
    int incx,
    const double beta,
    double *y,
    int incy)
{
    assert(0); // FIXME Implement
}

void
suzerain_bsmbsm_caPxpby(
    char trans,
    int S,
    int n,
    const float alpha[2],
    const float (*x)[2],
    int incx,
    const float beta[2],
    float (*y)[2],
    int incy)
{
    assert(0); // FIXME Implement
}

void
suzerain_bsmbsm_zaPxpby(
    char trans,
    int S,
    int n,
    const double alpha[2],
    const double (*x)[2],
    int incx,
    const double beta[2],
    double (*y)[2],
    int incy)
{
    assert(0); // FIXME Implement
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
