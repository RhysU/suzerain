/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
 * filterop.c: Discrete filtering operator routines
 * $Id$
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <suzerain/blas_et_al.h>
#include <suzerain/error.h>
#include <suzerain/filterop.h>

// TODO Account for non-banded nature of periodicity (Toeplitz)

static int suzerain_filterop_operator_bandwidths_cookcabot2005(
    const double *method_params,
    int *klat,
    int *kuat,
    int *klbt,
    int *kubt)
{
    SUZERAIN_UNUSED(method_params);
    *klat = *kuat = 2;
    *klbt = *kubt = 4;  // TODO Confirm 4 not 5
    return SUZERAIN_SUCCESS;
}

static int suzerain_filterop_operator_bandwidths(
    const enum suzerain_filterop_method method,
    const double *method_params,
    int *klat,
    int *kuat,
    int *ldat,
    int *klbt,
    int *kubt,
    int *ldbt)
{
    int retval;

    // Add new method-specific routines here, when necessary
    switch (method) {
        case SUZERAIN_FILTEROP_COOKCABOT2005:
            retval = suzerain_filterop_operator_bandwidths_cookcabot2005(
                        method_params, klat, kuat, klbt, kubt);
        default:
            SUZERAIN_ERROR_NULL("Unknown method", SUZERAIN_ESANITY);
    }

    assert(*klat >= 0);
    assert(*kuat >= 0);
    assert(*klbt >= 0);
    assert(*kubt >= 0);

    if (retval == SUZERAIN_SUCCESS) {
        *ldat = 2*(*klat) + *kuat + 1; // Preparing for factorization
        *ldbt =    *klbt  + *kubt + 1; // No factorization
    }

    return retval;
}

static int suzerain_filterop_operator_assemble_cookcabot2005(
    const double *method_params,
    suzerain_filterop_workspace *w)
{
    SUZERAIN_UNUSED(method_params);
    SUZERAIN_UNUSED(w);
    SUZERAIN_ERROR("Unimplemented", SUZERAIN_ESANITY); // FIXME
}

static int suzerain_filterop_operator_assemble(
    const enum suzerain_filterop_method method,
    const double *method_params,
    suzerain_filterop_workspace *w)
{
    int retval;

    // Add new method-specific routines here, when necessary
    switch (method) {
        case SUZERAIN_FILTEROP_COOKCABOT2005:
            retval = suzerain_filterop_operator_assemble_cookcabot2005(
                        method_params, w);
        default:
            SUZERAIN_ERROR_NULL("Unknown method", SUZERAIN_ESANITY);
    }
    return retval;
}

suzerain_filterop_workspace *
suzerain_filterop_alloc(
    const int n,
    const enum suzerain_filterop_method method,
    const double *method_params)
{
    /* Parameter sanity checks */
    if (n < 1) {
        SUZERAIN_ERROR_NULL("n must be strictly positive", SUZERAIN_EINVAL);
    }

    /* Compute bandwidth-related information before allocating workspace */
    int klat, kuat, ldat;
    int klbt, kubt, ldbt;
    if (SUZERAIN_SUCCESS != suzerain_filterop_operator_bandwidths(method,
                method_params, &klat, &kuat, &ldat, &klbt, &kubt, &ldbt)) {
        SUZERAIN_ERROR_NULL("Unable to compute bandwidths for method",
                            SUZERAIN_EFAILED);
    }

    /* Allocate space for the workspace struct and initialize to zero */
    suzerain_filterop_workspace * const w
        = suzerain_blas_calloc(1, sizeof(suzerain_filterop_workspace));
    if (w == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }

    /* Populate initial workspace values */
    w->method = method;
    w->n      = n;
    w->klat   = klat;
    w->kuat   = kuat;
    w->ldat   = ldat;
    w->klbt   = klbt;
    w->kubt   = kubt;
    w->ldbt   = ldbt;

    /* Allocate one contiguous block for B and A, in that order */
    w->B_T = suzerain_blas_calloc((w->n * w->ldbt) + (w->n * w->ldat),
                                  sizeof(w->B_T[0]));
    if (w->B_T == NULL) {
        suzerain_filterop_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for B_T and A_T",
                            SUZERAIN_ENOMEM);
    }
    w->A_T = w->B_T + (w->n * w->ldbt);

    /* Allocate space for the pivot vector for A_T's factorization */
    w->ipiva = suzerain_blas_calloc(w->n, sizeof(w->ipiva[0]));
    if (w->ipiva == NULL) {
        suzerain_filterop_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for ipiva",
                            SUZERAIN_ENOMEM);
    }

    /* Initialize A_T and B_T depending on the method */
    if (SUZERAIN_SUCCESS != suzerain_filterop_operator_assemble(
                method, method_params, w)) {
        suzerain_filterop_free(w);
        SUZERAIN_ERROR_NULL("failed to assemble A_T and B_T",
                            SUZERAIN_EFAILED);
    }

    return w;
}

void
suzerain_filterop_free(
    suzerain_filterop_workspace *w)
{
    if (w) {
        // B_T and A_T were allocated as a single block, so just one free
        w->A_T = NULL;
        suzerain_blas_free(w->B_T);
        w->B_T = NULL;
        suzerain_blas_free(w->ipiva);
        w->ipiva = NULL;
    }
}

int
suzerain_filterop_factorize(
    suzerain_filterop_workspace *w)
{
    /* Protect the user from incorrect API usage */
    if (w->ipiva[0] > 0) {
        SUZERAIN_ERROR("Unable to factor an already-factored operator",
                       SUZERAIN_EINVAL);
    }

    /* Compute in-place LU factorization of the operator */
    const int info = suzerain_lapack_dgbtrf(w->n, w->n, w->klat, w->kuat,
                                            w->A_T, w->ldat, w->ipiva);
    if (info) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_dgbtrf reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    /* Factorization overwrote w->ipiv[0] with something nonnegative;
     * This is our flag to check that factorization indeed occurred */

    return SUZERAIN_SUCCESS;
}

int
suzerain_filterop_apply(
    const double alpha,
    const double *x,
    const int incx,
    double *y,
    const int incy,
    const suzerain_filterop_workspace *w)
{
    return suzerain_blas_dgbmv('T', w->n, w->n, w->klbt, w->kubt,
                               alpha, w->B_T, w->ldbt, x, incx,
                               /* overwrite */ 0,      y, incy);
}

int
suzerain_filterop_solve(
    const int nrhs,
    double *X,
    const int ldx,
    const suzerain_filterop_workspace * w)
{
    if (w->ipiva[0] <= 0) {
        SUZERAIN_ERROR(
                "suzerain_filterop_factorize not called before solve",
                SUZERAIN_EINVAL);
    }

    const int info = suzerain_lapack_dgbtrs('T', w->n, w->klat, w->kuat,
                                            nrhs, w->A_T, w->ldat, w->ipiva,
                                            X, ldx);
    if (info) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_dgbtrs reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_filterop_filter(
    const double alpha,
    const double *x,
    const int incx,
    double *y,
    const suzerain_filterop_workspace *w)
{
    int info;

    /* Apply B operator out-of-place */
    info = suzerain_filterop_apply(alpha, x, incx, y, /* contiguous */ 1, w);
    if (info) return info;

    info = suzerain_filterop_solve(/* one vector */ 1, y, w->n, w);

    return info;
}
