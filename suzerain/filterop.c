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
#include <suzerain/gbmatrix.h>
#include <suzerain/filterop.h>

// TODO Account for non-banded nature of periodicity (Toeplitz)

// Indexing utilities for banded matrices
static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

static int suzerain_filterop_operator_bandwidths_cookcabot2005(
    const double *method_params,
    int *klat,
    int *kuat,
    int *klbt,
    int *kubt)
{
    SUZERAIN_UNUSED(method_params);
    *klat = *kuat = 2;
    *klbt = *kubt = 4;
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
            break;
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
    // Unpack method-specific parameters using defaults whenever NULL
    const double alpha = method_params ? method_params[0] : 66624./100000.;
    SUZERAIN_UNUSED(alpha);

    static const int inc = 1;   // Column vectors are contiguous in memory
    const int n = w->n;         // Filtering matrices are square

    // Compute coefficients for CookCabot2005 filter
    const double alpha_0 = 1.0;
    const double alpha_1 = alpha;
    const double alpha_2 = (  1.-     alpha)/  2.;
    const double a_0     = ( 58.-105.*alpha)/128.;
    const double a_1     = ( 14.+ 11.*alpha)/ 32.;
    const double a_2     = ( 18.- 11.*alpha)/ 64.;
    const double a_3     = (  2.-  3.*alpha)/ 32.;
    const double a_4     = (- 2.+  3.*alpha)/256.;


    // Place nonzero B^T values in w->B_T using w->klbt, w->kubt, and w->ldbt.
    {
        // Banded matrix access has form a[(ku + i)*inc + j*(lda - inc)].
        // Incorporate the ku offset and decrement ld to speed indexing.
        // Further, increment kl anticipating calls like imin(m, j + kl + 1).
        int ku = w->kubt, ld = w->ldbt;
        double * bt_j = w->B_T;
        bt_j += ku*inc;
        ld -= inc;

        // With the row elements of B_T coded explicitly
        // there is no need for kl+1
        // ++kl;

        assert(w->klbt == 4);
        assert(w->kubt == 4);

        for (int j = 0; j < n; bt_j += ld, ++j) {

            // Start at j - ku, go to j + kl, Bandwidth is 4
            bt_j[(j-4)*inc] = a_4;
            bt_j[(j-3)*inc] = a_3;
            bt_j[(j-2)*inc] = a_2;
            bt_j[(j-1)*inc] = a_1;
            bt_j[(j  )*inc] = a_0;
            bt_j[(j+1)*inc] = a_1;
            bt_j[(j+2)*inc] = a_2;
            bt_j[(j+3)*inc] = a_3;
            bt_j[(j+4)*inc] = a_4;
        }
    }

    // Place nonzero A^T values in w->A_T using w->klat, w->kuat, and w->ldat.
    {
        // Again, access has form a[(ku + i)*inc + j*(lda - inc)].
        int kl = w->klat, ku = w->kuat, ld = w->ldat;
        double * at_j = w->A_T + kl; // Accounts for factorization-ready data
        at_j += ku*inc;
        ld -= inc;

        // With the row elements of A_T coded explicitly
        // there is no need for kl+1
        // ++kl;

        assert(w->klat == 2);
        assert(w->kuat == 2);

        for (int j = 0; j < n; at_j += ld, ++j) {

            // Start at j - ku, go to j + kl, Bandwidth is 2
            at_j[(j-2)*inc] = alpha_2;
            at_j[(j-1)*inc] = alpha_1;
            at_j[(j  )*inc] = alpha_0;
            at_j[(j+1)*inc] = alpha_1;
            at_j[(j+2)*inc] = alpha_2;
        }
    }

    // Note: The matrices A_T and B_T were populated in
    // BLAS-General-Band storage mode. The following loop
    // would populate elements from a matrix
    // with elements A_T[i,j] = 100*i+j:
    // for (int j = 0; j < n; at_j += ld, ++j) {
    // const int il = imax(0, j - ku), iu = imin(m, j + kl);
    //   for (int i = il; i < iu; ++i) {
    //      at_j[i*inc] =  100*i + j;
    //   }
    // }
    // where kl = w->klat + 1
    // A sample matrix of bandwidth 2 with 7 elements
    // generated with the loop above:
    // const double good_A_T[] = {
    //     // ku2      ku1     diag      kl1      kl2
    //     -55555,  -55555,       0,     100,     200,
    //     -55555,       1,     101,     201,     301,
    //          2,     102,     202,     302,     402,
    //        103,     203,     303,     403,     503,
    //        204,     304,     404,     504,     604,
    //        305,     405,     505,     605,  -55555,
    //        406,     506,     606,  -55555,  -55555
    // };

    return SUZERAIN_SUCCESS;
}

static int suzerain_filterop_operator_assemble(
    const double *method_params,
    suzerain_filterop_workspace *w)
{
    int retval;

    // Add new method-specific routines here, when necessary
    switch (w->method) {
        case SUZERAIN_FILTEROP_COOKCABOT2005:
            retval = suzerain_filterop_operator_assemble_cookcabot2005(
                        method_params, w);
            break;
        default:
            SUZERAIN_ERROR_NULL("Unknown method", SUZERAIN_ESANITY);
    }

    return retval;
}

static int suzerain_filterop_operator_boundaries(
    suzerain_filterop_workspace *w)
{
    int retval = SUZERAIN_SUCCESS;
    static const int inc = 1;   // Column vectors are contiguous in memory
    const int n = w->n;         // Filtering matrices are square

    // Add new boundary-type routines here, when necessary
    switch (w->b_first) {
        case SUZERAIN_FILTEROP_BOUNDARY_IGNORE:
            /* Nothing to be done */
            break;

        case SUZERAIN_FILTEROP_BOUNDARY_NOFILTER:
            // Replace near-boundary schemes on B^T
            {
                // Banded matrix access has form a[(ku + i)*inc + j*(lda - inc)].
                // Incorporate the ku offset and decrement ld to speed indexing.
                int ku = w->kubt, ld = w->ldbt;
                int nbs = imax(w->kuat, w->kubt); // nbs: near-boundary stencils
                double * bt_j = w->B_T;
                bt_j += ku*inc;
                ld -= inc;

                for (int j = 0; j < nbs; bt_j += ld, ++j) {

                    // Start at j - ku, go down to 1
                    for (int joff = w->kubt; joff > 0; --joff) {
                        bt_j[(j-joff)*inc] = 0.;
                    }
                    // Diagonal element
                    bt_j[(j  )*inc] = 1.;
                    // Start at 1, go to j + kl
                    for (int joff = 1; joff < w->klbt+1; ++joff) {
                        bt_j[(j+joff)*inc] = 0.;
                    }
                }
            }

            // Replace near-boundary schemes on A^T
            {
                // Again, access has form a[(ku + i)*inc + j*(lda - inc)].
                int kl = w->klat, ku = w->kuat, ld = w->ldat;
                int nbs = imax(w->klat, w->klbt); // nbs: near-boundary stencils
                double * at_j = w->A_T + kl; // Accounts for factorization-ready data
                at_j += ku*inc;
                ld -= inc;

                // The number of near-boundary schemes that need to be replaced
                // goes by the width of the stencil (widest side), in this case,
                // that of B_T
                for (int j = 0; j < nbs; at_j += ld, ++j) {

                    // Start at j - ku, go down to 1
                    for (int joff = w->kuat; joff > 0; --joff) {
                        at_j[(j-joff)*inc] = 0.;
                    }
                    // Diagonal element
                    at_j[(j  )*inc] = 1.;
                    // Start at 1, go to j + kl
                    for (int joff = 1; joff < w->klat+1; ++joff) {
                        at_j[(j+joff)*inc] = 0.;
                    }
                }
            }
            break;

        case SUZERAIN_FILTEROP_BOUNDARY_SYMMETRY: // TODO
        case SUZERAIN_FILTEROP_BOUNDARY_PERIODIC: // TODO
        default:
            SUZERAIN_ERROR_NULL("Unknown first boundary treatment",
                                SUZERAIN_ESANITY);
    }

    switch (w->b_last) {
        case SUZERAIN_FILTEROP_BOUNDARY_IGNORE:
            /* Nothing to be done */
            break;

        case SUZERAIN_FILTEROP_BOUNDARY_NOFILTER:
            // Replace near-boundary schemes on B^T
            {
                // Banded matrix access has form a[(ku + i)*inc + j*(lda - inc)].
                // Incorporate the ku offset and decrement ld to speed indexing.
                int ku = w->kubt, ld = w->ldbt;
                int nbs = imax(w->klat, w->klbt); // nbs: near-boundary stencils
                double * bt_j = w->B_T;

                // Walk the upper end forwards:
                // Point initially to element (n-nbs) in the diagonal
                bt_j += (ku*inc + (n-nbs)*(w->ldbt)*inc);
                ld -= inc;

                // Walk the matrix forwards
                for (int j = 0; j < nbs; bt_j += ld, ++j) {

                    // Start at j - ku, go down to 1
                    for (int joff = w->kubt; joff > 0; --joff) {
                        bt_j[(j-joff)*inc] = 0.;
                    }
                    // Diagonal element
                    bt_j[(j  )*inc] = 1.;
                    // Start at 1, go to j + kl
                    for (int joff = 1; joff < w->klbt+1; ++joff) {
                        bt_j[(j+joff)*inc] = 0.;
                    }
                }

                // Note: To walk the upper end backwards instead, the initial bt_j offset
                // points initially to the last diagonal element, as
                // bt_j += (ku*inc + (n-1)*(w->ldbt)*inc);
                // and the outer loop changes to
                // for (int j = 0; j > -nbs; bt_j -= ld, --j) { ... }

            }

            // Replace near-boundary schemes on A^T
            {
                // Again, access has form a[(ku + i)*inc + j*(lda - inc)].
                int kl = w->klat, ku = w->kuat, ld = w->ldat;
                int nbs = imax(w->klat, w->klbt); // nbs: near-boundary stencils
                double * at_j = w->A_T + kl; // Accounts for factorization-ready data

                // Walk the upper end forwards:
                // Point initially to element (n-nbs) in the diagonal
                at_j += (ku*inc + (n-nbs)*(w->ldat)*inc);
                ld -= inc;

                // Walk the matrix forwards
                for (int j = 0; j < nbs; at_j += ld, ++j) {

                    // Start at j - ku, go down to 1
                    for (int joff = w->kuat; joff > 0; --joff) {
                        at_j[(j-joff)*inc] = 0.;
                    }
                    // Diagonal element
                    at_j[(j  )*inc] = 1.;
                    // Start at 1, go to j + kl
                    for (int joff = 1; joff < w->klat+1; ++joff) {
                        at_j[(j+joff)*inc] = 0.;
                    }
                }

                // Note: To walk the upper end backwards instead, the initial at_j offset
                // points initially to the last diagonal element, as
                // at_j += (ku*inc + (n-1)*(w->ldat)*inc);
                // and the outer loop changes to
                // for (int j = 0; j > -nbs; at_j -= ld, --j) { ... }
            }
            break;

        case SUZERAIN_FILTEROP_BOUNDARY_SYMMETRY: // TODO
        case SUZERAIN_FILTEROP_BOUNDARY_PERIODIC: // TODO
        default:
            SUZERAIN_ERROR_NULL("Unknown last boundary treatment",
                                SUZERAIN_ESANITY);
    }

    return retval;
}

suzerain_filterop_workspace *
suzerain_filterop_alloc(
    const int n,
    const enum suzerain_filterop_method method,
    const double *method_params,
    const enum suzerain_filterop_boundary_treatment b_first,
    const enum suzerain_filterop_boundary_treatment b_last)
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
            = calloc(1, sizeof(suzerain_filterop_workspace));
    if (w == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }

    /* Populate initial workspace values */
    w->method  = method;
    w->b_first = b_first;
    w->b_last  = b_last;
    w->n       = n;
    w->klat    = klat;
    w->kuat    = kuat;
    w->ldat    = ldat;
    w->klbt    = klbt;
    w->kubt    = kubt;
    w->ldbt    = ldbt;

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
                method_params, w)) {
        suzerain_filterop_free(w);
        SUZERAIN_ERROR_NULL("failed to assemble A_T and B_T",
                            SUZERAIN_EFAILED);
    }

    /* Apply chosen boundary treatment */
    if (SUZERAIN_SUCCESS != suzerain_filterop_operator_boundaries(w)) {
        suzerain_filterop_free(w);
        SUZERAIN_ERROR_NULL("failed to modify A_T and B_T for boundaries",
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
    free(w);
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
    // TODO Periodic boundary treatments require additional processing

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
    // TODO Periodic boundary treatments require additional processing

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
