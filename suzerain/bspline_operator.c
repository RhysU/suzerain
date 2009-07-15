/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
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
 *
 * bspline_operator.c: operator creation routines for a bspline basis
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline_operator.h>
#include <suzerain/error.h>

int
suzerain_bspline_operator_create(suzerain_bspline_operator_workspace *w);

int
suzerain_bspline_operator_bandwidths(suzerain_bspline_operator_workspace *w);

/* Compute the BLAS-compatible offset to a(i,j) for general banded matrices
 * a(i,j) -> storage(ku+i-j,j) where storage is column-major with LDA lda
 * Note missing constant one compared with Fortran because C 0-indexes arrays */
inline int gb_matrix_offset(int lda, int kl, int ku, int i, int j) {
    return ((j)*(lda)+(ku+i-j));
}

/* Determine if indices fall within the band of a general banded matrix */
inline int gb_matrix_in_band(int lda, int kl, int ku, int i, int j) {
    return ((j)-(ku) <= (i) && (i) <= (j)+(kl));
}

suzerain_bspline_operator_workspace *
suzerain_bspline_operator_alloc(int order,
                                int nderivatives,
                                int nbreakpoints,
                                const double * breakpoints,
                                enum suzerain_bspline_operator_method method)
{
    /* Parameter sanity checks */
    if (order < 1) {
        SUZERAIN_ERROR_NULL("order must be at least 1", SUZERAIN_EINVAL);
    }

    if (nbreakpoints < 2) {
        SUZERAIN_ERROR_NULL("nbreakpoints must be at least 2",
                            SUZERAIN_EINVAL);
    }

    if (nderivatives < 0) {
        SUZERAIN_ERROR_NULL("nderivatives must be at least 0",
                            SUZERAIN_EINVAL);
    }

    switch (method) {
    case SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available
         * For now, continue since only one method is implemented */
        break;
    default:
        SUZERAIN_ERROR_NULL("unknown method", SUZERAIN_ESANITY);
    }

    /* Allocate workspace */
    suzerain_bspline_operator_workspace * const w
        = malloc(sizeof(suzerain_bspline_operator_workspace));
    if (w == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }
    w->kl = malloc((nderivatives+1)*sizeof(w->kl[0]));
    if (w->kl == NULL) {
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace kl",
                            SUZERAIN_ENOMEM);
    }
    w->ku = malloc((nderivatives+1)*sizeof(w->ku[0]));
    if (w->ku == NULL) {
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace ku",
                            SUZERAIN_ENOMEM);
    }
    w->lda = malloc((nderivatives+1)*sizeof(w->lda[0]));
    if (w->lda == NULL) {
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace lda",
                            SUZERAIN_ENOMEM);
    }
    w->storagesize = malloc((nderivatives+1)*sizeof(w->storagesize[0]));
    if (w->storagesize == NULL) {
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace lda",
                            SUZERAIN_ENOMEM);
    }

    /* Setup workspaces to use GSL B-spline functionality */
    w->bw = gsl_bspline_alloc(order, nbreakpoints);
    if (w->bw == NULL) {
        free(w->storagesize);
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failure allocating bspline workspace",
                            SUZERAIN_ENOMEM);
    }
    gsl_vector_const_view breakpoints_view
        = gsl_vector_const_view_array(breakpoints, nbreakpoints);
    if (gsl_bspline_knots(&breakpoints_view.vector, w->bw)) {
        gsl_bspline_free(w->bw);
        free(w->storagesize);
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failure seting bspline breakpoints",
                            SUZERAIN_EFAILED);
    }
    w->dbw = gsl_bspline_deriv_alloc(order);
    if (w->dbw == NULL) {
        gsl_bspline_free(w->bw);
        free(w->storagesize);
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failure allocating bspline derivative workspace",
                            SUZERAIN_ENOMEM);
    }
    w->db = gsl_matrix_alloc(order, nderivatives + 1);
    if (w->db == NULL) {
        gsl_bspline_deriv_free(w->dbw);
        gsl_bspline_free(w->bw);
        free(w->storagesize);
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failure allocating db working matrix",
                            SUZERAIN_ENOMEM);
    }

    /* Save bspline operator parameters in workspace */
    w->order         = order;
    w->nbreakpoints  = nbreakpoints;
    w->nderivatives  = nderivatives;
    w->ncoefficients = gsl_bspline_ncoeffs(w->bw);
    w->method        = method;

    /* Storage parameters for BLAS/lapack-compatible general band matrix */
    if (suzerain_bspline_operator_bandwidths(w)) {
        gsl_matrix_free(w->db);
        gsl_bspline_deriv_free(w->dbw);
        gsl_bspline_free(w->bw);
        free(w->storagesize);
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failure determining operator bandwidths",
                            SUZERAIN_ESANITY);
    }
    for (int k = 0; k <= nderivatives; ++k) {
        w->lda[k]         = w->kl[k] + w->ku[k] + 1;
        w->storagesize[k] = w->lda[k] * w->ncoefficients;
    }

    /* Allocate space for pointers to matrices */
    w->D = malloc((w->nderivatives + 1) * sizeof(w->D[0]));
    if (w->D == NULL) {
        gsl_matrix_free(w->db);
        gsl_bspline_deriv_free(w->dbw);
        gsl_bspline_free(w->bw);
        free(w->storagesize);
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix pointers",
                            SUZERAIN_ENOMEM);
    }
    /* Allocate memory for each matrix separately */
    for (int k = 0; k <= w->nderivatives; ++k) {
        if (posix_memalign((void **) &(w->D[k]),
                           16 /* byte boundary */,
                           w->storagesize[k]*sizeof(w->D[0][0]))) {
            for (int i = k-1; i >= 0; ++i) free(w->D[i]);
            free(w->D);
            gsl_matrix_free(w->db);
            gsl_bspline_deriv_free(w->dbw);
            gsl_bspline_free(w->bw);
            free(w->storagesize);
            free(w->lda);
            free(w->ku);
            free(w->kl);
            free(w);
            SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                                SUZERAIN_ENOMEM);
        }
    }

    /* Workspace is fully assembled below here, can use ..._operator_free */

    /* Calculate operator matrices. */
    if (suzerain_bspline_operator_create(w)) {
        suzerain_bspline_operator_free(w);
        SUZERAIN_ERROR_NULL("Error creating operator matrices",
                            SUZERAIN_EFAILED);
    }

    return w;
}

void
suzerain_bspline_operator_free(suzerain_bspline_operator_workspace * w)
{
    if (w != NULL) {
        for (int k = 0; k <= w->nderivatives; ++k) free(w->D[k]);
        free(w->D);
        gsl_matrix_free(w->db);
        gsl_bspline_deriv_free(w->dbw);
        gsl_bspline_free(w->bw);
        free(w->storagesize);
        free(w->lda);
        free(w->ku);
        free(w->kl);
        free(w);
    }
}

int
suzerain_bspline_operator_ncoefficients(
    const suzerain_bspline_operator_workspace *w)
{
    return w->ncoefficients;
}

int
suzerain_bspline_operator_apply(
    int nderivative,
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_operator_workspace *w)
{
    if (nderivative < 0 || w->nderivatives < nderivative) {
        SUZERAIN_ERROR("nderivative out of range", SUZERAIN_EINVAL);
    }
    if (ldb < w->ncoefficients) {
        SUZERAIN_ERROR("ldb < w->ncoefficients", SUZERAIN_EINVAL);
    }

    /* Allocate scratch space */
    double * scratch;
    if (posix_memalign((void **) &(scratch),
                       16 /* byte boundary */,
                       w->ncoefficients*sizeof(scratch[0]))) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (int i = 0; i < nrhs; ++i) {
        double * const bi = b + i*ldb;
        /* Compute bi := w->D[nderivative]*bi */
        suzerain_blas_dcopy(w->ncoefficients, bi, 1, scratch, 1);
        suzerain_blas_dgbmv(
            'N', w->ncoefficients, w->ncoefficients,
            w->kl[nderivative], w->ku[nderivative],
            1.0, w->D[nderivative], w->lda[nderivative],
            scratch, 1,
            0.0, bi, 1);
    }

    free(scratch);
    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_operator_bandwidths(suzerain_bspline_operator_workspace *w) {
    /* Access bspline workspace */
    gsl_bspline_workspace * const bw
        = (gsl_bspline_workspace *) w->bw;
    gsl_bspline_deriv_workspace * const bdw
        = (gsl_bspline_deriv_workspace *) w->dbw;
    gsl_matrix * const db
        = (gsl_matrix *) w->db;


    /* Compute the operator bandwidths based on the supplied method */
    switch (w->method) {
    case SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available
         * For now, continue since only one method is implemented */
        break;
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    /* FIXME Incorrect */
    for (int k = 0; k <= w->nderivatives; ++k) {
        w->kl[k] = (GSL_IS_ODD(w->order) ? w->order : w->order + 1) / 2;
        w->ku[k] = w->kl[k];
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_operator_create(suzerain_bspline_operator_workspace *w)
{
    /* Clear operator storage; zeros out values not explicitly set below */
    for (int k = 0; k <= w->nderivatives; ++k) {
        memset(w->D[k], 0, w->storagesize[k]*sizeof(w->D[0][0]));
    }

    /* Compute the operator matrices based on the supplied method */
    switch (w->method) {
    case SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available
         * For now, continue since only one method is implemented */
        break;
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    /* Access bspline workspace */
    gsl_bspline_workspace * const bw
        = (gsl_bspline_workspace *) w->bw;
    gsl_bspline_deriv_workspace * const bdw
        = (gsl_bspline_deriv_workspace *) w->dbw;
    gsl_matrix * const db
        = (gsl_matrix *) w->db;

    /* Evaluate basis at the Greville abscissae: d^k/dx^k B_j(\xi_i) */

    /* Defensively dereference into local constants */
    const int n             = w->ncoefficients;
    const int nderivatives  = w->nderivatives;
    const int * const lda   = w->lda;
    const int * const ku    = w->ku;
    const int * const kl    = w->kl;
    double ** const D       = w->D;

    for (int i = 0; i < n; ++i) {
        const double xi = gsl_bspline_greville_abscissa(i, bw);
        size_t jstart, jend;
        gsl_bspline_deriv_eval_nonzero(xi, nderivatives, db,
                                       &jstart, &jend, bw, bdw);

        for (int k = 0; k <= nderivatives; ++k) {
            for (int j = jstart; j <= jend; ++j) {
                const double value = gsl_matrix_get(db, j - jstart, k);

                if (gb_matrix_in_band(lda[k], kl[k], ku[k], i, j)) {
                    const int offset
                        = gb_matrix_offset(lda[k], kl[k], ku[k], i, j);
                    D[k][offset] = value;
                } else if (value == 0.0) {
                    /* OK: value outside band is identically zero */
                } else {
                    /* NOT COOL: nonzero value outside bandwidth */
                    char buffer[255];
                    snprintf(buffer, sizeof(buffer)/sizeof(buffer[0]),
                             "encountered non-zero entry outside band"
                             " for basis spline order %d"
                             " (piecewise degree %d);"
                             " (d/dx)^%d B_%d(\\xi_%d=%g) = %g"
                             " corresponds to (row=%d, column=%d) of"
                             " general banded matrix defined by "
                             "[lda=%d, kl=%d, ku=%d, n=%d]",
                             w->order, w->order-1,
                             k, j, i, xi, value,
                             i, j,
                             lda[k], kl[k], ku[k], n);
                    SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
                }
            }
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_operator_functioncoefficient_rhs(
    const suzerain_function * function,
    double * coefficient_rhs,
    const suzerain_bspline_operator_workspace *w)
{
    switch (w->method) {
    case SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available */
        /* For now, continue since only one method is implemented */
        break;
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    /* Evaluate the function at the Greville abscissae */
    const int n = w->ncoefficients;
    for (int i = 0; i < n; ++i) {
        const double x = gsl_bspline_greville_abscissa(i, w->bw);
        coefficient_rhs[i] = SUZERAIN_FN_EVAL(function, x);
    }

    return SUZERAIN_SUCCESS;
}

suzerain_bspline_operator_lu_workspace *
suzerain_bspline_operator_lu_alloc(
    const suzerain_bspline_operator_workspace *w)
{
    suzerain_bspline_operator_lu_workspace * const luw
        = malloc(sizeof(suzerain_bspline_operator_lu_workspace));
    if (luw == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }

    /* Banded LU operator dimensions depend on largest derivative band */
    const int max_kl = gsl_stats_int_max(w->kl, 1, w->nderivatives+1);
    const int max_ku = gsl_stats_int_max(w->ku, 1, w->nderivatives+1);

    /* Determine general banded matrix shape parameters */
    luw->ncoefficients = w->ncoefficients;
    luw->kl            = max_kl;
    luw->ku            = max_kl + max_ku; /* Increase ku per DGBTRF, DGBTRS */
    luw->lda           = luw->kl + luw->ku + 1;
    luw->storagesize   = luw->lda * luw->ncoefficients;

    /* Allocate memory for LU factorization pivot storage */
    luw->ipiv = malloc(w->ncoefficients * sizeof(luw->ipiv[0]));
    if (luw->ipiv == NULL) {
        free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }

    /* Allocate memory for matrix */
    /* memory aligned per MKL user guide numerical stability suggestion */
    if (posix_memalign((void **) &(luw->A),
                       16 /* byte boundary */,
                       luw->storagesize*sizeof(luw->A[0]))) {
        free(luw->ipiv);
        free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }

    return luw;
}

void
suzerain_bspline_operator_lu_free(suzerain_bspline_operator_lu_workspace * luw)
{
    if (luw != NULL) {
        free(luw->A);
        free(luw->ipiv);
        free(luw);
    }
}

int
suzerain_bspline_operator_lu_form_general(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bspline_operator_workspace * w,
    suzerain_bspline_operator_lu_workspace *luw)
{
    if (ncoefficients < 0) {
        SUZERAIN_ERROR("Number of coefficients cannot be negative",
                       SUZERAIN_EINVAL);
    }
    if (ncoefficients > w->nderivatives + 1) {
        SUZERAIN_ERROR("More coefficients provided than derivatives available",
                       SUZERAIN_EINVAL);
    }
    if (luw->ncoefficients < w->ncoefficients) {
        SUZERAIN_ERROR("Incompatible workspaces:"
                       " luw->ncoefficients < w->ncoefficients",
                       SUZERAIN_EINVAL);
    }
    /* Banded LU operator dimensions depend on largest derivative band */
    const int max_kl = gsl_stats_int_max(w->kl, 1, w->nderivatives+1);
    const int max_ku = gsl_stats_int_max(w->ku, 1, w->nderivatives+1);
    if (luw->ku < max_kl + max_ku) {
        SUZERAIN_ERROR("Incompatible workspaces: luw->ku too small",
                       SUZERAIN_EINVAL);
    }

    /* Clear operator storage; zeros out values not explicitly set below */
    memset(luw->A, 0, luw->storagesize * sizeof(luw->A[0]));

    /* Accumulate coefficients times workspace w derivative operators */
    {
        double * A        = luw->A;
        double ** const D = w->D;

        const int luw_lda         = luw->lda;
        const int w_ncoefficients = w->ncoefficients;
        for (int k = 0; k < ncoefficients; ++k) {
            if (coefficients[k] == 0.0) {
                continue; /* Skip identically zero coefficient */
            }

            double * const D_k          = w->D[k];
            const double coefficients_k = coefficients[k];
            const int incr_ku_k         = luw->ku - w->ku[k];
            const int w_lda_k           = w->lda[k];

            for (int j = 0; j < w_ncoefficients; ++j) {
                const int joffset_D = j*w_lda_k;
                const int joffset_A = j*luw_lda + incr_ku_k;
                for (int i = 0; i < w_lda_k; ++i) {
                    A[i + joffset_A] += coefficients_k * D_k[i + joffset_D];
                }
            }
        }
    }

    /* Compute LU factorization of the just-formed operator */
    {
        const int info = suzerain_lapack_dgbtrf(luw->ncoefficients,
                                                luw->ncoefficients,
                                                luw->kl,
                                                luw->ku - luw->kl, /* NB */
                                                luw->A,
                                                luw->lda,
                                                luw->ipiv);
        if (info) {
            SUZERAIN_ERROR("suzerain_lapack_dgbtrf reported an error",
                           SUZERAIN_ESANITY);
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_operator_lu_form_mass(
    const suzerain_bspline_operator_workspace * w,
    suzerain_bspline_operator_lu_workspace *luw)
{
    const double d_one = 1.0;
    int i_one = 1;
    return suzerain_bspline_operator_lu_form_general(i_one, &d_one, w, luw);
}

int
suzerain_bspline_operator_lu_solve(
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_operator_lu_workspace *luw)
{
    const int info = suzerain_lapack_dgbtrs('N',
                                            luw->ncoefficients,
                                            luw->kl,
                                            luw->ku - luw->kl, /* NB */
                                            nrhs,
                                            luw->A,
                                            luw->lda,
                                            luw->ipiv,
                                            b,
                                            ldb);
    if (info) {
        SUZERAIN_ERROR("suzerain_lapack_dgbtrs reported an error",
                       SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}
