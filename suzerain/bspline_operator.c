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

#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline_operator.h>
#include <suzerain/error.h>

/* Compute the BLAS-compatible offset to a(i,j) for general banded matrices */
/* a(i,j) -> storage(ku+i-j,j) where storage is column-major with LDA lda */
/* Note missing constant one compared with Fortran because C 0-indexes arrays */
#define GB_OFFSET(lda, ku, i, j) ((j)*(lda)+(ku+i-j))

int
suzerain_bspline_operator_create(suzerain_bspline_operator_workspace *w);

suzerain_bspline_operator_workspace *
suzerain_bspline_operator_alloc(int order,
                                int nderivatives,
                                int nbreakpoints,
                                const double * breakpoints,
                                enum suzerain_bspline_operator_method method)
{
    gsl_vector_const_view breakpoints_view
    = gsl_vector_const_view_array(breakpoints, nbreakpoints);

    int i;
    int bandwidth = -1 /* uninitialized */;
    suzerain_bspline_operator_workspace * w = NULL;

    /* Parameter sanity checks */
    if (order < 1) {
        SUZERAIN_ERROR_NULL("order must be at least 1", SUZERAIN_EINVAL);
    }

    if (nbreakpoints < 2) {
        SUZERAIN_ERROR_NULL("nbreakpoints must be at least 2",
                            SUZERAIN_EINVAL);
    }

    if (nderivatives < 1) {
        SUZERAIN_ERROR_NULL("nderivatives must be at least 1",
                            SUZERAIN_EINVAL);
    }

    /* Compute the bandwidth based on the supplied method and order */
    switch (method) {
    case SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE:
        /* Compute bandwidth of resulting operator matrices:
         *   order = 1, degree = 0 (piecewise constants),  bandwidth = 1
         *   order = 2, degree = 1 (piecewise linears),    bandwidth = 3
         *   order = 3, degree = 2 (piecewise quadratics), bandwidth = 3
         *   order = 4, degree = 3 (piecewise cubics),     bandwidth = 5
         *   order = 5, degree = 4 (piecewise quadratics), bandwidth = 5
         *   order = 6, degree = 5 (piecewise quintics),   bandwidth = 7
         *   ...
         **/
        /* TODO Bandwidth may be overly large for two reasons...
         * 1) Right continuity of the rightmost basis function, which may
         *    increase kl by one near lower right corner of matrices.
         * 2) For even orders, kl may be one less than ku
         */
        bandwidth = GSL_IS_ODD(order) ? order : order + 1;
        break;
    default:
        SUZERAIN_ERROR_NULL("unknown method", SUZERAIN_EINVAL);
    }
    if (bandwidth < 1) {
        SUZERAIN_ERROR_NULL("bandwidth not computed", SUZERAIN_ESANITY);
    }

    w = malloc(sizeof(suzerain_bspline_operator_workspace));
    if (w == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }

    /* Save bspline operator parameters */
    w->order        = order;
    w->nbreakpoints = nbreakpoints;
    w->nderivatives = nderivatives;
    w->n            = nbreakpoints + order - 2; /* assumes max continuity */
    w->method       = method;

    /* Storage parameters for BLAS/lapack-compatible general band matrix */
    w->kl          = (bandwidth - 1) / 2;
    w->ku          = (bandwidth - 1) / 2;
    w->lda         = w->kl + w->ku + 1;
    w->storagesize = w->lda * w->n;

    /* Setup workspace to use GSL B-spline functionality */
    w->bw = gsl_bspline_alloc(w->order, w->nbreakpoints);
    if (w->bw == NULL) {
        free(w);
        SUZERAIN_ERROR_NULL("failure allocating bspline workspace",
                            SUZERAIN_ENOMEM);
    }
    if (gsl_bspline_knots(&breakpoints_view.vector, w->bw)) {
        gsl_bspline_free(w->bw);
        free(w);
        SUZERAIN_ERROR_NULL("failure seting bspline breakpoints",
                            SUZERAIN_EFAILED);
    }

    /* Allocate space for pointers to matrices */
    w->D = malloc((w->nderivatives + 1) * sizeof(double *));
    if (w->D == NULL) {
        gsl_bspline_free(w->bw);
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix pointers",
                            SUZERAIN_ENOMEM);
    }
    /* Allocate memory for all matrices in one contiguous block */
    /* Memory aligned per MKL user guide numerical stability suggestion */
    if (posix_memalign((void **) &(w->D[0]),
                       16 /* byte boundary */,
                       (w->nderivatives + 1)*w->storagesize*sizeof(double))) {
        free(w->D);
        gsl_bspline_free(w->bw);
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }
    /* w->D[0] now points to D[0] which will be the mass matrix */
    /* Establish pointers for D[1] = d/dx, d[2] = d^2/dx^2, ...*/
    for (i = 0; i < nderivatives; ++i) {
        w->D[i+1] = w->D[i] + w->storagesize;
    }

    /* Workspace is fully assembled below here */
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
    free(w->D[0]);
    /* D[1], ..., D[nderivatives-1] allocated through w->D[0]; no free() */
    free(w->D);
    gsl_bspline_free(w->bw);
    free(w);
}

int
suzerain_bspline_operator_apply(
    int nderivative,
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_operator_workspace *w)
{
    double *scratch;
    int i;

    if (nderivative < 0 || w->nderivatives < nderivative) {
        SUZERAIN_ERROR("nderivative out of range", SUZERAIN_EINVAL);
    }
    if (ldb < w->n) {
        SUZERAIN_ERROR("ldb < w->n", SUZERAIN_EINVAL);
    }

    /* Allocate scratch space */
    if (posix_memalign((void **) &(scratch),
                       16 /* byte boundary */,
                       w->n*sizeof(double))) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (i = 0; i < nrhs; ++i) {
        double * const bi = b + i*ldb;
        /* Compute bi := w->D[nderivative]*bi */
        suzerain_blas_dcopy(w->n, bi, 1, scratch, 1);
        suzerain_blas_dgbmv(
            'N', w->n, w->n, w->kl, w->ku,
            1.0, w->D[nderivative], w->lda,
            scratch, 1,
            0.0, bi, 1);
    }

    free(scratch);
    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_operator_create(suzerain_bspline_operator_workspace *w)
{
    gsl_matrix *db;
    gsl_bspline_deriv_workspace *bdw;

    /* Setup workspaces to use GSL B-spline functionality */
    db = gsl_matrix_alloc(w->order, w->nderivatives + 1);
    if (db == NULL) {
        SUZERAIN_ERROR("failure allocating dB working matrix",
                       SUZERAIN_ENOMEM);
    }
    bdw = gsl_bspline_deriv_alloc(gsl_bspline_order(w->bw));
    if (bdw == NULL) {
        gsl_matrix_free(db);
        SUZERAIN_ERROR("failure allocating bspline derivative workspace",
                       SUZERAIN_ENOMEM);
    }

    /* Clear operator storage; zeros out values not explicitly set below */
    memset(w->D[0], 0, (w->nderivatives + 1) * w->storagesize * sizeof(double));

    /* Compute the operator matrices based on the supplied method */
    switch (w->method) {
    case SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available */
        /* For now, continue since only one method is implemented */
        break;
    default:
        gsl_bspline_deriv_free(bdw);
        gsl_matrix_free(db);
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    /* Evaluate basis functions at the Greville abscissae: d^k/dx^k B_j(\xi_i) */
    {
        gsl_bspline_workspace * const bw = (gsl_bspline_workspace *) w->bw;

        const int n             = w->n;
        const int nderivatives  = w->nderivatives;
        const int lda           = w->lda;
        const int ku            = w->ku;
        double ** const D       = w->D;

        int i, j, k;

        for (i = 0; i < n; ++i) {
            const double xi = gsl_bspline_greville_abscissa(i, bw);
            size_t jstart, jend;
            gsl_bspline_deriv_eval_nonzero(xi, nderivatives, db,
                                           &jstart, &jend, bw, bdw);
            for (k = 0; k <= nderivatives; ++k) {
                for (j = jstart; j <= jend; ++j) {
                    D[k][GB_OFFSET(lda, ku, i, j)]
                        = gsl_matrix_get(db, j - jstart, k);
                }
            }
        }
    }

    /* Tear down calls for GSL B-spline functionality */
    gsl_bspline_deriv_free(bdw);
    gsl_matrix_free(db);

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
    {
        int i;

        for (i = 0; i < w->n; ++i) {
            const double x = gsl_bspline_greville_abscissa(i, w->bw);
            coefficient_rhs[i] = SUZERAIN_FN_EVAL(function, x);
        }
    }

    return SUZERAIN_SUCCESS;
}

suzerain_bspline_operator_lu_workspace *
suzerain_bspline_operator_lu_alloc(
    const suzerain_bspline_operator_workspace *w)
{
    int i;
    suzerain_bspline_operator_lu_workspace * luw = NULL;

    luw = malloc(sizeof(suzerain_bspline_operator_lu_workspace));
    if (luw == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }

    /* Determine general banded matrix shape parameters */
    luw->n           = w->n;
    luw->kl          = w->kl;
    luw->ku          = w->kl + w->ku; /* Increase ku per DGBTRF, DGBTRS */
    luw->lda         = luw->kl + luw->ku + 1;
    luw->storagesize = luw->lda * luw->n;

    /* Allocate memory for LU factorization pivot storage */
    luw->ipiv = malloc(w->n * sizeof(int));
    if (luw->ipiv == NULL) {
        free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }

    /* Allocate memory for matrix */
    /* memory aligned per MKL user guide numerical stability suggestion */
    if (posix_memalign((void **) &(luw->A),
                       16 /* byte boundary */,
                       luw->storagesize*sizeof(double))) {
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
    free(luw->A);
    free(luw->ipiv);
    free(luw);
}

int
suzerain_bspline_operator_lu_form(
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
    if (luw->n < w->n) {
        SUZERAIN_ERROR("Incompatible workspaces: luw->n < w->n",
                       SUZERAIN_EINVAL);
    }
    if (luw->ku < w->ku + w->kl) {
        SUZERAIN_ERROR("Incompatible workspaces: luw->ku too small",
                       SUZERAIN_EINVAL);
    }

    /* Clear operator storage; zeros out values not explicitly set below */
    memset(luw->A, 0, luw->storagesize * sizeof(double));

    /* Accumulate coefficients times workspace w derivative operators */
    {
        const int luw_lda       = luw->lda;
        const int w_lda         = w->lda;
        const int w_storagesize = w->storagesize;

        double * A       = luw->A + luw->ku - w->ku; /* location of A(0,0) */
        double * const D = w->D[0];

        int i, joffset, k;

        /* Nasty loops avoid unnecessary dereferencing */
        for (joffset = 0;
             joffset < w_storagesize;
             joffset += w_lda, A += luw_lda) {

            double * const Dj = D + joffset;
            double * Ai       = A;

            for (i = 0;
                 i < w_lda;
                 ++i, ++Ai) {

                double * Dji = Dj + i;

                for (k = 0;
                     k < ncoefficients;
                     ++k, Dji += w_storagesize) {

                    *Ai += coefficients[k] * *Dji;
                }
            }
        }

    }

    /* Compute LU factorization of the just-formed operator */
    {
        const int info = suzerain_lapack_dgbtrf(luw->n,
                                                luw->n,
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
suzerain_bspline_operator_lu_solve(
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_operator_lu_workspace *luw)
{
    const int info = suzerain_lapack_dgbtrs('N',
                                            luw->n,
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
