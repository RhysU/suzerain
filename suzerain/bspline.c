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
 * bspline.c: bspline basis manipulation and operator routines
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
#include <suzerain/bspline.h>
#include <suzerain/error.h>

int
suzerain_bspline_determine_operator_bandwidths(suzerain_bspline_workspace *w);

int
suzerain_bspline_create_operator(suzerain_bspline_workspace *w);

int
compute_banded_collocation_derivative_submatrix(
    const int ioffset,
    const int joffset,
    const int nderivatives,
    const int * const kl,
    const int * const ku,
    const int * const lda,
    const int npoints,
    const double * points,
    gsl_bspline_workspace * bw,
    gsl_bspline_deriv_workspace * dbw,
    gsl_matrix * db,
    double ** const D);

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

/* Sum an array of integers */
inline int sum_int(size_t n, const int * x, size_t incx) {
    const int  * const x_end = x + n*incx;
    int sum = 0;
    while (x < x_end) {
        sum += *x;
        x   += incx;
    }
    return sum;
}

suzerain_bspline_workspace *
suzerain_bspline_alloc(int order,
                                int nderivatives,
                                int nbreakpoints,
                                const double * breakpoints,
                                enum suzerain_bspline_method method)
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
    case SUZERAIN_BSPLINE_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available
         * For now, continue since only one method is implemented */
        break;
    default:
        SUZERAIN_ERROR_NULL("unknown method", SUZERAIN_ESANITY);
    }

    /* Allocate workspace */
    suzerain_bspline_workspace * const w
        = malloc(sizeof(suzerain_bspline_workspace));
    if (w == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }
    /* Make workspace pointers NULL */
    w->kl          = NULL;
    w->ku          = NULL;
    w->lda         = NULL;
    w->storagesize = NULL;
    w->bw          = NULL;
    w->dbw         = NULL;
    w->db          = NULL;
    w->D           = NULL;
    /* Prepare workspace */
    w->kl = malloc((nderivatives+1)*sizeof(w->kl[0]));
    if (w->kl == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace kl",
                            SUZERAIN_ENOMEM);
    }
    w->ku = malloc((nderivatives+1)*sizeof(w->ku[0]));
    if (w->ku == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace ku",
                            SUZERAIN_ENOMEM);
    }
    w->lda = malloc((nderivatives+1)*sizeof(w->lda[0]));
    if (w->lda == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace lda",
                            SUZERAIN_ENOMEM);
    }
    w->storagesize = malloc((nderivatives+1)*sizeof(w->storagesize[0]));
    if (w->storagesize == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate workspace lda",
                            SUZERAIN_ENOMEM);
    }

    /* Setup workspaces to use GSL B-spline functionality */
    w->bw = gsl_bspline_alloc(order, nbreakpoints);
    if (w->bw == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failure allocating bspline workspace",
                            SUZERAIN_ENOMEM);
    }
    gsl_vector_const_view breakpoints_view
        = gsl_vector_const_view_array(breakpoints, nbreakpoints);
    if (gsl_bspline_knots(&breakpoints_view.vector, w->bw)) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failure seting bspline breakpoints",
                            SUZERAIN_EFAILED);
    }
    w->dbw = gsl_bspline_deriv_alloc(order);
    if (w->dbw == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failure allocating bspline derivative workspace",
                            SUZERAIN_ENOMEM);
    }
    w->db = gsl_matrix_alloc(order, nderivatives + 1);
    if (w->db == NULL) {
        suzerain_bspline_free(w);
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
    if (suzerain_bspline_determine_operator_bandwidths(w)) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failure determining operator bandwidths",
                            SUZERAIN_ESANITY);
    }
    /* Precompute derived storage parameters */
    for (int k = 0; k <= nderivatives; ++k) {
        w->lda[k]         = w->kl[k] + w->ku[k] + 1;
        w->storagesize[k] = w->lda[k] * w->ncoefficients;
    }

    /* Allocate space for pointers to matrices */
    w->D = malloc((w->nderivatives + 1) * sizeof(w->D[0]));
    if (w->D == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix pointers",
                            SUZERAIN_ENOMEM);
    }
    /* Allocate one block for all derivative operator matrices */
    const size_t total_storage
        = sum_int(w->nderivatives+1, w->storagesize, 1);
    if (posix_memalign((void **) &(w->D[0]),
                       16 /* byte boundary */,
                       total_storage*sizeof(w->D[0][0]))) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }
    for (int k = 0; k < w->nderivatives; ++k) {
        w->D[k+1] = w->D[k] + w->storagesize[k];
    }

    /* Calculate operator matrices. */
    if (suzerain_bspline_create_operator(w)) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("Error creating operator matrices",
                            SUZERAIN_EFAILED);
    }

    return w;
}

void
suzerain_bspline_free(suzerain_bspline_workspace * w)
{
    if (w != NULL) {
        if (w->D != NULL) {
            free(w->D[0]);
            w->D[0] = NULL;

            free(w->D);
            w->D = NULL;
        }

        gsl_matrix_free(w->db);
        w->db = NULL;

        gsl_bspline_deriv_free(w->dbw);
        w->dbw = NULL;

        gsl_bspline_free(w->bw);
        w->bw = NULL;

        free(w->storagesize);
        w->storagesize = NULL;

        free(w->lda);
        w->lda = NULL;

        free(w->ku);
        w->ku = NULL;

        free(w->kl);
        w->kl = NULL;

        free(w);
    }
}

int
suzerain_bspline_ncoefficients(
    const suzerain_bspline_workspace *w)
{
    return w->ncoefficients;
}

int
suzerain_bspline_apply_operator(
    int nderivative,
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_workspace *w)
{
    /* Parameter sanity checks */
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
suzerain_bspline_evaluate(
    int nderivative,
    const double * coefficients,
    int npoints,
    const double * points,
    double * values,
    int ldvalues,
    const suzerain_bspline_workspace *w)
{
    /* Parameter sanity checks */
    if (nderivative < 0 || w->nderivatives < nderivative) {
        SUZERAIN_ERROR("nderivative out of range", SUZERAIN_EINVAL);
    }

    /* Dereference workspace pointers */
    gsl_matrix * const db = w->db;
    gsl_bspline_workspace * const bw = w->bw;
    gsl_bspline_deriv_workspace * const dbw = w->dbw;

    /* Dereference fixed db parameters; db is row-major per gsl_matrix */
    /* See GSL manual section 8.4.2 for details on gsl_matrix layout */
    double * const db_data = db->data;
    int db_tda = db->tda;

    /* bspline support is always the piecewise degree plus one, so */
    /* we can determine the number of nonzero basis functions outside loops */
    const int dotlength = w->order;

    /* ldvalues == 0 signals that we only want derivative nderivative */
    const int kstart = (ldvalues == 0) ? nderivative : 0;

    size_t jstart, jend;
    for (int i = 0; i < npoints; ++i) {

        gsl_bspline_deriv_eval_nonzero(points[i], nderivative,
                db, &jstart, &jend, bw, dbw);

        const double * coefficient_start = coefficients + jstart;

        for (int k = kstart; k <= nderivative; ++k) {
            const double * db_start = db_data + k;
            const double value = suzerain_blas_ddot(
                    dotlength, coefficient_start, 1, db_start, db_tda);
            const int storage_offset = i + k*ldvalues;
            values[storage_offset] = value;
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_determine_operator_bandwidths(suzerain_bspline_workspace *w)
{
    /* Compute the operator bandwidths based on the supplied method */
    switch (w->method) {
    case SUZERAIN_BSPLINE_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available
         * For now, continue since only one method is implemented */
        break;
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    /* Initially maximum collocation operator band storage parameters */
    for (int k = 0; k <= w->nderivatives; ++k) {
        w->kl[k]          = w->order - 1;
        w->ku[k]          = w->kl[k];
        w->lda[k]         = w->kl[k] + w->ku[k] + 1;
        w->storagesize[k] = w->lda[k] * w->order;
    }

    /* Compute collocation points at which we will check bandwidth */
    double * const points = malloc(2*w->order*sizeof(w->bw->knots->data[0]));
    if (points == NULL) {
        SUZERAIN_ERROR("Unable to allocate space for collocation points",
                       SUZERAIN_ENOMEM);
    }
    for (int i = 0; i < w->order; ++i) {
        points[i] = gsl_bspline_greville_abscissa(i, w->bw);
        points[i + w->order] = gsl_bspline_greville_abscissa(
                w->ncoefficients-w->order+i,w->bw);
    }

    /* Allocate space for collocation operator submatrices */
    double ** const scratch
        = malloc(2*(w->nderivatives+1) * sizeof(scratch[0]));
    if (scratch == NULL) {
        free(points);
        SUZERAIN_ERROR("Unable to allocate scratch pointers",
                SUZERAIN_ENOMEM);
    }
    const size_t total_storage
        = 2 * sum_int(w->nderivatives+1, w->storagesize, 1);
    if (posix_memalign((void **) &(scratch[0]),
                       16 /* byte boundary */,
                       total_storage*sizeof(w->D[0][0]))) {
        free(scratch);
        free(points);
        SUZERAIN_ERROR_NULL("failed to allocate scratch",
                            SUZERAIN_ENOMEM);
    }
    for (int k = 1; k < 2*(w->nderivatives+1); ++k) {
        const int storagesize_km1
            = w->storagesize[(k-1) % (w->nderivatives+1)];
        scratch[k] = scratch[k-1] + storagesize_km1;
    }

    /* Create convenience views of our working space */
    const double * const ul_points = &(points[0]);
    const double * const lr_points = &(points[w->order]);
    double ** const ul_D = &(scratch[0]);
    double ** const lr_D = &(scratch[w->nderivatives+1]);

    /* Compute the upper-left- and lower-right- most submatrices */
    if (compute_banded_collocation_derivative_submatrix(
             0, 0,
             w->nderivatives, w->kl, w->ku, w->lda,
             w->order, ul_points, w->bw, w->dbw, w->db, ul_D)) {
        free(scratch[0]);
        free(scratch);
        free(points);
        SUZERAIN_ERROR("Error computing operator UL submatrices",
                       SUZERAIN_ESANITY);
    }
    const int lr_offset = w->ncoefficients - w->order;
    if (compute_banded_collocation_derivative_submatrix(
             lr_offset, lr_offset,
             w->nderivatives, w->kl, w->ku, w->lda,
             w->order, lr_points, w->bw, w->dbw, w->db, lr_D)) {
        free(scratch[0]);
        free(scratch);
        free(points);
        SUZERAIN_ERROR("Error computing operator LR submatrices",
                       SUZERAIN_ESANITY);
    }

    /* Reduce kl/ku for each zero off-diagonal in ul_D and lr_D */
    for (int k = 0; k <= w->nderivatives; ++k) {
        const int fixed_ku_k  = w->ku[k];
        const int fixed_kl_k  = w->kl[k];
        const int fixed_lda_k = w->lda[k];

        for (int i=0; i < fixed_ku_k; ++i) {
            const double ul_asum
                = suzerain_blas_dasum(w->order, ul_D[k]+i, fixed_lda_k);
            const double lr_asum
                = suzerain_blas_dasum(w->order, lr_D[k]+i, fixed_lda_k);
            if (ul_asum + lr_asum == 0.0) {
                --w->ku[k];
            } else {
                break; /* Skip all after nonzero superdiagonal */
            }
        }
        for (int i=fixed_lda_k-1; i > fixed_ku_k; --i) {
            const double ul_asum
                = suzerain_blas_dasum(w->order, ul_D[k]+i, fixed_lda_k);
            const double lr_asum
                = suzerain_blas_dasum(w->order, lr_D[k]+i, fixed_lda_k);
            if (ul_asum + lr_asum == 0.0) {
                --w->kl[k];
            } else {
                break; /* Skip all after nonzero subdiagonal */
            }
        }
    }

    free(scratch[0]);
    free(scratch);
    free(points);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_create_operator(suzerain_bspline_workspace *w)
{
    /* Compute the operator matrices based on the supplied method */
    switch (w->method) {
    case SUZERAIN_BSPLINE_COLLOCATION_GREVILLE:
        /* Logic will need to change here once multiple methods available
         * For now, continue since only one method is implemented */
        break;
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    /* Evaluate basis at the Greville abscissae: d^k/dx^k B_j(\xi_i) */
    double * const points
        = malloc(w->ncoefficients*sizeof(w->bw->knots->data[0]));
    if (points == NULL) {
        SUZERAIN_ERROR("Unable to allocate space for Greville abscissae",
                       SUZERAIN_ENOMEM);
    }
    for (int i = 0; i < w->ncoefficients; ++i) {
        points[i] = gsl_bspline_greville_abscissa(i, w->bw);
    }

    /* Compute the full derivative operator matrices */
    if (compute_banded_collocation_derivative_submatrix(
             0, 0, w->nderivatives, w->kl, w->ku, w->lda,
             w->ncoefficients, points, w->bw, w->dbw, w->db, w->D)) {
        free(points);
        SUZERAIN_ERROR("Error computing operator matrices", SUZERAIN_ESANITY);
    }

    free(points);

    return SUZERAIN_SUCCESS;
}

int
compute_banded_collocation_derivative_submatrix(
    const int ioffset,
    const int joffset,
    const int nderivatives,
    const int * const kl,
    const int * const ku,
    const int * const lda,
    const int npoints,
    const double * points,
    gsl_bspline_workspace * bw,
    gsl_bspline_deriv_workspace * dbw,
    gsl_matrix * db,
    double ** const D)
{
    /* Clear operator storage; zeros out values not explicitly set below */
    for (int k = 0; k <= nderivatives; ++k) {
        memset(D[k], 0, lda[k]*npoints*sizeof(D[0][0]));
    }

    /* Fill all nonzero entries in the operator storage */
    for (int i = 0; i < npoints; ++i) {
        size_t dbjstart, dbjend;
        gsl_bspline_deriv_eval_nonzero(points[i], nderivatives, db,
                                       &dbjstart, &dbjend, bw, dbw);

        /* Coerce dbjstart/dbjend to stay within the submatrix of interest */
        const size_t jstart = GSL_MAX(dbjstart, joffset);
        const size_t jend   = GSL_MIN(dbjend, npoints-1+joffset);

        for (int k = 0; k <= nderivatives; ++k) {
            for (int j = jstart; j <= jend; ++j) {
                const double value = gsl_matrix_get(db, j - dbjstart, k);
                const int in_band  = gb_matrix_in_band(
                        lda[k], kl[k], ku[k], i+ioffset, j /* no joffset */);
                const int offset = gb_matrix_offset(
                        lda[k], kl[k], ku[k], i /* no ioffset */, j-joffset);

                if (in_band) {
                    D[k][offset] = value;
                } else if (value == 0.0) {
                    /* OK: value outside band is identically zero */
                } else {
                    /* NOT COOL: nonzero value outside bandwidth */
                    const int order = bw->k;
                    char buffer[384];
                    snprintf(buffer, sizeof(buffer)/sizeof(buffer[0]),
                             "non-zero outside band"
                             " of basis spline order %d"
                             " (piecewise degree %d);"
                             " (d/dx)^%d B_%d(\\points_%d=%g) = %g"
                             " in (row=%d, column=%d, offset=%d) of"
                             " general band matrix"
                             " [lda=%d, kl=%d, ku=%d, n=%d]"
                             " with offset (row=%d, column=%d)",
                             order, order-1,
                             k, j, i, points[i], value,
                             i, j, offset,
                             lda[k], kl[k], ku[k], npoints,
                             ioffset, joffset);
                    SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
                }
            }
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_find_coefficient_rhs(
    const suzerain_function * function,
    double * coefficient_rhs,
    const suzerain_bspline_workspace *w)
{
    switch (w->method) {
    case SUZERAIN_BSPLINE_COLLOCATION_GREVILLE:
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

suzerain_bspline_lu_workspace *
suzerain_bspline_lu_alloc(
    const suzerain_bspline_workspace *w)
{
    /* Allocate space for the luw workspace */
    suzerain_bspline_lu_workspace * const luw
        = malloc(sizeof(suzerain_bspline_lu_workspace));
    if (luw == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }
    /* Make workspace pointers NULL */
    luw->A    = NULL;
    luw->ipiv = NULL;

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
        suzerain_bspline_lu_free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }

    /* Allocate memory for the matrix formed within lu_form_general */
    if (posix_memalign((void **) &(luw->A),
                       16 /* byte boundary */,
                       luw->storagesize*sizeof(luw->A[0]))) {
        suzerain_bspline_lu_free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }

    return luw;
}

void
suzerain_bspline_lu_free(suzerain_bspline_lu_workspace * luw)
{
    if (luw != NULL) {
        free(luw->A);
        luw->A = NULL;

        free(luw->ipiv);
        luw->ipiv = NULL;

        free(luw);
    }
}

int
suzerain_bspline_lu_form_general(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bspline_workspace * w,
    suzerain_bspline_lu_workspace *luw)
{
    /* Parameter sanity checks */
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
suzerain_bspline_lu_form_mass(
    const suzerain_bspline_workspace * w,
    suzerain_bspline_lu_workspace *luw)
{
    const double d_one = 1.0;
    int i_one = 1;
    return suzerain_bspline_lu_form_general(i_one, &d_one, w, luw);
}

int
suzerain_bspline_lu_solve(
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_lu_workspace *luw)
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
