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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline.h>
#include <suzerain/error.h>
#include <suzerain/macros.h>

int
suzerain_bspline_determine_operator_bandwidths(suzerain_bspline_workspace *w);

int
suzerain_bspline_create_operators(suzerain_bspline_workspace *w);

int
compute_banded_collocation_derivative_submatrix(
    const int ioffset,
    const int joffset,
    const int nderivatives,
    const int * const kl,
    const int * const ku,
    const int ld,
    const int npoints,
    const double * points,
    gsl_bspline_workspace * bw,
    gsl_bspline_deriv_workspace * dbw,
    gsl_matrix * db,
    double ** const D);

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
    w->ndof          = gsl_bspline_ncoeffs(w->bw);
    w->method        = method;

    /* Storage parameters for BLAS/lapack-compatible general band matrix */
    if (suzerain_bspline_determine_operator_bandwidths(w)) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failure determining operator bandwidths",
                            SUZERAIN_ESANITY);
    }
    /* Compute derived storage parameters.
     *
     * We choose to have all operators stored using the same ld and
     * capable of being manipulated using w->max_kl and w->max_ku.  This
     * is not as space efficient as it could be in the pathological cases,
     * but it does allow using BLAS dgb_acc to form linear combinations
     * of operators.
     */
    w->max_kl = gsl_stats_int_max(w->kl, 1, w->nderivatives+1);
    w->max_ku = gsl_stats_int_max(w->ku, 1, w->nderivatives+1);
    w->ld     = w->max_ku + 1 + w->max_kl;

    /* Allocate space for pointers to matrices */
    w->D = malloc((w->nderivatives + 1) * sizeof(w->D[0]));
    if (w->D == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix pointers",
                            SUZERAIN_ENOMEM);
    }
    /* Allocate one block for all derivative operator matrices */
    const size_t an_operators_storage = w->ld * w->ndof;
    const size_t all_operators_storage
        = an_operators_storage * (w->nderivatives+1);
    w->D[0] = suzerain_blas_malloc(all_operators_storage*sizeof(w->D[0][0]));
    if (w->D[0] == NULL) {
        suzerain_bspline_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }
    /* Compute pointers to the beginning of each derivative operator */
    /* Determine the start of raw storage for (k+1)-th operator */
    for (int k = 1; k <= w->nderivatives; ++k) {
        w->D[k] = w->D[k-1] + an_operators_storage;
    }
    /* Step pointer past any unused superdiagonals in the k-th operator */
    for (int k = 0; k <= w->nderivatives; ++k) {
        w->D[k] += (w->max_ku - w->ku[k]);
    }

    /* Calculate operator matrices */
    if (suzerain_bspline_create_operators(w)) {
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
            if (w->D[0] != NULL) {
                // Must free originally allocated memory offset, not the
                // computed offset.  See suzerain_bspline_alloc.
                w->D[0] -= (w->max_ku - w->ku[0]);
                free(w->D[0]);
                w->D[0] = NULL;
            }

            free(w->D);
            w->D = NULL;
        }

        gsl_matrix_free(w->db);
        w->db = NULL;

        gsl_bspline_deriv_free(w->dbw);
        w->dbw = NULL;

        gsl_bspline_free(w->bw);
        w->bw = NULL;

        free(w->ku);
        w->ku = NULL;

        free(w->kl);
        w->kl = NULL;

        free(w);
    }
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
    if (ldb < w->ndof) {
        SUZERAIN_ERROR("ldb < w->ndof", SUZERAIN_EINVAL);
    }

    /* Allocate scratch space; Required because BLAS operations' behavior on
     * aliased pointers is undefined. */
    double * const scratch
        = suzerain_blas_malloc(w->ndof*sizeof(scratch[0]));
    if (scratch == NULL) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (int i = 0; i < nrhs; ++i) {
        double * const bi = b + i*ldb;
        /* Compute bi := w->D[nderivative]*bi */
        suzerain_blas_dcopy(w->ndof, bi, 1, scratch, 1);
        suzerain_blas_dgbmv(
            'N', w->ndof, w->ndof,
            w->kl[nderivative], w->ku[nderivative],
            1.0, w->D[nderivative], w->ld,
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
    if (ldvalues && 0 < ldvalues && ldvalues < npoints) {
        SUZERAIN_ERROR("ldvalues too small for npoints", SUZERAIN_EINVAL);
    }

    /* Dereference workspace pointers */
    gsl_matrix * const db = w->db;
    gsl_bspline_workspace * const bw = w->bw;
    gsl_bspline_deriv_workspace * const dbw = w->dbw;

    /* Dereference fixed db parameters; db is row-major per gsl_matrix */
    /* See GSL manual section 8.4.2 for details on gsl_matrix layout */
    double * const db_data = db->data;
    const int db_tda = db->tda;

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

    /* Initially set maximum collocation operator band storage parameters */
    for (int k = 0; k <= w->nderivatives; ++k) {
        w->kl[k] = w->order - 1;
        w->ku[k] = w->kl[k];
    }
    /* Determine storage needs, defensively assuming nonuniform parameters */
    w->max_kl = gsl_stats_int_max(w->kl, 1, w->nderivatives+1);
    w->max_ku = gsl_stats_int_max(w->ku, 1, w->nderivatives+1);
    w->ld     = w->max_kl + 1 + w->max_ku;

    /* Compute collocation points at which we will check bandwidth */
    double * const points = malloc(2*w->order*sizeof(w->bw->knots->data[0]));
    if (points == NULL) {
        SUZERAIN_ERROR("Unable to allocate space for collocation points",
                       SUZERAIN_ENOMEM);
    }
    for (int i = 0; i < w->order; ++i) {
        points[i] = gsl_bspline_greville_abscissa(i, w->bw);
        points[i + w->order] = gsl_bspline_greville_abscissa(
                w->ndof-w->order+i, w->bw);
    }

    /* Allocate space for collocation operator submatrices */
    double ** const scratch
        = malloc(2*(w->nderivatives+1) * sizeof(scratch[0]));
    if (scratch == NULL) {
        free(points);
        SUZERAIN_ERROR("Unable to allocate scratch pointers",
                SUZERAIN_ENOMEM);
    }
    const size_t an_operators_storage = w->ld * w->order;
    const size_t total_storage = 2*(w->nderivatives+1)*an_operators_storage;
    scratch[0] = suzerain_blas_calloc(total_storage, sizeof(w->D[0][0]));
    if (scratch[0] == NULL) {
        free(scratch);
        free(points);
        SUZERAIN_ERROR_NULL("failed to allocate scratch",
                            SUZERAIN_ENOMEM);
    }
    for (int k = 1; k < 2*(w->nderivatives+1); ++k) {
        scratch[k] = scratch[k-1] + an_operators_storage;
    }

    /* Create convenience views of our working space */
    const double * const ul_points = &(points[0]);
    const double * const lr_points = &(points[w->order]);
    double ** const ul_D = &(scratch[0]);
    double ** const lr_D = &(scratch[w->nderivatives+1]);

    /* Compute the upper-left- and lower-right- most submatrices */
    if (compute_banded_collocation_derivative_submatrix(
             0, 0,
             w->nderivatives, w->kl, w->ku, w->ld,
             w->order, ul_points, w->bw, w->dbw, w->db, ul_D)) {
        free(scratch[0]);
        free(scratch);
        free(points);
        SUZERAIN_ERROR("Error computing operator UL submatrices",
                       SUZERAIN_ESANITY);
    }
    const int lr_offset = w->ndof - w->order;
    if (compute_banded_collocation_derivative_submatrix(
             lr_offset, lr_offset,
             w->nderivatives, w->kl, w->ku, w->ld,
             w->order, lr_points, w->bw, w->dbw, w->db, lr_D)) {
        free(scratch[0]);
        free(scratch);
        free(points);
        SUZERAIN_ERROR("Error computing operator LR submatrices",
                       SUZERAIN_ESANITY);
    }

    /* Reduce kl/ku for each zero off-diagonal in ul_D and lr_D */
    for (int k = 0; k <= w->nderivatives; ++k) {
        const int fixed_ku_k = w->ku[k];
        const int fixed_kl_k = w->kl[k];
        const int fixed_ld   = w->ld;

        for (int i=0; i < fixed_ku_k; ++i) {
            const double ul_asum
                = suzerain_blas_dasum(w->order, ul_D[k]+i, fixed_ld);
            const double lr_asum
                = suzerain_blas_dasum(w->order, lr_D[k]+i, fixed_ld);
            if (ul_asum + lr_asum == 0.0) {
                --w->ku[k];
            } else {
                break; /* Skip all after nonzero superdiagonal */
            }
        }
        for (int i=fixed_ld-1; i > fixed_ku_k; --i) {
            const double ul_asum
                = suzerain_blas_dasum(w->order, ul_D[k]+i, fixed_ld);
            const double lr_asum
                = suzerain_blas_dasum(w->order, lr_D[k]+i, fixed_ld);
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
suzerain_bspline_create_operators(suzerain_bspline_workspace *w)
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
        = malloc(w->ndof*sizeof(w->bw->knots->data[0]));
    if (points == NULL) {
        SUZERAIN_ERROR("Unable to allocate space for Greville abscissae",
                       SUZERAIN_ENOMEM);
    }
    for (int i = 0; i < w->ndof; ++i) {
        points[i] = gsl_bspline_greville_abscissa(i, w->bw);
    }

    /* Zero the full derivative operator matrices */
    w->D[0] -= (w->max_ku - w->ku[0]); /* See suzerain_bspline_alloc */
    memset(w->D[0], 0, w->ld*w->ndof*(w->nderivatives+1)*sizeof(w->D[0][0]));
    w->D[0] += (w->max_ku - w->ku[0]);

    /* Compute the full derivative operator matrices */
    if (compute_banded_collocation_derivative_submatrix(
             0, 0, w->nderivatives, w->kl, w->ku, w->ld,
             w->ndof, points, w->bw, w->dbw, w->db, w->D)) {
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
    const int ld,
    const int npoints,
    const double * points,
    gsl_bspline_workspace * bw,
    gsl_bspline_deriv_workspace * dbw,
    gsl_matrix * db,
    double ** const D)
{
    /* PRECONDITION: D[0]...D[nderivatives] must have been filled with zeros */

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
                const int in_band  = suzerain_gbmatrix_in_band(
                        ld, kl[k], ku[k], i+ioffset, j /* no joffset */);
                const int offset = suzerain_gbmatrix_offset(
                        ld, kl[k], ku[k], i /* no ioffset */, j-joffset);

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
                             " [ld=%d, kl=%d, ku=%d, n=%d]"
                             " with offset (row=%d, column=%d)",
                             order, order-1,
                             k, j, i, points[i], value,
                             i, j, offset,
                             ld, kl[k], ku[k], npoints,
                             ioffset, joffset);
                    SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
                }
            }
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bspline_find_interpolation_problem_rhs(
    const suzerain_function * function,
    double * rhs,
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
    const int n = w->ndof;
    for (int i = 0; i < n; ++i) {
        const double x = gsl_bspline_greville_abscissa(i, w->bw);
        rhs[i] = SUZERAIN_FN_EVAL(function, x);
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

    /* Determine general banded matrix shape parameters */
    luw->ndof          = w->ndof;
    luw->kl            = w->max_kl;
    luw->ku            = w->max_kl + w->max_ku; /* Increase per GBTRF, GBTRS */
    luw->ld            = luw->kl + luw->ku + 1;

    /* Allocate memory for LU factorization pivot storage */
    luw->ipiv = malloc(w->ndof * sizeof(luw->ipiv[0]));
    if (luw->ipiv == NULL) {
        suzerain_bspline_lu_free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }
    /* Flag that decomposition has not occurred; Allows catching errors */
    luw->ipiv[0] = -1;

    /* Allocate memory for the matrix formed within lu_form_general */
    luw->A = suzerain_blas_malloc(luw->ld*luw->ndof*sizeof(luw->A[0]));
    if (luw->A == NULL) {
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
    if (luw->ndof < w->ndof) {
        SUZERAIN_ERROR("Incompatible workspaces:"
                       " luw->ndof < w->ndof",
                       SUZERAIN_EINVAL);
    }
    /* Banded LU operator dimensions depend on largest derivative band */
    if (luw->ku < w->max_kl + w->max_ku) {
        SUZERAIN_ERROR("Incompatible workspaces: luw->ku too small",
                       SUZERAIN_EINVAL);
    }

    /* Clear operator storage, including superdiagonals not touched below */
    memset(luw->A, 0, luw->ld*luw->ndof*sizeof(luw->A[0]));

    /* Accumulate scaled derivative operators into luw->A */
    for (int k = 0; k < ncoefficients; ++k) {
        suzerain_blas_dgb_acc(
            luw->ndof, luw->ndof, w->max_kl, w->max_ku,
            coefficients[k], w->D[k] - (w->max_ku - w->ku[k]), w->ld,
            1.0, luw->A + w->max_kl, luw->ld);
    }

    /* Compute LU factorization of the just-formed operator */
    const int info = suzerain_lapack_dgbtrf(luw->ndof,
                                            luw->ndof,
                                            luw->kl,
                                            luw->ku - luw->kl, /* NB */
                                            luw->A,
                                            luw->ld,
                                            luw->ipiv);
    if (info) {
        SUZERAIN_ERROR("suzerain_lapack_dgbtrf reported an error",
                       SUZERAIN_ESANITY);
    }

    /* Factorization overwrote luw->ipiv[0] == -1; This is our flag to
     * check that factorization occurred in suzerain_bspline_lu_solve */

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
    if (luw->ipiv[0] == -1) {
        SUZERAIN_ERROR("One of suzerain_bspline_lu_form_* not called before solve",
                       SUZERAIN_EINVAL);
    }

    const int info = suzerain_lapack_dgbtrs('N',
                                            luw->ndof,
                                            luw->kl,
                                            luw->ku - luw->kl, /* NB */
                                            nrhs,
                                            luw->A,
                                            luw->ld,
                                            luw->ipiv,
                                            b,
                                            ldb);
    if (info) {
        SUZERAIN_ERROR("suzerain_lapack_dgbtrs reported an error",
                       SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}
