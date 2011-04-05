/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * bsplineop.c: B-spline operator routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline.h>
#include <suzerain/bsplineop.h>
#include <suzerain/error.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/pre_gsl.h>

static
int
suzerain_bsplineop_determine_operator_bandwidths(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w);

static
int
suzerain_bsplineop_determine_collocation_bandwidths(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w);

static
int
suzerain_bsplineop_determine_galerkin_bandwidths(
        suzerain_bsplineop_workspace *w);

static
int
suzerain_bsplineop_create_operators(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w);

static
int
suzerain_bsplineop_create_collocation_operators(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w);

static
int
suzerain_bsplineop_create_galerkin_operators(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w);

static
int
compute_banded_collocation_derivative_submatrix(
    const int ioffset,
    const int joffset,
    const int nderiv,
    const int * const kl,
    const int * const ku,
    const int ld,
    const int npoints,
    const double * points,
    gsl_bspline_workspace * bw,
    gsl_bspline_deriv_workspace * dbw,
    gsl_matrix * db,
    double ** const D);

static
int
suzerain_bsplineop_lu_solve_contiguous(
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw);

static
int
suzerain_bsplineop_lu_solve_noncontiguous(
    int nrhs,
    double *b,
    int incb,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw);


suzerain_bsplineop_workspace *
suzerain_bsplineop_alloc(
    gsl_bspline_workspace *bw,
    gsl_bspline_deriv_workspace *dbw,
    int nderiv,
    enum suzerain_bsplineop_method method)
{
    /* Parameter sanity checks */
    if (nderiv < 0) {
        SUZERAIN_ERROR_NULL("nderiv must be at least 0", SUZERAIN_EINVAL);
    }

    switch (method) {
    case SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE:
    case SUZERAIN_BSPLINEOP_GALERKIN_L2:
        break;
    default:
        SUZERAIN_ERROR_NULL("unknown method", SUZERAIN_ESANITY);
    }

    /* Allocate a temporary for use with gsl_bspline_eval_deriv_nonzero */
    gsl_matrix *db = gsl_matrix_alloc(bw->k, nderiv + 1);
    if (db == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate scratch space db",
                            SUZERAIN_ENOMEM);
    }

    /* Calloc workspace and additional storage for known-length members */
    suzerain_bsplineop_workspace * const w = suzerain_blas_calloc(1,
                sizeof(suzerain_bsplineop_workspace)
                + (nderiv+1)*sizeof(w->kl[0])
                + (nderiv+1)*sizeof(w->ku[0])
                + (nderiv+1)*sizeof(w->D[0])
            );
    if (w == NULL) {
        gsl_matrix_free(db);
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }
    w->kl = ((void *) w    ) + sizeof(suzerain_bsplineop_workspace);
    w->ku = ((void *) w->kl) + (nderiv+1)*sizeof(w->kl[0]);
    w->D  = ((void *) w->ku) + (nderiv+1)*sizeof(w->ku[0]);

    /* Save bspline operator parameters in workspace */
    w->method = method;
    w->k      = bw->k;
    w->n      = gsl_bspline_ncoeffs(bw);
    w->nderiv = nderiv;

    /* Storage parameters for BLAS/lapack-compatible general band matrix */
    if (suzerain_bsplineop_determine_operator_bandwidths(db, bw, dbw, w)) {
        gsl_matrix_free(db);
        suzerain_bsplineop_free(w);
        SUZERAIN_ERROR_NULL("failure determining operator bandwidths",
                            SUZERAIN_ESANITY);
    }
    /* Compute derived storage parameters.
     *
     * We choose to have all operators stored using the same ld and capable of
     * being manipulated using w->max_kl and w->max_ku.  This is not as space
     * efficient as it could be in the pathological cases, but it does allow
     * using BLAS dgb_acc to form linear combinations of operators.
     */
    w->max_kl = gsl_stats_int_max(w->kl, 1, w->nderiv+1);
    w->max_ku = gsl_stats_int_max(w->ku, 1, w->nderiv+1);
    w->ld     = w->max_ku + 1 + w->max_kl;

    /* Allocate and clear one aligned block for all operator matrices */
    const size_t an_operators_storage = w->ld * w->n;
    const size_t all_operators_storage
        = an_operators_storage * (w->nderiv+1);
    w->D[0] = suzerain_blas_calloc(all_operators_storage, sizeof(w->D[0][0]));
    if (w->D[0] == NULL) {
        gsl_matrix_free(db);
        suzerain_bsplineop_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }
    /* Compute pointers to the beginning of each derivative operator */
    /* Determine the start of raw storage for (k+1)-th operator */
    for (int k = 1; k <= w->nderiv; ++k) {
        w->D[k] = w->D[k-1] + an_operators_storage;
    }
    /* Step pointer past any unused superdiagonals in the k-th operator */
    for (int k = 0; k <= w->nderiv; ++k) {
        w->D[k] += (w->max_ku - w->ku[k]);
    }

    /* Calculate operator matrices */
    if (suzerain_bsplineop_create_operators(db, bw, dbw, w)) {
        gsl_matrix_free(db);
        suzerain_bsplineop_free(w);
        SUZERAIN_ERROR_NULL("Error creating operator matrices",
                            SUZERAIN_EFAILED);
    }

    /* Free temporary storage */
    gsl_matrix_free(db);

    return w;
}

void
suzerain_bsplineop_free(suzerain_bsplineop_workspace * w)
{
    if (w != NULL) {
        if (w->D != NULL) {
            if (w->D[0] != NULL) {
                // Must free originally allocated memory offset, not the
                // computed offset.  See suzerain_bsplineop_alloc.
                w->D[0] -= (w->max_ku - w->ku[0]);
                suzerain_blas_free(w->D[0]);
                w->D[0] = NULL;
            }
        }

        // D, kl, ku, and w were all allocated together
        suzerain_blas_free(w);
    }
}

int
suzerain_bsplineop_accumulate(
    int nderiv,
    int nrhs,
    double alpha,
    const double *x,
    int incx,
    int ldx,
    double beta,
    double *y,
    int incy,
    int ldy,
    const suzerain_bsplineop_workspace *w)
{
    /* Parameter sanity checks */
    if (nderiv < 0 || w->nderiv < nderiv) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (nrhs > 1 && ldx < w->n) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }
    if (nrhs > 1 && ldy < w->n) {
        SUZERAIN_ERROR("nrhs > 1 && ldy < w->n", SUZERAIN_EINVAL);
    }
    if (x == y) {
        /* BLAS' behavior on aliased pointers is undefined */
        SUZERAIN_ERROR("x == y not allowed", SUZERAIN_EINVAL);
    }

    for (int j = 0; j < nrhs; ++j) {
        const double * const  x_j = x + j*ldx;
        double       * const  y_j = y + j*ldy;
        suzerain_blas_dgbmv('N', w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                            alpha, w->D[nderiv], w->ld, x_j, incx,
                            beta, y_j, incy);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_accumulate_complex(
    int nderiv,
    int nrhs,
    const double alpha[2],
    const double (*x)[2],
    int incx,
    int ldx,
    const double beta[2],
    double (*y)[2],
    int incy,
    int ldy,
    const suzerain_bsplineop_workspace *w)
{
    /* Parameter sanity checks */
    if (nderiv < 0 || w->nderiv < nderiv) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (nrhs > 1 && ldx < w->n) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }
    if (nrhs > 1 && ldy < w->n) {
        SUZERAIN_ERROR("nrhs > 1 && ldy < w->n", SUZERAIN_EINVAL);
    }
    if ((void*) x == (void*) y) {
        /* BLAS' behavior on aliased pointers is undefined */
        SUZERAIN_ERROR("x == y not allowed", SUZERAIN_EINVAL);
    }

    for (int j = 0; j < nrhs; ++j) {
        const double (*const x_j)[2] = x + j*ldx;
        double       (*const y_j)[2] = y + j*ldy;
        suzerain_blasext_dgbmzv('N', w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                                alpha, w->D[nderiv], w->ld, x_j, incx,
                                beta, y_j, incy);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_apply(
    int nderiv,
    int nrhs,
    double alpha,
    double *x,
    int incx,
    int ldx,
    const suzerain_bsplineop_workspace *w)
{
    /* Parameter sanity checks */
    if (nderiv < 0 || w->nderiv < nderiv) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (nrhs > 1 && ldx < w->n) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }

    /* Allocate scratch space; Required because BLAS operations' behavior on
     * aliased pointers is undefined. */
    double * const scratch
        = suzerain_blas_malloc(w->n*sizeof(scratch[0]));
    const int incscratch = 1;
    if (scratch == NULL) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double * const x_j = x + j*ldx;
        /* Compute x_j := w->D[nderiv]*x_j */
        suzerain_blas_dcopy(w->n, x_j, incx, scratch, incscratch);
        suzerain_blas_dgbmv('N', w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                            alpha, w->D[nderiv], w->ld, scratch, incscratch,
                            0.0, x_j, incx);
    }

    suzerain_blas_free(scratch);
    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_apply_complex(
    int nderiv,
    int nrhs,
    double alpha,
    double (*x)[2],
    int incx,
    int ldx,
    const suzerain_bsplineop_workspace *w)
{
    /* Parameter sanity checks */
    if (nderiv < 0 || w->nderiv < nderiv) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (nrhs > 1 && ldx < w->n) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }

    /* Allocate scratch space; Required because BLAS operations' behavior on
     * aliased pointers is undefined. */
    double * const scratch
        = suzerain_blas_malloc(w->n*sizeof(scratch[0]));
    const int incscratch = 1;
    if (scratch == NULL) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double (*const x_j)[2] = x + j*ldx;
        /* Compute x_j := w->D[nderiv]*x_j for real/imaginary parts */
        for (int i = 0; i < 2; ++i) {
            suzerain_blas_dcopy(
                    w->n, &(x_j[0][i]), 2*incx, scratch, incscratch);
            suzerain_blas_dgbmv(
                    'N', w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                    alpha, w->D[nderiv], w->ld, scratch, incscratch,
                    0.0, &(x_j[0][i]), 2*incx);
        }
    }

    suzerain_blas_free(scratch);
    return SUZERAIN_SUCCESS;
}

static
int
suzerain_bsplineop_determine_operator_bandwidths(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w)
{
    /* Compute the operator bandwidths based on the supplied method */
    switch (w->method) {
    case SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE:
        return suzerain_bsplineop_determine_collocation_bandwidths(
                db, bw, dbw, w);
    case SUZERAIN_BSPLINEOP_GALERKIN_L2:
        return suzerain_bsplineop_determine_galerkin_bandwidths(w);
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }
}

static
int
suzerain_bsplineop_determine_collocation_bandwidths(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w)
{
    /* Initially set maximum collocation operator band storage parameters */
    for (int k = 0; k <= w->nderiv; ++k) {
        w->kl[k] = w->k - 1;
        w->ku[k] = w->kl[k];
    }
    /* Determine storage needs, defensively assuming nonuniform parameters */
    w->max_kl = gsl_stats_int_max(w->kl, 1, w->nderiv+1);
    w->max_ku = gsl_stats_int_max(w->ku, 1, w->nderiv+1);
    w->ld     = w->max_kl + 1 + w->max_ku;

    /* Compute collocation points at which we will check bandwidth */
    double * const points = suzerain_blas_malloc(
            2*w->k*sizeof(bw->knots->data[0]));
    if (points == NULL) {
        SUZERAIN_ERROR("Unable to allocate space for collocation points",
                       SUZERAIN_ENOMEM);
    }
    for (int i = 0; i < w->k; ++i) {
        points[i]        = gsl_bspline_greville_abscissa(i, bw);
        points[i + w->k] = gsl_bspline_greville_abscissa(w->n - w->k + i, bw);
    }

    /* Allocate space for collocation operator submatrices */
    double ** const scratch
        = suzerain_blas_malloc(2*(w->nderiv+1) * sizeof(scratch[0]));
    if (scratch == NULL) {
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Unable to allocate scratch pointers",
                SUZERAIN_ENOMEM);
    }
    const size_t an_operators_storage = w->ld * w->k;
    const size_t total_storage = 2*(w->nderiv+1)*an_operators_storage;
    scratch[0] = suzerain_blas_calloc(total_storage, sizeof(w->D[0][0]));
    if (scratch[0] == NULL) {
        suzerain_blas_free(scratch);
        suzerain_blas_free(points);
        SUZERAIN_ERROR_NULL("failed to allocate scratch",
                            SUZERAIN_ENOMEM);
    }
    for (int k = 1; k < 2*(w->nderiv+1); ++k) {
        scratch[k] = scratch[k-1] + an_operators_storage;
    }

    /* Create convenience views of our working space */
    const double * const ul_points = &(points[0]);
    const double * const lr_points = &(points[w->k]);
    double ** const ul_D = &(scratch[0]);
    double ** const lr_D = &(scratch[w->nderiv+1]);

    /* Compute the upper-left- and lower-right- most submatrices */
    if (compute_banded_collocation_derivative_submatrix(
             0, 0,
             w->nderiv, w->kl, w->ku, w->ld,
             w->k, ul_points, bw, dbw, db, ul_D)) {
        suzerain_blas_free(scratch[0]);
        suzerain_blas_free(scratch);
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Error computing operator UL submatrices",
                       SUZERAIN_ESANITY);
    }
    const int lr_offset = w->n - w->k;
    if (compute_banded_collocation_derivative_submatrix(
             lr_offset, lr_offset,
             w->nderiv, w->kl, w->ku, w->ld,
             w->k, lr_points, bw, dbw, db, lr_D)) {
        suzerain_blas_free(scratch[0]);
        suzerain_blas_free(scratch);
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Error computing operator LR submatrices",
                       SUZERAIN_ESANITY);
    }

    /* Reduce kl/ku for each zero off-diagonal in ul_D and lr_D */
    for (int k = 0; k <= w->nderiv; ++k) {
        const int fixed_ku_k = w->ku[k];
        const int fixed_ld   = w->ld;

        for (int i=0; i < fixed_ku_k; ++i) {
            const double ul_asum
                = suzerain_blas_dasum(w->k, ul_D[k]+i, fixed_ld);
            const double lr_asum
                = suzerain_blas_dasum(w->k, lr_D[k]+i, fixed_ld);
#pragma warning(push,disable:1572)
            if (ul_asum + lr_asum == 0.0) {
#pragma warning(pop)
                --w->ku[k];
            } else {
                break; /* Skip all after nonzero superdiagonal */
            }
        }
        for (int i=fixed_ld-1; i > fixed_ku_k; --i) {
            const double ul_asum
                = suzerain_blas_dasum(w->k, ul_D[k]+i, fixed_ld);
            const double lr_asum
                = suzerain_blas_dasum(w->k, lr_D[k]+i, fixed_ld);
#pragma warning(push,disable:1572)
            if (ul_asum + lr_asum == 0.0) {
#pragma warning(pop)
                --w->kl[k];
            } else {
                break; /* Skip all after nonzero subdiagonal */
            }
        }
    }

    suzerain_blas_free(scratch[0]);
    suzerain_blas_free(scratch);
    suzerain_blas_free(points);

    return SUZERAIN_SUCCESS;
}

static
int
suzerain_bsplineop_determine_galerkin_bandwidths(
        suzerain_bsplineop_workspace *w)
{
    for (int k = 0; k <= w->nderiv; ++k) {
        w->kl[k] = w->k - 1;
        w->ku[k] = w->kl[k];
    }

    return SUZERAIN_SUCCESS;
}

static
int
suzerain_bsplineop_create_operators(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w)
{
    /* Compute the operator bandwidths based on the supplied method */
    switch (w->method) {
    case SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE:
        return suzerain_bsplineop_create_collocation_operators(db, bw, dbw, w);
    case SUZERAIN_BSPLINEOP_GALERKIN_L2:
        return suzerain_bsplineop_create_galerkin_operators(db, bw, dbw, w);
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }
}

static
int
suzerain_bsplineop_create_collocation_operators(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w)
{
    /* Evaluate basis at collocation points: d^k/dx^k B_j(\xi_i) */
    double * const points = suzerain_blas_malloc(w->n*sizeof(points[0]));
    if (points == NULL) {
        SUZERAIN_ERROR("Unable to allocate space for Greville abscissae",
                       SUZERAIN_ENOMEM);
    }
    for (size_t i = 0; i < bw->n; ++i) {
        points[i] = gsl_bspline_greville_abscissa(i, bw);
    }

    /* Compute the full derivative operator matrices */
    if (compute_banded_collocation_derivative_submatrix(
             0, 0, w->nderiv, w->kl, w->ku, w->ld,
             w->n, points, bw, dbw, db, w->D)) {
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Error computing operator matrices", SUZERAIN_ESANITY);
    }

    suzerain_blas_free(points);

    return SUZERAIN_SUCCESS;
}

static
int
compute_banded_collocation_derivative_submatrix(
    const int ioffset,
    const int joffset,
    const int nderiv,
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
    /* PRECONDITION: D[0]...D[nderiv] must have been filled with zeros */

    /* Fill all nonzero entries in the operator storage */
    for (int i = 0; i < npoints; ++i) {
        size_t dbjstart, dbjend;
        gsl_bspline_deriv_eval_nonzero(points[i], nderiv, db,
                                       &dbjstart, &dbjend, bw, dbw);

        /* Coerce dbjstart/dbjend to stay within the submatrix of interest */
        const size_t jstart = GSL_MAX(dbjstart, (size_t) joffset);
        const int    jend   = GSL_MIN(dbjend, (size_t)(npoints-1+joffset));

        for (int k = 0; k <= nderiv; ++k) {
            for (int j = jstart; j <= jend; ++j) {
                const double value = gsl_matrix_get(db, j - dbjstart, k);
                const int in_band  = suzerain_gbmatrix_in_band(
                        ld, kl[k], ku[k], i+ioffset, j /* no joffset */);
                const int offset = suzerain_gbmatrix_offset(
                        ld, kl[k], ku[k], i /* no ioffset */, j-joffset);

                if (in_band) {
                    D[k][offset] = value;
#pragma warning(push,disable:1572)
                } else if (value == 0.0) {
#pragma warning(pop)
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

static inline int int_div_ceil(int dividend, int divisor) {
    // Modified from http://www.justlinux.com/forum/showthread.php?t=126876
    return ((dividend + divisor) - 1) / divisor;
}

static
int
suzerain_bsplineop_create_galerkin_operators(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w)
{
    /* PRECONDITION: D[0]...D[nderiv] must have been filled with zeros */

    /* Determine the integration order and retrieve Gauss-Legendre rule */
    /* Maximum quadratic piecewise polynomial order is 2*(w->k-1)   */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc(int_div_ceil(2*w->k - 1, 2));
    if (tbl == NULL) {
        SUZERAIN_ERROR_NULL("failed to obtain Gauss-Legendre rule from GSL",
                            SUZERAIN_ESANITY);
    }

    /* For each DOF in the basis... */
    for (int i = 0; i < w->n; ++i) {

        /* ...loop over the nontrivial knot intervals in support... */
        for (int j = i; j < i + w->k; ++j) {
            const double a = gsl_vector_get(bw->knots, j);
            const double b = gsl_vector_get(bw->knots, j + 1);
            if (SUZERAIN_UNLIKELY(a == b)) continue;

            /* ...evaluate all derivatives at each Gauss point... */
            /* scale by the Gauss weight, and accumulate into   */
            /* the banded derivative matrix. */

            /* ...at each Gauss point... */
            for (size_t k = 0; k < tbl->n; ++k) {
                double xk, wk;
                gsl_integration_glfixed_point(a, b, k, &xk, &wk, tbl);

                /* ...evaluate all derivatives at the Gauss point... */
                size_t istart, iend;
                gsl_bspline_deriv_eval_nonzero(xk, w->nderiv,
                                               db, &istart, &iend,
                                               bw, dbw);

                /* ...scale by the corresponding Gauss weight...         */
                /* ...times the i-th DOF evaluated at the Gauss point... */
                gsl_matrix_scale(db,
                                 wk * gsl_matrix_get(db, i - istart, 0));

                /* ...and accumulate into the banded derivative matrices. */
                for (size_t l = istart; l <= iend; ++l) {
                    for (size_t m = 0; m <= (size_t) w->nderiv; ++m) {
                        const int offset = suzerain_gbmatrix_offset(
                                w->ld, w->kl[m], w->ku[m], i, l);
                        w->D[m][offset] += gsl_matrix_get(db, l - istart, m);
                    }
                }
            }
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    /* Mass matrices are analytically guaranteed to be symmetric, positive
     * definite and so the sub- and super- diagonals are identical.  Any
     * floating point noise differences spoil the discrete symmetry and will
     * cause problems if people use BLAS DGMBV vs DSBMV-with-upper-diagonals vs
     * DSBMV-with-lower-diagonals to compute matrix-vector products or
     * factorizations.  Explicitly enforce the mass matrix's discrete symmetry
     * by copying the subdiagonals to the superdiagonals.
     */
    if (w->nderiv >= 0) {
        assert(w->kl[0] == w->ku[0]);
        for (int diag = 1; diag <= w->kl[0]; ++diag) {
            double * const sup = w->D[0] + w->ku[0] - diag;
            double * const sub = w->D[0] + w->ku[0] + diag;
            for (int i = 0; i < w->n - diag; ++i) {
                sup[(i+diag)*w->ld] = sub[i*w->ld];
            }
        }
    }

    /* Odd derivative operators have zeros along the interior of their main
     * diagonals.  Any floating point noise here is infinitely huge compared to
     * the magnitude of zero and fundamentally changes the operator.
     * Explicitly set tiny diagonal values to zero in odd derivatives for a
     * reasonable (order- and derivative-dependent) definition of tiny.
     */
    for (int n = 1; n <= w->nderiv; n += 2) {
        const double tiny = (2*w->k - 1) * n * GSL_DBL_EPSILON;
        for (int i = 0; i < w->n; ++i) {
            const int offset = suzerain_gbmatrix_offset(
                    w->ld, w->kl[n], w->ku[n], i, i);
            if (fabs(w->D[n][offset]) < tiny) {
                w->D[n][offset] = 0.0;
            }
        }
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_interpolation_rhs(
    const suzerain_function * function,
    double * rhs,
    gsl_bspline_workspace *bw,
    const suzerain_bsplineop_workspace *w)
{
    switch (w->method) {
    case SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE:
        {
            /* Evaluate the function at the collocation points */
            for (size_t i = 0; i < bw->n; ++i) {
                const double x = gsl_bspline_greville_abscissa(i, bw);
                rhs[i] = SUZERAIN_FN_EVAL(function, x);
            }
        }

        break;
    case SUZERAIN_BSPLINEOP_GALERKIN_L2:
        SUZERAIN_ERROR("unimplemented method", SUZERAIN_ESANITY); // FIXME
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_interpolation_rhs_complex(
    const suzerain_zfunction * zfunction,
    double (*rhs)[2],
    gsl_bspline_workspace *bw,
    const suzerain_bsplineop_workspace *w)
{
    switch (w->method) {
    case SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE:
        {
            /* Evaluate the function at the collocation points */
            for (size_t i = 0; i < bw->n; ++i) {
                const double x = gsl_bspline_greville_abscissa(i, bw);
                SUZERAIN_ZFN_EVAL(zfunction, x, rhs[i]);
            }
        }

        break;
    case SUZERAIN_BSPLINEOP_GALERKIN_L2:
        SUZERAIN_ERROR("unimplemented method", SUZERAIN_ESANITY); // FIXME
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

/********************************/
/* Real-valued LU functionality */
/********************************/

suzerain_bsplineop_lu_workspace *
suzerain_bsplineop_lu_alloc(
    const suzerain_bsplineop_workspace *w)
{
    /* Allocate space for the luw workspace */
    suzerain_bsplineop_lu_workspace * const luw
        = suzerain_blas_malloc(sizeof(suzerain_bsplineop_lu_workspace));
    if (luw == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }
    /* Make workspace pointers NULL */
    luw->A    = NULL;
    luw->ipiv = NULL;

    /* Determine general banded matrix shape parameters */
    luw->n   = w->n;
    luw->kl  = w->max_kl;
    luw->ku  = w->max_kl + w->max_ku; /* Increase per GBTRF, GBTRS */
    luw->ld  = luw->kl + luw->ku + 1;

    /* Allocate memory for LU factorization pivot storage */
    luw->ipiv = suzerain_blas_malloc(w->n * sizeof(luw->ipiv[0]));
    if (luw->ipiv == NULL) {
        suzerain_bsplineop_lu_free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }
    /* Flag that decomposition has not occurred; Allows catching errors */
    luw->ipiv[0] = -1;

    /* Allocate memory for the matrix formed within lu_form */
    luw->A = suzerain_blas_malloc(luw->ld*luw->n*sizeof(luw->A[0]));
    if (luw->A == NULL) {
        suzerain_bsplineop_lu_free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }

    return luw;
}

void
suzerain_bsplineop_lu_free(suzerain_bsplineop_lu_workspace * luw)
{
    if (luw != NULL) {
        suzerain_blas_free(luw->A);
        luw->A = NULL;

        suzerain_blas_free(luw->ipiv);
        luw->ipiv = NULL;

        suzerain_blas_free(luw);
    }
}

int
suzerain_bsplineop_lu_form(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace *luw)
{
    /* Parameter sanity checks */
    if (ncoefficients < 0) {
        SUZERAIN_ERROR("Number of coefficients cannot be negative",
                       SUZERAIN_EINVAL);
    }
    if (ncoefficients > w->nderiv + 1) {
        SUZERAIN_ERROR("More coefficients provided than derivatives available",
                       SUZERAIN_EINVAL);
    }
    if (luw->n < w->n) {
        SUZERAIN_ERROR("Incompatible workspaces:  luw->n < w->n",
                       SUZERAIN_EINVAL);
    }
    /* Banded LU operator dimensions depend on largest derivative band */
    if (luw->ku < w->max_kl + w->max_ku) {
        SUZERAIN_ERROR("Incompatible workspaces: luw->ku too small",
                       SUZERAIN_EINVAL);
    }

    /* Clear operator storage, including superdiagonals not touched below */
    memset(luw->A, 0, luw->ld*luw->n*sizeof(luw->A[0]));

    /* Accumulate scaled derivative operators into luw->A */
    for (int k = 0; k < ncoefficients; ++k) {
        suzerain_blas_dgb_acc(
            luw->n, luw->n, w->max_kl, w->max_ku,
            coefficients[k], w->D[k] - (w->max_ku - w->ku[k]), w->ld,
            1.0, luw->A + w->max_kl, luw->ld);
    }

    /* Compute LU factorization of the just-formed operator */
    const int info = suzerain_lapack_dgbtrf(luw->n,
                                            luw->n,
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
     * check that factorization occurred in suzerain_bsplineop_lu_solve */

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_lu_form_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace *luw)
{
    const double d_one = 1.0;
    const int i_one = 1;
    return suzerain_bsplineop_lu_form(i_one, &d_one, w, luw);
}

static
int
suzerain_bsplineop_lu_solve_contiguous(
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw)
{
    if (luw->ipiv[0] == -1) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_lu_form_* not called before solve",
                SUZERAIN_EINVAL);
    }

    const int info = suzerain_lapack_dgbtrs('N',
                                            luw->n,
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

static
int
suzerain_bsplineop_lu_solve_noncontiguous(
    int nrhs,
    double *b,
    int incb,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw)
{
    if (luw->ipiv[0] == -1) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_lu_form_* not called before solve",
                SUZERAIN_EINVAL);
    }

    /* Allocate scratch space because GBTRS requires contiguous vectors */
    double * const scratch
        = suzerain_blas_malloc(luw->n*sizeof(scratch[0]));
    if (scratch == NULL) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double * const b_j = b + j*ldb;

        suzerain_blas_dcopy(luw->n, b_j, incb, scratch, 1);
        const int info = suzerain_lapack_dgbtrs('N',
                                                luw->n,
                                                luw->kl,
                                                luw->ku - luw->kl, /* NB */
                                                1, /* One RHS at a time */
                                                luw->A,
                                                luw->ld,
                                                luw->ipiv,
                                                scratch,
                                                luw->n);
        if (info) {
            suzerain_blas_free(scratch);
            SUZERAIN_ERROR("suzerain_lapack_dgbtrs reported an error",
                        SUZERAIN_ESANITY);
        }
        suzerain_blas_dcopy(luw->n, scratch, 1, b_j, incb);
    }

    suzerain_blas_free(scratch);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_lu_solve(
    int nrhs,
    double *b,
    int incb,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw)
{
    if (incb == 1) {
        return suzerain_bsplineop_lu_solve_contiguous(
                nrhs, b, ldb, luw);
    } else {
        return suzerain_bsplineop_lu_solve_noncontiguous(
                nrhs, b, incb, ldb, luw);
    }
}

/***********************************/
/* Complex-valued LU functionality */
/***********************************/

suzerain_bsplineop_luz_workspace *
suzerain_bsplineop_luz_alloc(
    const suzerain_bsplineop_workspace *w)
{
    /* Allocate space for the luzw workspace */
    suzerain_bsplineop_luz_workspace * const luzw
        = suzerain_blas_malloc(sizeof(suzerain_bsplineop_luz_workspace));
    if (luzw == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }
    /* Make workspace pointers NULL */
    luzw->A    = NULL;
    luzw->ipiv = NULL;

    /* Determine general banded matrix shape parameters */
    luzw->n  = w->n;
    luzw->kl    = w->max_kl;
    luzw->ku    = w->max_kl + w->max_ku; /* Increase per GBTRF, GBTRS */
    luzw->ld    = luzw->kl + luzw->ku + 1;

    /* Allocate memory for LU factorization pivot storage */
    luzw->ipiv = suzerain_blas_malloc(w->n * sizeof(luzw->ipiv[0]));
    if (luzw->ipiv == NULL) {
        suzerain_bsplineop_luz_free(luzw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }
    /* Flag that decomposition has not occurred; Allows catching errors */
    luzw->ipiv[0] = -1;

    /* Allocate memory for the matrix formed within luz_form */
    luzw->A = suzerain_blas_malloc(luzw->ld*luzw->n*sizeof(luzw->A[0]));
    if (luzw->A == NULL) {
        suzerain_bsplineop_luz_free(luzw);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }

    return luzw;
}

void
suzerain_bsplineop_luz_free(suzerain_bsplineop_luz_workspace * luzw)
{
    if (luzw != NULL) {
        suzerain_blas_free(luzw->A);
        luzw->A = NULL;

        suzerain_blas_free(luzw->ipiv);
        luzw->ipiv = NULL;

        suzerain_blas_free(luzw);
    }
}

int
suzerain_bsplineop_luz_form(
    int ncoefficients,
    const double (*coefficients)[2],
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace *luzw)
{
    /* Parameter sanity checks */
    if (ncoefficients < 0) {
        SUZERAIN_ERROR("Number of coefficients cannot be negative",
                       SUZERAIN_EINVAL);
    }
    if (ncoefficients > w->nderiv + 1) {
        SUZERAIN_ERROR("More coefficients provided than derivatives available",
                       SUZERAIN_EINVAL);
    }
    if (luzw->n < w->n) {
        SUZERAIN_ERROR("Incompatible workspaces: luzw->n < w->n",
                       SUZERAIN_EINVAL);
    }
    /* Banded LU operator dimensions depend on largest derivative band */
    if (luzw->ku < w->max_kl + w->max_ku) {
        SUZERAIN_ERROR("Incompatible workspaces: luzw->ku too small",
                       SUZERAIN_EINVAL);
    }

    /* Clear operator storage, including superdiagonals not touched below */
    memset(luzw->A, 0, luzw->ld*luzw->n*sizeof(luzw->A[0]));

    /* Accumulate scaled derivative operators into luzw->A */
    const double z_one[2] = { 1.0, 0.0 };
    for (int k = 0; k < ncoefficients; ++k) {
        suzerain_blasext_zgb_dacc(
            luzw->n, luzw->n, w->max_kl, w->max_ku,
            coefficients[k], w->D[k] - (w->max_ku - w->ku[k]), w->ld,
            z_one, luzw->A + w->max_kl, luzw->ld);
    }

    /* Compute LU factorization of the just-formed operator */
    const int info = suzerain_lapack_zgbtrf(luzw->n,
                                            luzw->n,
                                            luzw->kl,
                                            luzw->ku - luzw->kl, /* NB */
                                            luzw->A,
                                            luzw->ld,
                                            luzw->ipiv);
    if (info) {
        SUZERAIN_ERROR("suzerain_lapack_zgbtrf reported an error",
                       SUZERAIN_ESANITY);
    }

    /* Factorization overwrote luzw->ipiv[0] == -1; This is our flag to
     * check that factorization occurred in suzerain_bsplineop_luz_solve */

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_luz_form_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace *luzw)
{
    const double z_one[2] = { 1.0, 0.0 };
    const int i_one = 1;
    return suzerain_bsplineop_luz_form(i_one, &z_one, w, luzw);
}

static
int
suzerain_bsplineop_luz_solve_contiguous(
    int nrhs,
    double (*b)[2],
    int ldb,
    const suzerain_bsplineop_luz_workspace *luzw)
{
    if (luzw->ipiv[0] == -1) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_luz_form_* not called before solve",
                SUZERAIN_EINVAL);
    }

    const int info = suzerain_lapack_zgbtrs('N',
                                            luzw->n,
                                            luzw->kl,
                                            luzw->ku - luzw->kl, /* NB */
                                            nrhs,
                                            (const double (*)[2]) luzw->A,
                                            luzw->ld,
                                            luzw->ipiv,
                                            b,
                                            ldb);
    if (info) {
        SUZERAIN_ERROR("suzerain_lapack_zgbtrs reported an error",
                       SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

static
int
suzerain_bsplineop_luz_solve_noncontiguous(
    int nrhs,
    double (*b)[2],
    int incb,
    int ldb,
    const suzerain_bsplineop_luz_workspace *luzw)
{
    if (luzw->ipiv[0] == -1) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_luz_form_* not called before solve",
                SUZERAIN_EINVAL);
    }

    /* Allocate scratch space because ZGBTRS requires contiguous vectors */
    double (* const scratch)[2]
        = suzerain_blas_malloc(luzw->n*sizeof(scratch[0]));
    if (scratch == NULL) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double (* const b_j)[2] = b + j*ldb;

        suzerain_blas_zcopy(
                luzw->n, (const double (*)[2]) b_j, incb, scratch, 1);
        const int info = suzerain_lapack_zgbtrs('N',
                                                luzw->n,
                                                luzw->kl,
                                                luzw->ku - luzw->kl, /* NB */
                                                1, /* One RHS at a time */
                                                (const double (*)[2]) luzw->A,
                                                luzw->ld,
                                                luzw->ipiv,
                                                scratch,
                                                luzw->n);
        if (info) {
            suzerain_blas_free(scratch);
            SUZERAIN_ERROR("suzerain_lapack_zgbtrs reported an error",
                        SUZERAIN_ESANITY);
        }
        suzerain_blas_zcopy(
                luzw->n, (const double (*)[2]) scratch, 1, b_j, incb);
    }

    suzerain_blas_free(scratch);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_luz_solve(
    int nrhs,
    double (*b)[2],
    int incb,
    int ldb,
    const suzerain_bsplineop_luz_workspace *luzw)
{
    if (incb == 1) {
        return suzerain_bsplineop_luz_solve_contiguous(
                nrhs, b, ldb, luzw);
    } else {
        return suzerain_bsplineop_luz_solve_noncontiguous(
                nrhs, b, incb, ldb, luzw);
    }
}
