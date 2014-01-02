/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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
 */

/** @file
 * @copydoc bsplineop.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/bsplineop.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>

#include <suzerain/common.h>
#include <suzerain/blas_et_al.h>
#include <suzerain/bspline.h>
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
suzerain_bsplineop_create_collocation_operator_transposes(
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
compute_banded_collocation_derivative_submatrix_transpose(
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

suzerain_bsplineop_workspace *
suzerain_bsplineop_alloc(
    gsl_bspline_workspace *bw,
    gsl_bspline_deriv_workspace *dbw,
    int nderiv,
    enum suzerain_bsplineop_method method)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(nderiv < 0)) {
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
    if (SUZERAIN_UNLIKELY(db == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate scratch space db",
                            SUZERAIN_ENOMEM);
    }

    /* Calloc workspace and additional storage for known-length members */
    suzerain_bsplineop_workspace * const w = suzerain_blas_calloc(1,
                sizeof(suzerain_bsplineop_workspace)
                + (nderiv+1)*sizeof(w->kl[0])
                + (nderiv+1)*sizeof(w->ku[0])
                + (nderiv+1)*sizeof(w->D_T[0])
            );
    if (SUZERAIN_UNLIKELY(w == NULL)) {
        gsl_matrix_free(db);
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }
    w->kl  = (void *)(((char *) w) + sizeof(suzerain_bsplineop_workspace));
    w->ku  = w->kl + (nderiv+1);
    w->D_T = (void *)(w->ku + (nderiv+1));

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
    w->D_T[0] = suzerain_blas_calloc(all_operators_storage,
                                     sizeof(w->D_T[0][0]));
    if (SUZERAIN_UNLIKELY(w->D_T[0] == NULL)) {
        gsl_matrix_free(db);
        suzerain_bsplineop_free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }
    /* Compute pointers to the beginning of each derivative operator */
    /* Determine the start of raw storage for (k+1)-th operator */
    for (int k = 1; k <= w->nderiv; ++k) {
        w->D_T[k] = w->D_T[k-1] + an_operators_storage;
    }
    /* Step pointer past any unused superdiagonals in the k-th operator */
    for (int k = 0; k <= w->nderiv; ++k) {
        w->D_T[k] += (w->max_ku - w->ku[k]);
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
        if (w->D_T != NULL) {
            if (w->D_T[0] != NULL) {
                // Must free originally allocated memory offset, not the
                // computed offset.  See suzerain_bsplineop_alloc.
                w->D_T[0] -= (w->max_ku - w->ku[0]);
                suzerain_blas_free(w->D_T[0]);
                w->D_T[0] = NULL;
            }
        }

        // D_T, kl, ku, and w were all allocated together
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
    if (SUZERAIN_UNLIKELY(nderiv < 0 || w->nderiv < nderiv)) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(nrhs > 1 && ldx < w->n)) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(nrhs > 1 && ldy < w->n)) {
        SUZERAIN_ERROR("nrhs > 1 && ldy < w->n", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(x == y)) {
        /* BLAS' behavior on aliased pointers is undefined */
        SUZERAIN_ERROR("x == y not allowed", SUZERAIN_EINVAL);
    }

    for (int j = 0; j < nrhs; ++j) {
        suzerain_blas_dgbmv('T',   w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                            alpha, w->D_T[nderiv], w->ld, x + j*ldx, incx,
                            beta,                         y + j*ldy, incy);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_accumulate_complex(
    int nderiv,
    int nrhs,
    const complex_double alpha,
    const complex_double *x,
    int incx,
    int ldx,
    const complex_double beta,
    complex_double *y,
    int incy,
    int ldy,
    const suzerain_bsplineop_workspace *w)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(nderiv < 0 || w->nderiv < nderiv)) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(nrhs > 1 && ldx < w->n)) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(nrhs > 1 && ldy < w->n)) {
        SUZERAIN_ERROR("nrhs > 1 && ldy < w->n", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(x == y)) {
        /* BLAS' behavior on aliased pointers is undefined */
        SUZERAIN_ERROR("x == y not allowed", SUZERAIN_EINVAL);
    }

    for (int j = 0; j < nrhs; ++j) {
        suzerain_blas_zgbmv_d_z(
                'T', w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                alpha, w->D_T[nderiv], w->ld, x + j*ldx, incx,
                beta,                         y + j*ldy, incy);
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
    if (SUZERAIN_UNLIKELY(nderiv < 0 || w->nderiv < nderiv)) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(nrhs > 1 && ldx < w->n)) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }

    /* Allocate scratch space; Required because BLAS operations' behavior on
     * aliased pointers is undefined. */
    double * const scratch
        = suzerain_blas_malloc(w->n*sizeof(scratch[0]));
    const int incscratch = 1;
    if (SUZERAIN_UNLIKELY(scratch == NULL)) {
        SUZERAIN_ERROR("failed to allocate scratch space", SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double * const x_j = x + j*ldx;
        /* Compute x_j := w->D[nderiv]*x_j */
        suzerain_blas_dcopy(w->n, x_j, incx, scratch, incscratch);
        suzerain_blas_dgbmv('T',   w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                            alpha, w->D_T[nderiv], w->ld, scratch, incscratch,
                            0.0,                          x_j,     incx);
    }

    suzerain_blas_free(scratch);
    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_apply_complex(
    int nderiv,
    int nrhs,
    double alpha,
    complex_double *x,
    int incx,
    int ldx,
    const suzerain_bsplineop_workspace *w)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(nderiv < 0 || w->nderiv < nderiv)) {
        SUZERAIN_ERROR("nderiv out of range", SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(nrhs > 1 && ldx < w->n)) {
        SUZERAIN_ERROR("nrhs > 1 && ldx < w->n", SUZERAIN_EINVAL);
    }

    /* Allocate scratch space; Required because BLAS operations' behavior on
     * aliased pointers is undefined. */
    double * const scratch = suzerain_blas_malloc(w->n*sizeof(scratch[0]));
    const int incscratch = 1;
    if (SUZERAIN_UNLIKELY(scratch == NULL)) {
        SUZERAIN_ERROR("failed to allocate scratch space", SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double * const xreal_j = (double *)(x + j*ldx);
        /* Compute x_j := w->D[nderiv]*x_j for real/imaginary parts */
        for (int c = 0; c < 2; ++c) {
            suzerain_blas_dcopy(
                    w->n, xreal_j + c, 2*incx, scratch, incscratch);
            suzerain_blas_dgbmv(
                    'T',   w->n, w->n, w->kl[nderiv], w->ku[nderiv],
                    alpha, w->D_T[nderiv], w->ld, scratch,     incscratch,
                    0.0,                          xreal_j + c, 2*incx);
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
    if (SUZERAIN_UNLIKELY(points == NULL)) {
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
    if (SUZERAIN_UNLIKELY(scratch == NULL)) {
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Unable to allocate scratch pointers",
                SUZERAIN_ENOMEM);
    }
    const size_t an_operators_storage = w->ld * w->k;
    const size_t total_storage = 2*(w->nderiv+1)*an_operators_storage;
    scratch[0] = suzerain_blas_calloc(total_storage, sizeof(w->D_T[0][0]));
    if (SUZERAIN_UNLIKELY(scratch[0] == NULL)) {
        suzerain_blas_free(scratch);
        suzerain_blas_free(points);
        SUZERAIN_ERROR_NULL("failed to allocate scratch", SUZERAIN_ENOMEM);
    }
    for (int k = 1; k < 2*(w->nderiv+1); ++k) {
        scratch[k] = scratch[k-1] + an_operators_storage;
    }

    /* Create convenience views of our working space */
    const double * const ul_points = &(points[0]);
    const double * const lr_points = &(points[w->k]);
    double ** const ul_D_T = &(scratch[0]);
    double ** const lr_D_T = &(scratch[w->nderiv+1]);

    /* Compute the upper-left- and lower-right- most submatrices */
    if (compute_banded_collocation_derivative_submatrix_transpose(
             0, 0,
             w->nderiv, w->kl, w->ku, w->ld,
             w->k, ul_points, bw, dbw, db, ul_D_T)) {
        suzerain_blas_free(scratch[0]);
        suzerain_blas_free(scratch);
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Error computing operator UL submatrices",
                       SUZERAIN_ESANITY);
    }
    const int lr_offset = w->n - w->k;
    if (compute_banded_collocation_derivative_submatrix_transpose(
             lr_offset, lr_offset,
             w->nderiv, w->kl, w->ku, w->ld,
             w->k, lr_points, bw, dbw, db, lr_D_T)) {
        suzerain_blas_free(scratch[0]);
        suzerain_blas_free(scratch);
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Error computing operator LR submatrices",
                       SUZERAIN_ESANITY);
    }

    /* Reduce kl/ku for each zero off-diagonal in ul_D_T and lr_D_T */
    for (int k = 0; k <= w->nderiv; ++k) {
        const int fixed_ku_k = w->ku[k];
        const int fixed_ld   = w->ld;

        for (int i=0; i < fixed_ku_k; ++i) {
            const double ul_asum
                = suzerain_blas_dasum(w->k, ul_D_T[k]+i, fixed_ld);
            const double lr_asum
                = suzerain_blas_dasum(w->k, lr_D_T[k]+i, fixed_ld);
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
                = suzerain_blas_dasum(w->k, ul_D_T[k]+i, fixed_ld);
            const double lr_asum
                = suzerain_blas_dasum(w->k, lr_D_T[k]+i, fixed_ld);
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
        return suzerain_bsplineop_create_collocation_operator_transposes(
                db, bw, dbw, w);
    case SUZERAIN_BSPLINEOP_GALERKIN_L2:
        return suzerain_bsplineop_create_galerkin_operators(db, bw, dbw, w);
    default:
        SUZERAIN_ERROR("unknown method", SUZERAIN_ESANITY);
    }
}

static
int
suzerain_bsplineop_create_collocation_operator_transposes(
        gsl_matrix *db,
        gsl_bspline_workspace *bw,
        gsl_bspline_deriv_workspace *dbw,
        suzerain_bsplineop_workspace *w)
{
    /* Evaluate basis at collocation points: d^k/dx^k B_i(\xi_j) */
    double * const points = suzerain_blas_malloc(w->n*sizeof(points[0]));
    if (SUZERAIN_UNLIKELY(points == NULL)) {
        SUZERAIN_ERROR("Unable to allocate space for Greville abscissae",
                       SUZERAIN_ENOMEM);
    }
    for (size_t i = 0; i < bw->n; ++i) {
        points[i] = gsl_bspline_greville_abscissa(i, bw);
    }

    /* Compute the full derivative operator matrix's transpose */
    if (compute_banded_collocation_derivative_submatrix_transpose(
             0, 0, w->nderiv, w->kl, w->ku, w->ld,
             w->n, points, bw, dbw, db, w->D_T)) {
        suzerain_blas_free(points);
        SUZERAIN_ERROR("Error computing operator matrices", SUZERAIN_ESANITY);
    }

    suzerain_blas_free(points);

    return SUZERAIN_SUCCESS;
}

static
int
compute_banded_collocation_derivative_submatrix_transpose(
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
    double ** const D_T)
{
    /* PRECONDITION: D_T[0]...D_T[nderiv] must have been filled with zeros */

    /*
     * Fill all nonzero entries in the transpose of the operator's storage.
     * B-spline basis evaluation routines make d^k/dx^k B_j(\xi_i) easier
     * to compute so we loop over B_j(\xi_i) while storing B_i(\xi_j).
     */

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
                const int in_band  = suzerain_gbmatrix_in_band( /*TRANS*/
                            kl[k], ku[k], j /* no joffset */, i+ioffset);
                const int offset = suzerain_gbmatrix_offset(    /*TRANS*/
                        ld, kl[k], ku[k], j-joffset, i /* no ioffset */);

                if (in_band) {
                    D_T[k][offset] = value;
#pragma warning(push,disable:1572)
                } else if (value == 0.0) {
#pragma warning(pop)
                    /* OK: value outside band is identically zero */
                } else {
                    /* NOT COOL: nonzero value outside bandwidth */
                    const int order = bw->k;
                    char buffer[384];
                    snprintf(buffer, sizeof(buffer),
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
    /* Modified from http://www.justlinux.com/forum/showthread.php?t=126876 */
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

    /*
     * B-spline basis evaluation routines make d^k/dx^k B_j(\xi_i) easier
     * to compute so we loop over B_j(\xi_i) while storing B_i(\xi_j).
     */

    /* Determine the integration order and retrieve Gauss-Legendre rule */
    /* Maximum quadratic piecewise polynomial order is 2*(w->k-1)   */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc(int_div_ceil(2*w->k - 1, 2));
    if (SUZERAIN_UNLIKELY(tbl == NULL)) {
        SUZERAIN_ERROR_NULL("failed to obtain Gauss-Legendre rule from GSL",
                            SUZERAIN_ESANITY);
    }

    /* For each DOF in the basis... */
    for (int i = 0; i < w->n; ++i) {

        /* ...loop over the nontrivial knot intervals in support... */
        for (int j = i; j < i + w->k; ++j) {
            const double a = gsl_vector_get(bw->knots, j);
            const double b = gsl_vector_get(bw->knots, j + 1);
#pragma warning(push,disable:1572)
            if (SUZERAIN_UNLIKELY(a == b)) continue;
#pragma warning(pop)

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

                /* ...and accumulate into derivative matrix's transpose. */
                for (size_t l = istart; l <= iend; ++l) {
                    for (size_t m = 0; m <= (size_t) w->nderiv; ++m) {
                        const int offset = suzerain_gbmatrix_offset( /*TRANS!*/
                                w->ld, w->kl[m], w->ku[m], l, i);
                        w->D_T[m][offset] += gsl_matrix_get(db, l - istart, m);
                    }
                }
            }
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    // TODO See Redmine ticket #2253
    /* Mass matrices are analytically guaranteed to be symmetric, positive
     * definite and so the sub- and super- diagonals are identical.  Any
     * floating point noise differences spoil the discrete symmetry and will
     * cause problems if people use BLAS DGMBV vs DSBMV-with-upper-diagonals vs
     * DSBMV-with-lower-diagonals to compute matrix-vector products or
     * factorizations.  Explicitly enforce the mass matrix's discrete symmetry
     * by copying the superdiagonals to the subdiagonals.
     */
    if (w->nderiv >= 0) {
        assert(w->kl[0] == w->ku[0]);
        for (int diag = 1; diag <= w->kl[0]; ++diag) {
            double * const sub = w->D_T[0] + w->ku[0] - diag;
            double * const sup = w->D_T[0] + w->ku[0] + diag;
            for (int i = 0; i < w->n - diag; ++i) {
                sub[(i+diag)*w->ld] = sup[i*w->ld];
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
            if (fabs(w->D_T[n][offset]) < tiny) {
                w->D_T[n][offset] = 0.0;
            }
        }
    }

    return SUZERAIN_SUCCESS;
}

/** Helper for  multiply_function_against_basis_function() */
struct multiply_function_against_basis_function_params {
    const suzerain_function * function;
    gsl_bspline_workspace *bw;
    size_t i;
    gsl_vector *Bk;
};

static double multiply_function_against_basis_function(double x, void *params)
{
    struct multiply_function_against_basis_function_params *p
        = (struct multiply_function_against_basis_function_params *) params;

    // Evaluate nonzero basis functions at x
    size_t istart, iend;
    const int error = gsl_bspline_eval_nonzero(x, p->Bk, &istart, &iend, p->bw);
    assert(!error); SUZERAIN_UNUSED(error);

    // If x is in the support of basis function i...
    if (istart <= p->i && p->i <= iend) {
        // ...return the product of the basis function and p->function...
#pragma warning(push,disable:981)
        return   gsl_vector_get(p->Bk, p->i - istart)
               * SUZERAIN_FN_EVAL(p->function, x);
#pragma warning(pop)
    } else {
        // ...otherwise return zero
        return 0.0;
    }
}

static int integrate_function_against_basis_functions(
        const suzerain_function * function,
        double * results,
        gsl_bspline_workspace *bw)
{
    // Prepare suzerain_function instance for use with integration scheme
    struct multiply_function_against_basis_function_params params = {
        function, bw, -1, gsl_vector_alloc(bw->k)
    };
    if (SUZERAIN_UNLIKELY(!params.Bk)) {
        SUZERAIN_ERROR("failed to allocate scratch space Bk", SUZERAIN_ENOMEM);
    }
    gsl_function product_f = {
        &multiply_function_against_basis_function, &params
    };

    // Prepare integration scheme using at most limit subintervals.
    // Limit chosen as no one should be expecting good B-spline interpolation
    // performance on functions requiring so many small subintervals.
    // Failure due to this choice will manifest as GSL_MAXITERs.
    const size_t limit = 3 * bw->knots->size;
    gsl_integration_workspace *iw = gsl_integration_workspace_alloc(limit);
    if (SUZERAIN_UNLIKELY(!iw)) {
        gsl_vector_free(params.Bk);
        SUZERAIN_ERROR("failed to allocate integration workspace",
                       SUZERAIN_ENOMEM);
    }

    // Defensively fill the results with NaNs until values are computed
    for (size_t i = 0; i < bw->n; ++i) results[i] = GSL_NAN;

    // Integrate the function against each basis function separately.
    // Precision is aggressive but GSL quits when it discovers roundoff error.
    // GSL_INTEG_GAUSS21 chosen empirically from a small sample of tests.
    // Bail as soon as we run into any trouble whatsoever...
    int stat = GSL_SUCCESS;
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    for (params.i = 0; params.i < bw->n && stat == GSL_SUCCESS; ++(params.i)) {
        double abserr;
        const double a = gsl_vector_get(bw->knots, params.i);
        const double b = gsl_vector_get(bw->knots, params.i + bw->k);
        stat = gsl_integration_qag(&product_f, a, b,
                                   GSL_DBL_EPSILON, 0.0,
                                   limit, GSL_INTEG_GAUSS21, iw,
                                   &results[params.i], &abserr);

        if (stat == GSL_EROUND) {  // ...except for roundoff/extrapolation error
            stat = GSL_SUCCESS;    // ...as we're asking for a strict tolerance
        }
    }
    gsl_set_error_handler(old_handler);

    // Free resources
    gsl_integration_workspace_free(iw);
    gsl_vector_free(params.Bk);

    if (SUZERAIN_UNLIKELY(stat != GSL_SUCCESS)) {
        char buffer[256];
        snprintf(buffer, sizeof(buffer),
                "Error(s) integrating function against basis: %s",
                gsl_strerror(stat));
        SUZERAIN_ERROR(buffer, stat);
    }
    return stat;
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

        return SUZERAIN_SUCCESS;

    case SUZERAIN_BSPLINEOP_GALERKIN_L2:
        return integrate_function_against_basis_functions(function, rhs, bw);
    }

    SUZERAIN_ERROR("unknown method in suzerain_bsplineop_interpolation_rhs",
                   SUZERAIN_ESANITY);
}

int
suzerain_bsplineop_interpolation_rhs_complex(
    const suzerain_zfunction * zfunction,
    complex_double *rhs,
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

        return SUZERAIN_SUCCESS;

    case SUZERAIN_BSPLINEOP_GALERKIN_L2: // FIXME Unimplemented
        break;
    }

    SUZERAIN_ERROR("unknown method in suzerain_bsplineop_interpolation_rhs_complex",
                   SUZERAIN_ESANITY);
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
    if (SUZERAIN_UNLIKELY(luw == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }

    /* Determine general banded matrix shape parameters */
    luw->n  = w->n;
    luw->kl = w->max_kl;
    luw->ku = w->max_ku;
    luw->ld = 2*luw->kl + luw->ku + 1; /* Increase per GBTRF, GBTRS */

    /* Allocate memory for LU factorization pivot storage */
    luw->ipiv = suzerain_blas_malloc(w->n * sizeof(luw->ipiv[0]));
    if (SUZERAIN_UNLIKELY(luw->ipiv == NULL)) {
        suzerain_bsplineop_lu_free(luw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }
    /* Flag that factorization has not occurred to help catch errors */
    luw->ipiv[0] = -1;

    /* Allocate memory for the matrix formed within lu_opform */
    luw->A_T = suzerain_blas_malloc(luw->ld*luw->n*sizeof(luw->A_T[0]));
    if (SUZERAIN_UNLIKELY(luw->A_T == NULL)) {
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
        suzerain_blas_free(luw->A_T);
        luw->A_T = NULL;

        suzerain_blas_free(luw->ipiv);
        luw->ipiv = NULL;

        suzerain_blas_free(luw);
    }
}

int
suzerain_bsplineop_lu_opaccumulate(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bsplineop_workspace * w,
    const double scale,
    suzerain_bsplineop_lu_workspace *luw)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(ncoefficients < 0)) {
        SUZERAIN_ERROR("Number of coefficients cannot be negative",
                       SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(ncoefficients > w->nderiv + 1)) {
        SUZERAIN_ERROR("More coefficients provided than derivatives available",
                       SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(luw->n < w->n)) {
        SUZERAIN_ERROR("Incompatible workspaces:  luw->n < w->n",
                       SUZERAIN_EINVAL);
    }
    /* Banded LU operator dimensions depend on largest derivative band */
    if (SUZERAIN_UNLIKELY(luw->ld < 2*w->max_kl + w->max_ku + 1)) {
        SUZERAIN_ERROR("Incompatible workspaces: luw->ld too small",
                       SUZERAIN_EINVAL);
    }
    /* Protect the user from incorrect API usage */
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(scale != 0.0 && luw->ipiv[0] >= 0)) {
#pragma warning(pop)
        SUZERAIN_ERROR("Unable to accumulate into a factored operator",
                       SUZERAIN_EINVAL);
    }

    /* Specialized handling for the first step: mass matrix accumulation */
    if (ncoefficients == 0) {

        /* No accumulation implies a scal operation on the operator */
        /* Tramples on the superdiagonals but it should not matter */
        suzerain_blas_dscal(luw->ld*luw->n, scale, luw->A_T, 1);

#pragma warning(push,disable:1572)
    } else if (scale == 0.0) {
#pragma warning(pop)

        /* Zero all operator storage including the superdiagonals */
        memset(luw->A_T, 0, luw->ld*luw->n*sizeof(luw->A_T[0]));

        /* Accumulate the mass matrix into operator storage, allowing */
        /* lower level routines to optimize for the factor of 1.      */
        suzerain_blas_dgb_acc(
            luw->n, luw->n, w->max_kl, w->max_ku,
            coefficients[0], w->D_T[0] - (w->max_ku - w->ku[0]), w->ld,
            1.0, luw->A_T + w->max_kl, luw->ld);

    } else {

        /* Accumulate the mass matrix into scaled operator storage */
        suzerain_blas_dgb_acc(
            luw->n, luw->n, w->max_kl, w->max_ku,
            coefficients[0], w->D_T[0] - (w->max_ku - w->ku[0]), w->ld,
            scale, luw->A_T + w->max_kl, luw->ld);

    }

    /* Accumulate any necessary higher derivative operators into luw->A_T */
    for (int k = 1; k < ncoefficients; ++k) {
        suzerain_blas_dgb_acc(
            luw->n, luw->n, w->max_kl, w->max_ku,
            coefficients[k], w->D_T[k] - (w->max_ku - w->ku[k]), w->ld,
            1.0, luw->A_T + w->max_kl, luw->ld);
    }

    /* Make -1*ipiv[0] indicate the number of opaccumulate calls performed */
    luw->ipiv[0] = GSL_MIN_INT(-1, luw->ipiv[0]-1);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_lu_opnorm(
    const suzerain_bsplineop_lu_workspace * luw,
    double * norm)
{
    /* Protect the user from incorrect API usage */
    if (SUZERAIN_UNLIKELY(luw->ipiv[0] >= 0)) {
        SUZERAIN_ERROR("Unable to obtain norm from a factored operator",
                       SUZERAIN_EINVAL);
    }

    /* As operator is not yet factored, it starts at A_T + kl per GBTRF */
    /* We want the infinity norm of A which is the one norm of A^T. */
    const int info = suzerain_blasext_dgbnorm1(luw->n, luw->n, luw->kl,
                                               luw->ku, luw->A_T + luw->kl,
                                               luw->ld, /* modified */ norm);
    if (info) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_blasext_dgbnorm1 reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_lu_factor(
    suzerain_bsplineop_lu_workspace * luw)
{
    /* Protect the user from incorrect API usage */
    if (SUZERAIN_UNLIKELY(luw->ipiv[0] >= 0)) {
        SUZERAIN_ERROR("Unable to factor an already-factored operator",
                       SUZERAIN_EINVAL);
    }

    /* Compute in-place LU factorization of the operator */
    const int info = suzerain_lapack_dgbtrf(luw->n, luw->n, luw->kl, luw->ku,
                                            luw->A_T, luw->ld, luw->ipiv);
    if (SUZERAIN_UNLIKELY(info)) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_dgbtrf reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    /* Factorization overwrote luw->ipiv[0] with something nonnegative;
     * This is our flag to check that factorization indeed occurred */

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_lu_rcond(
    const double norm,
    double * rcond,
    const suzerain_bsplineop_lu_workspace *luw)
{
    /* Protect the user from incorrect API usage */
    if (SUZERAIN_UNLIKELY(luw->ipiv[0] < 0)) {
        SUZERAIN_ERROR(
            "One of suzerain_bsplineop_lu_factor* not called before lu_rcond",
            SUZERAIN_EINVAL);
    }

    /* Allocate one chunk of memory for both work and iwork usage */
    double *work = suzerain_blas_malloc(   (3*luw->n)*sizeof(luw->A_T[0])
                                         + (  luw->n)*sizeof(int)       );
    if (SUZERAIN_UNLIKELY(work == NULL)) {
        SUZERAIN_ERROR("failed to allocate scratch space", SUZERAIN_ENOMEM);
    }

    /* The infinity condition number of A is one-based result for A^T. */
    const int info = suzerain_lapack_dgbcon('1', luw->n, luw->kl, luw->ku,
                                            luw->A_T, luw->ld, luw->ipiv,
                                            norm, rcond, work,
                                            (void *)(work + 3*luw->n));
    suzerain_blas_free(work);

    if (info) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_dgbcon reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_lu_opform(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace * luw)
{
    return suzerain_bsplineop_lu_opaccumulate(
            ncoefficients, coefficients, w, 0, luw);
}

int
suzerain_bsplineop_lu_opform_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace *luw)
{
    const double d_one = 1.0;
    return suzerain_bsplineop_lu_opform(1, &d_one, w, luw);
}

int
suzerain_bsplineop_lu_factor_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace *luw)
{
    const int info = suzerain_bsplineop_lu_opform_mass(w, luw);
    if (info) return info;
    return suzerain_bsplineop_lu_factor(luw);
}

static
int
suzerain_bsplineop_lu_solve_contiguous(
    int nrhs,
    double *B,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw)
{
    if (SUZERAIN_UNLIKELY(luw->ipiv[0] < 0)) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_lu_factor* not called before solve",
                SUZERAIN_EINVAL);
    }

    const int info = suzerain_lapack_dgbtrs('T', luw->n, luw->kl, luw->ku,
                                            nrhs, luw->A_T, luw->ld, luw->ipiv,
                                            B, ldb);
    if (SUZERAIN_UNLIKELY(info)) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_dgbtrs reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

static
int
suzerain_bsplineop_lu_solve_noncontiguous(
    int nrhs,
    double *B,
    int incb,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw)
{
    if (SUZERAIN_UNLIKELY(luw->ipiv[0] < -1)) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_lu_factor* not called before solve",
                SUZERAIN_EINVAL);
    }

    /* Allocate scratch space because GBTRS requires contiguous vectors */
    double * const scratch
        = suzerain_blas_malloc(luw->n*sizeof(scratch[0]));
    if (SUZERAIN_UNLIKELY(scratch == NULL)) {
        SUZERAIN_ERROR("failed to allocate scratch space", SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        double * const b_j = B + j*ldb;

        suzerain_blas_dcopy(luw->n, b_j, incb, scratch, 1);
        const int info = suzerain_lapack_dgbtrs('T', luw->n, luw->kl, luw->ku,
                                                1, /* One RHS at a time */
                                                luw->A_T, luw->ld, luw->ipiv,
                                                scratch, luw->n);
        if (info) {
            suzerain_blas_free(scratch);
            char buffer[80];
            snprintf(buffer, sizeof(buffer),
                    "suzerain_lapack_dgbtrs reported error %d", info);
            SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
        }
        suzerain_blas_dcopy(luw->n, scratch, 1, b_j, incb);
    }

    suzerain_blas_free(scratch);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_lu_solve(
    int nrhs,
    double *B,
    int incb,
    int ldb,
    const suzerain_bsplineop_lu_workspace *luw)
{
    if (incb == 1) {
        return suzerain_bsplineop_lu_solve_contiguous(
                nrhs, B, ldb, luw);
    } else {
        return suzerain_bsplineop_lu_solve_noncontiguous(
                nrhs, B, incb, ldb, luw);
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
    if (SUZERAIN_UNLIKELY(luzw == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate space for workspace",
                            SUZERAIN_ENOMEM);
    }

    /* Determine general banded matrix shape parameters */
    luzw->n  = w->n;
    luzw->kl = w->max_kl;
    luzw->ku = w->max_ku;
    luzw->ld = 2*luzw->kl + luzw->ku + 1; /* Increase per GBTRF, GBTRS */

    /* Allocate memory for LU factorization pivot storage */
    luzw->ipiv = suzerain_blas_malloc(w->n * sizeof(luzw->ipiv[0]));
    if (SUZERAIN_UNLIKELY(luzw->ipiv == NULL)) {
        suzerain_bsplineop_luz_free(luzw);
        SUZERAIN_ERROR_NULL("failed to allocate space for pivot storage",
                            SUZERAIN_ENOMEM);
    }
    /* Flag that decomposition has not occurred; Allows catching errors */
    luzw->ipiv[0] = -1;

    /* Allocate memory for the matrix formed within luz_form */
    luzw->A_T = suzerain_blas_malloc(luzw->ld*luzw->n*sizeof(luzw->A_T[0]));
    if (SUZERAIN_UNLIKELY(luzw->A_T == NULL)) {
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
        suzerain_blas_free(luzw->A_T);
        luzw->A_T = NULL;

        suzerain_blas_free(luzw->ipiv);
        luzw->ipiv = NULL;

        suzerain_blas_free(luzw);
    }
}

int
suzerain_bsplineop_luz_opaccumulate(
    int ncoefficients,
    const complex_double *coefficients,
    const suzerain_bsplineop_workspace *w,
    const complex_double scale,
    suzerain_bsplineop_luz_workspace *luzw)
{
    /* Parameter sanity checks */
    if (SUZERAIN_UNLIKELY(ncoefficients < 0)) {
        SUZERAIN_ERROR("Number of coefficients cannot be negative",
                       SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(ncoefficients > w->nderiv + 1)) {
        SUZERAIN_ERROR("More coefficients provided than derivatives available",
                       SUZERAIN_EINVAL);
    }
    if (SUZERAIN_UNLIKELY(luzw->n < w->n)) {
        SUZERAIN_ERROR("Incompatible workspaces: luzw->n < w->n",
                       SUZERAIN_EINVAL);
    }
    /* Banded LU operator dimensions depend on largest derivative band */
    if (SUZERAIN_UNLIKELY(luzw->ld < 2*w->max_kl + w->max_ku + 1)) {
        SUZERAIN_ERROR("Incompatible workspaces: luzw->ld too small",
                       SUZERAIN_EINVAL);
    }
    /* Protect the user from incorrect API usage */
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(scale != 0.0 && luzw->ipiv[0] >= 0)) {
#pragma warning(pop)
        SUZERAIN_ERROR("Unable to accumulate into a factored operator",
                       SUZERAIN_EINVAL);
    }

    /* Specialized handling for the first step: mass matrix accumulation */
    if (ncoefficients == 0) {

        /* No accumulation implies a scal operation on the operator */
        /* Tramples on the superdiagonals but it should not matter */
        suzerain_blas_zscal(luzw->ld*luzw->n, scale, luzw->A_T, 1);

#pragma warning(push,disable:1572)
    } else if (scale == 0.0) {
#pragma warning(pop)

        /* Zero all operator storage including the superdiagonals */
        memset(luzw->A_T, 0, luzw->ld*luzw->n*sizeof(luzw->A_T[0]));

        /* Accumulate the mass matrix into operator storage, allowing */
        /* lower level routines to optimize for the factor of 1.      */
        suzerain_blas_zgb_acc_d(
            luzw->n, luzw->n, w->max_kl, w->max_ku,
            coefficients[0], w->D_T[0] - (w->max_ku - w->ku[0]), w->ld,
            1, luzw->A_T + w->max_kl, luzw->ld);

    } else {

        /* Accumulate the mass matrix into scaled operator storage */
        suzerain_blas_zgb_acc_d(
            luzw->n, luzw->n, w->max_kl, w->max_ku,
            coefficients[0], w->D_T[0] - (w->max_ku - w->ku[0]), w->ld,
            scale, luzw->A_T + w->max_kl, luzw->ld);

    }

    /* Accumulate necessary any higher derivative operators into luzw->A_T */
    for (int k = 1; k < ncoefficients; ++k) {
        suzerain_blas_zgb_acc_d(
            luzw->n, luzw->n, w->max_kl, w->max_ku,
            coefficients[k], w->D_T[k] - (w->max_ku - w->ku[k]), w->ld,
            1, luzw->A_T + w->max_kl, luzw->ld);
    }

    /* Make -1*ipiv[0] indicate the number of opaccumulate calls performed */
    luzw->ipiv[0] = GSL_MIN_INT(-1, luzw->ipiv[0]-1);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_luz_opnorm(
    const suzerain_bsplineop_luz_workspace * luzw,
    double * norm)
{
    /* Protect the user from incorrect API usage */
    if (SUZERAIN_UNLIKELY(luzw->ipiv[0] >= 0)) {
        SUZERAIN_ERROR("Unable to obtain norm from a factored operator",
                       SUZERAIN_EINVAL);
    }

    /* As operator is not yet factored, it starts at A_T + kl per GBTRF */
    /* We want the infinity norm of A which is the one norm of A^T. */
    const int info = suzerain_blasext_zgbnorm1(
                                luzw->n, luzw->n, luzw->kl, luzw->ku,
                                luzw->A_T + luzw->kl, luzw->ld,
                                /* modified */ norm);
    if (SUZERAIN_UNLIKELY(info)) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_blasext_zgbnorm1 reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_luz_factor(
    suzerain_bsplineop_luz_workspace * luzw)
{
    /* Protect the user from incorrect API usage */
    if (SUZERAIN_UNLIKELY(luzw->ipiv[0] >= 0)) {
        SUZERAIN_ERROR("Unable to factor an already-factored operator",
                       SUZERAIN_EINVAL);
    }

    /* Compute in-place LU factorization of the operator */
    const int info = suzerain_lapack_zgbtrf(luzw->n, luzw->n, luzw->kl,
                                            luzw->ku, luzw->A_T, luzw->ld,
                                            luzw->ipiv);
    if (SUZERAIN_UNLIKELY(info)) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_zgbtrf reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    /* Factorization overwrote luzw->ipiv[0] == -1; This is our flag to
     * check that factorization occurred in suzerain_bsplineop_luz_solve */

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_luz_rcond(
    const double norm,
    double * rcond,
    const suzerain_bsplineop_luz_workspace *luzw)
{
    if (SUZERAIN_UNLIKELY(luzw->ipiv[0] < 0)) {
        SUZERAIN_ERROR(
            "One of suzerain_bsplineop_luz_factor* not called before luz_rcond",
            SUZERAIN_EINVAL);
    }

    /* Allocate one chunk of memory for both work and rwork usage */
    double *work = suzerain_blas_malloc(   (2*(2*luzw->n))*sizeof(luzw->A_T[0])
                                         + (     luzw->n )*sizeof(luzw->A_T[0]));
    if (SUZERAIN_UNLIKELY(work == NULL)) {
        SUZERAIN_ERROR("failed to allocate scratch space", SUZERAIN_ENOMEM);
    }

    /* The infinity condition number of A is one-based result for A^T. */
    const int info = suzerain_lapack_zgbcon('1', luzw->n, luzw->kl, luzw->ku,
                                            luzw->A_T, luzw->ld, luzw->ipiv,
                                            norm, /* modified */ rcond,
                                            (complex_double *) work,
                                            work + 2*(2*luzw->n));
    suzerain_blas_free(work);

    if (SUZERAIN_UNLIKELY(info)) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_zgbcon reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

static
int
suzerain_bsplineop_luz_solve_contiguous(
    int nrhs,
    complex_double *B,
    int ldb,
    const suzerain_bsplineop_luz_workspace *luzw)
{
    if (SUZERAIN_UNLIKELY(luzw->ipiv[0] < 0)) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_luz_factor* not called before solve",
                SUZERAIN_EINVAL);
    }

    const int info = suzerain_lapack_zgbtrs('T', luzw->n, luzw->kl, luzw->ku,
                                            nrhs, luzw->A_T, luzw->ld,
                                            luzw->ipiv, B, ldb);
    if (SUZERAIN_UNLIKELY(info)) {
        char buffer[80];
        snprintf(buffer, sizeof(buffer),
                 "suzerain_lapack_zgbtrs reported error %d", info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }

    return SUZERAIN_SUCCESS;
}

static
int
suzerain_bsplineop_luz_solve_noncontiguous(
    int nrhs,
    complex_double *B,
    int incb,
    int ldb,
    const suzerain_bsplineop_luz_workspace *luzw)
{
    if (SUZERAIN_UNLIKELY(luzw->ipiv[0] < 0)) {
        SUZERAIN_ERROR(
                "One of suzerain_bsplineop_luz_factor* not called before solve",
                SUZERAIN_EINVAL);
    }

    /* Allocate scratch space because ZGBTRS requires contiguous vectors */
    complex_double * const scratch
            = suzerain_blas_malloc(luzw->n*sizeof(scratch[0]));
    if (SUZERAIN_UNLIKELY(scratch == NULL)) {
        SUZERAIN_ERROR("failed to allocate scratch space",
                       SUZERAIN_ENOMEM);
    }

    for (int j = 0; j < nrhs; ++j) {
        complex_double * const b_j = B + j*ldb;

        suzerain_blas_zcopy(luzw->n, b_j, incb, scratch, 1);
        const int info = suzerain_lapack_zgbtrs('T', luzw->n, luzw->kl,
                                                luzw->ku,
                                                1, /* One RHS at a time */
                                                luzw->A_T, luzw->ld, luzw->ipiv,
                                                scratch, luzw->n);
        if (info) {
            suzerain_blas_free(scratch);
            char buffer[80];
            snprintf(buffer, sizeof(buffer),
                    "suzerain_lapack_zgbtrs reported error %d", info);
            SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
        }
        suzerain_blas_zcopy(luzw->n, scratch, 1, b_j, incb);
    }

    suzerain_blas_free(scratch);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bsplineop_luz_solve(
    int nrhs,
    complex_double *B,
    int incb,
    int ldb,
    const suzerain_bsplineop_luz_workspace *luzw)
{
    if (incb == 1) {
        return suzerain_bsplineop_luz_solve_contiguous(
                nrhs, B, ldb, luzw);
    } else {
        return suzerain_bsplineop_luz_solve_noncontiguous(
                nrhs, B, incb, ldb, luzw);
    }
}

int
suzerain_bsplineop_luz_opform(
    int ncoefficients,
    const complex_double *coefficients,
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace * luzw)
{
    return suzerain_bsplineop_luz_opaccumulate(
            ncoefficients, coefficients, w, 0, luzw);
}

int
suzerain_bsplineop_luz_opform_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace * luzw)
{
    const complex_double z_one = 1;
    return suzerain_bsplineop_luz_opform(1, &z_one, w, luzw);
}

int
suzerain_bsplineop_luz_factor_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace * luzw)
{
    const int info = suzerain_bsplineop_luz_opform_mass(w, luzw);
    if (SUZERAIN_UNLIKELY(info)) return info;
    return suzerain_bsplineop_luz_factor(luzw);
}
