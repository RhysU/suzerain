/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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

#include <stdlib.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <suzerain/bspline_operator.h>
#include <suzerain/suzerain_error.h>

/* Compute the BLAS-compatible offset to a(i,j) for general banded matrices */
/* a(i,j) -> storage(ku+i-j,j) where storage is column-major with LDA lda */
/* Note missing constant one compared with Fortran because C 0-indexes arrays */
#define GB_OFFSET(lda, ku, i, j) ((j)*(lda)+(ku+i-j))

suzerain_bspline_operator_workspace *
suzerain_bspline_operator_alloc(int order,
                                int nbreakpoints,
                                int nderivatives,
                                enum suzerain_bspline_operator_method method)
{

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

    /* Save bspline parameters */
    w->order        = order;
    w->nbreakpoints = nbreakpoints;
    w->nderivatives = nderivatives;
    w->n            = nbreakpoints + order - 2; /* assumes max continuity */

    /* Storage parameters for BLAS/lapack-compatible general band matrix */
    w->kl          = (bandwidth - 1) / 2;
    w->ku          = (bandwidth - 1) / 2;
    w->lda         = w->kl + w->ku + 1;
    w->storagesize = w->lda * w->n;

    /* Allocate space for pointers to matrices */
    w->D = malloc((w->nderivatives+1) * sizeof(double *));
    if (w->D == NULL) {
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix pointers",
                            SUZERAIN_ENOMEM);
    }
    /* allocate memory for all matrices in one contiguous block */
    w->D[0] = malloc((w->nderivatives+1) * w->storagesize * sizeof(double));
    if (w->D[0] == NULL) {
        free(w->D);
        free(w);
        SUZERAIN_ERROR_NULL("failed to allocate space for matrix storage",
                            SUZERAIN_ENOMEM);
    }
    /* w->D[0] now points to D[0] which will be the mass matrix */
    /* Establish pointers for D[1] = d/dx, d[2] = d^2/dx^2, ...*/
    for (i = 0; i < nderivatives; ++i) {
        w->D[i+1] = w->D[i] + w->storagesize;
    }

    return w;
}

void
suzerain_bspline_operator_free(suzerain_bspline_operator_workspace * w)
{
    free(w->D[0]);
    /* D[1], ..., D[nderivatives-1] allocated through w->D[0]; no free() */
    free(w->D);
    free(w);
}

int
suzerain_bspline_operator_create(const double * breakpoints,
                                 suzerain_bspline_operator_workspace *w)
{
    gsl_vector_const_view breakpoints_view
        = gsl_vector_const_view_array(breakpoints, w->nbreakpoints);
    gsl_bspline_workspace *bw;
    gsl_bspline_deriv_workspace *bdw;
    gsl_matrix *db;

    /* Setup calls to use GSL B-spline functionality */
    bw = gsl_bspline_alloc(w->order, w->nbreakpoints);
    if (bw == NULL) {
        SUZERAIN_ERROR("failure allocating bspline workspace",
                       SUZERAIN_ENOMEM);
    }
    if (w->n != gsl_bspline_ncoeffs(bw)) {
        gsl_bspline_free(bw);
        SUZERAIN_ERROR("bspline coefficient count does not match workspace",
                       SUZERAIN_EINVAL);
    }
    bdw = gsl_bspline_deriv_alloc(gsl_bspline_order(bw));
    if (bdw == NULL) {
        gsl_bspline_free(bw);
        SUZERAIN_ERROR("failure allocating bspline derivative workspace",
                       SUZERAIN_ENOMEM);
    }
    if (gsl_bspline_knots(&breakpoints_view.vector, bw)) {
        gsl_bspline_deriv_free(bdw);
        gsl_bspline_free(bw);
        SUZERAIN_ERROR("failure seting bspline breakpoints",
                       SUZERAIN_EFAILED);
    }


    /* Tear down calls for GSL B-spline functionality */
    gsl_bspline_deriv_free(bdw);
    gsl_bspline_free(bw);

    return SUZERAIN_SUCCESS;
}
