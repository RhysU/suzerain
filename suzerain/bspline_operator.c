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

    /* Per gsl_bspline_ncoeffs; assumes highest continuity possible used */
    const int ncoeff = nbreakpoints + order - 2;

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
        /* Compute bandwidth of resulting operator matrices
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

    /* Save bspline details */
    w->order  = order;
    w->nbreakpoints = nbreakpoints;
    w->nderivatives = nderivatives;
    /* per gsl_bspline_ncoeffs */

    /* Storage parameters per http://www.intel.com/software/products/mkl/docs/WebHelp/appendices/mkl_appB_MA.html */
    w->D_kl          = (bandwidth - 1) / 2;
    w->D_ku          = (bandwidth - 1) / 2;
    w->D_lda         = w->D_kl + w->D_ku + 1;
    w->D_storagesize = w->D_lda * ncoeff;
    w->M_kl          = w->D_kl;
    w->M_ku          = w->D_kl + w->D_ku;  /* Extra for LU by DGBTRF */
    w->M_lda         = w->M_kl + w->M_ku + 1;
    w->M_storagesize = w->M_lda * ncoeff;

    /* Allocate space for mass matrix storage */
    /* calloc ensures zeros for all matrix indices */
    w->M = calloc(w->M_storagesize, sizeof(double));
    if (w->M == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for mass matrix",
                            SUZERAIN_ENOMEM);
        free(w);
    }

    /* Allocate space for pointers to operator matrices */
    w->D = malloc(w->nderivatives * sizeof(double *));
    if (w->D == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for derivative pointers",
                            SUZERAIN_ENOMEM);
        free(w->M);
        free(w);
    }
    /* simultaneously allocate memory for all derivative matrices */
    /* calloc ensures zeros in all matrix indices */
    w->D[0] = calloc(w->nderivatives * w->D_storagesize, sizeof(double));
    if (w->D[0] == NULL) {
        SUZERAIN_ERROR_NULL("failed to allocate space for derivative matrices",
                            SUZERAIN_ENOMEM);
        free(w->D);
        free(w->M);
        free(w);
    }
    /* w->D[0] now points to D[0] = d/dx */
    /* Establish pointers for higher derivatives, e.g. D[1] = d^2/dx^2, ...*/
    for (i = 1; i < nderivatives; ++i) {
        w->D[i] = w->D[i-1] + w->D_storagesize;
    }

    return w;
}

void
suzerain_bspline_operator_free(suzerain_bspline_operator_workspace * w)
{
    free(w->D[0]);
    /* D[1], ..., D[nderivatives-1] allocated through w->D[0]; no free() */
    free(w->D);
    free(w->M);
    free(w);
}

int
suzerain_bspline_operator_create(const double * breakpoints,
                                 suzerain_bspline_operator_workspace *w)
{
    /* IMPLEMENT ME */

    return SUZERAIN_SUCCESS;
}
