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
 * bspline_operator.h: operator creation routines for a bspline basis
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_BSPLINE_OPERATOR_H
#define PECOS_SUZERAIN_BSPLINE_OPERATOR_H

#include <suzerain/function.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS   /* empty */
#endif

__BEGIN_DECLS

/* Conditional forward declarations of implementation-related structs */
#ifndef __GSL_MATRIX_DOUBLE_H__
typedef struct gsl_matrix gsl_matrix;
#endif
#ifndef __GSL_BSPLINE_H__
typedef struct gsl_bspline_workspace gsl_bspline_workspace;
typedef struct gsl_bspline_deriv_workspace gsl_bspline_deriv_workspace;
#endif

enum suzerain_bspline_method {
    SUZERAIN_BSPLINE_COLLOCATION_GREVILLE = 1
};

typedef struct {
    int order;          /* spline order, piecewise degree is order - 1 */
    int nbreakpoints;
    int nderivatives;
    int ncoefficients;
    enum suzerain_bspline_method method;
    int * kl;
    int * ku;
    int * lda;
    int * storagesize;
    gsl_bspline_workspace * bw;
    gsl_bspline_deriv_workspace * dbw;
    gsl_matrix * db;
    double **D;
} suzerain_bspline_workspace;

suzerain_bspline_workspace *
suzerain_bspline_alloc(
    int order,
    int nderivatives,
    int nbreakpoints,
    const double * breakpoints,
    enum suzerain_bspline_method method);

void
suzerain_bspline_free(
    suzerain_bspline_workspace *w);

int
suzerain_bspline_ncoefficients(
    const suzerain_bspline_workspace *w);

int
suzerain_bspline_apply_operator(
    int nderivative,
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_workspace *w);

int
suzerain_bspline_find_coefficient_rhs(
    const suzerain_function * function,
    double * coefficient_rhs,
    const suzerain_bspline_workspace *w);

typedef struct {
    int ncoefficients;
    int kl;
    int ku;
    int lda;
    int storagesize;
    int *ipiv;
    double *A;
} suzerain_bspline_lu_workspace;

suzerain_bspline_lu_workspace *
suzerain_bspline_lu_alloc(
    const suzerain_bspline_workspace *w);

void
suzerain_bspline_lu_free(
    suzerain_bspline_lu_workspace *luw);

int
suzerain_bspline_lu_form_general(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bspline_workspace * w,
    suzerain_bspline_lu_workspace *luw);

int
suzerain_bspline_lu_form_mass(
    const suzerain_bspline_workspace * w,
    suzerain_bspline_lu_workspace *luw);

int
suzerain_bspline_lu_solve(
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_lu_workspace *luw);

__END_DECLS

#endif // PECOS_SUZERAIN_BSPLINE_OPERATOR_H
