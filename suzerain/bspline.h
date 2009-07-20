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
 * bspline.h: bspline basis manipulation and operator routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_BSPLINE_OPERATOR_H
#define PECOS_SUZERAIN_BSPLINE_OPERATOR_H

#include <suzerain/function.h>

/** @file
 * Provides bspline basis evaluation and derivative operator routines.  These
 * functions provide bspline basis function evaluation, differentiation,
 * derivative linear operator construction, and operator application routines.
 *
 * Internally, the code builds upon the GNU Scientific Library (GSL) bspline
 * routines and uses banded matrix BLAS and LAPACK functionality where
 * possible.
 */

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

/** Indicates the method chosen to compute derivative operators.  */
enum suzerain_bspline_method {
    /** Form derivative operators using the Greville abscissae.  */
    SUZERAIN_BSPLINE_COLLOCATION_GREVILLE = 1
};

/**
 * Encapsulates basis function and derivative operator information.
 * Callers obtain a workspace using suzerain_bspline_alloc() and
 * release it using suzerain_bspline_free().
 */
typedef struct {

    /**
     * Spline order per GSL/PPPACK conventions
     * For example, 3 denotes piecewise quadratics and 4 denotes piecewise
     * cubics.
     **/
    int order;

    /** Number of breakpoints in the basis */
    int nbreakpoints;

    /** Maximum derivative requested */
    int nderivatives;

    /** Number of coefficients or degrees of freedom in the basis */
    int ncoefficients;

    /** Method chosen to form derivative operators */
    enum suzerain_bspline_method method;

    /** @name Banded derivative matrix storage details
     * Each of the \c 0 through \c nderivative (inclusive) derivative operators
     * is stored using BLAS banded matrix conventions.
     * @{
     */

    /** Number of subdiagonals in each derivative */
    int * kl;

    /** Number of superdiagonals in each derivative */
    int * ku;

    /** Leading dimension in each derivative */
    int * lda;

    /** Size of each banded derivative storage measured using <tt>double</tt>s */
    int * storagesize;

    /**
     * Raw data storage for each banded derivative operator matrix
     * \c D[0] is the storage for the 0th derivative, D[1] is the storage for
     * the 1st derivative, etc..
     **/
    double **D;

    /* @} */

    /**
     * GSL workspace for basis function evaluation
     * \internal
     **/
    gsl_bspline_workspace * bw;

    /**
     * GSL workspace for basis function derivative evaluation
     * \internal
     **/
    gsl_bspline_deriv_workspace * dbw;

    /**
     * Used when computing basis function derivatives
     * \internal
     **/
    gsl_matrix * db;

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
suzerain_bspline_evaluate(
    int nderivative,
    const double * coefficients,
    int npoints,
    const double * points,
    double * values,
    int ldvalues,
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
