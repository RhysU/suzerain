/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
 * filterop.h: Discrete filtering operator routines
 * $Id$
 */

#ifndef SUZERAIN_FILTEROP_H
#define SUZERAIN_FILTEROP_H

#include <suzerain/common.h>
#include <suzerain/complex.h>
#include <suzerain/function.h>

/** @file
 * Provides discrete, banded filtering operators.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Indicates the method chosen to compute filter operators.  */
enum suzerain_filterop_method {

    /**
     * Form derivative operators using hyperviscosity per <a
     * href="http://dx.doi.org/10.1016/j.jcp.2004.09.011">Cook and Cabot JCP
     * 2005</a>.
     */
    SUZERAIN_FILTEROP_COOKCABOT2005 = 1
};

/** Indicates the chosen filtering behavior at the boundaries.  */
enum suzerain_filterop_boundary_treatment {

    /**
     * Apply the interior scheme at all points, including the points near the
     * wall.  Out-of-domain values are treated as if they are zero.
     */
    SUZERAIN_FILTEROP_BOUNDARY_IGNORE = 1,

    /**
     * No filtering is performed at the boundary points nearest the
     * edge of the domain.
     */
    SUZERAIN_FILTEROP_BOUNDARY_NOFILTER,

    /**
     * Reflect near-boundary points outside the domain and apply
     * the interior scheme.
     */
    SUZERAIN_FILTEROP_BOUNDARY_SYMMETRY,

    /**
     * Enforce periodic conditions across the domain boundaries.
     */
    SUZERAIN_FILTEROP_BOUNDARY_PERIODIC
};

/**
 * Encapsulates filter operator information.  Callers obtain a workspace using
 * suzerain_filterop_alloc() and release it using suzerain_filterop_free().
 */
typedef struct suzerain_filterop_workspace {

    /** Method chosen to form derivative operators */
    enum suzerain_filterop_method method;

    /**
     * Boundary treatment chosen for first vector index.  That is,
     * <tt>y[0]</tt>.
     */
    enum suzerain_filterop_boundary_treatment b_first;

    /**
     * Boundary treatment chosen for last vector index.  That is,
     * <tt>y[n - 1]</tt>.
     */
    enum suzerain_filterop_boundary_treatment b_last;

    /** Number of degrees of freedom in each vector */
    int n;

    /** Number of subdiagonals in operator \f$A^\mathrm{T}\f$. */
    int klat;

    /** Number of superdiagonals in operator \f$A^\mathrm{T}\f$. */
    int kuat;

    /** Leading dimension of storage for operator \f$A^{\mathrm{T}}\f$. */
    int ldat;

    /**
     * Raw data storage for the unfactored and factored \f$A^{\mathrm{T}}\f$.
     * The transpose of the operator is stored according to LAPACK \c GBTRF
     * conventions regarding \c kl, \c ku, and \c n.  In particular, this means
     * <tt>A + kl</tt> is the starting general band storage location for the
     * transpose of the unfactored operator.
     */
    double *A_T;

    /**
     *  Pivot matrix \c P from the \c LUP decomposition of \f$A^{\mathrm{T}}\f$.
     */
    int *ipiva;

    /** Number of subdiagonals in operator \f$B^\mathrm{T}\f$. */
    int klbt;

    /** Number of superdiagonals in operator \f$B^\mathrm{T}\f$. */
    int kubt;

    /** Leading dimension of storage for operator \f$B^{\mathrm{T}}\f$. */
    int ldbt;

    /**
     * Raw data storage for the unfactored \f$B^{\mathrm{T}}\f$.  The operator
     * is stored according to BLAS \c GBMV conventions regarding \c klbt, \c
     * kubt, and \c n.
     */
    double *B_T;

} suzerain_filterop_workspace;

/** @name Allocation and deallocation */
/**@{*/

/**
 * Allocate a filtering operator workspace.
 *
 * @param[in] n             Length of the state vector to be filtered.
 * @param[in] method        Filtering method or scheme to perform.
 * @param[in] method_params Method-specific adjustable parameters.
 *                          Set \c NULL to use published defaults.  Length of
 *                          \c method_params required depends upon the method.
 *                          For example, SUZERAIN_FILTEROP_COOKCABOT2005 takes
 *                          only <tt>method_params[0]</tt>.
 * @param[in] b_first       Boundary treatment for the first vector index.
 * @param[in] b_last        Boundary treatment for the last vector index.
 *
 * @return a workspace instance on success.
 *      On failure calls suzerain_error() and returns NULL.
 *
 * \memberof suzerain_filterop_workspace
 */
suzerain_filterop_workspace *
suzerain_filterop_alloc(
    const int n,
    const enum suzerain_filterop_method method,
    const double *method_params,
    const enum suzerain_filterop_boundary_treatment b_first,
    const enum suzerain_filterop_boundary_treatment b_last);

/**
 * Factorize the \f$A^\mathrm{T}\f$ operator within the workspace.
 *
 * @param[in,out] w Workspace to modify
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_filterop_workspace
 */
int
suzerain_filterop_factorize(
    suzerain_filterop_workspace *w);

/**
 * Free a previously allocated workspace.
 *
 * @param[in] w Workspace to free.
 *
 * \memberof suzerain_filterop_workspace
 */
void
suzerain_filterop_free(
    suzerain_filterop_workspace *w);

/**@}*/

/** @name Primitive operations */
/**@{*/

/**
 * Apply the scaled \f$\alpha{}B\f$ filtering operator to real coefficients \c
 * x storing the result in \c y.
 *
 * @param[in]  alpha Real scaling factor \f$\alpha\f$ to apply.
 * @param[in]  x    Coefficients to be multiplied.
 * @param[in]  incx Stride between elements stored in \c x.
 * @param[out] y    Storage for the result.
 * @param[in]  incy Stride between elements stored in \c x.
 * @param[in]  w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_filterop_workspace
 */
int
suzerain_filterop_apply(
    const double alpha,
    const double *x,
    const int incx,
    double *y,
    const int incy,
    const suzerain_filterop_workspace *w);

/**
 * Apply \f$A^{-1}\f$ in place to matrix \c X using the previously factored
 * operator \f$A^{\mathrm{T}}\f$.
 *
 * @param[in]     nrhs Number of right hand sides within \c X.
 * @param[in,out] X    Right hand sides to solve and the resulting solution.
 * @param[in]     ldx  Leading dimension between columns in matrix \c X.
 * @param[in]     w    Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_filterop_workspace
 */
int
suzerain_filterop_solve(
    const int nrhs,
    double *X,
    const int ldx,
    const suzerain_filterop_workspace * w);

/**@}*/

/** @name High level operations */
/**@{*/

/**
 * Apply the filtering operator \f$F = \alpha{}A^{-1}B\f$ to the vector \c x
 * overwriting vector \c y with the result.
 *
 * @param[in]  alpha Real scaling factor \f$\alpha\f$ to apply.
 * @param[in]  x    Coefficients to be multiplied.
 * @param[in]  incx Stride between elements stored in \c x.
 * @param[out] y    Contiguous storage for the result.  Will be overwritten.
 * @param[in]  w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_filterop_workspace
 */
int
suzerain_filterop_filter(
    const double alpha,
    const double *x,
    const int incx,
    double *y,
    const suzerain_filterop_workspace *w);

/**@}*/

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SUZERAIN_FILTEROP_H
