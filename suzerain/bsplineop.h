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

#ifndef SUZERAIN_BSPLINEOP_H
#define SUZERAIN_BSPLINEOP_H

/** @file
 * Logic for working with B-spline operators.
 */

#include <gsl/gsl_bspline.h>

#include <suzerain/complex.h>
#include <suzerain/function.h>

/* TODO Add usage example */

/** @file
 * Provides B-spline-based banded operator construction and application
 * routines.  Both real- and complex- valued operator application and inversion
 * are supported.  The derivative operators map a function's spline
 * coefficients to an approximation of the derivative's spline coefficients.
 *
 * For collocation-based operators, the approximation \f$\tilde{\gamma}(x) =
 * \sum_{i} \beta_{i} B_{i}(x)\f$ to the continuous operator
 * \f$\mathcal{D}\f$ acting on
 * \f$\tilde{\phi}(x) = \sum_{i} \alpha_{i} B_{i}(x)\f$ must satisfy
 * \f$\tilde{\gamma}(x_j) = \left(\mathcal{D}\tilde{\phi}\right)(x_j)\f$
 * at each collocation point \f$x_j\f$.  This implies
 * \f[
 *        \sum_{i} \beta_{i} B_{i}(x_j)
 *      = \sum_{i} \alpha_{i} \mathcal{D}(B_i)(x_j)
 *      \quad \forall{} j.
 * \f]
 * The \f$k\f$-th collocation derivative operator matrix is given by
 * \f$D^k_{i,j}=\left(B^{(k)}_j(x_i)\right)\f$.  Conveniently, the
 * matrix-vector product of a collocation derivative matrix with a B-spline
 * coefficient vector is equivalent to evaluating the linear combination of the
 * basis functions at the collocation points.
 *
 * For Galerkin-based operators, the relationship \f$ \left( B_j,
 * \tilde{\gamma}\right) = \left(B_j, \mathcal{D}\tilde{\phi}\right)\f$ must
 * hold for each spline basis index \f$j\f$ where \f$(f,g)\f$ denotes the
 * \f$L_2\f$ inner product.  This implies
 * \f[
 *          \sum_{i} \beta_{i} \left( B_j, B_i \right)
 *      =   \sum_{i} \alpha_{j} \left( B_j, \mathcal{D}B_i \right).
 * \f]
 * The \f$k\f$-th Galerkin derivative operator matrix is given by
 * \f$D^k_{i,j}=\left(B_j,B^{(k)}_i\right)\f$.  For either method, applying
 * the \f$k\f$-th derivative operator requires solving the linear problem
 * \f$D^{0} \beta = D^{k} \alpha\f$.  The method to use is chosen
 * through the ::suzerain_bsplineop_method value provided to
 * suzerain_bsplineop_alloc().
 *
 * All operator matrices are stored \e transposed in column-major (Fortran)
 * storage.  The operator transposes are stored so that discrete operator
 * application (which is nothing but a matrix-vector product) has an optimal
 * memory access pattern.  Multiple threads may call the routines
 * simultaneously provided that each thread has its own workspace instances.
 * The behavior is undefined if two threads simultaneously share a workspace.
 * Internally, the code builds upon the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Basis-Splines.html">
 * B-spline routines</a> and uses banded matrix BLAS and LAPACK functionality
 * where possible.
 *
 * @see <a href="http://dx.doi.org/10.1006/jcph.2001.6919">A Critical
 * Evaluation of the Resolution Properties of B-Spline and Compact Finite
 * Difference Methods</a> by Kwok, Moser, and Jimenez published in the
 * <em>Journal of Computational Physics</em>, volume 174, number 2, pages
 * 510-551 (December 2001) for more details on B-spline derivative operators
 * and their properties.
 * @see <a href="http://dx.doi.org/10.1016/j.apnum.2004.04.002"> Higher order
 * B-spline collocation at the Greville abscissae</a> by Johnson published in
 * <em>Applied Numerical Mathematics</em>, volume 52, number 1, pages 63-75
 * (January 2005) for information about Greville abscissae.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Indicates the method chosen to compute derivative operators.  */
enum suzerain_bsplineop_method {

    /**
     * Form derivative operators using collocation at the Greville abscissae.
     */
    SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE = 1,

    /**
     * Form derivative operators using Galerkin's method with the \f$L^2\f$
     * inner product.
     */
    SUZERAIN_BSPLINEOP_GALERKIN_L2 
};

/**
 * Encapsulates basis function and derivative operator information.  Callers
 * obtain a workspace using suzerain_bsplineop_alloc() and release it using
 * suzerain_bsplineop_free().
 */
typedef struct suzerain_bsplineop_workspace {

    /** Method chosen to form derivative operators */
    enum suzerain_bsplineop_method method;

    /**
     * Spline order per GSL/PPPACK conventions.  For example, 3 denotes
     * piecewise quadratics and 4 denotes piecewise cubics.
     **/
    int k;

    /** Number of degrees of freedom in the basis */
    int n;

    /** Maximum derivative operator within this workspace, inclusive. */
    int nderiv;

    /** @name Banded derivative matrix storage details
     * The \e transpose of each of the \c 0 through \c nderiv (inclusive)
     * derivative operators is stored using BLAS banded matrix conventions.
     * @{
     */

    /** Number of subdiagonals in each derivative operator's transpose */
    int * kl;

    /** Number of superdiagonals in each derivative operator's transpose */
    int * ku;

    /** Maximum of all values in \c kl */
    int max_kl;

    /** Maximum of all values in \c ku */
    int max_ku;

    /** Leading dimension for all derivative operators */
    int ld;

    /**
     * Storage for each banded derivative operator matrix.
     *
     * This general band storage can be accessed in two different ways:
     * \li <tt>D_T[k]</tt> is the band storage for the transpose of the
     *     <tt>k</tt>-th derivative using <tt>kl[k]</tt>, <tt>ku[k]</tt>,
     *     and <tt>ld</tt>.  This view is optimal for BLAS \c gbmv operations
     *     using <tt>TRANS == 'T'</tt>.
     * \li <tt>D_T[k] - (max_ku - ku[k])</tt> is the same band storage
     *     for the transpose of the <tt>k</tt>-th derivative viewed using
     *     <tt>max_kl</tt>, <tt>max_ku</tt>, and <tt>ld</tt>.  This view is
     *     optimal for BLAS \c gb_acc and \c gb_add operations.
     */
    double **D_T;

    /* @} */

} suzerain_bsplineop_workspace;

/** @name Allocation and deallocation */
/**@{*/

/**
 * Allocate a B-spline operator workspace.
 *
 * @param[in] bw  A \c gsl_bspline_workspace instance which fixes the B-spline
 *                basis function order and breakpoints.
 * @param[in] dbw A \c gsl_bspline_deriv_workspace with an order
 *                at least that of \c bw.
 * @param[in] nderiv Highest derivative operator requested.  The zeroth
 *                   through \c nderiv-th operator will be available in
 *                   banded matrix form.
 * @param[in] method Specifies the method used to compute derivative operators.
 *
 * @return a workspace instance with operators ready for use on success.
 *      On failure calls suzerain_error() and returns NULL.
 *
 * \memberof suzerain_bsplineop_workspace
 */
suzerain_bsplineop_workspace *
suzerain_bsplineop_alloc(
    gsl_bspline_workspace *bw,
    gsl_bspline_deriv_workspace *dbw,
    int nderiv,
    enum suzerain_bsplineop_method method);

/**
 * Free a previously allocated workspace.
 *
 * @param[in] w Workspace to free.
 *
 * \memberof suzerain_bsplineop_workspace
 */
void
suzerain_bsplineop_free(
    suzerain_bsplineop_workspace *w);

/**@}*/

/** @name Real-valued operations */
/**@{*/

/**
 * Apply the <tt>nderiv</tt>-th derivative operator to real coefficients \c x.
 * Multiplies the precomputed banded derivative operator scaled by \c alpha
 * against one or more coefficient vectors stored in \c x.  Results overwrite
 * \c x.  Each coefficient vector is of length <tt>w->n</tt>.
 *
 * @param[in] nderiv Derivative operator to apply.  May be zero.
 * @param[in] nrhs Number of coefficient vectors stored in \c x.
 * @param[in] alpha Real scaling factor to apply.
 * @param[in,out] x Coefficients to be multiplied.
 * @param[in] incx Stride between elements stored in \c x.
 * @param[in] ldx Leading dimension of the data stored in \c x.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bsplineop_apply_complex() for a way to apply an operator
 *      to complex-valued coefficients.
 * @see suzerain_bsplineop_accumulate() for a way to accumulate the
 *      effect of an operator against other storage.
 *
 * \memberof suzerain_bsplineop_workspace
 */
int
suzerain_bsplineop_apply(
    int nderiv,
    int nrhs,
    double alpha,
    double *x,
    int incx,
    int ldx,
    const suzerain_bsplineop_workspace *w);

/**
 * Apply the <tt>nderiv</tt>-th derivative operator to real coefficients \c x
 * accumulating the results in \c y.  Multiplies \c alpha times the precomputed
 * banded derivative operator against one or more coefficient vectors stored in
 * \c x.  Results are added to \c beta times \c y.  Each coefficient vector is
 * of length <tt>w->n</tt>.  \c x and \c y cannot be aliased.
 *
 * @param[in] nderiv Derivative operator to apply.  May be zero.
 * @param[in] nrhs Number of vectors stored in \c x and \c y.
 * @param[in] alpha Multiplicative factor to use on \c x.
 * @param[in] x Coefficients to be multiplied.
 * @param[in] incx Stride between elements stored in \c x.
 * @param[in] ldx Leading dimension of the data stored in \c x.
 * @param[in] beta Multiplicative factor to use on \c y.
 * @param[in,out] y Locations in which to accumulate the result.
 * @param[in] incy Stride between elements stored in \c y.
 * @param[in] ldy Leading dimension of the data stored in \c y.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bsplineop_accumulate_complex() for a way to apply an operator
 *      to complex-valued coefficients.
 * @see suzerain_bsplineop_apply() for a way to apply an operator in
 *      place.
 *
 * \memberof suzerain_bsplineop_workspace
 */
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
    const suzerain_bsplineop_workspace *w);

/**
 * Determine the right hand side of the interpolation problem <tt>D[0] x =
 * rhs</tt>.  Here <tt>D[0]</tt> is the zeroth derivative operator (i.e. mass
 * matrix), \c x are the basis function coefficients that will best represent
 * \c function for the given method, and \c rhs is the vector computed by this
 * routine.
 *
 * @param[in]  function Real-valued function to use when computing \c rhs.
 * @param[out] rhs Output vector containing computed right hand side
 * @param[in]  bw Workspace to use which must have been the one used
 *             to create \c w.
 * @param[in]  w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see suzerain_bsplineop_method() for details on available methods.
 *      Different methods will necessarily compute the right hand side
 *      differently.
 * @see suzerain_bsplineop_lu_opaccumulate(), suzerain_bsplineop_lu_factor(),
 *      and suzerain_bsplineop_lu_solve() for how to solve the linear
 *      equation for \c x.
 * @see suzerain_bsplineop_interpolation_rhs_complex() for a way to
 *      perform this operation for a complex-valued function.
 *
 * \memberof suzerain_bsplineop_workspace
 */
int
suzerain_bsplineop_interpolation_rhs(
    const suzerain_function * function,
    double * rhs,
    gsl_bspline_workspace *bw,
    const suzerain_bsplineop_workspace *w);

/**@}*/

/** @name Complex-valued operations */
/**@{*/

/**
 * Apply the <tt>nderiv</tt>-th derivative operator to complex coefficients \c
 * x.  Multiplies the precomputed banded derivative operator scaled by
 * real-valued \c alpha against one or more coefficient vectors stored in \c x.
 * Results overwrite \c x.  Each coefficient vector is of length <tt>w->n</tt>.
 *
 * @param[in] nderiv Derivative operator to apply.  May be zero.
 * @param[in] nrhs Number of coefficient vectors stored in \c x.
 * @param[in] alpha Real scaling factor to apply.
 * @param[in,out] x Coefficients to be multiplied.
 * @param[in] incx Stride between elements stored in \c x.
 * @param[in] ldx Leading dimension of the data stored in \c x.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bsplineop_apply() for a way to apply an operator
 *      to real-valued coefficients.
 * @see suzerain_bsplineop_accumulate() for a way to accumulate the
 *      effect of an operator against other storage.
 *
 * \memberof suzerain_bsplineop_workspace
 */
int
suzerain_bsplineop_apply_complex(
    int nderiv,
    int nrhs,
    double alpha,
    complex_double *x,
    int incx,
    int ldx,
    const suzerain_bsplineop_workspace *w);

/**
 * Apply the <tt>nderiv</tt>-th derivative operator to complex coefficients \c
 * x accumulating the complex result in \c y.  Multiplies \c alpha times the
 * precomputed banded derivative operator against one or more coefficient
 * vectors stored in \c x.  Results are added to \c beta times \c y.  Each
 * coefficient vector is of length <tt>w->n</tt>.  \c x and \c y cannot be
 * aliased.
 *
 * @param[in] nderiv Derivative operator to apply.  May be zero.
 * @param[in] nrhs Number of vectors stored in \c x and \c y.
 * @param[in] alpha Multiplicative factor to use on \c x.
 * @param[in] x Coefficients to be multiplied.
 * @param[in] incx Stride between elements stored in \c x.
 * @param[in] ldx Leading dimension of the data stored in \c x.
 * @param[in] beta Multiplicative factor to use on \c y.
 * @param[in,out] y Locations in which to accumulate the result.
 * @param[in] incy Stride between elements stored in \c y.
 * @param[in] ldy Leading dimension of the data stored in \c y.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bsplineop_accumulate() for a way to apply an operator
 *      to real-valued coefficients.
 * @see suzerain_bsplineop_apply() for a way to apply an operator in
 *      place.
 *
 * \memberof suzerain_bsplineop_workspace
 */
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
    const suzerain_bsplineop_workspace *w);

/**
 * Determine the right hand side of the interpolation problem <tt>D[0] x
 * = rhs</tt>.  Here <tt>D[0]</tt> is the zeroth derivative operator (i.e. mass
 * matrix), \c x are the basis function coefficients that will best represent
 * \c function for the given method, and \c rhs is the vector computed by this
 * routine.
 *
 * @param[in]  zfunction Complex-valued function to use when computing \c rhs.
 * @param[out] rhs Output vector containing computed right hand side
 * @param[in]  bw Workspace to use which must have been the one used
 *             to create \c w.
 * @param[in]  w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see suzerain_bsplineop_method() for details on available methods.
 *      Different methods will necessarily compute the right hand side
 *      differently.
 * @see suzerain_bsplineop_luz_accumulate(), suzerain_bsplineop_luz_factor(),
 *      and suzerain_bsplineop_luz_solve() for how to solve the linear
 *      equation for \c x.
 * @see suzerain_bsplineop_interpolation_rhs() for a way to
 *      perform this operation for a real-valued function.
 *
 * \memberof suzerain_bsplineop_workspace
 */
int
suzerain_bsplineop_interpolation_rhs_complex(
    const suzerain_zfunction * zfunction,
    complex_double *rhs,
    gsl_bspline_workspace *bw,
    const suzerain_bsplineop_workspace *w);

/**@}*/

/**
 * Encapsulates an operator and a LU decomposition from a real-valued linear
 * combination of derivative operators.  Callers obtain a workspace using
 * suzerain_bsplineop_lu_alloc() and release it using
 * suzerain_bsplineop_lu_free().
 *
 * As manipulating an unfactored operator is easiest using the unfactored \c kl
 * and \c ku, and as LAPACK's \c GBTRF and \c GBTRS are written using the
 * unfactored values of \c kl and \c ku, we do too.  The unfactored \e
 * transpose of the operator uses column-major band storage per
 * <tt>GBTRF</tt>'s \c AB "on entry" documentation.  After factorization, the
 * factorization of the transposed operator is stored per <tt>GBTRS</tt>'s "on
 * entry" documentation.
 *
 * @see suzerain_bsplineop_luz_workspace for manipulating the complex-valued
 * case.
 */
typedef struct suzerain_bsplineop_lu_workspace {

    /** Number of degrees of freedom in the basis */
    int n;

    /** Number of subdiagonals in the transpose of the unfactored operator. */
    int kl;

    /**
     * Number of superdiagonals in the transpose of the unfactored operator.
     * The number of superdiagonals in the upper triangular portion of the
     * factored operator transpose is <tt>kl + ku</tt>.
     */
    int ku;

    /** Leading dimension of the unfactored and factored operator storage */
    int ld;

    /**
     *  Pivot matrix \c P from the \c LUP decomposition of the transposed
     *  operator.
     */
    int *ipiv;

    /**
     * Raw data storage for the unfactored and factored operator transposes.
     * The transpose of the operator is stored according to LAPACK \c GBTRF
     * conventions regarding \c kl, \c ku, \c ld, and \c n.  In particular,
     * this means <tt>A + kl</tt> is the starting general band storage location
     * for the transpose of the unfactored operator.
     */
    double *A_T;

} suzerain_bsplineop_lu_workspace;

/** @name Allocation and deallocation */
/**@{*/

/**
 * Allocate a B-spline operator real-valued LU workspace.
 *
 * @param[in] w B-spline workspace on which to base the new workspace's
 *      dimensions.
 *
 * @return a workspace instance ready to be used with
 *      suzerain_bsplineop_lu_opform() or suzerain_bsplineop_lu_opform_mass().
 *      On failure calls suzerain_error() and returns NULL.
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
suzerain_bsplineop_lu_workspace *
suzerain_bsplineop_lu_alloc(
    const suzerain_bsplineop_workspace *w);

/**
 * Free a previously allocated workspace.
 *
 * @param[in] luw Workspace to free.
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
void
suzerain_bsplineop_lu_free(
    suzerain_bsplineop_lu_workspace * luw);

/**@}*/

/** @name Operations */
/**@{*/

/**
 * Build an operator by accumulating a linear combination of derivative
 * operators from the workspace \c w.  That is, \f[
 *  \mbox{luw} \leftarrow
 *  \sum_{j=0}^{\mbox{ncoefficients}} \mbox{coefficients}_{j} \, D[j]
 *  + \mbox{scale} \mbox{luw}
 * \f].
 *
 * @param[in] ncoefficients Number of coefficients stored in \c coefficients.
 * @param[in] coefficients Coefficients to use in the linear combination of
 *      derivative operators.
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in] scale Scale factor used on the current operator content.
 * @param[in,out] luw Workspace in which to accumulate the operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_opaccumulate(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bsplineop_workspace * w,
    const double scale,
    suzerain_bsplineop_lu_workspace * luw);

/**
 * Compute the \f$p = \infty\f$ norm of the unfactored operator stored in
 * \c luw.  That is, find \f$ \left|\left|A\right|\right|_{\infty} \f$.
 *
 * \param[in]  luw  Workspace containing the operator to investigate.
 * \param[out] norm The infinity norm of the operator.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_opnorm(
    const suzerain_bsplineop_lu_workspace * luw,
    double * norm);

/**
 * Factorize the operator in the current workspace.  That is, \f[ \mbox{luw}
 * \leftarrow \mbox{LU}\left(\mbox{luw}\right) \f].
 *
 * @param[in,out] luw Workspace containing the desired operator and
 *      in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_factor(
    suzerain_bsplineop_lu_workspace * luw);

/**
 * Estimate the reciprocal condition number of the factored operator stored in
 * \c luw.  That is, \f$ \left( \left|\left|A\right|\right|_{\infty}
 * \left|\left|A^{-1}\right|\right|_{\infty} \right)^{-1} \f$.
 *
 * \param[in]  norm  The infinity norm of the unfactored operator
 *                   obtained from suzerain_bsplineop_lu_opnorm().
 * \param[out] rcond The reciprocal of the condition number of the operator.
 * \param[in]  luw   Workspace containing the factored operator to investigate.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_rcond(
    const double norm,
    double * rcond,
    const suzerain_bsplineop_lu_workspace * luw);

/**
 * Solves the equation <tt>A X = B</tt> using the factored operator
 * stored in \c luw.
 *
 * @param[in] nrhs Number of right hand sides, i.e. columns, stored in
 *      matrix \c B.
 * @param[in,out] B Matrix of right hand sides to solve and the resulting
 *      solutions.  On input, contains data \c B.  On output, contains
 *      solutions \c X.
 * @param[in] incb Stride between elements in matrix \c B.
 * @param[in] ldb Leading dimension of matrix \c B.
 * @param[in] luw Workspace containing the factored operator to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_solve(
    int nrhs,
    double *B,
    int incb,
    int ldb,
    const suzerain_bsplineop_lu_workspace * luw);

/**@}*/

/** @name Convenience operations */
/**@{*/

/**
 * Form an operator from a linear combination of derivative operators
 * from the workspace \c w.  That is, \f[
 *  \mbox{luw} \leftarrow
 *  \sum_{j=0}^{\mbox{ncoefficients}} \mbox{coefficients}_{j} \, D[j]
 * \f].
 *
 * @param[in] ncoefficients Number of coefficients stored in \c coefficients.
 * @param[in] coefficients Coefficients to use in the linear combination of
 *      derivative operators.
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in,out] luw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bsplineop_lu_opaccumulate().
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_opform(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace * luw);

/**
 * Form the zeroth derivative operator from the workspace \c w.
 * That is, \f[ \mbox{luw} \leftarrow D[0] \f].
 *
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in,out] luw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bsplineop_lu_opaccumulate().
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_opform_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace * luw);

/**
 * Form and factorize the zeroth derivative operator from the workspace \c w.
 * That is, \f[ \mbox{luw} \leftarrow \mbox{LU}\left(D[0]\right) \f].
 *
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in,out] luw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bsplineop_lu_opaccumulate() and suzerain_bsplineop_lu_factor().
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_lu_factor_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_lu_workspace * luw);

/**@}*/

/**
 * Encapsulates an operator and a LU decomposition from a complex-valued linear
 * combination of derivative operators.  Callers obtain a workspace using
 * suzerain_bsplineop_luz_alloc() and release it using
 * suzerain_bsplineop_luz_free().
 *
 * As manipulating an unfactored operator is easiest using the unfactored \c kl
 * and \c ku, and as LAPACK's \c GBTRF and \c GBTRS are written using the
 * unfactored values of \c kl and \c ku, we do too.  The unfactored \e
 * transpose of the operator uses column-major band storage per
 * <tt>GBTRF</tt>'s \c AB "on entry" documentation.  After factorization, the
 * factorization of the transposed operator is stored per <tt>GBTRS</tt>'s "on
 * entry" documentation.
 *
 * @see suzerain_bsplineop_lu_workspace for manipulating the real-valued case.
 */
typedef struct suzerain_bsplineop_luz_workspace {

    /** Number of degrees of freedom in the basis */
    int n;

    /** Number of subdiagonals in the transpose of the unfactored operator. */
    int kl;

    /**
     * Number of superdiagonals in the transpose of the unfactored operator.
     * The number of superdiagonals in the upper triangular portion of the
     * factored operator transpose is <tt>kl + ku</tt>.
     */
    int ku;

    /** Leading dimension of the unfactored and factored operator storage */
    int ld;

    /** Pivot matrix \c P from the \c LUP decomposition of the operator. */
    int *ipiv;

    /**
     * Raw data storage for the unfactored and factored operator transposes.
     * The transpose of the operator is stored according to LAPACK \c GBTRF
     * conventions regarding \c kl, \c ku, \c ld, and \c n.  In particular,
     * this means <tt>A + kl</tt> is the starting general band storage location
     * for the transpose of the unfactored operator.
     */
    complex_double *A_T;

} suzerain_bsplineop_luz_workspace;

/** @name Allocation and deallocation */
/**@{*/

/**
 * Allocate a B-spline operator LU complex-valued workspace.
 *
 * @param[in] w B-spline workspace on which to base the new workspace's
 *      dimensions.
 *
 * @return a workspace instance ready to be used with
 *      suzerain_bsplineop_luz_opform() or suzerain_bsplineop_luz_opform_mass().
 *      On failure calls suzerain_error() and returns NULL.
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
suzerain_bsplineop_luz_workspace *
suzerain_bsplineop_luz_alloc(
    const suzerain_bsplineop_workspace *w);

/**
 * Free a previously allocated workspace.
 *
 * @param[in] luzw Workspace to free.
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
void
suzerain_bsplineop_luz_free(
    suzerain_bsplineop_luz_workspace * luzw);

/**@}*/

/** @name Operations */
/**@{*/

/**
 * Build an operator by accumulating a linear combination of derivative
 * operators from the workspace \c w.  That is, \f[
 *  \mbox{luzw} \leftarrow
 *  \sum_{j=0}^{\mbox{ncoefficients}} \mbox{coefficients}_{j} \, D[j]
 *  + \mbox{scale} \mbox{luzw}
 * \f].
 *
 * @param[in] ncoefficients Number of coefficients stored in \c coefficients.
 * @param[in] coefficients Coefficients to use in the linear combination of
 *      derivative operators.
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in] scale Scale factor used on the current operator content.
 * @param[in,out] luzw Workspace in which to accumulate the operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
int
suzerain_bsplineop_luz_opaccumulate(
    int ncoefficients,
    const complex_double *coefficients,
    const suzerain_bsplineop_workspace * w,
    const complex_double scale,
    suzerain_bsplineop_luz_workspace * luzw);

/**
 * Compute the \f$p = \infty\f$ norm of the unfactored operator stored in
 * \c luzw.  That is, find \f$ \left|\left|A\right|\right|_{\infty} \f$.
 *
 * \param[in]  luzw Workspace containing the operator to investigate.
 * \param[out] norm The infinity norm of the operator.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
int
suzerain_bsplineop_luz_opnorm(
    const suzerain_bsplineop_luz_workspace * luzw,
    double * norm);

/**
 * Factorize the operator in the current workspace.  That is, \f[ \mbox{luzw}
 * \leftarrow \mbox{LU}\left(\mbox{luzw}\right) \f].
 *
 * @param[in,out] luzw Workspace containing the desired operator and
 *      in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
int
suzerain_bsplineop_luz_factor(
    suzerain_bsplineop_luz_workspace * luzw);

/**
 * Estimate the reciprocal condition number of the factored operator stored in
 * \c luzw.  That is, \f$ \left( \left|\left|A\right|\right|_{\infty}
 * \left|\left|A^{-1}\right|\right|_{\infty} \right)^{-1} \f$.
 *
 * \param[in]  norm  The infinity norm of the unfactored operator
 *                   obtained from suzerain_bsplineop_luz_norm().
 * \param[out] rcond The reciprocal of the condition number of the operator.
 * \param[in]  luzw Workspace containing the factored operator to investigate.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
int
suzerain_bsplineop_luz_rcond(
    const double norm,
    double * rcond,
    const suzerain_bsplineop_luz_workspace * luzw);

/**
 * Solves the equation <tt>A X = B</tt> using the factored operator
 * stored in \c luzw.
 *
 * @param[in] nrhs Number of right hand sides, i.e. columns, stored in
 *      matrix \c B.
 * @param[in,out] B Matrix of right hand sides to solve and the resulting
 *      solutions.  On input, contains data \c B.  On output, contains
 *      solutions \c X.
 * @param[in] incb Stride between elements in matrix \c B.
 * @param[in] ldb  Leading dimension of matrix \c B.
 * @param[in] luzw Workspace containing the factored operator to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
int
suzerain_bsplineop_luz_solve(
    int nrhs,
    complex_double *B,
    int incb,
    int ldb,
    const suzerain_bsplineop_luz_workspace * luzw);

/**@}*/

/** @name Convenience operations */
/**@{*/

/**
 * Form an operator from a linear combination of derivative operators
 * from the workspace \c w.  That is, \f[
 *  \mbox{luzw} \leftarrow
 *  \sum_{j=0}^{\mbox{ncoefficients}} \mbox{coefficients}_{j} \, D[j]
 * \f].
 *
 * @param[in] ncoefficients Number of coefficients stored in \c coefficients.
 * @param[in] coefficients Coefficients to use in the linear combination of
 *      derivative operators.
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in,out] luzw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bsplineop_luz_opaccumulate().
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
int
suzerain_bsplineop_luz_opform(
    int ncoefficients,
    const complex_double *coefficients,
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace * luzw);

/**
 * Form the zeroth derivative operator from the workspace \c w.
 * That is, \f[ \mbox{luzw} \leftarrow D[0] \f].
 *
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in,out] luzw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bsplineop_luz_opaccumulate().
 *
 * \memberof suzerain_bsplineop_luz_workspace
 */
int
suzerain_bsplineop_luz_opform_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace * luzw);

/**
 * Form and factorize the zeroth derivative operator from the workspace \c w.
 * That is, \f[ \mbox{luw} \leftarrow \mbox{LU}\left(D[0]\right) \f].
 *
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in,out] luzw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bsplineop_luz_opaccumulate()
 *      and suzerain_bsplineop_luz_factor().
 *
 * \memberof suzerain_bsplineop_lu_workspace
 */
int
suzerain_bsplineop_luz_factor_mass(
    const suzerain_bsplineop_workspace * w,
    suzerain_bsplineop_luz_workspace * luzw);

/**@}*/

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // SUZERAIN_BSPLINEOP_H
