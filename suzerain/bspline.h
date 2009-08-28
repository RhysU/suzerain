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
 * bspline.h: B-spline basis manipulation and operator routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef PECOS_SUZERAIN_BSPLINE_OPERATOR_H
#define PECOS_SUZERAIN_BSPLINE_OPERATOR_H

#include <suzerain/function.h>

/** @file
 * Provides B-spline basis evaluation and derivative operator routines.  These
 * functions provide B-spline basis function evaluation, differentiation,
 * derivative linear operator construction, and operator application routines.
 * Unless otherwise noted, all matrices are stored in column-major (Fortran)
 * storage.  Multiple threads may call the routines simultaneously provided
 * that each thread has its own workspace instances.  The behavior is undefined
 * if two threads simultaneously share a workspace.  Internally, the code
 * builds upon the <a href="http://www.gnu.org/software/gsl/">GNU Scientific
 * Library</a> (GSL) <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Basis-Splines.html">
 * B-spline routines</a> and uses banded matrix BLAS and LAPACK functionality
 * where possible.
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
__BEGIN_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/* Conditional forward declarations of implementation-related structs */
#ifndef __GSL_MATRIX_DOUBLE_H__
/** Forward declaration of GSL's \c gsl_matrix. */
typedef struct gsl_matrix gsl_matrix;
#endif
#ifndef __GSL_BSPLINE_H__
/** Forward declaration of GSL's \c gsl_bspline_workspace. */
typedef struct gsl_bspline_workspace gsl_bspline_workspace;
/** Forward declaration of GSL's \c gsl_bspline_deriv_workspace. */
typedef struct gsl_bspline_deriv_workspace gsl_bspline_deriv_workspace;
#endif

/** Indicates the method chosen to compute derivative operators.  */
enum suzerain_bspline_method {

    /**
     * Form derivative operators using collocation at the Greville abscissae.
     **/
    SUZERAIN_BSPLINE_COLLOCATION_GREVILLE = 1
};

/**
 * Encapsulates basis function and derivative operator information.
 * Callers obtain a workspace using suzerain_bspline_alloc() and
 * release it using suzerain_bspline_free().
 */
typedef struct suzerain_bspline_workspace {

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

    /** Number of degrees of freedom in the basis */
    int ndof;

    /** Method chosen to form derivative operators */
    enum suzerain_bspline_method method;

    /** @name Banded derivative matrix storage details
     * Each of the \c 0 through \c nderivative (inclusive) derivative operators
     * is stored using BLAS banded matrix conventions.
     * @{
     */

    /** Number of subdiagonals in each derivative operator */
    int * kl;

    /** Number of superdiagonals in each derivative operator */
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
     * \li <tt>D[k]</tt> is the band storage for the <tt>k</tt>-th derivative
     *     using <tt>kl[k]</tt>, <tt>ku[k]</tt>, and <tt>ld</tt>.  This view
     *     is optimal for BLAS \c gbmv operations.
     * \li <tt>D[k] - (max_ku - ku[k])</tt> is the same band storage
     *     for the <tt>k</tt>-th derivative viewed using <tt>max_kl</tt>,
     *     <tt>max_ku</tt>, and <tt>ld</tt>.  This view is optimal for
     *     BLAS \c gb_acc and \c gb_add operations.
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

/**
 * Allocate a B-spline workspace.
 *
 * @param[in] order Bspline order desired.  This is the piecewise degree plus
 *          one.  For example, 4 represents piecewise cubics.
 * @param[in] nderivatives Highest derivative operator requested.  The zeroth
 *          through \c nderivatives-th operator will be available in
 *          banded matrix form.
 * @param[in] nbreakpoints Number of breakpoints to use.  This can be thought
 *          of as the number of piecewise intervals plus one.
 * @param[in] breakpoints Breakpoint locations in a length \c nbreakpoints
 *          array.  Breakpoint locations may be nonuniform but must be
 *          strictly increasing.  The values are copied for internal usage.
 * @param[in] method Specifies the method used to compute derivative operators.
 *
 * @return a workspace instance with operators ready for use on success.
 *      On failure calls suzerain_error() and returns NULL.
 */
suzerain_bspline_workspace *
suzerain_bspline_alloc(
    int order,
    int nderivatives,
    int nbreakpoints,
    const double * breakpoints,
    enum suzerain_bspline_method method);

/**
 * Frees a previously allocated workspace.
 *
 * @param[in] w Workspace to free.
 */
void
suzerain_bspline_free(
    suzerain_bspline_workspace *w);

/**
 * Apply the <tt>nderivative</tt>-th derivative operator to coefficients \c b.
 * Multiplies the precomputed banded derivative operator against one or more
 * coefficient vectors stored in \c b.  Results overwrite \c b.  Each
 * coefficient vector is of length suzerain_bspline_ndof().
 *
 * @param[in] nderivative Derivative operator to apply.  May be zero.
 * @param[in] nrhs Number of coefficient vectors stored in b.
 * @param[in,out] b Coefficients to be multiplied.
 * @param[in] ldb Leading dimension of the data stored in \c b.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_apply_operator(
    int nderivative,
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_workspace *w);

/**
 * Evaluate a function and its derivatives based upon a linear combination of
 * the basis functions.  At each \c point in \c points, evaluate the function
 * and its derivatives determined from a linear combination of the supplied \c
 * coefficients times the B-spline basis functions at \c point.  Derivatives \c
 * 0 through \c nderivative (inclusive) are computed and stored as columns
 * within \c values.
 * That is, \f[
 *  \mbox{values}\left[i + k\,\mbox{ldvalues}\right] = \frac{d^k}{dx^k}
 *  \sum_{j=0}^{\mbox{ndof}} \mbox{coefficients}_{j}\,B_{j}(\mbox{points}_{i})
 * \f]
 * for \f$ i\in\left\{0,\ldots,\mbox{npoints}\right\} \f$,
 *\f$ k\in\left\{0,\ldots,\mbox{nderivative}\right\} \f$,
 * \f$\mbox{ndof} = \f$ suzerain_bspline_ndof(w), and B-spline basis
 * functions \f$ B_{j} \f$.  If only a single derivative is desired, passing \c
 * 0 to \c ldvalues will cause only that single derivative to be written in the
 * first column of \c values.
 *
 * \note It is more efficient to compute a function and its derivatives
 * simultaneously than to request each derivative separately.  This is due to
 * the recurrence relationship used to compute B-spline derivatives.
 *
 * @param[in] nderivative Maximum requested derivative.  This may be higher
 *      than the number of derivatives requested in suzerain_bspline_alloc().
 * @param[in] coefficients Expansion coefficients for a function in terms
 *      of the B-spline basis.  Must be of length suzerain_bspline_ndof().
 * @param[in] npoints Number of evaluation points.
 * @param[in] points Points at which to evaluate the function.
 * @param[out] values Matrix of values resulting from evaluating the function
 *      and its derivatives.  Matrix dimensions are suzerain_bspline_ndof()
 *      by \c nderivative if \c ldvalues >= \c suzerain_bspline_ndof().
 *      If \c ldvalues is zero, only a single column is returned in \c values.
 * @param[in] ldvalues Leading dimension of the output matrix \c values.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_evaluate(
    int nderivative,
    const double * coefficients,
    int npoints,
    const double * points,
    double * values,
    int ldvalues,
    const suzerain_bspline_workspace *w);

/**
 * Determine the right hand side of the interpolation problem <tt>D[0] x
 * =rhs</tt>.  Here <tt>D[0]</tt> is the zeroth derivative operator (i.e. mass
 * matrix), \c x are the basis function coefficients that will best represent
 * \c function for the given method, and \c rhs is the vector computed by this
 * routine.
 *
 * @param[in] function Function to use when computing \c rhs.
 * @param[out] rhs Output vector containing computed right hand side
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see suzerain_bspline_method() for details on available methods.
 *      Different methods will necessarily compute the right hand side
 *      differently.
 * @see suzerain_bspline_lu_form_mass() and suzerain_bspline_lu_solve()
 *      for how to solve the linear equation for \c x.
 */
int
suzerain_bspline_find_interpolation_problem_rhs(
    const suzerain_function * function,
    double * rhs,
    const suzerain_bspline_workspace *w);

/**
 * Encapsulates the LU decomposition of a linear combination of derivative
 * operators.  Callers obtain a workspace using suzerain_bspline_lu_alloc() and
 * release it using suzerain_bspline_lu_free().
 */
typedef struct suzerain_bspline_lu_workspace {

    /** Number of degrees of freedom in the basis */
    int ndof;

    /** Number of subdiagonals in the factored operator. */
    int kl;

    /**
     * Number of superdiagonals in the factored operator.
     * Note that, for a given B-spline basis, suzerain_bspline_lu_workspace::ku
     * is larger than the corresponding suzerain_bspline_workspace::ku
     * according to the requirements of LAPACK's \c dgbtrf.
     */
    int ku;

    /** Leading dimension in the factored operator derivative */
    int ld;

    /** Pivot matrix \c P from the \c LUP decomposition of the operator. */
    int *ipiv;

    /**
     * Raw data storage for the factored operator.  Operator is stored
     * according to BLAS banded matrix conventions according to \c kl, \c ku,
     * \c ld, and \c ndof
     **/
    double *A;

} suzerain_bspline_lu_workspace;


/**
 * Allocates a B-spline operator LU workspace.
 *
 * @param[in] w B-spline workspace on which to base the new workspace's
 *      dimensions.
 *
 * @return a workspace instance ready to be used with
 *      suzerain_bspline_lu_form_general() or suzerain_bspline_lu_form_mass().
 *      On failure calls suzerain_error() and returns NULL.
 */
suzerain_bspline_lu_workspace *
suzerain_bspline_lu_alloc(
    const suzerain_bspline_workspace *w);

/**
 * Frees a previously allocated workspace.
 *
 * @param[in] luw Workspace to free.
 */
void
suzerain_bspline_lu_free(
    suzerain_bspline_lu_workspace *luw);

/**
 * Forms the LU decomposition of a linear combination of derivative operators
 * from the workspace \c w.  That is, \f[
 *  \mbox{luw} \leftarrow \mbox{LU}\left(
 *      \sum_{j=0}^{\mbox{ncoefficients}} \mbox{coefficients}_{j} \, D[j]
 *  \right)
 * \f]
 *
 * @param[in] ncoefficients Number of coefficients stored in \c coefficients.
 * @param[in] coefficients Coefficients to use in the linear combination of
 *      derivative operators.
 * @param[in] w Workspace containing desired derivative operators.
 * @param[in,out] luw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_lu_form_general(
    int ncoefficients,
    const double * coefficients,
    const suzerain_bspline_workspace * w,
    suzerain_bspline_lu_workspace *luw);

/**
 * Forms the LU decomposition of the zeroth derivative operator, also
 * called the mass matrix, from workspace \c w.
 *
 * @param[in] w Workspace containing the desired mass matrix.
 * @param[in,out] luw  Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bspline_lu_form_general().
 */
int
suzerain_bspline_lu_form_mass(
    const suzerain_bspline_workspace * w,
    suzerain_bspline_lu_workspace *luw);

/**
 * Solves the equation <tt>A x = b</tt> using the factored operator
 * stored in \c luw.
 *
 * @param[in] nrhs Number of right hand sides, i.e. columns, stored in
 *      matrix \c b.
 * @param[in,out] b Matrix of right hand sides to solve and the resulting
 *      solutions.  On input, contains data \c b.  On output, contains
 *      solutions \c x.
 * @param[in] ldb Leading dimension of matrix b.
 * @param[in] luw Workspace containing the factored operator to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_lu_solve(
    int nrhs,
    double *b,
    int ldb,
    const suzerain_bspline_lu_workspace *luw);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // PECOS_SUZERAIN_BSPLINE_OPERATOR_H
