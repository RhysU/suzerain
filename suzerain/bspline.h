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
#ifndef __SUZERAIN_BSPLINE_OPERATOR_H
#define __SUZERAIN_BSPLINE_OPERATOR_H

#include <suzerain/common.h>
#include <suzerain/function.h>

/** @file
 * Provides B-spline basis evaluation and derivative operator routines.  These
 * functions provide B-spline basis function evaluation, differentiation,
 * derivative operator construction, and operator application routines.
 * Both real- and complex- valued operator application and inversion
 * are supported.
 *
 * The derivative operators map a function's spline coefficients to an
 * approximation of the derivative's spline coefficients.  For
 * collocation-based operators, the approximation \f$\tilde{\gamma}(x) =
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
 * through the ::suzerain_bspline_method value provided to
 * suzerain_bspline_alloc().
 *
 * Assuming you already have a suzerain_function instance \c f, this snippet
 * interpolates a general function onto a cubic spline basis and then takes a
 * derivative:
 * \code
 *  // Declare the breakpoints used in our basis
 *  const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
 *  const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
 *
 *  // Create a cubic (k=4) workspace with 1 derivative available
 *  suzerain_bspline_workspace *w
 *      = suzerain_bspline_alloc(4, 1, nbreak, breakpoints,
 *          SUZERAIN_BSPLINE_COLLOCATION_GREVILLE);
 *
 *  // Factor the mass matrix D[0] to allow solving D[0] x = b
 *  suzerain_bspline_lu_workspace *mass
 *      = suzerain_bspline_lu_alloc(w);
 *  suzerain_bspline_lu_form_mass(w, mass);
 *
 *  // Set up and solve for coefficients that interpolate f
 *  double x[w->ndof]; // C99
 *  suzerain_bspline_find_interpolation_problem_rhs(&f, x, w);
 *  suzerain_bspline_lu_solve(1, x, 1, mass);
 *
 *  // Solve D[0] x' = D[1] x by forming right hand side and solving
 *  suzerain_bspline_apply_operator(1, 1, 1, x, 1, w->ndof, w);
 *  suzerain_bspline_lu_solve(1, x, 1, mass);
 *  // x now contains an approximation to the derivative of f
 *
 *  // Deallocate workspaces
 *  suzerain_bspline_lu_free(luw);
 *  suzerain_bspline_free(w);
 * \endcode
 *
 * All matrices are stored in column-major (Fortran) storage.  Multiple threads
 * may call the routines simultaneously provided that each thread has its own
 * workspace instances.  The behavior is undefined if two threads
 * simultaneously share a workspace.  Internally, the code builds upon the <a
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
     * Spline order per GSL/PPPACK conventions.  For example, 3 denotes
     * piecewise quadratics and 4 denotes piecewise cubics.
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
 * Obtain the number of degrees of freedom for the basis in the given workspace.
 * This information is also available via <tt>w->ndof</tt>.
 *
 * @param[in] w Workspace to use.
 *
 * @return the number of degrees of freedom.
 */
inline
int
suzerain_bspline_ndof(
    const suzerain_bspline_workspace *w)
{
    return w->ndof;
}

/**
 * Apply the <tt>nderivative</tt>-th derivative operator to real coefficients
 * \c x accumulating the results in \c y.  Multiplies \c alpha times the
 * precomputed banded derivative operator against one or more coefficient
 * vectors stored in \c x.  Results are added to \c beta times \c y.  Each
 * coefficient vector is of length suzerain_bspline_ndof().  \c x and \c y
 * cannot be aliased.  Increments and leading dimensions are specified in
 * <tt>double</tt>-sized units.
 *
 * @param[in] nderivative Derivative operator to apply.  May be zero.
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
 * @see suzerain_bspline_zaccumulate_operator() for a way to apply an operator
 *      to complex-valued coefficients.
 * @see suzerain_bspline_apply_operator() for a way to apply an operator in
 *      place.
 */
int
suzerain_bspline_accumulate_operator(
    int nderivative,
    int nrhs,
    double alpha,
    const double *x,
    int incx,
    int ldx,
    double beta,
    double *y,
    int incy,
    int ldy,
    const suzerain_bspline_workspace *w);

/**
 * Apply the <tt>nderivative</tt>-th derivative operator to complex
 * coefficients \c x accumulating the complex result in \c y.  Multiplies \c
 * alpha times the precomputed banded derivative operator against one or more
 * coefficient vectors stored in \c x.  Results are added to \c beta times \c
 * y.  Each coefficient vector is of length suzerain_bspline_ndof().  \c x and
 * \c y cannot be aliased.  Increments and leading dimensions are specified in
 * <tt>double[2]</tt>-sized units.
 *
 * @param[in] nderivative Derivative operator to apply.  May be zero.
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
 * @see suzerain_bspline_accumulate_operator() for a way to apply an operator
 *      to real-valued coefficients.
 * @see suzerain_bspline_apply_operator() for a way to apply an operator in
 *      place.
 */
int
suzerain_bspline_zaccumulate_operator(
    int nderivative,
    int nrhs,
    const double alpha[2],
    const double (*x)[2],
    int incx,
    int ldx,
    const double beta[2],
    double (*y)[2],
    int incy,
    int ldy,
    const suzerain_bspline_workspace *w);

/**
 * Apply the <tt>nderivative</tt>-th derivative operator to real coefficients
 * \c x.  Multiplies the precomputed banded derivative operator scaled by \c
 * alpha against one or more coefficient vectors stored in \c x.  Results
 * overwrite \c x.  Each coefficient vector is of length
 * suzerain_bspline_ndof().  Increments and leading dimensions are specified in
 * <tt>double</tt>-sized units.
 *
 * @param[in] nderivative Derivative operator to apply.  May be zero.
 * @param[in] nrhs Number of coefficient vectors stored in \c x.
 * @param[in] alpha Real scaling factor to apply.
 * @param[in,out] x Coefficients to be multiplied.
 * @param[in] incx Stride between elements stored in \c x.
 * @param[in] ldx Leading dimension of the data stored in \c x.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bspline_zapply_operator() for a way to apply an operator
 *      to complex-valued coefficients.
 * @see suzerain_bspline_accumulate_operator() for a way to accumulate the
 *      effect of an operator against other storage.
 */
int
suzerain_bspline_apply_operator(
    int nderivative,
    int nrhs,
    double alpha,
    double *x,
    int incx,
    int ldx,
    const suzerain_bspline_workspace *w);

/**
 * Apply the <tt>nderivative</tt>-th derivative operator to complex
 * coefficients \c x.  Multiplies the precomputed banded derivative operator
 * scaled by real-valued \c alpha against one or more coefficient vectors
 * stored in \c x.  Results overwrite \c x.  Each coefficient vector is of
 * length suzerain_bspline_ndof().  Increments and leading dimensions are
 * specified in <tt>double[2]</tt>-sized units.
 *
 * @param[in] nderivative Derivative operator to apply.  May be zero.
 * @param[in] nrhs Number of coefficient vectors stored in \c x.
 * @param[in] alpha Real scaling factor to apply.
 * @param[in,out] x Coefficients to be multiplied.
 * @param[in] incx Stride between elements stored in \c x.
 * @param[in] ldx Leading dimension of the data stored in \c x.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 * @see suzerain_bspline_apply_operator() for a way to apply an operator
 *      to real-valued coefficients.
 * @see suzerain_bspline_accumulate_operator() for a way to accumulate the
 *      effect of an operator against other storage.
 */
int
suzerain_bspline_zapply_operator(
    int nderivative,
    int nrhs,
    double alpha,
    double (*x)[2],
    int incx,
    int ldx,
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
 * For collocation methods like ::SUZERAIN_BSPLINE_COLLOCATION_GREVILLE, obtain
 * the <tt>j</tt>-th collocation point \f$x_j\f$.
 *
 * @param[in] j      Index of the desired collocation point.
 *                   Must be in the range <tt>[0,w->ndof)</tt>.
 * @param[out] x_j Location of collocation point \f$x_j\f$.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_collocation_point(
    int j,
    double *x_j,
    const suzerain_bspline_workspace *w);

/**
 * For collocation methods like ::SUZERAIN_BSPLINE_COLLOCATION_GREVILLE obtain
 * all collocation points \f$x_j\f$ for x in <tt>[0,w->ndof)</tt>.
 *
 * @param[out] x   Memory in which to store collocation points.
 * @param[in] incx Increment between collocation points, measured in
 *                 <tt>sizeof(double)</tt>.
 * @param[in] w Workspace to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_collocation_points(
    double *x,
    int incx,
    const suzerain_bspline_workspace *w);

/**
 * Encapsulates the LU decomposition of a real-valued linear combination of
 * derivative operators.  Callers obtain a workspace using
 * suzerain_bspline_lu_alloc() and release it using suzerain_bspline_lu_free().
 *
 * @see suzerain_bspline_luz_workspace for manipulating the complex-valued case.
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
 * @param[in] incb Stride between elements in matrix \c b.
 * @param[in] ldb Leading dimension of matrix \c b.
 * @param[in] luw Workspace containing the factored operator to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_lu_solve(
    int nrhs,
    double *b,
    int incb,
    int ldb,
    const suzerain_bspline_lu_workspace *luw);

/**
 * Encapsulates the LU decomposition of a complex-valued linear combination of
 * derivative operators.  Callers obtain a workspace using
 * suzerain_bspline_luz_alloc() and release it using suzerain_bspline_luz_free().
 *
 * @see suzerain_bspline_lu_workspace for manipulating the real-valued case.
 */
typedef struct suzerain_bspline_luz_workspace {

    /** Number of degrees of freedom in the basis */
    int ndof;

    /** Number of subdiagonals in the factored operator. */
    int kl;

    /**
     * Number of superdiagonals in the factored operator.
     * Note that, for a given B-spline basis, suzerain_bspline_luz_workspace::ku
     * is larger than the corresponding suzerain_bspline_workspace::ku
     * according to the requirements of LAPACK's \c zgbtrf.
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
    double (*A)[2];

} suzerain_bspline_luz_workspace;


/**
 * Allocates a B-spline operator LU workspace.
 *
 * @param[in] w B-spline workspace on which to base the new workspace's
 *      dimensions.
 *
 * @return a workspace instance ready to be used with
 *      suzerain_bspline_luz_form_general() or suzerain_bspline_luz_form_mass().
 *      On failure calls suzerain_error() and returns NULL.
 */
suzerain_bspline_luz_workspace *
suzerain_bspline_luz_alloc(
    const suzerain_bspline_workspace *w);

/**
 * Frees a previously allocated workspace.
 *
 * @param[in] luzw Workspace to free.
 */
void
suzerain_bspline_luz_free(
    suzerain_bspline_luz_workspace *luzw);

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
 * @param[in,out] luzw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_luz_form_general(
    int ncoefficients,
    const double (*coefficients)[2],
    const suzerain_bspline_workspace * w,
    suzerain_bspline_luz_workspace *luzw);

/**
 * Forms the LU decomposition of the zeroth derivative operator, also
 * called the mass matrix, from workspace \c w.
 *
 * @param[in] w Workspace containing the desired mass matrix.
 * @param[in,out] luzw Workspace in which to store the factored operator.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 *
 * @see This is a convenience wrapper function built atop
 *      suzerain_bspline_luz_form_general().
 */
int
suzerain_bspline_luz_form_mass(
    const suzerain_bspline_workspace * w,
    suzerain_bspline_luz_workspace *luzw);

/* TODO Add noncontiguous luz_solve support to match lu_solve */

/**
 * Solves the equation <tt>A x = b</tt> using the factored operator
 * stored in \c luzw.
 *
 * @param[in] nrhs Number of right hand sides, i.e. columns, stored in
 *      matrix \c b.
 * @param[in,out] b Matrix of right hand sides to solve and the resulting
 *      solutions.  On input, contains data \c b.  On output, contains
 *      solutions \c x.
 * @param[in] ldb Leading dimension of matrix \c b.
 * @param[in] luzw Workspace containing the factored operator to use.
 *
 * @return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bspline_luz_solve(
    int nrhs,
    double (*b)[2],
    int ldb,
    const suzerain_bspline_luz_workspace *luzw);

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif // __SUZERAIN_BSPLINE_OPERATOR_H
