/*--------------------------------------------------------------------------
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 * Adapted from the GNU Scientific Library
 *
 * Extensions Copyright (C) 2011-2014 Rhys Ulerich
 * Extensions Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_FUNCTION_H
#define SUZERAIN_FUNCTION_H

/** @file
 * Provides a function evaluation interface compatible with the <a
 * href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a>'s <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Providing-the-function-to-solve.html">
 * <tt>gsl_function</tt></a>.
 *
 * For example, the following snippet establishes a suzerain_function instance
 * which performs an arbitrary order polynomial evaluation via GSL's <a
 * href="http://www.gnu.org/software/gsl/manual/html_node/Polynomial-Evaluation.html">
 * <tt>gsl_poly_eval</tt></a>:
 * \code
 * // Declarations
 * typedef struct { int n; double c[]; } poly_params; // Flexible array
 *
 * double poly_f(double x, void *params)
 * {
 *     poly_params *p = (poly_params *) params;
 *     return gsl_poly_eval(p->c, p->n, x);
 * }
 *
 * // Within some method
 * poly_params *p = (poly_params *)
 *                   malloc(sizeof(poly_params) + 3*sizeof(double));
 * p->n    = 3;
 * p->c[0] = 1.1; // Constant
 * p->c[1] = 2.2; // Linear
 * p->c[2] = 3.3; // Quadratic
 * suzerain_function f = {poly_f, p};
 *
 * const double result = SUZERAIN_FN_EVAL(&f,1.0) // equals 6.6
 * \endcode
 */

#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Represents an arbitrary real-valued function of two parameters with the
 * second parameter fixed.  Used to fix some parameters of a multivariate
 * function so that it can be passed to a routine expecting a univariate
 * function.  This structure must be binary compatible with GSL's \c
 * gsl_function.
 *
 * @see The function evaluation macro SUZERAIN_FN_EVAL.
 * @see suzerain_function for the complex-valued analogue.
 */
typedef struct suzerain_function {

    /**
     * Underlying function to be evaluated
     *
     * @param x Single parameter that will be varied when invoked.
     * @param params Arbitrary data to be passed to the function.
     */
    double (* function) (double x, void * params);

    /**
     * Fixed data to be provided to \c function on each invocation.
     * May be safely mutated between invocations.
     */
    void * params;

} suzerain_function ;

/**
 * \def SUZERAIN_FN_EVAL(F,x)
 * Evaluate \c F(x) supplying <tt>F->params</tt> to <tt>F->function</tt>.
 *
 * @param F Function to evaluate given as a suzerain_function.
 * @param x Point at which \c F should be evaluated.
 * @return The value \c F(x)
 */
#define SUZERAIN_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/**
 * Represents an arbitrary complex-valued function of two parameters with the
 * second parameter fixed.
 *
 * @see The function evaluation macro SUZERAIN_ZFN_EVAL.
 * @see suzerain_zfunction for the real-valued analogue.
 */
typedef struct suzerain_zfunction {

    /**
     * Underlying function to be evaluated
     *
     * @param[in]  x Single parameter that will be varied when invoked.
     * @param[in]  params Arbitrary data to be passed to the function.
     * @param[out] z Complex-valued result of evaluation.
     */
    void (* function) (double x, void * params, complex_double z);

    /**
     * Fixed data to be provided to \c function on each invocation.
     * May be safely mutated between invocations.
     */
    void * params;

} suzerain_zfunction ;

/**
 * \def SUZERAIN_ZFN_EVAL(F,x,z)
 * Evaluate <tt>z == F(x)</tt> supplying <tt>F->params</tt> to
 * <tt>F->function</tt>.  Unlike SUZERAIN_FN_EVAL, this macro is a statement
 * unto itself.
 *
 * @param F Function to evaluate given as a suzerain_zfunction.
 * @param x Point at which \c F should be evaluated.
 * @param z Complex-valued result of evaluation.
 */
#define SUZERAIN_ZFN_EVAL(F,x,z) \
    do { (*((F)->function))(x,(F)->params,z); } while (0)

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_FUNCTION_H */
