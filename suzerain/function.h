/*
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 * Adapted from the GNU Scientific Library
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __SUZERAIN_FUNCTION_H__
#define __SUZERAIN_FUNCTION_H__

/** @file
 * Provides a function evaluation interface compatible with the GNU Scientific
 * Library's \c gsl_function.
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

/**
 * Represents an arbitrary function of two parameters with the second parameter
 * fixed.  Used to fix some parameters of a multivariate function so that it
 * can be passed to a routine expecting a univariate function.  This structure
 * must be binary compatible with GSL's \c gsl_function.
 *
 * @see The function evaluation macro SUZERAIN_FN_EVAL.
 */
typedef struct {

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
 * Evalute \c F(x) supplying <tt>F->params</tt> to it under the covers.
 * @param F Function to evaluate.
 * @param x Point at which \c F should be evaluated.
 */
#define SUZERAIN_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
__END_DECLS
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif /* __SUZERAIN_FUNCTION_H__ */
