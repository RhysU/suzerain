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
 * common.h: C common definitions, utility macros, and inline functions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_COMMON_H
#define __SUZERAIN_COMMON_H

#include <suzerain/suzerain-config.h>

// Release-mode specific #defines whenever NDEBUG in effect.  We do not use
// SUZERAIN_NDEBUG to determine this behavior as it would impact all source
// code that needs to include our header due to <suzerain/suzerain-config.h>
#ifdef NDEBUG
/** Disable GNU Scientific Library range checking */
#define GSL_RANGE_CHECK_OFF
/** Disable assertions within Boost */
#define BOOST_DISABLE_ASSERTS
/** Enable release mode within Boost uBLAS */
#define BOOST_UBLAS_NDEBUG
#endif

/** Enable boost::numeric::ublas::shallow_array_adaptor<T> */
#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR 1

#if __GNUC__ >= 3
/** Create an inline assembly comment, if supported */
#define SUZERAIN_ASM_INLINE_COMMENT(comment) __asm__("# "__FILE__##":"##SUZERAIN_ASM_INLINE_COMMENT_HELPER_TOSTRING(__LINE__)##" "## #comment)
/** Helps ASM_INLINE_COMMENT convert __LINE__ to a (char *) */
#define SUZERAIN_ASM_INLINE_COMMENT_HELPER_STRINGIFY(x) #x
/** Helps ASM_INLINE_COMMENT convert __LINE__ to a (char *) */
#define SUZERAIN_ASM_INLINE_COMMENT_HELPER_TOSTRING(x) SUZERAIN_ASM_INLINE_COMMENT_HELPER_STRINGIFY(x)
#else
/** Create an inline assembly comment, if supported */
#define SUZERAIN_ASM_INLINE_COMMENT(comment)
#endif

#if __GNUC__ >= 3
/**
Provides hint to the compiler to optimize for the expression being false.
@param expr boolean expression.
@returns value of expression.
*/
#define SUZERAIN_UNLIKELY(expr) __builtin_expect(expr, 0)
#else
/**
Provides hint to the compiler to optimize for the expression being false.
@param expr boolean expression.
@returns value of expression.
**/
#define SUZERAIN_UNLIKELY(expr) expr
#endif

#if __GNUC__ >= 3
/** Strongly recommend that a function be inlined */
#define SUZERAIN_FORCEINLINE __forceinline
#else
/** Strongly recommend that a function be inlined */
#define SUZERAIN_FORCEINLINE inline
#endif

/**
 * Explicitly marks a variable as unused to avoid compiler warnings.
 *
 * @param variable Variable to be marked.
 */
#define SUZERAIN_UNUSED(variable) do {(void)(variable);} while (0)

// Required standard C functionality appearing through Suzerain
// C++-specific functionality is included in common.hpp
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/**
 * Compute the BLAS-compatible offset to <tt>a(i,j)</tt> for general
 * banded matrices <tt>a(i,j) -> storage(ku+i-j,j)</tt> where storage
 * is column-major with LDA lda.  When comparing with BLAS documentation, note
 * missing constant one compared to Fortran due to C's 0-indexing.
 *
 * @param lda Leading dimension for column-major ordering
 * @param kl Number of superdiagonals
 * @param ku Number of subdiagonals
 * @param i Row index desired
 * @param j Column index desired
 *
 * @return Column-major offset where entry <tt>(i,j)</tt> is stored.
 */
inline
int suzerain_gbmatrix_offset(int lda, int kl, int ku, int i, int j) {
    /* Unused parameters present to be consistent with gbmatrix_in_band */
    return j*lda+(ku+i-j);
}

/** Determine if indices fall within the band of a general banded matrix.
 *
 * @param lda Leading dimension for column-major ordering
 * @param kl Number of superdiagonals
 * @param ku Number of subdiagonals
 * @param i Row index desired
 * @param j Column index desired
 *
 * @return true if entry <tt>(i,j)</tt> falls within the matrix's band.
 */
inline
int suzerain_gbmatrix_in_band(int lda, int kl, int ku, int i, int j) {
    /* Unused parameters present to be consistent with gbmatrix_offset */
    return j-ku <= i && i <= j+kl;
}

#endif // __SUZERAIN_COMMON_H
