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

#ifndef SUZERAIN_COMMON_H
#define SUZERAIN_COMMON_H

/** @file
 * Common includes, definitions, utility macros, and inline functions for C.
 */

#include <suzerain/suzerain-config.h>

/*
 * Release-mode specific #defines whenever NDEBUG in effect.  We do not use
 * SUZERAIN_NDEBUG to determine this behavior as it would impact all source
 * code that needs to include our header due to <suzerain/suzerain-config.h>
 */
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
 * Provides hint to the compiler to optimize for the expression being true.
 * @param cond expression to be evaluted in boolean context.
 * @returns value of expression.
 */
#define SUZERAIN_LIKELY(cond) __builtin_expect(!!(cond), 1)
/**
 * Provides hint to the compiler to optimize for the expression being false.
 * @param cond expression to be evaluted in boolean context.
 * @returns value of expression.
 */
#define SUZERAIN_UNLIKELY(cond) __builtin_expect(!!(cond), 0)
#else
/**
 * Provides hint to the compiler to optimize for the expression being true.
 * @param cond expression to be evaluted in boolean context.
 * @returns value of expression.
 */
#define SUZERAIN_LIKELY(cond) !!(cond)
/**
 * Provides hint to the compiler to optimize for the expression being false.
 * @param cond expression to be evaluted in boolean context.
 * @returns value of expression.
 */
#define SUZERAIN_UNLIKELY(cond) !!(cond)
#endif

#ifdef __INTEL_COMPILER
/** Strongly recommend that a function be inlined */
#define SUZERAIN_FORCEINLINE __forceinline
#elif __GNUC__ >= 3
/** Strongly recommend that a function be inlined */
#define SUZERAIN_FORCEINLINE __attribute__((always_inline))
#else
/** Strongly recommend that a function be inlined */
#define SUZERAIN_FORCEINLINE inline
#endif

#ifdef __INTEL_COMPILER
/** Prevent inlining a given function */
#define SUZERAIN_DONTINLINE __attribute__((noinline))
#elif __GNUC__ >= 3
/** Prevent inlining a given function */
#define SUZERAIN_DONTINLINE __attribute__((noinline))
#else
/** Prevent inlining a given function */
#define SUZERAIN_DONTINLINE
#endif

#ifndef __cplusplus
/** Define synonym for C99 restrict keyword */
#define SUZERAIN_RESTRICT restrict
#else
/** Define synonym for C++ __restrict__ keyword */
#define SUZERAIN_RESTRICT __restrict__
#endif

// Works around -Wstrict-aliasing warnings similarly to
// http://www.xemacs.org/Documentation/sources/xemacs/src/compiler.h.
// Seemingly requires using typeof for both GCC/Intel correctness.
#ifdef SUZERAIN_typeof
/** Cast to <tt>const void *</tt> without tripping strict aliasing warnings */
#define SUZERAIN_CONST_CAST_VOIDP(l)                          \
    (((union { const void *p; SUZERAIN_typeof(l) v; })(l)).p)
/** Cast to <tt>      void *</tt> without tripping strict aliasing warnings */
#define SUZERAIN_CAST_VOIDP(l)                                \
    (((union {       void *p; SUZERAIN_typeof(l) v; })(l)).p)
#else
/** Cast to <tt>const void *</tt> without tripping strict aliasing warnings */
#define SUZERAIN_CONST_CAST_VOIDP(l)                          \
    (((union { const void *p; typeof(l)          v; })(l)).p)
/** Cast to <tt>      void *</tt> without tripping strict aliasing warnings */
#define SUZERAIN_CAST_VOIDP(l)                                \
    (((union {       void *p; typeof(l)          v; })(l)).p)
#endif

/* Push/pop-like macros for suppressing GCC warnings lifted from   */
/* http://dbp-consulting.com/tutorials/SuppressingGCCWarnings.html */
#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 402
# define SUZERAIN_GCC_DIAG_STR(s) #s
# define SUZERAIN_GCC_DIAG_JOINSTR(x,y) SUZERAIN_GCC_DIAG_STR(x ## y)
# define SUZERAIN_GCC_DIAG_DO_PRAGMA(x) _Pragma (#x)
# define SUZERAIN_GCC_DIAG_PRAGMA(x) SUZERAIN_GCC_DIAG_DO_PRAGMA(GCC diagnostic x)
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#  define SUZERAIN_GCC_DIAG_OFF(x) SUZERAIN_GCC_DIAG_PRAGMA(push) \
        SUZERAIN_GCC_DIAG_PRAGMA(ignored SUZERAIN_GCC_DIAG_JOINSTR(-W,x))
#  define SUZERAIN_GCC_DIAG_ON(x) SUZERAIN_GCC_DIAG_PRAGMA(pop)
# else
#  define SUZERAIN_GCC_DIAG_OFF(x) \
    SUZERAIN_GCC_DIAG_PRAGMA(ignored SUZERAIN_GCC_DIAG_JOINSTR(-W,x))
#  define SUZERAIN_GCC_DIAG_ON(x)  \
    SUZERAIN_GCC_DIAG_PRAGMA(warning SUZERAIN_GCC_DIAG_JOINSTR(-W,x))
# endif
#else
# define SUZERAIN_GCC_DIAG_OFF(x)
# define SUZERAIN_GCC_DIAG_ON(x)
#endif

/**
 * Explicitly marks a variable as unused to avoid compiler warnings.
 *
 * @param variable Variable to be marked.
 */
#define SUZERAIN_UNUSED(variable) do {(void)(variable);} while (0)

/*
 * Required standard C functionality appearing through Suzerain
 * Care taken to include C++ header versions for C++ compilation
 * C++-specific functionality is included in common.hpp
 */
#ifndef __cplusplus
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <time.h>
#else
#include <cassert>
#include <cctype>
#include <cerrno>
/* cinttypes not commonly available pre-C++0x */
#include <climits>
#include <cmath>
#include <cstddef>
/* cstdint not commonly available pre-C++0x */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#endif
#include <fcntl.h>
#include <signal.h>
#include <strings.h>
#include <unistd.h>

/* Boost.Preprocessor is useful in either C or C++ mode */
#include <boost/preprocessor.hpp>

#endif // SUZERAIN_COMMON_H
