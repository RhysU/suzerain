/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2011-2014 Rhys Ulerich
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

#ifndef SUZERAIN_COMPLEX_H
#define SUZERAIN_COMPLEX_H

/** @file
 * Compatibility file for C99 \c _Complex and C++ <tt>std::complex</tt> types.
 *
 * Modified from "Complex Arithmetic: In the Intersection of C and C++"
 * by Randy Meyers and Dr. Thomas Plum (http://drdobbs.com/cpp/184401628).
 * Intent is for C++ source to be able to seamlessly call C99-based functions
 * without needlessly polluting the C++ global namespace.
 */

#ifdef __cplusplus

#include <cmath>
#include <complex>

/** Allows C99/C++ <tt>float</tt>-based complex type interoperability */
typedef std::complex<float>       complex_float;

/** Allows C99/C++ <tt>double</tt>-based complex type interoperability */
typedef std::complex<double>      complex_double;

/** Allows C99/C++ <tt>long double</tt>-based complex type interoperability */
typedef std::complex<long double> complex_long_double;

#else   // #ifdef __cplusplus

// Note that <tgmath.h> includes <math.h> and <complex.h>
#include <tgmath.h>

/** Allows C99/C++ <tt>float</tt>-based complex type interoperability */
typedef _Complex float       complex_float;

/** Allows C99/C++ <tt>double</tt>-based complex type interoperability */
typedef _Complex double      complex_double;

/** Allows C99/C++ <tt>long double</tt>-based complex type interoperability */
typedef _Complex long double complex_long_double;

/** Construct a <tt>float</tt>-based complex with value <tt>r+I*i</tt> */
#define complex_float(r,i)       ((float)      (r) + ((float)      (i))*_Complex_I)

/** Construct a <tt>double</tt>-based complex with value <tt>r+I*i</tt> */
#define complex_double(r,i)      ((double)     (r) + ((double)     (i))*_Complex_I)

/** Construct a <tt>long double</tt>-based complex with value <tt>r+I*i</tt> */
#define complex_long_double(r,i) ((long double)(r) + ((long double)(i))*_Complex_I)

#endif  // #ifdef __cplusplus

#endif  // #ifndef SUZERAIN_COMPLEX_H
