/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012-2014 Rhys Ulerich
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

#ifndef SUZERAIN_KAHAN_H
#define SUZERAIN_KAHAN_H

/** @file
 * Provides routines performing <a
 * href="http://en.wikipedia.org/wiki/Kahan_summation_algorithm">Kahan
 * summation</a>.
 */

#include <suzerain/complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Compute the Kahan summation of <code>a[0],...,a[n-1]</code>. */
double suzerain_kahan(const double * a, int n);

/** @copydoc suzerain_kahan(const double*, int) */
float suzerain_kahanf(const float * a, int n);

/** @copydoc suzerain_kahan(const double*, int) */
complex_float suzerain_kahanc(const complex_float * a, int n);

/** @copydoc suzerain_kahan(const double*, int) */
complex_double suzerain_kahanz(const complex_double * a, int n);

/** Compute the Kahan summation of <code>a0,...,a2</code>. */
double suzerain_kahan3(const double a0,
                       const double a1,
                       const double a2);

/** @copydoc suzerain_kahan3(const double, const double, const double) */
float suzerain_kahanf3(const float a0,
                       const float a1,
                       const float a2);

/** @copydoc suzerain_kahan3(const double, const double, const double) */
complex_float suzerain_kahanc3(const complex_float a0,
                               const complex_float a1,
                               const complex_float a2);

/** @copydoc suzerain_kahan3(const double, const double, const double) */
complex_double suzerain_kahanz3(const complex_double a0,
                                const complex_double a1,
                                const complex_double a2);

/** Compute the Kahan summation of <code>a0,...,a3</code>. */
double suzerain_kahan4(const double a0,
                       const double a1,
                       const double a2,
                       const double a3);

/** @copydoc suzerain_kahan4(const double, const double,
 *                           const double, const double) */
float suzerain_kahanf4(const float a0,
                       const float a1,
                       const float a2,
                       const float a3);

/** @copydoc suzerain_kahan4(const double, const double,
 *                           const double, const double) */
complex_float suzerain_kahanc4(const complex_float a0,
                               const complex_float a1,
                               const complex_float a2,
                               const complex_float a3);

/** @copydoc suzerain_kahan4(const double, const double,
 *                           const double, const double) */
complex_double suzerain_kahanz4(const complex_double a0,
                                const complex_double a1,
                                const complex_double a2,
                                const complex_double a3);

/** Compute the Kahan summation of <code>a0,...,a4</code>. */
double suzerain_kahan5(const double a0,
                       const double a1,
                       const double a2,
                       const double a3,
                       const double a4);

/** @copydoc suzerain_kahan5(const double, const double, const double
 *                           const double, const double) */
float suzerain_kahanf5(const float a0,
                       const float a1,
                       const float a2,
                       const float a3,
                       const float a4);

/** @copydoc suzerain_kahan5(const double, const double, const double
 *                           const double, const double) */
complex_float suzerain_kahanc5(const complex_float a0,
                               const complex_float a1,
                               const complex_float a2,
                               const complex_float a3,
                               const complex_float a4);

/** @copydoc suzerain_kahan5(const double, const double, const double
 *                           const double, const double) */
complex_double suzerain_kahanz5(const complex_double a0,
                                const complex_double a1,
                                const complex_double a2,
                                const complex_double a3,
                                const complex_double a4);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_KAHAN_H */
