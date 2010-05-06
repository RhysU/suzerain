/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * htstretch.c: Routines providing hyperbolic tangent grid stretching
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_HTSTRETCH_H__
#define __SUZERAIN_HTSTRETCH_H__

/** @file
 * Provides utilities atop POSIX methods.
 */

#ifdef __cplusplus
extern "C" {
#endif

double
suzerain_htstretch_onesided_continuous(const double I,
                                       const double delta,
                                       const double xi);

double
suzerain_htstretch_onesided_discrete(const double I,
                                     const double delta,
                                     const int N,
                                     const int j);

double
suzerain_htstretch_twosided_continuous(const double I,
                                       const double delta,
                                       const double xi);

double
suzerain_htstretch_twosided_discrete(const double I,
                                     const double delta,
                                     const int N,
                                     const int j);

double
suzerain_htstretch_twosided_discrete_delta(const double I,
                                           const int N,
                                           const int k,
                                           const double u_crit);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_HTSTRETCH_H__ */
