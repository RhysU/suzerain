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

#ifdef __cplusplus
extern "C" {
#endif

double
suzerain_htstretch1(const double delta,
                    const double L,
                    const double x);

double
suzerain_htstretch1_ddelta(const double delta,
                           const double L,
                           const double x);

double
suzerain_htstretch1_dL(const double delta,
                       const double L,
                       const double x);

double
suzerain_htstretch1_dx(const double delta,
                       const double L,
                       const double x);

double
suzerain_htstretch2(const double delta,
                    const double L,
                    const double x);

double
suzerain_htstretch2_ddelta(const double delta,
                           const double L,
                           const double x);

double
suzerain_htstretch2_dL(const double delta,
                       const double L,
                       const double x);

double
suzerain_htstretch2_dx(const double delta,
                       const double L,
                       const double x);

double
suzerain_htstretch_twosided_delta(const double L,
                                  const double xi_crit,
                                  const double u_crit);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_HTSTRETCH_H__ */
