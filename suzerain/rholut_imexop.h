/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
 * rholut_imexop.h: hybrid implicit/explicit operator apply and solve
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_RHOLUT_IMEXOP_H__
#define __SUZERAIN_RHOLUT_IMEXOP_H__

#include <suzerain/bsplineop.h>

/** @file
 * Provides implicit operator apply and solve routines for a single wall-normal
 * pencil of state information.  Meant to be used in conjunction with compute
 * kernels found in rholut.hpp.  Coded in C99 to facilitate profiling.
 */

#ifdef __cplusplus
extern "C" {
#endif

void suzerain_rholut_imexop_apply(const double         km,
                                  const double         kn,
                                  const double         Re,
                                  const double         Pr,
                                  const double         Ma,
                                  const double         alpha,
                                  const double         gamma,
                                  const double * const nu,
                                  const int            inc_nu,
                                  const double * const ux,
                                  const int            inc_ux,
                                  const double * const uy,
                                  const int            inc_uy,
                                  const double * const uz,
                                  const int            inc_uz,
                                  const double * const nuux,
                                  const int            inc_nuux,
                                  const double * const nuuy,
                                  const int            inc_nuuy,
                                  const double * const nuuz,
                                  const int            inc_nuuz,
                                  const double * const m_gradrho,
                                  const int            inc_m_gradrho,
                                  const double * const ex_gradrho,
                                  const int            inc_ex_gradrho,
                                  const double * const ey_gradrho,
                                  const int            inc_ey_gradrho,
                                  const double * const ez_gradrho,
                                  const int            inc_ez_gradrho,
                                  const double * const e_divm,
                                  const int            inc_e_divm,
                                  const double * const e_nablarho,
                                  const int            inc_e_nablarho,
                                  const double phi,
                                  const double (*in_rho )[2],
                                  const double (*in_rhou)[2],
                                  const double (*in_rhov)[2],
                                  const double (*in_rhow)[2],
                                  const double (*in_rhoe)[2],
                                  double (*out_rho )[2],
                                  double (*out_rhou)[2],
                                  double (*out_rhov)[2],
                                  double (*out_rhow)[2],
                                  double (*out_rhoe)[2],
                                  const suzerain_bsplineop_workspace * w);

void suzerain_rholut_imexop_build();



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SUZERAIN_RHOLUT_IMEXOP_H__ */
