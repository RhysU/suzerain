//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_BL_H
#define SUZERAIN_BL_H

/** @file
 * Compute boundary layer quantities of interest.
 */

#include <suzerain/bspline.h>
#include <suzerain/bsplineop.h>

#ifdef __cplusplus
extern "C" {
#endif

// TODO Document
typedef struct {
    double a;
    double gamma;
    double mu;
    double Pr;
    double p__x;
    double rho;
    double T;
    double u;
    double u__x;
    double u__y;
    double v;
} suzerain_bl_local;

// TODO Document
typedef struct {
    double delta;
    double deltastar;
    double theta;
} suzerain_bl_thick;

// TODO Document
typedef struct {
    double beta;
    double Cf;
    double gamma_e;
    double K_e;
    double K_s;
    double K_w;
    double Lambda_n;
    double Ma_e;
    double p_exi;
    double Pr_w;
    double Re_delta;
    double Re_deltastar;
    double Re_theta;
    double shapefactor;
    double T_ratio;
    double v_wallplus;
} suzerain_bl_qoi;

// FIXME Account for any nondimensionalization
// TODO Document
int
suzerain_bl_compute_qoi(
//      double code_Ma,
//      double code_Re,
        const suzerain_bl_local * wall,
        const suzerain_bl_local * edge,
        const suzerain_bl_thick * thick,
              suzerain_bl_qoi   * qoi);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BL_H */
