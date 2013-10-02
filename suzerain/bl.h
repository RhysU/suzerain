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

/**
 * Information characterizing local state in a boundary layer.
 */
typedef struct {
    double a;      /**< Sound speed with units \f$a_0\f$.                    */
    double gamma;  /**< Nondimensional Ratio of specific heats \f$\gamma\f$. */
    double mu;     /**< Dynamic viscosity with units \f$\mu_0\f$.            */
    double Pr;     /**< Nondimensional Prandtl number \f$C_p \mu / \kappa.   */
    double p__x;   /**< Streamwise pressure gradient with units
                        \f$p_0 / l_0 = \rho_0 a_0^2 / l_0\f$.                */
    double rho;    /**< Density with units \f$\rho_0\f$.                     */
    double T;      /**< Temperature with units \f$T_0\f$.                    */
    double u;      /**< Streamwise velocity with units \f$u_0\f$.            */
    double u__x;   /**< Streamwise derivative of streamwise velocity
                        with units \f$u_0 / l_0\f$.                          */
    double u__y;   /**< Wall-normal derivative of streamwise velocity
                        with units \f$u_0 / l_0\f$.                          */
    double v;      /**< Wall-normal velocity with units \f$u_0\f$.           */
} suzerain_bl_local;

/**
 * Information characterizing boundary layer thickness in various ways.
 * Each member has units of \f$l_0\f$.
 */
typedef struct {
    double delta;     /**< Boundary layer thickness.  */
    double deltastar; /**< */
    double theta;     /**< */
} suzerain_bl_thick;

/**
 * Nondimensional boundary layer quantities of interest.
 */
typedef struct {
    double  beta;         /**<  */
    double  Cf;           /**<  */
    double  gamma_e;      /**<  */
    double  K_e;          /**<  */
    double  K_s;          /**<  */
    double  K_w;          /**<  */
    double  Lambda_n;     /**<  */
    double  Ma_e;         /**<  */
    double  p_exi;        /**<  */
    double  Pr_w;         /**<  */
    double  Re_delta;     /**<  */
    double  Re_deltastar; /**<  */
    double  Re_theta;     /**<  */
    double  shapefactor;  /**<  */
    double  T_ratio;      /**<  */
    double  v_wallplus;   /**<  */
} suzerain_bl_qoi;

/**
 * FIXME Document
 * \memberof suzerain_bsplineop_workspace
 */
int
suzerain_bl_compute_qoi(
        double code_Ma,
        double code_Re,
        const suzerain_bl_local * wall,
        const suzerain_bl_local * edge,
        const suzerain_bl_thick * thick,
              suzerain_bl_qoi   * qoi);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BL_H */
