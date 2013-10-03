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
 * Quantities possess dimensional units or nondimensional scaling as
 * detailed for each member.
 */
typedef struct {
    double a;      /**< Sound speed with units \f$a_0\f$.                     */
    double gamma;  /**< Nondimensional Ratio of specific heats \f$\gamma\f$.  */
    double mu;     /**< Dynamic viscosity with units \f$\mu_0\f$.             */
    double Pr;     /**< Nondimensional Prandtl number \f$C_p \mu / \kappa\f$. */
    double p__x;   /**< Streamwise pressure gradient with units
                        \f$p_0 / l_0 = \rho_0 a_0^2 / l_0\f$.                 */
    double rho;    /**< Density with units \f$\rho_0\f$.                      */
    double T;      /**< Temperature with units \f$T_0\f$.                     */
    double u;      /**< Streamwise velocity with units \f$u_0\f$.             */
    double u__x;   /**< Streamwise derivative of streamwise velocity
                        with units \f$u_0 / l_0\f$.                           */
    double u__y;   /**< Wall-normal derivative of streamwise velocity
                        with units \f$u_0 / l_0\f$.                           */
    double v;      /**< Wall-normal velocity with units \f$u_0\f$.            */
} suzerain_bl_local;

/**
 * Information characterizing wall-related viscous scales.
 * Quantities possess dimensional units or nondimensional scaling as
 * detailed for each member.
 */
typedef struct {
    double tau_w;    /**< Wall shear stress \f$\tau_w\f$
                          with units \f$\mu_0 u_0 / l_0\f$.                */
    double u_tau;    /**< Wall velocity \f$u_\tau\f$ with units \f$u_0\f$. */
    double delta_nu; /**< Wall length scale with units \f$l_0\f$.          */
} suzerain_bl_viscous;

/**
 * Compute viscous-related wall scalings.
 *
 * \param[in ] wall    Local state information from the wall.
 * \param[out] viscous Populated on success.
 *                     See type documentation for contents.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bl_compute_viscous(
        const suzerain_bl_local   * wall,
              suzerain_bl_viscous * viscous);

/**
 * Information characterizing boundary layer thickness in various ways.
 * Each member has units of \f$l_0\f$.
 */
typedef struct {
    double delta;     /**< Boundary layer thickness \f$\delta\f$. */
    double deltastar; /**< Displacement thickness \f$\delta^\ast\f$. */
    double theta;     /**< Momentum thickness \f$\theta\f$. */
} suzerain_bl_thick;

/**
 * Nondimensional boundary layer quantities of interest.
 *
 * Many of the pressure gradient parameters are defined within Cal and
 * Castillo's <a href="http://dx.doi.org/10.1063/1.2991433">"Similarity
 * analysis of favorable pressure gradient turbulent boundary layers with
 * eventual quasilaminarization."</a> in Physics of Fluids 20 (2008).
 */
typedef struct {
    double beta;         /**< The Clauser parameter \f$\beta =
                              \frac{\delta^\ast}{\tau_w}
                              \frac{\partial p}{\partial \xi}\f$. */
    double Cf;           /**< The skin friction coefficient \f$C_f =
                              \frac{2 \tau_w}{\rho_e u_e^2}\f$. */
    double gamma_e;      /**< The ratio of specific heats at the edge. */
    double K_e;          /**< Launder's acceleration parameter \f$K =
                              \frac{\mu}{\rho_e u_e^2} \,
                              \frac{\partial{}u_e}{\partial\xi}\f$
                              computed taking \f$\mu = \mu_e\f$. */
    double K_s;          /**< The Pohlhausen parameter \f$K_s =
                              \frac{\delta^2}{\nu_e}
                              \frac{\partial u_e}{\partial \xi}\f$. */
    double K_w;          /**< Launder's acceleration parameter \f$K =
                              \frac{\mu}{\rho_e u_e^2} \,
                              \frac{\partial{}u_e}{\partial\xi}\f$
                              computed taking \f$\mu = \mu_w\f$. */
    double Lambda_n;     /**< Pressure parameter \f$\Lambda_n =
                              -\frac{\delta}{\tau_w}
                              \frac{\partial p}{\partial \xi}\f$. */
    double Ma_e;         /**< The local Mach number at the edge. */
    double p_exi;        /**< The inviscid-friendly pressure parameter
                              \f$p_{e,\xi}^\ast = \frac{\delta}{\rho_e u_e^2}
                              \frac{\partial p}{\partial \xi}\f$. */
    double Pr_w;         /**< The Prandtl number at the wall. */
    double Re_delta;     /**< Reynolds number based on boundary layer
                              thickness \f$\delta\f$ and \f$\nu_e\f$. */
    double Re_deltastar; /**< Reynolds number based on displacement
                              thickness \f$\delta^\ast\f$ and \f$\nu_e\f$. */
    double Re_theta;     /**< Reynolds number based on momentum
                              thickness \f$\theta\f$ and \f$\nu_e\f$. */
    double ratio_rho;    /**< The ratio of edge to wall density. */
    double ratio_nu;     /**< The ratio of edge to wall kinematic viscosity. */
    double ratio_T;      /**< The ratio of edge to wall temperature. */
    double shapefactor;  /**< The shape factor \f$\delta^\ast / \theta. */
    double v_wallplus;   /**< The wall transpiration rate in plus units.*/
} suzerain_bl_qoi;

/**
 * Compute nondimensional boundary layer parameters.
 *
 * \param[in ] code_Ma Mach number \f$u_0/a_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code.
 * \param[in ] code_Re Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code.
 * \param[in ] wall    Local state information from the wall.
 * \param[in ] viscous Viscous-related wall scaling information
 * \param[in ] edge    Local state information from the boundary layer edge.
 * \param[in ] thick   Thickness information for the boundary layer.
 * \param[out] qoi     Populated on success.
 *                     See type documentation for contents.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bl_compute_qoi(
        double code_Ma,
        double code_Re,
        const suzerain_bl_local   * wall,
        const suzerain_bl_viscous * viscous,
        const suzerain_bl_local   * edge,
        const suzerain_bl_thick   * thick,
              suzerain_bl_qoi     * qoi);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BL_H */
