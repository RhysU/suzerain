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

#ifndef SUZERAIN_CHANNEL_H
#define SUZERAIN_CHANNEL_H

/** @file
 * Compute channel quantities of interest.
 *
 * Beware that many of these computations neglect covariance-related
 * correction terms that should be present when working with functions
 * of expected values.  For example, when computed from mean velocity
 * \f$\bar{u}\f$ and speed of sound \f$\bar{a}\f$,
 * \f[
 *     \operatorname{E}[\mbox{Ma}] \approx \operatorname{E}[u]/E[a]
 *                  +       \operatorname{Var}(a)
 *                          \operatorname{E}[u]
 *                         /\operatorname{E}[a]^3
 *                  -       \operatorname{Cov}[a,u]
 *                          \operatorname{E}[a]^2
 * \f]
 * are the appropriate corrections to be second-order accurate.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Information characterizing local state in a channel flow.
 * Each nondimensional quantity is scaled as documented.
 */
typedef struct suzerain_channel_local {
    double y;      /**< Distance from the wall divided by \f$l_0\f$.        */
    double a;      /**< Sound speed divided by \f$a_0\f$.                   */
    double gamma;  /**< Nondimensional Ratio of specific heats \f$\gamma\f$.*/
    double mu;     /**< Dynamic viscosity divided by \f$\mu_0\f$.           */
    double Pr;     /**< Nondimensional Prandtl number \f$C_p \mu/\kappa\f$. */
    double rho;    /**< Density divided by \f$\rho_0\f$.                    */
    double T;      /**< Temperature divided by \f$T_0\f$.                   */
    double T__y;   /**< Wall-normal derivative of temperature
                        divided by \f$T_0 / l_0\f$.                         */
    double u;      /**< Streamwise velocity divided by \f$u_0\f$.           */
    double u__y;   /**< Wall-normal derivative of streamwise velocity
                        divided by \f$u_0 / l_0\f$.                         */
    double v;      /**< Wall-normal velocity divided by \f$u_0\f$.          */
} suzerain_channel_local;

/**
 * Information characterizing wall-related viscous scales.
 * Each nondimensional quantity is scaled as documented.
 */
typedef struct suzerain_channel_viscous {
    double tau_w;    /**< Wall shear stress \f$\tau_w\f$
                          divided by \f$\rho_0 u_0^2\f$.                   */
    double u_tau;    /**< Wall velocity \f$u_\tau\f$ divided by \f$u_0\f$. */
    double delta_nu; /**< Wall length scale divided by \f$l_0\f$.          */
} suzerain_channel_viscous;

/**
 * Compute viscous-related wall scalings.
 *
 * \param[in ] code_Re Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code.
 * \param[in ] wall    Local state information from the wall.
 * \param[out] viscous Populated on success.
 *                     See type documentation for contents.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_channel_compute_viscous(
        double code_Re,
        const suzerain_channel_local   * wall,
              suzerain_channel_viscous * viscous);

/**
 * Nondimensional channel flow quantities of general interest.
 */
typedef struct suzerain_channel_qoi {
    double cf;           /**< The skin friction coefficient \f$c_f =
                              \frac{2 \tau_w}{\rho_c u_c^2}\f$ based on
                              centerline quantities. */
    double gamma_c;      /**< The ratio of specific heats at the centerline. */
    double Ma_c;         /**< The local Mach number at the centerline. */
    double Ma_tau;       /**< The friction Mach number
                              \f$\mbox{Ma}_\tau = u_\tau / a_w\f$. */
    double Pr_w;         /**< The Prandtl number at the wall. */
    double neg_Bq;       /**< The negated nondimensional heat flux
                              \f$-B_q = \frac{k \nabla T}{\rho_w C_p
                              u_\tau T_w} = \frac{\mu \nabla T}{\mbox{Pr}
                              \rho_w C_p u_\tau T_w} \f$. */
    double ratio_rho;    /**< The ratio of center to wall density. */
    double ratio_nu;     /**< The ratio of center to wall kinematic viscosity. */
    double ratio_T;      /**< The ratio of center to wall temperature. */
    double Re_c;         /**< Reynolds number based on centerline state
                              and channel half width. */
    double Re_tau;       /**< Reynolds number based on viscous scales
                              and channel half width. */
    double v_wallplus;   /**< The wall transpiration rate in plus units.*/
} suzerain_channel_qoi;

/**
 * Compute nondimensional channel flow parameters.
 *
 * \param[in ] code_Ma Mach number \f$u_0/a_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code.
 * \param[in ] code_Re Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code.
 * \param[in ] wall    Local state information from the wall.
 * \param[in ] viscous Viscous-related wall scaling information
 * \param[in ] center  Local state information from the channel centerline.
 * \param[out] qoi     Populated on success.
 *                     See type documentation for contents.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_channel_compute_qoi(
        double code_Ma,
        double code_Re,
        const suzerain_channel_local       * wall,
        const suzerain_channel_viscous     * viscous,
        const suzerain_channel_local       * center,
              suzerain_channel_qoi         * qoi);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_CHANNEL_H */
