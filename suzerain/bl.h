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

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Information characterizing local state in a boundary layer.
 * Each nondimensional quantity is scaled as documented.
 */
typedef struct suzerain_bl_local {
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
    double y;      /**< Distance from the wall divided by \f$l_0\f$.        */
} suzerain_bl_local;

/**
 * Information characterizing wall-related viscous scales.
 * Each nondimensional quantity is scaled as documented.
 */
typedef struct suzerain_bl_viscous {
    double tau_w;    /**< Wall shear stress \f$\tau_w\f$
                          divided by \f$\rho_0 u_0^2\f$.                   */
    double u_tau;    /**< Wall velocity \f$u_\tau\f$ divided by \f$u_0\f$. */
    double delta_nu; /**< Wall length scale divided by \f$l_0\f$.          */
} suzerain_bl_viscous;

/**
 * Compute viscous-related wall scalings.
 *
 * \param[in ] code_Re Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code>.
 * \param[in ] wall    Local state information from the wall.
 * \param[out] viscous Populated on success.
 *                     See type documentation for contents.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bl_compute_viscous(
        double code_Re,
        const suzerain_bl_local   * wall,
              suzerain_bl_viscous * viscous);

// TODO Find reference discussing NASA's technique mentioned below.
/**
 * Find the boundary layer edge within <tt>[lowerbnd, upperbnd]</tt> given a
 * B-spline coefficient representation of the specific total enthalpy \f$H_0 =
 * \frac{\rho E + p}{\rho}\f$ in \c coeffs_H0.  The procedure looks for the
 * second derivative of \f$H_0 \approx{} 0\f$ based on discussions with T.
 * Oliver about what NASA has done in practice to get robust results.  In an
 * incompressible context, other quantities might be substituted for \f$H_0\f$
 * in this computation.
 *
 * \param[in ] coeffs_H0 Coefficient representation of \f$H_0\f$
 *                       using the basis provided in \c w and \c dw.
 * \param[in]  lowerbnd  Lower bound (inclusive) for the search.
 * \param[in]  upperbnd  Upper bound (inclusive) for the search.
 * \param[out] location  Location at which edge is detected.
 * \param[in]  dB        Temporary storage to use of size <tt>w->k</tt> by
 *                       no less than <tt>3</tt>.
 * \param[in]  w         Workspace possessing non-constant second derivatives.
 * \param[in]  dw        Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*location</code>.  On recoverable error (e.g., no edge detected) sets
 * <code>*location</code> to be <tt>NaN</tt> and returns one of
 * #suzerain_error_status.  On unrecoverable error, additionally calls
 * suzerain_error().
 */
int
suzerain_bl_find_edge(
    const double * coeffs_H0,
    const double lowerbnd,
    const double upperbnd,
    double * location,
    gsl_matrix * dB,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw);

/**
 * Compute the displacement thickness \f$\delta_1\f$ (sometimes written
 * \f$\delta^\ast\f$) given the edge location and a B-spline coefficient
 * representation of streamwise momentum \f$\rho u\f$.  The method computes
 * \f[
 *  \delta_1 = \int_0^\infty
 *  \left(1 - \frac{\rho u}{\rho_\infty u_\infty}\right)
 *  \, \mathrm{d}y
 * \f]
 * where the value from the final B-spline collocation point is taken to
 * be "at infinity".
 *
 * \param[in ] coeffs_rhou B-spline coefficients for \f$\rho u\f$.
 * \param[out] delta1      The computed displacement thickness.
 * \param[in]  dB          Temporary storage to use of size <tt>w->k</tt> by
 *                         no less than <tt>1</tt>.
 * \param[in ] w           Workspace to use.
 * \param[in ] dw          Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*delta1</code>.  On recoverable error sets <code>*delta1</code> to be
 * <tt>NaN</tt> and returns one of #suzerain_error_status.  On unrecoverable
 * error, additionally calls suzerain_error().
 */
int
suzerain_bl_displacement_thickness(
    const double * coeffs_rhou,
    double * delta1,
    gsl_matrix * dB,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw);

/**
 * Compute the momentum thickness \f$\delta_2\f$ (sometimes written
 * \f$\theta\f$) given the edge location and a B-spline coefficient
 * representation of streamwise momentum \f$\rho u\f$ and velocity \f$u\f$.
 * The method computes
 * \f[
 *  \delta_2 = \int_0^\infty
 *  \frac{\rho u}{\rho_\infty u_\infty} \left(1 - \frac{u}{u_\infty}\right)
 *  \, \mathrm{d}y
 * \f]
 * where the value from the final B-spline collocation point is taken to be "at
 * infinity".  Among many other places, this definition appears in equation
 * (13.47) on page 324 of Leipmann and Roshko's <a
 * href="http://www.worldcat.org/title/elements-of-gasdynamics/oclc/636935705">
 * Elements of Gasdynamics</a>.
 *
 * \param[in ] coeffs_rhou B-spline coefficients for \f$\rho u\f$.
 * \param[in ] coeffs_u    B-spline coefficients for \f$u\f$.
 * \param[out] delta2      The computed momentum thickness.
 * \param[in]  dB          Temporary storage to use of size <tt>w->k</tt> by
 *                         no less than <tt>1</tt>.
 * \param[in ] w           Workspace to use.
 * \param[in ] dw          Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*delta2</code>.  On recoverable error sets <code>*delta2</code> to be
 * <tt>NaN</tt> and returns one of #suzerain_error_status.  On unrecoverable
 * error, additionally calls suzerain_error().
 */
int
suzerain_bl_momentum_thickness(
    const double * coeffs_rhou,
    const double * coeffs_u,
    double * delta2,
    gsl_matrix * dB,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw);

/**
 * Compute the energy thickness \f$\delta_3\f$ (not to be confused with
 * the enthalpy thickness) given the edge location and a B-spline coefficient
 * representation of streamwise momentum \f$\rho u\f$ and velocity \f$u\f$.
 * The method computes
 * \f[
 *  \delta_3 = \int_0^\infty
 *  \frac{\rho u}{\rho_\infty u_\infty}
 *  \left(1 - \left(\frac{u}{u_\infty}\right)^2\right)
 *  \, \mathrm{d}y
 * \f]
 * where the value from the final B-spline collocation point is taken to be "at
 * infinity".  Among many other places, this definition appears in equation
 * (10.97) on page 258 of Schlichting and Gersten's <a
 * href="http://www.worldcat.org/title/boundary-layer-theory-with-22-tables/oclc/615466700">
 * Boundary Layer Theory</a>.
 *
 * \param[in ] coeffs_rhou B-spline coefficients for \f$\rho u\f$.
 * \param[in ] coeffs_u    B-spline coefficients for \f$u\f$.
 * \param[out] delta3      The computed momentum thickness.
 * \param[in]  dB          Temporary storage to use of size <tt>w->k</tt> by
 *                         no less than <tt>1</tt>.
 * \param[in ] w           Workspace to use.
 * \param[in ] dw          Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*delta3</code>.  On recoverable error sets <code>*delta3</code> to be
 * <tt>NaN</tt> and returns one of #suzerain_error_status.  On unrecoverable
 * error, additionally calls suzerain_error().
 */
int
suzerain_bl_energy_thickness(
    const double * coeffs_rhou,
    const double * coeffs_u,
    double * delta3,
    gsl_matrix * dB,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw);

/**
 * Compute the enthalpy thickness \f$\delta_H\f$ (sometimes written
 * \f$\delta_h\f$ and not to be confused with the
 * energy thickness) given the edge location and a B-spline coefficient
 * representation of streamwise momentum \f$\rho u\f$ and specific total
 * enthalpy \f$H_0 = \frac{\rho E + p}{\rho}\f$.
 * The method computes
 * \f[
 *  \delta_H = \int_0^\infty
 *  \frac{\rho u}{\rho_\infty u_\infty}
 *  \left(1 - \frac{H_0}{H_{0,\infty}}\right)
 *  \, \mathrm{d}y
 * \f]
 * where the value from the final B-spline collocation point is taken to be "at
 * infinity".  Among many other places, this definition appears in equation
 * This definition appears in equation (13.48) on page 324 of Leipmann and
 * Roshko's <a
 * href="http://www.worldcat.org/title/elements-of-gasdynamics/oclc/636935705">
 * Elements of Gasdynamics</a> but is there called "energy thickness".
 *
 * \param[in ] coeffs_rhou B-spline coefficients for \f$\rho u\f$.
 * \param[in ] coeffs_H0   B-spline coefficients for \f$H_0\f$.
 * \param[out] deltaH      The computed enthalpy thickness.
 * \param[in]  dB          Temporary storage to use of size <tt>w->k</tt> by
 *                         no less than <tt>1</tt>.
 * \param[in ] w           Workspace to use.
 * \param[in ] dw          Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*deltaH</code>.  On recoverable error sets <code>*deltaH</code> to be
 * <tt>NaN</tt> and returns one of #suzerain_error_status.  On unrecoverable
 * error, additionally calls suzerain_error().
 */
int
suzerain_bl_enthalpy_thickness(
    const double * coeffs_rhou,
    const double * coeffs_H0,
    double * deltaH,
    gsl_matrix * dB,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw);

/**
 * Information characterizing boundary layer thickness in various ways.  Each
 * nondimensional quantity has been divided by \f$l_0\f$.  Using \ref
 * suzerain_bl_compute_thicknesses is the recommended way to populate this
 * data.
 */
typedef struct suzerain_bl_thicknesses {
    double delta;     /**< Boundary layer thickness \f$\delta\f$. */
    double delta1;    /**< Displacement thickness \f$\delta_1\f$
                           (sometimes written \f$\delta^\ast\f$). */
    double delta2;    /**< Momentum thickness \f$\delta_2\f$
                           (sometimes written \f$\theta\f$).      */
    double delta3;    /**< Energy thickness \f$\delta_3\f$ .      */
    double deltaH;    /**< Energy thickness \f$\delta_H\f$ .      */
} suzerain_bl_thicknesses;

/**
 * Compute boundary layer thickness parameters.  Requires a B-spline
 * coefficient representation of specific total enthalpy \f$H_0 = \frac{\rho E
 * + p}{\rho}\f$, streamwise momentum \f$\rho u\f$, and velocity \f$u\f$.  This
 * is a convenience method around \ref suzerain_bl_find_edge, \ref
 * suzerain_bl_displacement_thickness, \ref suzerain_bl_momentum_thickness,
 * \ref suzerain_bl_energy_thickness, and \ref suzerain_bl_enthalpy_thickness
 * packing the results into a \ref suzerain_bl_thicknesses structure.
 *
 * \param[in ] coeffs_H0   Coefficient representation of \f$H_0\f$.
 * \param[in ] coeffs_rhou Coefficient representation of \f$\rho u\f$.
 * \param[in ] coeffs_u    Coefficient representation of \f$u\f$.
 * \param[out] thicknesses Populated on success.
 *                         See type documentation for contents.
 * \param[in ] w           Workspace to use.
 * \param[in ] dw          Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bl_compute_thicknesses(
    const double * coeffs_H0,
    const double * coeffs_rhou,
    const double * coeffs_u,
    suzerain_bl_thicknesses * thick,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw);

/**
 * Boundary layer Reynolds numbers of interest.  These are computed from the
 * indicated length scales and some dynamic viscosity \f$\mu\f$.
 */
typedef struct suzerain_bl_reynolds {
    double delta;     /**< Reynolds number based on boundary layer
                           thickness \f$\delta\f$. */
    double delta1;    /**< Reynolds number based on displacement
                           thickness \f$\delta_1\f$. */
    double delta2;    /**< Reynolds number based on momentum
                           thickness \f$\delta_2\f$. */
    double delta3;    /**< Reynolds number based on energy
                           thickness \f$\delta_3\f$. */
    double deltaH;    /**< Reynolds number based on enthalpy
                           thickness \f$\delta_H\f$. */
} suzerain_bl_reynolds;

/**
 * Compute boundary layer Reynolds numbers.
 *
 * \param[in ] code_Re  Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used to scale
 *                      nondimensional quantities.  For dimensional
 *                      calculations, use <code>1</code>.
 * \param[in ] edge     Local state information from the boundary layer edge.
 * \param[in ] thick    Thickness information for the boundary layer.
 * \param[out] reynolds Populated on success.
 *                      See type documentation for contents.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bl_compute_reynolds(
    double code_Re,
    const suzerain_bl_local       * edge,
    const suzerain_bl_thicknesses * thick,
          suzerain_bl_reynolds    * reynolds);

/**
 * Nondimensional boundary layer quantities of general interest.
 */
typedef struct suzerain_bl_qoi {
    double cf;           /**< The skin friction coefficient \f$c_f =
                              \frac{2 \tau_w}{\rho_e u_e^2}\f$. */
    double gamma_e;      /**< The ratio of specific heats at the edge. */
    double Ma_e;         /**< The local Mach number at the edge. */
    double Ma_tau;       /**< The friction Mach number
                              \f$\mbox{Ma}_\tau = u_\tau / a_w\f$. */
    double Pr_w;         /**< The Prandtl number at the wall. */
    double Bq;           /**< The nondimensional heat flux
                              \f$ B_q = - \frac{k \nabla T}{\rho_w C_p
                              u_\tau T_w} = \frac{\mu \nabla T}{\mbox{Pr}
                              \rho_w C_p u_\tau T_w} \f$. */
    double ratio_rho;    /**< The ratio of edge to wall density. */
    double ratio_nu;     /**< The ratio of edge to wall kinematic viscosity. */
    double ratio_T;      /**< The ratio of edge to wall temperature. */
    double shapefactor;  /**< The shape factor \f$\delta_1 / \delta_2. */
    double v_wallplus;   /**< The wall transpiration rate in plus units.*/
} suzerain_bl_qoi;

/**
 * Compute nondimensional boundary layer parameters.
 *
 * \param[in ] code_Ma Mach number \f$u_0/a_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code>.
 * \param[in ] code_Re Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used to scale
 *                     nondimensional quantities.  For dimensional
 *                     calculations, use <code>1</code>.
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
        const suzerain_bl_local       * wall,
        const suzerain_bl_viscous     * viscous,
        const suzerain_bl_local       * edge,
        const suzerain_bl_thicknesses * thick,
              suzerain_bl_qoi         * qoi);

/**
 * Nondimensional boundary layer pressure gradient parameters.
 *
 * Most of these quantities are defined and discussed within Cal and Castillo's
 * <a href="http://dx.doi.org/10.1063/1.2991433">"Similarity analysis of
 * favorable pressure gradient turbulent boundary layers with eventual
 * quasilaminarization."</a> in Physics of Fluids 20 (2008).
 *
 * One notable exception is \f$p_{e,x}^\ast\f$ which is motivated by the radial
 * nozzle problem implemented in \ref radial_nozzle.h and discussed at length
 * in <tt>writeups/baseflow.tex</tt>.
 */
typedef struct suzerain_bl_pg {
    double Clauser;      /**< The Clauser parameter \f$\beta =
                              \frac{\delta^\ast}{\tau_w}
                              \frac{\partial p}{\partial x}\f$. */
    double Lambda_n;     /**< Pressure parameter \f$\Lambda_n =
                              -\frac{\delta}{\tau_w}
                              \frac{\partial p}{\partial x}\f$
                              which Cal and Castillo attribute to
                              Narasimha and Sreenivasan. */
    double Launder_e;    /**< Launder's acceleration parameter \f$K =
                              \frac{\mu}{\rho_e u_e^2} \,
                              \frac{\partial{}u_e}{\partial x}\f$
                              computed taking \f$\mu = \mu_e\f$. */
    double Launder_w;    /**< Launder's acceleration parameter \f$K =
                              \frac{\mu}{\rho_e u_e^2} \,
                              \frac{\partial{}u_e}{\partial x}\f$
                              computed taking \f$\mu = \mu_w\f$. */
    double p_ex;         /**< The inviscid-friendly pressure parameter
                              \f$p_{e,x}^\ast = \frac{\delta}{\rho_e u_e^2}
                              \frac{\partial p}{\partial x}\f$. */
    double Pohlhausen;   /**< The Pohlhausen parameter \f$K_s =
                              \frac{\delta^2}{\nu_e}
                              \frac{\partial u_e}{\partial x}\f$. */
} suzerain_bl_pg;

/**
 * Compute nondimensional boundary layer pressure gradient parameters.
 * Each nondimensional quantity is scaled as documented.
 *
 * \param[in ] code_Ma   Mach number \f$u_0/a_0\f$ used to scale
 *                       nondimensional quantities.  For dimensional
 *                       calculations, use <code>1</code>.
 * \param[in ] code_Re   Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used to
 *                       scale nondimensional quantities.  For dimensional
 *                       calculations, use <code>1</code>.
 * \param[in ] wall      Local state information from the wall.
 * \param[in ] viscous   Viscous-related wall scaling information
 * \param[in ] edge      Local state information from the boundary layer edge.
 * \param[in ] edge_p__x Streamwise derivative of pressure
 *                       at the boundary layer edge divided by
                         \f$p_0 / l_0 = \rho_0 a_0^2 / l_0\f$.
 * \param[in ] edge_u__x Streamwise derivative of streamwise velocity at the
 *                       boundary layer edge divided by \f$u_0 / l_0\f$.
 * \param[in ] thick     Thickness information for the boundary layer.
 * \param[out] pg        Populated on success.
 *                       See type documentation for contents.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bl_compute_pg(
        double code_Ma,
        double code_Re,
        const suzerain_bl_local       * wall,
        const suzerain_bl_viscous     * viscous,
        const suzerain_bl_local       * edge,
        const double                    edge_p__x,
        const double                    edge_u__x,
        const suzerain_bl_thicknesses * thick,
              suzerain_bl_pg          * pg);

/**
 * Compute Reynolds numbers accounting for an inviscid baseflow.  The results
 * are analogous to the nondimensional parameters computed by \ref
 * suzerain_bl_compute_reynolds.  However, unlike that method, the quantities
 * are correct for integral quantities numbers in the presence of a known,
 * potentially nonuniform, inviscid background flow.
 *
 * The following definitions, derived in writeups/thicknesses.pdf, are used:
 * \f{align}{
 *   \mbox{Re}_{\delta}   &= \frac{\rho_e u_e \delta}{\mu_e}
 * \\
 *   \mbox{Re}_{\delta_1} &= \mu_e^{-1} \int_0^\infty
 *     \left(\rho u\right)_\mbox{inv}
 *   - \left(\rho u\right)_\mbox{vis}
 *   \, \mathrm{d}y
 * \\
 *   \mbox{Re}_{\delta_2} &= \mu_e^{-1} \int_0^\infty
 *     \left(\rho u\right)_\mbox{vis}
 *   - \left(\rho u^2\right)_\mbox{vis} u_\mbox{inv}^{-1}
 *   \, \mathrm{d}y
 * \\
 *   \mbox{Re}_{\delta_3} &= \mu_e^{-1} \int_0^\infty
 *     \left(\rho u\right)_\mbox{vis}
 *   - \left(\rho u^3\right)_\mbox{vis} u^{-2}_\mbox{inv}
 *   \, \mathrm{d}y
 * \\
 *   \mbox{Re}_{\delta_H} &= \mu_e^{-1} \int_0^\infty
 *     \left(\rho u\right)_\mbox{vis}
 *   - \left(\rho H_0 u\right)_\mbox{vis} {H}^{-1}_{0,\mbox{inv}}
 *   \, \mathrm{d}y
 * \f}
 * The method requires a B-spline coefficient representation of both the
 * viscous and inviscid flow specific total enthalpy \f$H_0 = \frac{\rho E +
 * p}{\rho}\f$, streamwise momentum \f$\rho u\f$, and velocity \f$u\f$.
 *
 * \param[in ] code_Re         Reynolds number \f$\rho_0 u_0 l_0/\mu_0\f$ used
 *                             to scale nondimensional quantities.  For
 *                             dimensional calculations, use <code>1</code>.
 * \param[in ] coeffs_vis_H0   Coefficients for viscous \f$H_0\f$ profile.
 *                             These value must be contiguous in memory.
 * \param[in ] coeffs_vis_rhou Coefficients for viscous \f$\rho u\f$ profile.
 *                             These value must be contiguous in memory.
 * \param[in ] coeffs_vis_u    Coefficients for viscous \f$u\f$ profile.
 *                             These value must be contiguous in memory.
 * \param[in ] inv_stride      Stride between inviscid profile coefficients.
 *                             May be zero to indicate uniform inviscid data.
 * \param[in ] coeffs_inv_H0   Coefficients for inviscid \f$H_0\f$ profile.
 *                             These are strided per \c inv_stride.
 * \param[in ] coeffs_inv_rhou Coefficients for inviscid \f$\rho u\f$ profile.
 *                             These are strided per \c inv_stride.
 * \param[in ] coeffs_inv_u    Coefficients for inviscid \f$u\f$ profile.
 *                             These are strided per \c inv_stride.
 * \param[in ] edge            Local information from the boundary layer edge.
 * \param[out] reynolds        Populated on success.
 *                             See type documentation for contents.
 * \param[in ] w               Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error()
 *      and returns one of #suzerain_error_status.  A best effort attempt
 *      is made to compute as many of the Reynolds numbers as possible.
 *      In this circumstance, the #suzerain_error_status returned will
 *      reflect the first problematic integral thickness computation.
 */
int
suzerain_bl_compute_reynolds_baseflow(
    const double                          code_Re,
    const double                  * const coeffs_vis_H0,
    const double                  * const coeffs_vis_rhou,
    const double                  * const coeffs_vis_u,
    const int                             inv_stride,
    const double                  * const coeffs_inv_H0,
    const double                  * const coeffs_inv_rhou,
    const double                  * const coeffs_inv_u,
    const suzerain_bl_local       * const edge,
    suzerain_bl_reynolds          * const reynolds,
    gsl_bspline_workspace         * const w);

/**
 * Compute integral length scales accounting for an inviscid baseflow.  The
 * results are analogous to the nondimensional parameters computed by \ref
 * suzerain_bl_compute_reynolds.  However, unlike that method, the quantities
 * are correct for integral quantities numbers in the presence of a known,
 * potentially nonuniform, inviscid background flow.
 *
 * The integral thicknesses are defined to be the values such that
 * the following equations, derived in writeups/thicknesses.pdf, are satisfied:
 * \f{align}{
 *     \int_{\delta_1}^\infty
 *         \left(\rho u\right)_\mbox{inv}
 *       - \left(\rho u\right)_\mbox{vis}
 *     \, \mathrm{d}y
 *   &=
 *     \int_0^{\delta_1}
 *         \left(\rho u\right)_\mbox{vis}
 *     \, \mathrm{d}y
 * \\
 *     \int_{\delta_1 + \delta_2}^\infty
 *         \left(\rho u^2\right)_\mbox{inv}
 *       - \left(\rho u^2\right)_\mbox{vis}
 *     \, \mathrm{d}y
 *   &=
 *     \int_0^{\delta_1 + \delta_2}
 *         \left(\rho u^2\right)_\mbox{vis}
 *     \, \mathrm{d}y
 * \\
 *     \int_{\delta_1 + \delta_3}^\infty
 *         \left(\rho u^3\right)_\mbox{inv}
 *       - \left(\rho u^3\right)_\mbox{vis}
 *     \, \mathrm{d}y
 *   &=
 *     \int_0^{\delta_1 + \delta_3}
 *         \left(\rho u^3\right)_\mbox{vis}
 *     \, \mathrm{d}y
 * \\
 *     \int_{\delta_1 + \delta_H}^\infty
 *         \left(\rho H_0 u\right)_\mbox{inv}
 *       - \left(\rho H_0 u\right)_\mbox{vis}
 *     \, \mathrm{d}y
 *   &=
 *     \int_0^{\delta_1 + \delta_H}
 *         \left(\rho H_0 u\right)_\mbox{vis}
 *     \, \mathrm{d}y
 * \\
 * \f}
 * The method requires a B-spline coefficient representation of both the
 * viscous and inviscid flow specific total enthalpy \f$H_0 = \frac{\rho E +
 * p}{\rho}\f$, streamwise momentum \f$\rho u\f$, and velocity \f$u\f$.
 *
 * \param[in ] coeffs_vis_H0   Coefficients for viscous \f$H_0\f$ profile.
 *                             These value must be contiguous in memory.
 * \param[in ] coeffs_vis_rhou Coefficients for viscous \f$\rho u\f$ profile.
 *                             These value must be contiguous in memory.
 * \param[in ] coeffs_vis_u    Coefficients for viscous \f$u\f$ profile.
 *                             These value must be contiguous in memory.
 * \param[in ] inv_stride      Stride between inviscid profile coefficients.
 *                             May be zero to indicate uniform inviscid data.
 * \param[in ] coeffs_inv_H0   Coefficients for inviscid \f$H_0\f$ profile.
 *                             These are strided per \c inv_stride.
 * \param[in ] coeffs_inv_rhou Coefficients for inviscid \f$\rho u\f$ profile.
 *                             These are strided per \c inv_stride.
 * \param[in ] coeffs_inv_u    Coefficients for inviscid \f$u\f$ profile.
 *                             These are strided per \c inv_stride.
 * \param[out] thick           Populated on success.
 *                             See type documentation for contents.
 * \param[in ] w               Workspace to use.
 * \param[in ] dw              Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error()
 *      and returns one of #suzerain_error_status.  A best effort attempt
 *      is made to compute as many of the integral thicknesses as possible.
 *      In this circumstance, the #suzerain_error_status returned will
 *      reflect the first problematic integral thickness computation.
 */
int
suzerain_bl_compute_thicknesses_baseflow(
    const double                  * const coeffs_vis_H0,
    const double                  * const coeffs_vis_rhou,
    const double                  * const coeffs_vis_u,
    const int                             inv_stride,
    const double                  * const coeffs_inv_H0,
    const double                  * const coeffs_inv_rhou,
    const double                  * const coeffs_inv_u,
    suzerain_bl_thicknesses       * const thick,
    gsl_bspline_workspace         * const w,
    gsl_bspline_deriv_workspace   * const dw);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SUZERAIN_BL_H */
