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

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

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

// TODO Find reference discussing NASA's technique mentioned below.
/**
 * Find the boundary layer edge within <tt>[lower, upper]</tt> given a B-spline
 * coefficient representation of the specific total enthalpy \f$H_0 =
 * \frac{\rho E + p}{\rho}\f$ in \c coeffs_H0.  The procedure looks for the
 * second derivative of \f$H_0 \approx{} 0\f$ based on discussions with T.
 * Oliver about what NASA has done in practice to get robust results.
 *
 * \param[in ] coeffs_H0 Coefficient representation of \f$H_0\f$
 *                       using the basis provided in \c w and \c dw.
 * \param[out] location  Location at which edge is detected.
 * \param[in]  dB        Temporary storage to use of size <tt>w->k</tt> by
 *                       no less than <tt>3</tt>.
 * \param[in]  w         Workspace to use.
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
    double * location,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw);

/**
 * Compute the displacement thickness \f$\delta^\ast\f$ given the edge
 * location and a B-spline coefficient representation of streamwise
 * momentum \f$\rho u\f$.  The method computes
 * \f[
 *  \delta^\ast = \int_0^\infty
 *  \left(1 - \frac{\rho u}{\rho_e \u_e}\right)
 *  \, \mathrm{d}y
 *  .
 * \f]
 *
 * The method propagates a <tt>NaN</tt> \c edge_location into a <tt>NaN</tt>
 * \c deltastar result considering the computation to be successful
 * when this happens.
 *
 * \param[in ] edge_location Location of the boundary layer edge possibly
 *                           computed by suzerain_bl_find_edge().
 * \param[in ] coeffs_rho_u  B-spline coefficients for \f$rho u\f$.
 * \param[out] deltastar     The computed displacement thickness.
 * \param[in]  dB            Temporary storage to use of size <tt>w->k</tt> by
 *                           no less than <tt>1</tt>.
 * \param[in ] w             Workspace to use.
 * \param[in ] dw            Workspace to use.
 * \param[in ] iw            Workspace to use.  Result precision may be limited
 *                           by the \c n parameter passed to
 *                           \c gsl_integration_workspace_alloc.
 * \param[in ] epsabs        Absolute error limit for adaptive integration.
 * \param[in ] epsrel        Relative error limit for adaptive integration.
 * \param[out] abserr        Estimate of the absolute error in the result.
 *
 * \return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*deltastar</code>.  On recoverable error sets <code>*deltastar</code>
 * to be <tt>NaN</tt> <code>*deltastar</code> to be <tt>NaN</tt> and returns
 * one of #suzerain_error_status.  On unrecoverable error, additionally calls
 * suzerain_error().
 */
int
suzerain_bl_compute_deltastar(
    const double edge_location,
    const double * coeffs_rho_u,
    double * deltastar,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw,
    gsl_integration_workspace *iw,
    const double epsabs,
    const double epsrel,
    double *abserr);

/**
 * Compute the momentum thickness \f$\theta\f$ given the edge location and a
 * B-spline coefficient representation of streamwise momentum \f$\rho u\f$ and
 * velocity \f$u\f$.  The method computes
 * \f[
 *  \theta = \int_0^\infty
 *  \frac{\rho u}{\rho_e u_e} \left(1 - \frac{u}{\u_e}\right)
 *  \, \mathrm{d}y
 *  .
 * \f]
 *
 * The method propagates a <tt>NaN</tt> \c edge_location into a <tt>NaN</tt>
 * \c theta result considering the computation to be successful
 * when this happens.
 *
 * \param[in ] edge_location Location of the boundary layer edge possibly
 *                           computed by suzerain_bl_find_edge().
 * \param[in ] coeffs_rho_u  B-spline coefficients for \f$rho u\f$.
 * \param[in ] coeffs_u      B-spline coefficients for \f$u\f$.
 * \param[out] theta         The computed momentum thickness.
 * \param[in]  dB            Temporary storage to use of size <tt>w->k</tt> by
 *                           no less than <tt>1</tt>.
 * \param[in ] w             Workspace to use.
 * \param[in ] dw            Workspace to use.
 * \param[in ] iw            Workspace to use.  Result precision may be limited
 *                           by the \c n parameter passed to
 *                           \c gsl_integration_workspace_alloc.
 * \param[in ] epsabs        Absolute error limit for adaptive integration.
 * \param[in ] epsrel        Relative error limit for adaptive integration.
 * \param[out] abserr        Estimate of the absolute error in the result.
 *
 * \return ::SUZERAIN_SUCCESS on success and returns the answer in
 * <code>*theta</code>.  On recoverable error sets <code>*theta</code> to be
 * <tt>NaN</tt> <code>*theta</code> to be <tt>NaN</tt> and returns one of
 * #suzerain_error_status.  On unrecoverable error, additionally calls
 * suzerain_error().
 */
int
suzerain_bl_compute_theta(
    const double edge_location,
    const double * coeffs_rho_u,
    const double * coeffs_u,
    double * theta,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw,
    gsl_integration_workspace *iw,
    const double epsabs,
    const double epsrel,
    double *abserr);

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
 * Compute boundary layer thickness parameters.  Requires a B-spline
 * coefficient representation of specific total enthalpy \f$H_0 = \frac{\rho E
 * + p}{\rho}\f$, streamwise momentum \f$\rho u\f$, and velocity \f$u\f$.  This
 * is a convenience method around \ref suzerain_bl_find_edge \ref
 * suzerain_bl_compute_deltastar and \ref suzerain_bl_compute_theta packing the
 * results into a \ref suzerain_bl_thick structure.
 *
 * \param[in ] coeffs_H0    Coefficient representation of \f$H_0\f$.
 * \param[in ] coeffs_rho_u Coefficient representation of \f$\rho u\f$.
 * \param[in ] coeffs_u     Coefficient representation of \f$u\f$.
 * \param[out] thick        Populated on success.
 *                          See type documentation for contents.
 * \param[in ] w            Workspace to use.
 * \param[in ] dw           Workspace to use.
 *
 * \return ::SUZERAIN_SUCCESS on success.  On error calls suzerain_error() and
 *      returns one of #suzerain_error_status.
 */
int
suzerain_bl_compute_thick(
    const double * coeffs_H0,
    const double * coeffs_rho_u,
    const double * coeffs_u,
    suzerain_bl_thick * thick,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw);

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
                              \frac{\partial p}{\partial x}\f$. */
    double Cf;           /**< The skin friction coefficient \f$C_f =
                              \frac{2 \tau_w}{\rho_e u_e^2}\f$. */
    double gamma_e;      /**< The ratio of specific heats at the edge. */
    double K_e;          /**< Launder's acceleration parameter \f$K =
                              \frac{\mu}{\rho_e u_e^2} \,
                              \frac{\partial{}u_e}{\partial x}\f$
                              computed taking \f$\mu = \mu_e\f$. */
    double K_s;          /**< The Pohlhausen parameter \f$K_s =
                              \frac{\delta^2}{\nu_e}
                              \frac{\partial u_e}{\partial x}\f$. */
    double K_w;          /**< Launder's acceleration parameter \f$K =
                              \frac{\mu}{\rho_e u_e^2} \,
                              \frac{\partial{}u_e}{\partial x}\f$
                              computed taking \f$\mu = \mu_w\f$. */
    double Lambda_n;     /**< Pressure parameter \f$\Lambda_n =
                              -\frac{\delta}{\tau_w}
                              \frac{\partial p}{\partial x}\f$. */
    double Ma_e;         /**< The local Mach number at the edge. */
    double p_ex;         /**< The inviscid-friendly pressure parameter
                              \f$p_{e,x}^\ast = \frac{\delta}{\rho_e u_e^2}
                              \frac{\partial p}{\partial x}\f$. */
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
