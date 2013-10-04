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

/** @file
 * @copydoc bl.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/bl.h>

#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>

#include <suzerain/common.h>
#include <suzerain/bspline.h>
#include <suzerain/error.h>

static inline double square(double x) { return x*x; }

// Defensively NaN assuming p points to a type containing only doubles
#define FILL_WITH_NANS(p)                                    \
    for (size_t i = 0; i < sizeof(*p)/sizeof(double); ++i) { \
        ((double *) p)[i] = GSL_NAN;                         \
    }

int
suzerain_bl_compute_viscous(
        const suzerain_bl_local   * wall,
              suzerain_bl_viscous * viscous)
{
    FILL_WITH_NANS(viscous);

    // Compute dimensional quantities in "code units" each having [units    ]
    viscous->tau_w    = wall->mu * wall->u__y;                 // [\mu u / l]
    viscous->u_tau    = sqrt(viscous->tau_w / wall->rho);      // [u        ]
    viscous->delta_nu = wall->mu / wall->rho / viscous->u_tau; // [l        ]

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_find_edge(
    const double * coeffs_H0,
    double * location,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Everything hinges on the second derivative of H0 crossing 0.0 */
    enum { nderiv = 2, threshold = 0 };

    /* Assume failure until proven otherwise. */
    *location = GSL_NAN;
    int status = SUZERAIN_EFAILED;

    /* Start by evaluating function on 0th breakpoint... */
    double flower, lower = gsl_bspline_breakpoint(0, w);
    suzerain_bspline_linear_combination(
            nderiv, coeffs_H0, 1, &lower, &flower, 0, dB, w, dw);

    /* ...look breakpoint-by-breakpoint for upwards crossing of threshold... */
    double fupper, upper;
    for (size_t i = 1; i < w->nbreak; ++i, lower = upper, flower = fupper) {
        upper = gsl_bspline_breakpoint(i, w);
        suzerain_bspline_linear_combination(
                nderiv, coeffs_H0, 1, &upper, &fupper, 0, dB, w, dw);

        /* ...when found, polish the crossing bracket into a location... */
        if (flower < threshold && fupper >= threshold) {
            status = suzerain_bspline_crossing(
                    nderiv, coeffs_H0, threshold, &lower, &upper, 100,
                    GSL_DBL_EPSILON, GSL_DBL_EPSILON, location, dB, w, dw);
            break;
        }

    }
    /* ...if never found, location remains NAN and status reflects failure. */

    return status;
}

// Parameters necessary for deltastar_function.
typedef struct {
    const double                *coeffs_rho_u;
    gsl_matrix                  *dB;
    gsl_bspline_workspace       *w;
    gsl_bspline_deriv_workspace *dw;
    double                       rho_u_edge;
} deltastar_params;

// A GSL-ready way to evaluate the displacement thickness integrand
static
double deltastar_function(double x, void * params)
{
    deltastar_params * p = (deltastar_params *) params;
    double rho_u = GSL_NAN;
    suzerain_bspline_linear_combination(
            0, p->coeffs_rho_u, 1U, &x, &rho_u, 0U, p->dB, p->w, p->dw);
    return (1 - rho_u / p->rho_u_edge);
}

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
    double *abserr)
{
    *deltastar = GSL_NAN;           // Be defensive
    if (gsl_isnan(edge_location)) { // Propagate NaN but succeed
        return SUZERAIN_SUCCESS;
    }

    // Compute edge momentum given edge_location and coeffs_rho_u
    double rho_u_edge = GSL_NAN;
    int status = suzerain_bspline_linear_combination(
            0, coeffs_rho_u, 1, &edge_location, &rho_u_edge, 0, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) {
        SUZERAIN_ERROR_VAL("failed to compute rho_u_edge",
                           SUZERAIN_EFAILED, status);
    }

    // Integrate to obtain displacement thickness
    deltastar_params params = { coeffs_rho_u, dB, w, dw, rho_u_edge };
    gsl_function f          = { deltastar_function, &params };
    status = gsl_integration_qag(&f, gsl_bspline_breakpoint(0, w),
            edge_location, epsabs, epsrel, iw->limit, GSL_INTEG_GAUSS61, iw,
            deltastar, abserr);
    if (SUZERAIN_UNLIKELY(status != GSL_SUCCESS)) {
        SUZERAIN_ERROR("integration of deltastar failed", status);
    }
    return status;
}

// Parameters necessary for theta_function.
typedef struct {
    const double                *coeffs_rho_u;
    const double                *coeffs_u;
    gsl_matrix                  *dB;
    gsl_bspline_workspace       *w;
    gsl_bspline_deriv_workspace *dw;
    double                       rho_u_edge;
    double                       u_edge;
} theta_params;

// A GSL-ready way to evaluate the momentum thickness integrand
static
double theta_function(double x, void * params)
{
    theta_params * p = (theta_params *) params;
    double rho_u = GSL_NAN;
    suzerain_bspline_linear_combination(
            0, p->coeffs_rho_u, 1U, &x, &rho_u, 0U, p->dB, p->w, p->dw);
    double u = GSL_NAN;
    suzerain_bspline_linear_combination(
            0, p->coeffs_u, 1U, &x, &u, 0U, p->dB, p->w, p->dw);
    return (rho_u / p->rho_u_edge) * (1 - u / p->u_edge);
}

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
    double *abserr)
{
    *theta = GSL_NAN;               // Be defensive
    if (gsl_isnan(edge_location)) { // Propagate NaN but succeed
        return SUZERAIN_SUCCESS;
    }

    // Compute edge momentum given edge_location and coeffs_rho_u
    double rho_u_edge = GSL_NAN;
    int status = suzerain_bspline_linear_combination(
            0, coeffs_rho_u, 1, &edge_location, &rho_u_edge, 0, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) {
        SUZERAIN_ERROR_VAL("failed to compute rho_u_edge",
                           SUZERAIN_EFAILED, status);
    }

    // Compute edge velocity given edge_location and coeffs_rho_u
    double u_edge = GSL_NAN;
    status = suzerain_bspline_linear_combination(
            0, coeffs_u, 1, &edge_location, &u_edge, 0, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) {
        SUZERAIN_ERROR_VAL("failed to compute u_edge",
                           SUZERAIN_EFAILED, status);
    }

    // Integrate to obtain momentum thickness
    theta_params params = { coeffs_rho_u, coeffs_u,
                            dB, w, dw, rho_u_edge, u_edge };
    gsl_function f      = { theta_function, &params };
    status = gsl_integration_qag(&f, gsl_bspline_breakpoint(0, w),
            edge_location, epsabs, epsrel, iw->limit, GSL_INTEG_GAUSS61, iw,
            theta, abserr);
    if (SUZERAIN_UNLIKELY(status != GSL_SUCCESS)) {
        SUZERAIN_ERROR("integration of deltastar failed", status);
    }
    return status;
}

int
suzerain_bl_compute_thick(
    const double * coeffs_H0,
    const double * coeffs_rho_u,
    const double * coeffs_u,
    suzerain_bl_thick * thick,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    FILL_WITH_NANS(thick);

    return SUZERAIN_EUNIMPL; // FIXME Ticket #2963
}

int
suzerain_bl_compute_qoi(
        const double code_Ma,
        const double code_Re,
        const suzerain_bl_local   * const wall,
        const suzerain_bl_viscous * const viscous,
        const suzerain_bl_local   * const edge,
        const suzerain_bl_thick   * const thick,
              suzerain_bl_qoi     * const qoi)
{
    FILL_WITH_NANS(qoi);

    // Nondimensional quantities are computed with the first line being the
    // quantity and the second line being any needed "code unit" correction.
    qoi->beta         = thick->deltastar / viscous->tau_w * edge->p__x
                      * code_Re / square(code_Ma);
    qoi->Cf           = 2 * viscous->tau_w / edge->rho / square(edge->u)
                      / code_Re;
    qoi->gamma_e      = edge->gamma
                      * 1;
    qoi->K_e          = edge->mu * edge->u__x / edge->rho / square(edge->u)
                      / code_Re;
    qoi->K_s          = square(thick->delta) * edge->rho / edge->mu * edge->u__x
                      * code_Re;
    qoi->K_w          = wall->mu * edge->u__x / edge->rho / square(edge->u)
                      / code_Re;
    qoi->Lambda_n     = - thick->delta / viscous->tau_w * edge->p__x
                      * code_Re / square(code_Ma);
    qoi->Ma_e         = edge->u / edge->a
                      * code_Ma;
    qoi->p_ex         = thick->delta / edge->rho / square(edge->u) * edge->p__x
                      / square(code_Ma);
    qoi->Pr_w         = wall->Pr;
    qoi->Re_delta     = edge->rho * edge->u * thick->delta     / edge->mu
                      * code_Re;
    qoi->Re_deltastar = edge->rho * edge->u * thick->deltastar / edge->mu
                      * code_Re;
    qoi->Re_theta     = edge->rho * edge->u * thick->theta     / edge->mu
                      * code_Re;
    qoi->ratio_nu     = (edge->mu / edge->rho) / (wall->mu / wall->rho)
                      * 1;
    qoi->ratio_rho    = edge->rho / wall->rho
                      * 1;
    qoi->ratio_T      = edge->T / wall->T
                      * 1;
    qoi->shapefactor  = thick->deltastar / thick->theta
                      * 1;
    qoi->v_wallplus   = wall->v / viscous->u_tau
                      * 1;

    return SUZERAIN_SUCCESS;
}
