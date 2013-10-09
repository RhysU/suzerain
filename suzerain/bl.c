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

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>

#include <suzerain/common.h>
#include <suzerain/bspline.h>
#include <suzerain/error.h>
#include <suzerain/pre_gsl.h>

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

    /* Ensure the basis has at least non-constant derivatives of interest */
    if (SUZERAIN_UNLIKELY(w->k < nderiv + 2)) {
        SUZERAIN_ERROR("Basis must have non-constant second derivatives",
                       SUZERAIN_EINVAL);
    }

    /* Assume failure until proven otherwise. */
    *location  = GSL_NAN;
    int status = SUZERAIN_EFAILED;

    /* Start by evaluating function on 0th breakpoint... */
    double flower = GSL_NAN, lower = gsl_bspline_breakpoint(0, w);
    suzerain_bspline_linear_combination(
            nderiv, coeffs_H0, 1, &lower, &flower, 0, dB, w, dw);

    /* ...look breakpoint-by-breakpoint for upwards crossing of threshold... */
    double fupper = GSL_NAN, upper = GSL_NAN;
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

int
suzerain_bl_compute_deltastar(
    const double edge_location,
    const double * coeffs_rho_u,
    double * deltastar,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    *deltastar = GSL_NAN;           // Be defensive
    if (gsl_isnan(edge_location)) { // Propagate NaN but succeed
        return SUZERAIN_SUCCESS;
    }

    /* Compute edge momentum given edge_location and coeffs_rho_u */
    double rho_u_edge = GSL_NAN;
    int status = suzerain_bspline_linear_combination(
            0, coeffs_rho_u, 1, &edge_location, &rho_u_edge, 0, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) {
        SUZERAIN_ERROR_VAL("failed to compute rho_u_edge",
                           SUZERAIN_EFAILED, status);
    }

    /* Obtain an appropriate order Gauss-Legendre integration rule  */
    /* Integrand (1 - \rho{}u/{\rho{}u}_oo) has same order as basis */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc((w->k + 1)/2);
    if (SUZERAIN_UNLIKELY(tbl == NULL)) {
        SUZERAIN_ERROR("failed to obtain Gauss-Legendre rule from GSL",
                       SUZERAIN_ESANITY);
    }

    /* Notice any failure below will be catastrophic   */
    /* so accumulate integral into *deltastar directly */
    *deltastar = 0;

    /* Accumulate the breakpoint-by-breakpoint contributions to result */
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval with edge being upper limit */
        const double a = gsl_bspline_breakpoint(i,   w);
        if (a > edge_location) break;
        const double b = GSL_MIN_DBL(gsl_bspline_breakpoint(i+1, w),
                                     edge_location);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_deriv_eval_nonzero(xj, 0,
                    dB, &kstart, &kend, w, dw);

            /* Accumulate basis linear combinations to evaluate rho_u */
            double rho_u = 0;
            for (size_t k = kstart; k <= kend; ++k) {
                rho_u += coeffs_rho_u[k] * gsl_matrix_get(dB, k - kstart, 0);
            }

            /* Then use integrand scaled by weight to accumulate result */
            *deltastar += wj * (1 - rho_u / rho_u_edge);
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_compute_theta(
    const double edge_location,
    const double * coeffs_rho_u,
    const double * coeffs_u,
    double * theta,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    *theta = GSL_NAN;               // Be defensive
    if (gsl_isnan(edge_location)) { // Propagate NaN but succeed
        return SUZERAIN_SUCCESS;
    }

    /* Compute edge momentum given edge_location and coeffs_rho_u */
    double rho_u_edge = GSL_NAN;
    int status = suzerain_bspline_linear_combination(
            0, coeffs_rho_u, 1, &edge_location, &rho_u_edge, 0, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) {
        SUZERAIN_ERROR_VAL("failed to compute rho_u_edge",
                           SUZERAIN_EFAILED, status);
    }

    /* Compute edge velocity given edge_location and coeffs_u */
    double u_edge = GSL_NAN;
    status = suzerain_bspline_linear_combination(
            0, coeffs_u, 1, &edge_location, &u_edge, 0, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) {
        SUZERAIN_ERROR_VAL("failed to compute u_edge",
                           SUZERAIN_EFAILED, status);
    }

    /* Obtain an appropriate order Gauss-Legendre integration rule     */
    /* Product (u / u_oo)(1 - \rho{}u/{\rho{}u}_oo) must be integrated */
    /* That is, solve 2*(k - 1) = 2*n - 1 for number of Gauss points n */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc(w->k);
    if (SUZERAIN_UNLIKELY(tbl == NULL)) {
        SUZERAIN_ERROR("failed to obtain Gauss-Legendre rule from GSL",
                       SUZERAIN_ESANITY);
    }

    /* Notice any failure below will be catastrophic */
    /* so accumulate integral into *theta directly   */
    *theta = 0;

    /* Accumulate the breakpoint-by-breakpoint contributions to result */
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval with edge being upper limit */
        const double a = gsl_bspline_breakpoint(i,   w);
        if (a > edge_location) break;
        const double b = GSL_MIN_DBL(gsl_bspline_breakpoint(i+1, w),
                                     edge_location);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_deriv_eval_nonzero(xj, 0,
                    dB, &kstart, &kend, w, dw);

            /* Accumulate basis linear combinations to evaluate rho_u */
            double rho_u = 0;
            double u     = 0;
            for (size_t k = kstart; k <= kend; ++k) {
                const double Bk = gsl_matrix_get(dB, k - kstart, 0);
                rho_u += coeffs_rho_u[k] * Bk;
                u     += coeffs_u[k]     * Bk;
            }

            /* Then use integrand scaled by weight to accumulate result */
            *theta += wj * (u / u_edge) * (1 - rho_u / rho_u_edge);
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    return SUZERAIN_SUCCESS;
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

    int status     = SUZERAIN_SUCCESS;
    gsl_matrix *dB = NULL;

    /* Allocate working storage */
    if (NULL == (dB = gsl_matrix_alloc(w->k, 3))) {
        SUZERAIN_ERROR_REPORT("failed to allocate dB",
                              (status = SUZERAIN_ENOMEM));
    }
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    /* Compute edge location */
    status = suzerain_bl_find_edge(coeffs_H0, &thick->delta, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    /* Compute displacement thickness */
    status = suzerain_bl_compute_deltastar(
            thick->delta, coeffs_rho_u, &thick->deltastar, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    /* Compute momentum thickness */
    status = suzerain_bl_compute_theta(
            thick->delta, coeffs_rho_u, coeffs_u, &thick->theta, dB, w, dw);
    /* Done regardless of status */

done:

    gsl_matrix_free(dB);
    return status;
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
    qoi->Cf           = 2 * viscous->tau_w / edge->rho / square(edge->u)
                      / code_Re;
    qoi->Clauser      = thick->deltastar / viscous->tau_w * edge->p__x
                      * code_Re / square(code_Ma);
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
