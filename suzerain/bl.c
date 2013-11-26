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
        const double code_Re,
        const suzerain_bl_local   * wall,
              suzerain_bl_viscous * viscous)
{
    FILL_WITH_NANS(viscous);

    // Compute quantities of interest in "code units" divided by  [units    ]
    viscous->tau_w    = wall->mu * wall->u__y;                 // [\mu u / l]
    viscous->u_tau    = sqrt(viscous->tau_w / wall->rho);      // [u        ]
    viscous->delta_nu = wall->mu / wall->rho / viscous->u_tau; // [l        ]

    // Adjust for code Reynolds number per writeups/viscous_nondim.pdf which
    // details why and how how these scaling factors appear.  Provided that
    // they are addressed here, they don't "leak" into other computations.
    const double sqrt_code_Re = sqrt(code_Re);
    viscous->tau_w    /= code_Re;
    viscous->u_tau    /= sqrt_code_Re;
    viscous->delta_nu /= sqrt_code_Re;

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
            const double epsabs = 10*GSL_DBL_EPSILON;
            status = suzerain_bspline_crossing(
                    nderiv, coeffs_H0, threshold, &lower, &upper, 100,
                    epsabs, /*not relative*/ 0, location, dB, w, dw);
            break;
        }

    }
    /* ...if never found, location remains NAN and status reflects failure. */

    return status;
}

int
suzerain_bl_displacement_thickness(
    const double * coeffs_rho_u,
    double * delta1,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Obtain an appropriate order Gauss-Legendre integration rule */
    /* when integrand has same order as basis                      */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc((w->k + 1)/2);
    if (SUZERAIN_UNLIKELY(tbl == NULL)) {
        *delta1 = GSL_NAN;
        SUZERAIN_ERROR("failed to obtain Gauss-Legendre rule from GSL",
                       SUZERAIN_ESANITY);
    }

    /* Momentum at infinity taken from final B_spline collocation point */
    /* which happens to be the value of the final coefficient */
    const double rho_u_edge = coeffs_rho_u[w->n - 1];

    /* Accumulate the breakpoint-by-breakpoint contributions into *delta1 */
    *delta1 = 0;
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval for integration */
        const double a = gsl_bspline_breakpoint(i,   w);
        const double b = gsl_bspline_breakpoint(i+1, w);

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
            *delta1 += wj * (1 - rho_u / rho_u_edge);
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_momentum_thickness(
    const double * coeffs_rho_u,
    const double * coeffs_u,
    double * delta2,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Obtain an appropriate order Gauss-Legendre integration rule     */
    /* when integrand has twice the order of the basis.                */
    /* That is, solve 2*(k - 1) = 2*n - 1 for number of Gauss points n */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc(w->k);
    if (SUZERAIN_UNLIKELY(tbl == NULL)) {
        *delta2 = GSL_NAN;
        SUZERAIN_ERROR("failed to obtain Gauss-Legendre rule from GSL",
                       SUZERAIN_ESANITY);
    }

    /* State at infinity taken from final B_spline collocation point */
    /* which happens to be the value of the final coefficient */
    const double rho_u_edge = coeffs_rho_u[w->n - 1];
    const double     u_edge = coeffs_u    [w->n - 1];

    /* Accumulate the breakpoint-by-breakpoint contributions into *delta2 */
    *delta2 = 0;
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval for integration */
        const double a = gsl_bspline_breakpoint(i,   w);
        const double b = gsl_bspline_breakpoint(i+1, w);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_deriv_eval_nonzero(xj, 0,
                    dB, &kstart, &kend, w, dw);

            /* Accumulate basis linear combinations to evaluate rho_u, u */
            double rho_u = 0;
            double u     = 0;
            for (size_t k = kstart; k <= kend; ++k) {
                const double Bk = gsl_matrix_get(dB, k - kstart, 0);
                rho_u += coeffs_rho_u[k] * Bk;
                u     += coeffs_u[k]     * Bk;
            }

            /* Then use integrand scaled by weight to accumulate result */
            *delta2 += wj * (rho_u / rho_u_edge) * (1 - u / u_edge);
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_energy_thickness(
    const double * coeffs_rho_u,
    const double * coeffs_u,
    double * delta3,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Obtain an appropriate order Gauss-Legendre integration rule     */
    /* when integrand has three times the order of the basis.          */
    /* That is, solve 3*(k - 1) = 2*n - 1 for number of Gauss points n */
    gsl_integration_glfixed_table * const tbl
        = gsl_integration_glfixed_table_alloc((3*w->k+1)/2 - 1);
    if (SUZERAIN_UNLIKELY(tbl == NULL)) {
        *delta3 = GSL_NAN;
        SUZERAIN_ERROR("failed to obtain Gauss-Legendre rule from GSL",
                       SUZERAIN_ESANITY);
    }

    /* State at infinity taken from final B_spline collocation point */
    /* which happens to be the value of the final coefficient */
    const double rho_u_edge = coeffs_rho_u[w->n - 1];
    const double     u_edge = coeffs_u    [w->n - 1];

    /* Accumulate the breakpoint-by-breakpoint contributions into *delta3 */
    *delta3 = 0;
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

        /* Determine i-th breakpoint interval for integration */
        const double a = gsl_bspline_breakpoint(i,   w);
        const double b = gsl_bspline_breakpoint(i+1, w);

        for (size_t j = 0; j < tbl->n; ++j) {

            /* Get j-th Gauss point xj and weight wj */
            gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

            /* Evaluate basis functions at point xj */
            size_t kstart, kend;
            gsl_bspline_deriv_eval_nonzero(xj, 0,
                    dB, &kstart, &kend, w, dw);

            /* Accumulate basis linear combinations to evaluate rho_u, u */
            double rho_u = 0;
            double u     = 0;
            for (size_t k = kstart; k <= kend; ++k) {
                const double Bk = gsl_matrix_get(dB, k - kstart, 0);
                rho_u += coeffs_rho_u[k] * Bk;
                u     += coeffs_u[k]     * Bk;
            }

            /* Then use integrand scaled by weight to accumulate result */
            *delta3 += wj * (rho_u / rho_u_edge) * (1 - square(u / u_edge));
        }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_enthalpy_thickness(
    const double * coeffs_rho_u,
    const double * coeffs_H0,
    double * deltaH,
    gsl_matrix *dB,
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    // The form of the enthalpy thickness equation is nothing but
    // the displacement thickness with H_0 replacing u.
    return suzerain_bl_momentum_thickness(
            coeffs_rho_u, coeffs_H0, deltaH, dB, w, dw);
}

int
suzerain_bl_compute_thicknesses(
    const double * coeffs_H0,
    const double * coeffs_rho_u,
    const double * coeffs_u,
    suzerain_bl_thicknesses * thick,
    gsl_bspline_workspace * w,
    gsl_bspline_deriv_workspace * dw)
{

    FILL_WITH_NANS(thick);

    int status     = SUZERAIN_SUCCESS;
    gsl_matrix *dB = NULL;

    if (NULL == (dB = gsl_matrix_alloc(w->k, 3))) {
        SUZERAIN_ERROR_REPORT("failed to allocate dB",
                              (status = SUZERAIN_ENOMEM));
    }
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_find_edge(
            coeffs_H0, &thick->delta, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_displacement_thickness(
            coeffs_rho_u, &thick->delta1, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_momentum_thickness(
            coeffs_rho_u, coeffs_u, &thick->delta2, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_energy_thickness(
            coeffs_rho_u, coeffs_u, &thick->delta3, dB, w, dw);
    if (SUZERAIN_UNLIKELY(status != SUZERAIN_SUCCESS)) goto done;

    status = suzerain_bl_enthalpy_thickness(
            coeffs_rho_u, coeffs_H0, &thick->deltaH, dB, w, dw);
    /* Done regardless of status */

done:

    gsl_matrix_free(dB);
    return status;
}

int
suzerain_bl_compute_qoi(
        const double code_Ma,
        const double code_Re,
        const suzerain_bl_local       * const wall,
        const suzerain_bl_viscous     * const viscous,
        const suzerain_bl_local       * const edge,
        const suzerain_bl_thicknesses * const thick,
              suzerain_bl_qoi         * const qoi)
{
    FILL_WITH_NANS(qoi);

    // Nondimensional quantities are computed with the first line being the
    // quantity and the final line being any needed "code unit" correction.
    // Notice viscous->tau_w and viscous->u_tau already account for code_Re;
    // see the suzerain_bl_viscous struct declaration to check their scaling.
    qoi->cf          = 2 * viscous->tau_w / edge->rho / square(edge->u)
                     * 1;
    qoi->gamma_e     = edge->gamma
                     * 1;
    qoi->Ma_e        = edge->u / edge->a
                     * code_Ma;
    qoi->Pr_w        = wall->Pr;
    qoi->neg_Bq      = (wall->mu * wall->T__y)
                     / (qoi->Pr_w * wall->rho * viscous->u_tau * wall->T)
                     / code_Re;
    qoi->ratio_nu    = (edge->mu / edge->rho) / (wall->mu / wall->rho)
                     * 1;
    qoi->ratio_rho   = edge->rho / wall->rho
                     * 1;
    qoi->ratio_T     = edge->T / wall->T
                     * 1;
    qoi->Re_delta    = edge->rho * edge->u * thick->delta  / edge->mu
                     * code_Re;
    qoi->Re_delta1   = edge->rho * edge->u * thick->delta1 / edge->mu
                     * code_Re;
    qoi->Re_delta2   = edge->rho * edge->u * thick->delta2 / edge->mu
                     * code_Re;
    qoi->shapefactor = thick->delta1 / thick->delta2
                     * 1;
    qoi->v_wallplus  = wall->v / viscous->u_tau
                     * 1;

    return SUZERAIN_SUCCESS;
}

int
suzerain_bl_compute_pg(
        const double code_Ma,
        const double code_Re,
        const suzerain_bl_local       * const wall,
        const suzerain_bl_viscous     * const viscous,
        const suzerain_bl_local       * const edge,
        const double                          edge_p__x,
        const double                          edge_u__x,
        const suzerain_bl_thicknesses * const thick,
              suzerain_bl_pg          * const pg)
{
    FILL_WITH_NANS(pg);

    // Nondimensional quantities are computed with the first line being the
    // quantity and the second line being any needed "code unit" correction.
    // Notice viscous->tau_w and viscous->u_tau already account for code_Re;
    // see the suzerain_bl_viscous struct declaration to check their scaling.
    pg->Clauser      = thick->delta1 / viscous->tau_w * edge_p__x
                     / square(code_Ma);
    pg->Launder_e    = edge->mu * edge_u__x / edge->rho / square(edge->u)
                     / code_Re;
    pg->Launder_w    = wall->mu * edge_u__x / edge->rho / square(edge->u)
                     / code_Re;
    pg->Lambda_n     = - thick->delta / viscous->tau_w * edge_p__x
                     / square(code_Ma);
    pg->p_ex         = thick->delta / edge->rho / square(edge->u) * edge_p__x
                     / square(code_Ma);
    pg->Pohlhausen   = square(thick->delta) * edge->rho / edge->mu * edge_u__x
                     * code_Re;

    return SUZERAIN_SUCCESS;
}
