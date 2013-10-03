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

#include <gsl/gsl_matrix.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

static inline double square(double x) { return x*x; }

int
suzerain_bl_compute_viscous(
        const suzerain_bl_local   * wall,
              suzerain_bl_viscous * viscous)
{
    {   // Defensively NaN output assuming suzerain_bl_viscous is all doubles
        double * const p = (double *) viscous;
        const size_t N = sizeof(*viscous)/sizeof(double);
        for (size_t i = 0; i < N; ++i) p[i] = INFINITY / INFINITY;
    }

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
    gsl_bspline_workspace *w,
    gsl_bspline_deriv_workspace *dw)
{
    /* Everything hinges on the second derivative of H0 crossing 0.0 */
    enum { nderiv = 2, threshold = 0 };

    /* Allocate working storage for function evaluation. */
    gsl_matrix *dB = gsl_matrix_alloc(w->k, nderiv + 1);
    if (SUZERAIN_UNLIKELY(dB == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate scratch space dB",
                            SUZERAIN_ENOMEM);
    }

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

    /* Free working storage and return status */
    gsl_matrix_free(dB);
    return status;
}

int
suzerain_bl_compute_deltastar(
    const double edge_location,
    const double * coeffs_rho_u,
    double * deltastar,
    gsl_bspline_workspace *w)
{
    return SUZERAIN_EUNIMPL; // FIXME
}

int
suzerain_bl_compute_theta(
    const double edge_location,
    const double * coeffs_rho_u,
    const double * coeffs_u,
    double * theta,
    gsl_bspline_workspace *w)
{
    return SUZERAIN_EUNIMPL; // FIXME
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
    {   // Defensively NaN output assuming suzerain_bl_qoi is all doubles
        double * const p = (double *) qoi;
        const size_t N = sizeof(*qoi)/sizeof(double);
        for (size_t i = 0; i < N; ++i) p[i] = INFINITY / INFINITY;
    }

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
