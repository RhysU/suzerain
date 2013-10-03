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
    /* Allocate working storage for function evaluation. */
    gsl_matrix *dB = gsl_matrix_alloc(w->k, 3); /* nderiv=2 + 1 */
    if (SUZERAIN_UNLIKELY(dB == NULL)) {
        SUZERAIN_ERROR_NULL("failed to allocate scratch space dB",
                            SUZERAIN_ENOMEM);
    }

    /* Use somewhat high-level crossing first function to find edge */
    double lower = gsl_bspline_breakpoint(0,             w);
    double upper = gsl_bspline_breakpoint(w->nbreak - 1, w);
    const int status = suzerain_bspline_crossing_first(
            /* lower-towards-upper */ 1, 2, coeffs_H0, 0.0, &lower, &upper,
            100, GSL_DBL_EPSILON, GSL_DBL_EPSILON, location, dB, w, dw);

    /* Free working storage and return status */
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
    qoi->p_exi        = thick->delta / edge->rho / square(edge->u) * edge->p__x
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
