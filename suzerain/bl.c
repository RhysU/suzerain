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

#include <suzerain/common.h>
#include <suzerain/error.h>

static inline double square(double x) { return x*x; }

int
suzerain_bl_compute_qoi(
//      const double code_Ma,
//      const double code_Re,
        const suzerain_bl_local * const wall,
        const suzerain_bl_local * const edge,
        const suzerain_bl_thick * const thick,
              suzerain_bl_qoi   * const qoi)
{
    {   // Defensively NaN output assuming suzerain_bl_qoi is all doubles
        double * const p = (double *) qoi;
        const size_t N = sizeof(*qoi)/sizeof(double);
        for (size_t i = 0; i < N; ++i) p[i] = INFINITY / INFINITY;
    }

    // Compute any necessary information redundantly contained in wall or edge.
    const double edge_nu = edge->mu / edge->rho;

    // Compute dimensional quantities in "code units" each having [units]
    const double tau_w    = wall->mu * wall->u__y;        // [\mu_0 u_0 / l_0]
    const double u_tau    = sqrt(tau_w / wall->rho);      // [u_0]
    const double delta_nu = wall->mu / wall->rho / u_tau; // [l_0]

    qoi->beta         = thick->deltastar / tau_w * edge->p__x;
    qoi->Cf           = 2 * tau_w / edge->rho / square(edge->u);
    qoi->gamma_e      = edge->gamma;
    qoi->K_e          = edge->mu * edge->u__x / edge->rho / square(edge->u);
    qoi->K_s          = square(thick->delta) / edge_nu * edge->u__x;
    qoi->K_w          = wall->mu * edge->u__x / edge->rho / square(edge->u);
    qoi->Lambda_n     = - thick->delta / tau_w * edge->p__x;
    qoi->Ma_e         = edge->u / edge->a;
    qoi->p_exi        = thick->delta / edge->rho / square(edge->u) * edge->p__x;
    qoi->Pr_w         = wall->Pr;
    qoi->Re_delta     = edge->rho * edge->u * thick->delta     / edge->mu;
    qoi->Re_deltastar = edge->rho * edge->u * thick->deltastar / edge->mu;
    qoi->Re_theta     = edge->rho * edge->u * thick->theta     / edge->mu;
    qoi->shapefactor  = thick->deltastar / thick->theta;
    qoi->T_ratio      = edge->T / wall->T;
    qoi->v_wallplus   = wall->v / u_tau;

    return SUZERAIN_SUCCESS;
}
