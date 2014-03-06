//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
 * @copydoc channel.h
 */

#include <suzerain/channel.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

static inline double square(double x) { return x*x; }

// Defensively NaN assuming p points to a type containing only doubles
#define FILL_WITH_NANS(p)                                    \
    for (size_t i = 0; i < sizeof(*p)/sizeof(double); ++i) { \
        ((double *) p)[i] = GSL_NAN;                         \
    }

int
suzerain_channel_compute_viscous(
        const double code_Re,
        const suzerain_channel_local   * wall,
              suzerain_channel_viscous * viscous)
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
suzerain_channel_compute_qoi(
        const double code_Ma,
        const double code_Re,
        const suzerain_channel_local       * const wall,
        const suzerain_channel_viscous     * const viscous,
        const suzerain_channel_local       * const center,
              suzerain_channel_qoi         * const qoi)
{
    FILL_WITH_NANS(qoi);

    // Nondimensional quantities are computed with the first line being the
    // quantity and the final line being any needed "code unit" correction.
    // Notice viscous->tau_w and viscous->u_tau already account for code_Re;
    // see the suzerain_channel_viscous struct declaration to check their scaling.
    qoi->cf          = 2 * viscous->tau_w / center->rho / square(center->u)
                     * 1;
    qoi->gamma_c     = center->gamma
                     * 1;
    qoi->Ma_c        = center->u / center->a
                     * code_Ma;
    qoi->Ma_tau      = viscous->u_tau / wall->a
                     * code_Ma;
    qoi->Pr_w        = wall->Pr;
    qoi->Bq          = - (wall->mu * wall->T__y)
                     / (qoi->Pr_w * wall->rho * viscous->u_tau * wall->T)
                     / code_Re;
    qoi->ratio_nu    = (center->mu / center->rho) / (wall->mu / wall->rho)
                     * 1;
    qoi->ratio_rho   = center->rho / wall->rho
                     * 1;
    qoi->ratio_T     = center->T / wall->T
                     * 1;
    qoi->Re_c        = center->rho * center->u * center->y / center->mu
                     * code_Re;
    qoi->Re_tau      = center->y / viscous->delta_nu
                     * 1;
    qoi->v_wallplus  = wall->v / viscous->u_tau
                     * 1;

    return SUZERAIN_SUCCESS;
}
