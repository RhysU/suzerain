/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013 Rhys Ulerich
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

/** @file
 * @copydoc radial_nozzle.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/radial_nozzle.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>

// Parameters for \ref nozzle per conventions of gsl_odeiv2_system->function
typedef struct params_type {
    double Ma02;   // Ma0^2
    double gam0m1; // gamma_0 - 1
} params_type;

// Compute u' and a^2 given r, u, Ma02 = Ma0**2, and gam0m1 = gam0 - 1
static inline
void nozzle_upa2(const double R,
                 const double u,
                 const double Ma02,
                 const double gam0m1,
                 double *up,
                 double *a2)
{
    *up = -(u/R) * (2 + Ma02*gam0m1 - Ma02*(gam0m1  )*(u*u))
                 / (2 + Ma02*gam0m1 - Ma02*(gam0m1+2)*(u*u));
    *a2 = 1 + 0.5*Ma02*gam0m1*(1 - u*u);
}

// Find [u; log rho; p]' given r, x=[u; log rho; p], Ma02=Ma0**2, gam0m1=gam0-1
static
int
nozzle_f(double R,
         const double y[],
         double dydt[],
         void * params)
{
    // Unpack
    const double Ma02   = ((params_type *)params)->Ma02;
    const double gam0m1 = ((params_type *)params)->gam0m1;
    const double u      =            y[0];
    const double rho    = gsl_sf_exp(y[1]);
    const double p      =            y[2];
    SUZERAIN_UNUSED(p);

    // Compute
    double up, a2;
    nozzle_upa2(R, u, Ma02, gam0m1, &up, &a2);
    const double logrhop  = -Ma02*R*u*up / a2;
    const double pp       = -Ma02*R*rho*u*up;

    // Pack
    dydt[0] = up;
    dydt[1] = logrhop;
    dydt[2] = pp;

    return GSL_SUCCESS;
}

suzerain_radial_nozzle_solution *
suzerain_radial_nozzle_solver(
    const double         Ma0,
    const double         gam0,
    const double         rho1,
    const double         u1,
    const double         p1,
    const double * const R,
    const size_t         size)
{
    return NULL;
}
