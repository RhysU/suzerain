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
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>

// Compute u' and a^2 given r, u, Ma02 = Ma0**2, and gam0m1 = gam0 - 1
// Compare nozzle_upa2 function source within writeups/notebooks/nozzle.m
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

// Parameters for \ref nozzle per conventions of gsl_odeiv2_system->function
typedef struct params_type {
    double Ma02;   // Ma0^2
    double gam0m1; // gamma_0 - 1
} params_type;

// Find [u; log rho; p]' given r, x=[u; log rho; p], Ma02=Ma0**2, gam0m1=gam0-1
// Compare nozzle_f function source within writeups/notebooks/nozzle.m
static
int
nozzle_f(double R,
         const double y[],
         double dydt[],
         void *params)
{
    // Unpack
    const double Ma02   = ((params_type *)params)->Ma02;
    const double gam0m1 = ((params_type *)params)->gam0m1;
    const double u      =     y[0];
    const double rho    = exp(y[1]);
    const double p      =     y[2];
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

// Compare nozzle function source within writeups/notebooks/nozzle.m
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
    // Sanity check incoming arguments and realizability
    if (Ma0  <= 0) SUZERAIN_ERROR_NULL("Ma0 <= 0",  SUZERAIN_EDOM);
    if (gam0 <= 1) SUZERAIN_ERROR_NULL("gam0 <= 1", SUZERAIN_EDOM);
    if (rho1 <= 0) SUZERAIN_ERROR_NULL("rho1 <= 0", SUZERAIN_EDOM);
    if (p1   <= 0) SUZERAIN_ERROR_NULL("p1 <= 0",   SUZERAIN_EDOM);
    if (size == 0) SUZERAIN_ERROR_NULL("size == 0", SUZERAIN_EINVAL);
    if (!(u1*u1 < 2 / (Ma0*Ma0 * (gam0-1)) + 1)) {
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "Ma0=%g, gam0=%g, u1=%g imply a*a <= 0", Ma0, gam0, u1);
        SUZERAIN_ERROR_NULL(msg, SUZERAIN_EDOM);
    }

    // Allocate solution storage into which integrator results will be saved
    suzerain_radial_nozzle_solution * const s = malloc(
        sizeof(suzerain_radial_nozzle_solution)
        + size * sizeof(suzerain_radial_nozzle_state));
    if (!s) SUZERAIN_ERROR_NULL("Solution allocation failed", SUZERAIN_ENOMEM);
    s->state[0].R = s->state[0].u = s->state[0].rho = s->state[0].p = GSL_NAN;

    // Use GNU Scientific Library ODE integrator on [u; log rho; p]' system
    double current_R = R[0];
    double y[3] = { u1, log(rho1), p1 };
    params_type params = { Ma0*Ma0, gam0 - 1 };
    gsl_odeiv2_system sys = { &nozzle_f, NULL, sizeof(y)/sizeof(y[0]), &params};
    const double abstol = GSL_SQRT_DBL_EPSILON;
    const double reltol = GSL_SQRT_DBL_EPSILON;
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(
            &sys, gsl_odeiv2_step_rkf45, sqrt(abstol), abstol, reltol);
    if (!driver) {
        free(s);
        SUZERAIN_ERROR_NULL("Driver allocation failed", GSL_EFAILED);
    }
    int error = 0;
    for (size_t i = 1; i < size && !error; ++i) { // Advance state to R[i]
        error = gsl_odeiv2_driver_apply(driver, &current_R, R[i], y);
        s->state[i].R   = current_R;
        s->state[i].u   =     y[0];
        s->state[i].rho = exp(y[1]); // (log rho) -> (rho)
        s->state[i].p   =     y[2];
    }
    gsl_odeiv2_driver_free(driver);
    if (error) {
        char msg[128];
        snprintf(msg, sizeof(msg), "gsl_odeiv2_driver_apply reports %d: %s",
                error, gsl_strerror(error));
        free(s);
        SUZERAIN_ERROR_NULL(msg, GSL_EFAILED);
    }

    // Save solution parameters and initial conditions
    s->Ma0          = Ma0;
    s->gam0         = gam0;
    s->size         = size;
    s->state[0].R   = R[0];
    s->state[0].u   = u1;
    s->state[0].rho = rho1;
    s->state[0].p   = p1;

    // Compute sound speed squared and derivative information from state
    for (size_t i = 0; i < size; ++i) {
        nozzle_upa2(s->state[i].R,
                    s->state[i].u,
                    params.Ma02,
                    params.gam0m1,
                    &s->state[i].up,
                    &s->state[i].a2);
        s->state[i].pp   = -params.Ma02
                         * s->state[i].rho * s->state[i].u * s->state[i].up;
        s->state[i].rhop = s->state[i].pp / s->state[i].a2;
    }

    return s;
}

inline // Permits inlining by functions appearing later in this file
double
suzerain_radial_nozzle_delta(
    const suzerain_radial_nozzle_solution * s,
    const size_t i)
{
    assert(i < s->size);
    return sqrt(gsl_pow_2(s->state[i].R) - gsl_pow_2(s->state[0].R));
}

double
suzerain_radial_nozzle_qoi_Mae(
    const suzerain_radial_nozzle_solution * s,
    const size_t i)
{
    assert(i < s->size);
    return  (s->Ma0 * s->state[0].R * fabs(s->state[i].u) )
          / (         s->state[i].R * sqrt(s->state[i].a2));
}

double
suzerain_radial_nozzle_qoi_pexi(
    const suzerain_radial_nozzle_solution * s,
    const size_t i)
{
    assert(i < s->size);
    const double sgn_u = s->state[i].u >= 0 ? 1 : -1;
    const double delta = suzerain_radial_nozzle_delta(s, i);
    return (sgn_u * s->state[i].R * delta * s->state[i].pp)
         / (   gsl_pow_2(s->Ma0) * s->state[0].R
             * s->state[i].rho * gsl_pow_2(s->state[i].u));
}

inline // Suggests inlining within suzerain_radial_nozzle_cartesian_conserved
void
suzerain_radial_nozzle_cartesian_primitive(
    const suzerain_radial_nozzle_solution * s,
    const size_t i,
    const double Ma,
    double *rho,
    double *u,
    double *v,
    double *p,
    double *rhop,
    double *up,
    double *vp,
    double *pp)
{
    assert(i < s->size);
    const double x       = s->state[0].R;
    const double y       = suzerain_radial_nozzle_delta(s, i);
    const double inv_R   = 1 / s->state[i].R;
    const double sgn_u   = s->state[i].u >= 0 ? 1 : -1;
    const double Ma02Ma2 = gsl_pow_2(s->Ma0) / gsl_pow_2(Ma);

    const suzerain_radial_nozzle_state * const t = &s->state[i];
    *rho  = t->rho;
    *u    = fabs(t->u) * x * inv_R;
    *v    =      t->u  * y * inv_R;
    *p    = Ma02Ma2 * t->p;
    *rhop = x * sgn_u * t->rhop * inv_R;
    *up   = sgn_u * gsl_pow_2(inv_R) * (   gsl_pow_2(x) * t->up
                                         + gsl_pow_2(y) * t->u  * inv_R);
    *vp   = x * y * sgn_u * gsl_pow_2(inv_R) * (   t->up
                                                 - t->u  * inv_R);
    *pp   = x * sgn_u * Ma02Ma2 * t->pp * inv_R;
}

void
suzerain_radial_nozzle_cartesian_conserved(
    const suzerain_radial_nozzle_solution * s,
    const size_t i,
    const double Ma,
    double *r,
    double *ru,
    double *rv,
    double *rE,
    double *rp,
    double *rup,
    double *rvp,
    double *rEp)
{
    // Delegate to compute local primitive state
    double rho, u, v, p, rhop, up, vp, pp;
    suzerain_radial_nozzle_cartesian_primitive(
        s, i, Ma, &rho, &u, &v, &p, &rhop, &up, &vp, &pp);

    // Convert primitive to conserved state
    *r  = rho;
    *ru = rho*u;
    *rv = rho*v;
    *rE = p / (s->gam0 - 1) + Ma*Ma / 2 * rho * (u*u + v*v);

    // Convert primitive derivatives to conserved derivatives
    *rp  = rhop;
    *rup = rho*up + rhop*u;
    *rvp = rho*vp + rhop*v;
    *rEp = pp / (s->gam0 - 1)
         + Ma*Ma * (
                 rhop / 2 * (u*u  + v*v )
              +  rho      * (u*up + v*vp)
           );
}
