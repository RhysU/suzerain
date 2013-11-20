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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

// Find pointwise solution details given state, Ma2=Ma**2, gamm1=gam-1
// Compare nozzle_details function source within writeups/notebooks/nozzle.m
static inline
void nozzle_details(const double R,
                    const double u,
                    const double rho,
                    const double p,
                    const double Ma02,
                    const double gam0m1,
                    double *up,
                    double *rhop,
                    double *pp,
                    double *a2)
{
    const double u2 = u*u;
    const double C  = (2/Ma02 + gam0m1*(1 - u2));
    *up   = (u*C) / (R*(2*u2 - C));
    *pp   = -Ma02*rho*u*(*up);
    *a2   = 1 + Ma02/2*gam0m1*(1 - u2);
    *rhop = *pp / *a2;
    SUZERAIN_UNUSED(p);
}

// Parameters for \ref nozzle per conventions of gsl_odeiv2_system->function
typedef struct params_type {
    double Ma02;   // Ma0^2
    double gam0m1; // gamma_0 - 1
} params_type;

// Find [u; rho; p]' given r, x=[u; rho; p], Ma02=Ma0**2, gam0m1=gam0-1
// Compare nozzle_rhs function source within writeups/notebooks/nozzle.m
static int nozzle_rhs(double R,
                      const double y[],
                      double dydt[],
                      void *params)
{
    double a2;
    nozzle_details(R, y[0], y[1], y[2],
                   ((params_type *)params)->Ma02,
                   ((params_type *)params)->gam0m1,
                   &dydt[0], &dydt[1], &dydt[2], &a2);
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

    // Use GNU Scientific Library ODE integrator on [u; rho; p]' system
    double current_R = R[0];
    double y[3] = { u1, rho1, p1 };
    params_type params = { Ma0*Ma0, gam0 - 1 };
    gsl_odeiv2_system sys = { &nozzle_rhs, 0, sizeof(y)/sizeof(y[0]), &params};
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
        s->state[i].u   = y[0];
        s->state[i].rho = y[1];
        s->state[i].p   = y[2];
    }
    gsl_odeiv2_driver_free(driver);
    if (error) {
        char msg[128];
        snprintf(msg, sizeof(msg), "gsl_odeiv2_driver_apply reports %d: %s",
                error, gsl_strerror(error));
        free(s);
        SUZERAIN_ERROR_NULL(msg, GSL_EFAILED);
    }

    // Save scenario parameters and initial conditions (marking solution valid)
    s->Ma0          = Ma0;
    s->gam0         = gam0;
    s->size         = size;
    s->state[0].R   = R[0];
    s->state[0].u   = u1;
    s->state[0].rho = rho1;
    s->state[0].p   = p1;

    // Compute derivatives and sound speed squared directly from state
    for (size_t i = 0; i < size; ++i) {
        nozzle_details(s->state[i].R,
                       s->state[i].u,
                       s->state[i].rho,
                       s->state[i].p,
                       s->Ma0 * s->Ma0,
                       s->gam0 - 1,
                       &s->state[i].up,
                       &s->state[i].rhop,
                       &s->state[i].pp,
                       &s->state[i].a2);
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
    double *rho_xi,
    double *u_xi,
    double *v_xi,
    double *p_xi,
    double *rho_y,
    double *u_y,
    double *v_y,
    double *p_y)
{
    assert(i < s->size);
    const double x       = s->state[0].R;
    const double x2      = gsl_pow_2(x);
    const double y       = suzerain_radial_nozzle_delta(s, i);
    const double y2      = gsl_pow_2(y);
    const double inv_R   = 1 / s->state[i].R;
    const double inv_R2  = gsl_pow_2(inv_R);
    const double Ma2Ma02 = gsl_pow_2(Ma / s->Ma0);

    const suzerain_radial_nozzle_state * const t = &s->state[i];
    *rho    = t->rho;
    *u      = t->u * x * inv_R;
    *v      = t->u * y * inv_R;
    *p      = t->p * Ma2Ma02;
    *rho_xi = t->rhop * x * inv_R;
    *u_xi   =         inv_R2 * (x2 * t->up + y2 * inv_R * t->u);
    *v_xi   = x * y * inv_R2 * (     t->up -      inv_R * t->u);
    *p_xi   = t->pp   * Ma2Ma02 * x * inv_R;
    *rho_y  = t->rhop * y * inv_R;
    *u_y    = x * y * inv_R2 * (     t->up -      inv_R * t->u);
    *v_y    =         inv_R2 * (y2 * t->up + x2 * inv_R * t->u);
    *p_y    = t->pp * Ma2Ma02 * y * inv_R;
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
    double *p,
    double *r_xi,
    double *ru_xi,
    double *rv_xi,
    double *rE_xi,
    double *p_xi,
    double *r_y,
    double *ru_y,
    double *rv_y,
    double *rE_y,
    double *p_y)
{
    // Delegate to compute local primitive state
    double rho, u, v, rho_xi, u_xi, v_xi, rho_y, u_y, v_y;
    suzerain_radial_nozzle_cartesian_primitive(
        s, i, Ma, &rho,    &u,    &v,    p,
                  &rho_xi, &u_xi, &v_xi, p_xi,
                  &rho_y,  &u_y,  &v_y,  p_y);

    // Convert primitive to conserved state
    const double invgam0m1 = 1/(s->gam0 - 1), Ma2 = Ma*Ma, u2 = u*u, v2 = v*v;
    *r  = rho;
    *ru = rho*u;
    *rv = rho*v;
    *rE = *p * invgam0m1 + Ma2 / 2 * rho * (u2 + v2);

    // Convert streamwise primitive derivatives to conserved derivatives
    *r_xi  = rho_xi;
    *ru_xi = rho*u_xi + rho_xi*u;
    *rv_xi = rho*v_xi + rho_xi*v;
    *rE_xi = *p_xi * invgam0m1
           + Ma2 * (
                   rho_xi / 2 * (u2     + v2    )
                +  rho        * (u*u_xi + v*v_xi)
             );

    // Convert wall-normal primitive derivatives to conserved derivatives
    *r_y   = rho_y ;
    *ru_y  = rho*u_y  + rho_y *u;
    *rv_y  = rho*v_y  + rho_y *v;
    *rE_y  = *p_y * invgam0m1
           + Ma2 * (
                   rho_y  / 2 * (u2     + v2    )
                +  rho        * (u*u_y  + v*v_y )
             );
}
