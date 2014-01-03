/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2013-2014 Rhys Ulerich
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
 * @copydoc radialflow.h
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/radialflow.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_poly.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

// Find pointwise solution details given state, Ma2=Ma**2, gamm1=gam-1
// Compare nozzle_details function source within notebooks/nozzle.m
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
// Compare nozzle_rhs function source within notebooks/nozzle.m
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

// Compare nozzle function source within notebooks/nozzle.m
suzerain_radialflow_solution *
suzerain_radialflow_solver(
    const double         Ma0,
    const double         gam0,
    const double         u1,
    const double         rho1,
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
    suzerain_radialflow_solution * const s = malloc(
        sizeof(suzerain_radialflow_solution)
        + size * sizeof(suzerain_radialflow_state));
    if (!s) SUZERAIN_ERROR_NULL("Solution allocation failed", SUZERAIN_ENOMEM);
    s->state[0].R = s->state[0].u = s->state[0].rho = s->state[0].p = GSL_NAN;

    // Use GNU Scientific Library ODE integrator on [u; rho; p]' system
    double current_R = R[0];
    double y[3] = { u1, rho1, p1 };
    params_type params = { Ma0*Ma0, gam0 - 1 };
    gsl_odeiv2_system sys = { &nozzle_rhs, 0, sizeof(y)/sizeof(y[0]), &params};
    const double abstol = 0;
    const double reltol = GSL_SQRT_DBL_EPSILON;
    gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new(
            &sys, gsl_odeiv2_step_rkf45, GSL_SQRT_DBL_EPSILON, abstol, reltol);
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
suzerain_radialflow_delta(
    const suzerain_radialflow_solution * s,
    const size_t i)
{
    assert(i < s->size);
    return sqrt(gsl_pow_2(s->state[i].R) - gsl_pow_2(s->state[0].R));
}

double
suzerain_radialflow_qoi_Mae(
    const suzerain_radialflow_solution * s,
    const size_t i)
{
    assert(i < s->size);
    return  (s->Ma0 * s->state[0].R * fabs(s->state[i].u) )
          / (         s->state[i].R * sqrt(s->state[i].a2));
}

double
suzerain_radialflow_qoi_pexi(
    const suzerain_radialflow_solution * s,
    const size_t i)
{
    assert(i < s->size);
    const double delta = suzerain_radialflow_delta(s, i);
    return - (s->state[i].R * delta * s->state[i].up)
           / (s->state[0].R * fabs(s->state[i].u));
}

double
suzerain_radialflow_qoi_Te(
    const suzerain_radialflow_solution * s,
    const size_t i,
    const double Ma)
{
    assert(i < s->size);
    return s->gam0 * gsl_pow_2(Ma / s->Ma0) * (s->state[i].p / s->state[i].rho);
}

// FIXME Improve output names to reflect returned edge state!
int
suzerain_radialflow_qoi_match(
    const double delta,
    const double gam0,
    const double Ma_e,
    const double p_exi,
    const double T_e,
    double *Ma0,
    double *R0,
    double *R,
    double *uR,
    double *rhoR,
    double *pR)
{

    // Solve quadratic for R0
    double x[2];
    const double Ma_e2  = gsl_pow_2(Ma_e );
    const double delta2 = gsl_pow_2(delta);
    const int nreal = gsl_poly_solve_quadratic(
                          - fabs(Ma_e2 - 1) * fabs(p_exi),
                          delta,
                          - Ma_e2 * delta2 * fabs(p_exi) * GSL_SIGN(Ma_e2 - 1),
                          x+0, x+1);

    // If real roots exist...
    if (nreal > 0) {
        // ...take the largest one possible...
        *R0   = x[nreal - 1];
        *Ma0  = 1 / sqrt(1/Ma_e2 + (gam0 - 1)*delta2/gsl_pow_2(*R0)/2);
        *R    = sqrt(gsl_pow_2(*R0) + delta2);
        *uR   = - (*R / *R0) * GSL_SIGN(p_exi*(Ma_e2 - 1));
        *rhoR = 1;
        *pR   = T_e == 0
              ? 1
              : T_e * gsl_pow_2(*Ma0) * *rhoR/gam0/Ma_e2;
    } else {
        // ...otherwise defensive failure...
        *R0   = GSL_NAN;
        *Ma0  = GSL_NAN;
        *R    = GSL_NAN;
        *uR   = GSL_NAN;
        *rhoR = GSL_NAN;
        *pR   = GSL_NAN;
    }
    // ...notice success requires a positive R0!
    return *R0 >= 0 ? SUZERAIN_SUCCESS : SUZERAIN_EFAILED;
}

inline // Suggests inlining within suzerain_radialflow_cartesian_conserved
void
suzerain_radialflow_cartesian_primitive(
    const suzerain_radialflow_solution * s,
    const size_t i,
    const double Ma,
    double * const rho,
    double * const u,
    double * const v,
    double * const p,
    double * const a,
    double * const rho_xi,
    double * const u_xi,
    double * const v_xi,
    double * const p_xi,
    double * const a_xi,
    double * const rho_y,
    double * const u_y,
    double * const v_y,
    double * const p_y,
    double * const a_y)
{
    assert(i < s->size);
    // If u(R0) >= 0, compute at coords (x,delta).  Otherwise (-x,delta).
    // Possibly flipping from the first to the second quadrant permits
    // flow to always be in positive x direction regardless of sonic-ness.
    const double x       = (s->state[0].u >= 0 ? 1 : -1) * s->state[0].R;
    const double x_inv_R = x / s->state[i].R;
    const double x2      = gsl_pow_2(x);
    const double y       = suzerain_radialflow_delta(s, i);
    const double y2      = gsl_pow_2(y);
    const double y_inv_R = y / s->state[i].R;
    const double inv_R   = 1 / s->state[i].R;
    const double inv_R2  = gsl_pow_2(inv_R);
    const double Ma2Ma02 = gsl_pow_2(Ma / s->Ma0);

    const suzerain_radialflow_state * const t = &s->state[i];
    *rho    = t->rho;
    *u      = t->u * x_inv_R;
    *v      = t->u * y_inv_R;
    *p      = t->p * Ma2Ma02;
    *rho_xi = t->rhop * x_inv_R;
    *u_xi   =          inv_R2 * (x2 * t->up + y2 * inv_R * t->u);
    *v_xi   = x_inv_R*y_inv_R * (     t->up -      inv_R * t->u);
    *p_xi   = t->pp   * Ma2Ma02 * x_inv_R;
    *rho_y  = t->rhop * y_inv_R;
    *u_y    = x_inv_R*y_inv_R * (     t->up -      inv_R * t->u);
    *v_y    =          inv_R2 * (y2 * t->up + x2 * inv_R * t->u);
    *p_y    = t->pp   * Ma2Ma02 * y_inv_R;

    // Beware (*a) contains Ma/Ma0*a when computing ap, *a_xi, and *a_y
    *a              = sqrt(t->a2 * Ma2Ma02);
    const double ap = (1 - s->gam0) / 2 * Ma * Ma * t->u * t->up / (*a);
    *a_xi           = x_inv_R * ap;
    *a_y            = y_inv_R * ap;
}

void
suzerain_radialflow_cartesian_conserved(
    const suzerain_radialflow_solution * s,
    const size_t i,
    const double Ma,
    double * const r,
    double * const ru,
    double * const rv,
    double * const rE,
    double * const p,
    double * const r_xi,
    double * const ru_xi,
    double * const rv_xi,
    double * const rE_xi,
    double * const p_xi,
    double * const r_y,
    double * const ru_y,
    double * const rv_y,
    double * const rE_y,
    double * const p_y)
{
    // Delegate to compute local primitive state and {p,a}_{,_xi,_y}
    double rho, u, v, a, rho_xi, u_xi, v_xi, a_xi, rho_y, u_y, v_y, a_y;
    suzerain_radialflow_cartesian_primitive(
            s, i, Ma, &rho,    &u,    &v,    p,    &a,
                      &rho_xi, &u_xi, &v_xi, p_xi, &a_xi,
                      &rho_y,  &u_y,  &v_y,  p_y,  &a_y);

    // Though expected from formulation, method does not assume H = constant
    // This permits both detecting formulation and coding errors as well as
    // correctly propagating any roundoff-ish drift in H from solution.
    const double invgam0m1 = 1/(s->gam0 - 1);
    const double Ma22      = Ma*Ma/2;
    const double H         = a*a*invgam0m1 + Ma22*(u*u + v*v);
    const double H_xi      = 2*(a*a_xi*invgam0m1 + Ma22*(u*u_xi + v*v_xi));
    const double H_y       = 2*(a*a_y *invgam0m1 + Ma22*(u*u_y  + v*v_y ));

    // Convert to conserved state using \rho H = \rho E + p
    *r  = rho;
    *ru = rho*u;
    *rv = rho*v;
    *rE = rho*H - *p;

    // Convert streamwise primitive derivatives to conserved derivatives
    *r_xi  = rho_xi;
    *ru_xi = rho*u_xi + rho_xi*u;
    *rv_xi = rho*v_xi + rho_xi*v;
    *rE_xi = rho*H_xi + rho_xi*H - *p_xi;

    // Convert wall-normal primitive derivatives to conserved derivatives
    *r_y  = rho_y ;
    *ru_y = rho*u_y + rho_y*u;
    *rv_y = rho*v_y + rho_y*v;
    *rE_y = rho*H_y + rho_y*H - *p_y;
}
