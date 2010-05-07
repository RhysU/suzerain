/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * htstretch.c: Routines providing hyperbolic tangent grid stretching
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.h>
#pragma hdrstop
#include <gsl/gsl_roots.h>
#include <suzerain/error.h>
#include <suzerain/htstretch.h>

#pragma warning(disable:981)

double
suzerain_htstretch1(const double delta,
                    const double L,
                    const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    return 1 + tanh(delta*(x/L - 1))/tanh(delta);
}

double
suzerain_htstretch1_ddelta(const double delta,
                           const double L,
                           const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    const double xoverLlessone = x/L-1;
    const double cosh_expr     = cosh(delta*xoverLlessone);
    const double tanh_expr     = tanh(delta*xoverLlessone);
    const double sinh_expr     = sinh(delta);
    return xoverLlessone/(tanh(delta)*cosh_expr*cosh_expr)
        -  tanh_expr/(sinh_expr*sinh_expr);
}

double
suzerain_htstretch1_dL(const double delta,
                       const double L,
                       const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    const double xoverLlessone = x/L-1;
    const double cosh_expr     = cosh(delta*xoverLlessone);
    return -x*delta/(L*L*tanh(delta)*cosh_expr*cosh_expr);
}

double
suzerain_htstretch1_dx(const double delta,
                       const double L,
                       const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    const double cosh_expr = cosh(delta*(x/L-1));
    return delta/(L*tanh(delta)*cosh_expr*cosh_expr);
}


double
suzerain_htstretch2(const double delta,
                    const double L,
                    const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    return 0.5*(1 + tanh(delta*(x/L - 0.5))/tanh(delta/2.));
}

double
suzerain_htstretch2_ddelta(const double delta,
                           const double L,
                           const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    const double xoverLlesshalf = x/L - 0.5;
    const double sinh_halfdelta = sinh(delta/2);
    const double cosh_term      = cosh(delta*xoverLlesshalf);
    return -  ((L-2*x)*sinh(delta)+L*sinh((2*x/L-1)*delta))
            / (8*L*cosh_term*cosh_term*sinh_halfdelta*sinh_halfdelta);
}

double
suzerain_htstretch2_dL(const double delta,
                       const double L,
                       const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    const double xoverLlesshalf = x/L - 0.5;
    const double cosh_term      = cosh(delta*xoverLlesshalf);
    return -x*delta/(tanh(delta/2)*2*L*L*cosh_term*cosh_term);
}

double
suzerain_htstretch2_dx(const double delta,
                       const double L,
                       const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
    const double xoverLlesshalf = x/L - 0.5;
    const double cosh_term      = cosh(delta*xoverLlesshalf);
    return delta/(tanh(delta/2)*2*L*cosh_term*cosh_term);
}

typedef struct delta_problem_params {
    double L;
    double crit_x;
    double crit_val;
} delta_problem_params;

static
double htstretch1_delta_problem_f(double delta, void *params)
{
    const delta_problem_params * const dpp = (delta_problem_params *) params;

    return   suzerain_htstretch1(delta, dpp->L, dpp->crit_x)
           - 0.0 /* == suzerain_htstretch1(  0.0, dpp->L, dpp->crit_x) */
           - dpp->crit_val;
}

static
double htstretch2_delta_problem_f(double delta, void *params)
{
    const delta_problem_params * const dpp = (delta_problem_params *) params;

    return   suzerain_htstretch2(delta, dpp->L, dpp->crit_x)
           - 0.0 /* == suzerain_htstretch2(  0.0, dpp->L, dpp->crit_x) */
           - dpp->crit_val;
}

int
suzerain_htstretch1_find_delta(double *delta,
                               const double L,
                               const double crit_x,
                               const double crit_val)
{
    // Initial bounds guess and sanity check
    const double delta_lower = 1e-8;
    const double crit_val_lower = suzerain_htstretch1(delta_lower, L, crit_x);
    const double delta_upper = 1e3;
    const double crit_val_upper = suzerain_htstretch1(delta_upper, L, crit_x);
    if (crit_val >= crit_val_lower) {
        SUZERAIN_ERROR("(crit_x, crit_val) incompatible with delta_lower",
                       SUZERAIN_FAILURE);
    }
    if (crit_val <= crit_val_upper) {
        SUZERAIN_ERROR("(crit_x, crit_val) incompatible with delta_upper",
                       SUZERAIN_FAILURE);
    }

    // Initialize nonlinear function to solve
    delta_problem_params dpp;
    dpp.L          = L;
    dpp.crit_x     = crit_x;
    dpp.crit_val   = crit_val;
    gsl_function F = { &htstretch1_delta_problem_f, &dpp };

    // Initialize the solver
    gsl_root_fsolver *fsolver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    if (!fsolver) {
        SUZERAIN_ERROR("Unable to allocate solver", SUZERAIN_ENOMEM);
    }
    if (gsl_root_fsolver_set(fsolver, &F, delta_lower, delta_upper)) {
        gsl_root_fsolver_free(fsolver);
        SUZERAIN_ERROR("Unable to initialize solver", SUZERAIN_ESANITY);
    }

    // Iterate until convergence or tolerance reached
    for (int i = 0; i < 100; ++i)
    {
        const int err = gsl_root_fsolver_iterate(fsolver);
        if (err) SUZERAIN_ERROR("Error during solver iteration", err);
        const double bracket =   gsl_root_fsolver_x_upper(fsolver)
                               - gsl_root_fsolver_x_lower(fsolver);
        if (fabs(bracket) < 1000*GSL_DBL_EPSILON) break;
    }

    // Save the root to delta and tear down the solver
    *delta = gsl_root_fsolver_root(fsolver);
    gsl_root_fsolver_free(fsolver);

    return SUZERAIN_SUCCESS;
}
