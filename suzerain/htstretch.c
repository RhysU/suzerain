/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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
 * @copydoc htstretch.h
 */

#include <suzerain/htstretch.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include <suzerain/common.h>
#include <suzerain/error.h>

#pragma warning(disable:981)

double
suzerain_htstretch1(const double delta,
                    const double L,
                    const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return x / L;
    } else {
        return 1 + tanh(delta*(x/L - 1))/tanh(delta);
    }
}

double
suzerain_htstretch1_ddelta(const double delta,
                           const double L,
                           const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return 0.0;
    } else {
        const double xoverLlessone = x/L-1;
        const double cosh_expr     = cosh(delta*xoverLlessone);
        const double tanh_expr     = tanh(delta*xoverLlessone);
        const double sinh_expr     = sinh(delta);
        return xoverLlessone/(tanh(delta)*cosh_expr*cosh_expr)
                -  tanh_expr/(sinh_expr*sinh_expr);
    }
}

double
suzerain_htstretch1_dL(const double delta,
                       const double L,
                       const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return -x / (L * L);
    } else {
        const double xoverLlessone = x/L-1;
        const double cosh_expr     = cosh(delta*xoverLlessone);
        return -x*delta/(L*L*tanh(delta)*cosh_expr*cosh_expr);
    }
}

double
suzerain_htstretch1_dx(const double delta,
                       const double L,
                       const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return 1.0 / L;
    } else {
        const double cosh_expr = cosh(delta*(x/L-1));
        return delta/(L*tanh(delta)*cosh_expr*cosh_expr);
    }
}

double
suzerain_htstretch2(const double delta,
                    const double L,
                    const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return x / L;
    } else {
        return 0.5*(1 + tanh(delta*(x/L - 0.5))/tanh(delta/2.));
    }
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
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return 0.0;
    } else {
        return -  ((L-2*x)*sinh(delta)+L*sinh((2*x/L-1)*delta))
                / (8*L*cosh_term*cosh_term*sinh_halfdelta*sinh_halfdelta);
    }
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
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return -x / (L * L);
    } else {
        return -x*delta/(tanh(delta/2)*2*L*L*cosh_term*cosh_term);
    }
}

double
suzerain_htstretch2_dx(const double delta,
                       const double L,
                       const double x)
{
    assert(L > 0);
    assert(0 <= x && x <= L);
#pragma warning(push,disable:1572)
    if (SUZERAIN_UNLIKELY(delta == 0.0)) {
#pragma warning(pop)
        return 1.0 / L;
    } else {
        const double xoverLlesshalf = x/L - 0.5;
        const double cosh_term      = cosh(delta*xoverLlesshalf);
        return delta/(tanh(delta/2)*2*L*cosh_term*cosh_term);
    }
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
double htstretch1_delta_problem_df(double delta, void *params)
{
    const delta_problem_params * const dpp = (delta_problem_params *) params;
    return suzerain_htstretch1_ddelta(delta, dpp->L, dpp->crit_x);
}

static
void htstretch1_delta_problem_fdf(double delta,
                                  void *params,
                                  double *f,
                                  double *df)
{
    *f  = htstretch1_delta_problem_f( delta, params);
    *df = htstretch1_delta_problem_df(delta, params);
}

static
double htstretch2_delta_problem_f(double delta, void *params)
{
    const delta_problem_params * const dpp = (delta_problem_params *) params;
    return   suzerain_htstretch2(delta, dpp->L, dpp->crit_x)
           - 0.0 /* == suzerain_htstretch2(  0.0, dpp->L, dpp->crit_x) */
           - dpp->crit_val;
}

static
double htstretch2_delta_problem_df(double delta, void *params)
{
    const delta_problem_params * const dpp = (delta_problem_params *) params;
    return suzerain_htstretch2_ddelta(delta, dpp->L, dpp->crit_x);
}

static
void htstretch2_delta_problem_fdf(double delta,
                                  void *params,
                                  double *f,
                                  double *df)
{
    *f  = htstretch2_delta_problem_f( delta, params);
    *df = htstretch2_delta_problem_df(delta, params);
}

static
int find_delta(gsl_function_fdf *problem_fdf,
               const double delta_lower,
               const double delta_upper,
               const double epsabs,
               const int maxiter,
               double *delta)
{
    // Initial bounds sanity check.
    // Hopefully ensures root existence for hyperbolic tangent functions,
    // At least for the class of things we care about here.
    const double crit_val
        = ((delta_problem_params*) problem_fdf->params)->crit_val;
    const double crit_val_lower = GSL_FN_FDF_EVAL_F(problem_fdf, delta_lower);
    const double crit_val_upper = GSL_FN_FDF_EVAL_F(problem_fdf, delta_upper);
    if (crit_val >= crit_val_lower) {
        SUZERAIN_ERROR("(crit_x, crit_val) incompatible with delta_lower",
                       SUZERAIN_FAILURE);
    }
    if (crit_val <= crit_val_upper) {
        SUZERAIN_ERROR("(crit_x, crit_val) incompatible with delta_upper",
                       SUZERAIN_FAILURE);
    }

    // Initialize a bracketing solver
    gsl_function problem_f = { problem_fdf->f, problem_fdf->params };
    gsl_root_fsolver *fsolver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    if (!fsolver) {
        SUZERAIN_ERROR("Unable to allocate solver", SUZERAIN_ENOMEM);
    }
    if (gsl_root_fsolver_set(fsolver, &problem_f, delta_lower, delta_upper)) {
        gsl_root_fsolver_free(fsolver);
        SUZERAIN_ERROR("Unable to initialize solver", SUZERAIN_ESANITY);
    }

    // Iterate until either iteration or tolerance criteria is met
    int iter = 0;
    do {
        *delta = gsl_root_fsolver_root(fsolver);
        const double residual = GSL_FN_EVAL(&problem_f, *delta);
        if (GSL_SUCCESS == gsl_root_test_residual(residual, epsabs))
            break;

        const int err = gsl_root_fsolver_iterate(fsolver);
        if (err) {
            gsl_root_fsolver_free(fsolver);
            SUZERAIN_ERROR("Error during solver iteration", err);
        }
    } while (iter++ < maxiter);

    // Tear down solver and finish up
    gsl_root_fsolver_free(fsolver);
    return (iter < maxiter) ? SUZERAIN_SUCCESS : SUZERAIN_EMAXITER;
}

int
suzerain_htstretch1_find_delta(const double L,
                               const double x_crit,
                               const double s_crit,
                               const double epsabs,
                               const int maxiter,
                               double *delta)
{
    assert(L > 0);
    assert(0.0 <= x_crit && x_crit <= L);
    assert(0.0 <= s_crit && s_crit <= 1.0);
    assert(epsabs > 0);
    assert(maxiter > 0);
    assert(delta != NULL);

    delta_problem_params dpp = { L, x_crit, s_crit };
    gsl_function_fdf fdf     = { &htstretch1_delta_problem_f,
                                 &htstretch1_delta_problem_df,
                                 &htstretch1_delta_problem_fdf,
                                 &dpp };

    return find_delta(&fdf, 1e-8, 1e3, epsabs, maxiter, delta);
}

int
suzerain_htstretch2_find_delta(const double L,
                               const double x_crit,
                               const double u_crit,
                               const double epsabs,
                               const int maxiter,
                               double *delta)
{
    assert(L > 0);
    assert(0.0 <= x_crit && x_crit <= L);
    assert(0.0 <= u_crit && u_crit <= 1.0);
    assert(epsabs > 0);
    assert(maxiter > 0);
    assert(delta != NULL);

    delta_problem_params dpp = { L, x_crit, u_crit };
    gsl_function_fdf fdf     = { &htstretch2_delta_problem_f,
                                 &htstretch2_delta_problem_df,
                                 &htstretch2_delta_problem_fdf,
                                 &dpp };

    return find_delta(&fdf, 1e-8, 1e3, epsabs, maxiter, delta);
}
