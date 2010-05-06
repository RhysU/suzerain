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

typedef struct vinokur_delta_problem_params {
    double L;
    int    k;
    int    N;
    double u_crit;
} vinokur_delta_problem_params;

static
double doublesided_f(double delta, void *params)
{
    const vinokur_delta_problem_params * const vp
        = (vinokur_delta_problem_params *) params;

    const double koverNlesshalf = ((double) vp->k)/((double) vp->N) - 1./2.;
    double retval;
    retval  = delta*koverNlesshalf;
    retval  = tanh(retval);
    retval /= tanh(delta/2.);
    retval += 1.;
    retval *= 1./2.;
    retval -= vp->u_crit;

    return retval;
}

static
double doublesided_df(double delta, void *params)
{
    const vinokur_delta_problem_params * const vp
        = (vinokur_delta_problem_params *) params;

    const double koverNlesshalf = ((double) vp->k)/((double) vp->N) - 1./2.;
    const double csch_expr = 1./sinh(delta/2.);
    const double sech_expr = 1./cosh(koverNlesshalf*delta);
    const double twokoverNlessone = (2.*vp->k)/vp->N - 1.;
    const double retval = -(csch_expr*csch_expr*sech_expr*sech_expr)*(
                (vp->N-2*vp->k)*sinh(delta)+vp->N*sinh(twokoverNlessone*delta)
            )/(8*vp->N);

    return retval;
}

static
void doublesided_fdf(double delta, void *params, double *f, double *df)
{
    const vinokur_delta_problem_params * const vp
        = (vinokur_delta_problem_params *) params;
    const double koverNlesshalf = ((double) vp->k)/((double) vp->N) - 1./2.;

    // Function evaluation
    *f  = delta*koverNlesshalf;
    *f  = tanh(*f);
    *f /= tanh(delta/2.);
    *f += 1.;
    *f *= 1./2.;
    *f -= vp->u_crit;

    // Derivative evaluation
    const double csch_expr = 1./sinh(delta/2.);
    const double sech_expr = 1./cosh(koverNlesshalf*delta);
    const double twokoverNlessone = (2.*vp->k)/vp->N - 1.;
    *df = -(csch_expr*csch_expr*sech_expr*sech_expr)*(
                (vp->N-2*vp->k)*sinh(delta)+vp->N*sinh(twokoverNlessone*delta)
            )/(8*vp->N);
}
