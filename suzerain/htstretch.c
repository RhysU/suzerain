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

double
suzerain_htstretch_onesided_continuous(const double I,
                                       const double delta,
                                       const double xi)
{
#pragma warning(push,disable:981)
    return 1 + tanh(delta*(xi/I - 1))/tanh(delta);
#pragma warning(pop)
}

double
suzerain_htstretch_onesided_discrete(const double I,
                                     const double delta,
                                     const int N,
                                     const int j)
{
    return suzerain_htstretch_onesided_continuous(I, delta, j*I/N);
}

double
suzerain_htstretch_twosided_continuous(const double I,
                                       const double delta,
                                       const double xi)
{
#pragma warning(push,disable:981)
    return 0.5*(1 + tanh(delta*(xi/I - 0.5))/tanh(delta/2.));
#pragma warning(pop)
}

double
suzerain_htstretch_twosided_discrete(const double I,
                                     const double delta,
                                     const int N,
                                     const int j)
{
    return suzerain_htstretch_twosided_continuous(I, delta, j*I/N);
}

typedef struct vinokur_delta_problem_params {
    double I;
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
#pragma warning(push,disable:981)
    const double retval = -(csch_expr*csch_expr*sech_expr*sech_expr)*(
                (vp->N-2*vp->k)*sinh(delta)+vp->N*sinh(twokoverNlessone*delta)
            )/(8*vp->N);
#pragma warning(pop)

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
#pragma warning(push,disable:981)
    *df = -(csch_expr*csch_expr*sech_expr*sech_expr)*(
                (vp->N-2*vp->k)*sinh(delta)+vp->N*sinh(twokoverNlessone*delta)
            )/(8*vp->N);
#pragma warning(pop)
}
