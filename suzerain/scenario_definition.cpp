//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// scenario_definition.cpp: classes handling problem scenario parameters
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling problem scenario parameters which are either
 * reference quantities or nondimensional parameters describing a particular
 * problem setup.
 */

namespace suzerain {

static void parse_positive(const std::string& s, real_t *t, const char *n)
{
    const real_t v = suzerain::exprparse<real_t>(s, n);
    suzerain::validation::ensure_positive(v, n);
    *t = v;
}

static void parse_nonnegative(const std::string& s, real_t *t, const char *n)
{
    const real_t v = suzerain::exprparse<real_t>(s, n);
    suzerain::validation::ensure_nonnegative(v, n);
    *t = v;
}

namespace problem {

ScenarioDefinition::ScenarioDefinition()
    : IDefinition("Nondimensional scenario parameters"),
      Re(std::numeric_limits<real_t>::quiet_NaN()),
      Ma(std::numeric_limits<real_t>::quiet_NaN()),
      Pr(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rho(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rhou(std::numeric_limits<real_t>::quiet_NaN()),
      alpha(std::numeric_limits<real_t>::quiet_NaN()),
      beta(std::numeric_limits<real_t>::quiet_NaN()),
      gamma(std::numeric_limits<real_t>::quiet_NaN())
{
    this->initialize_options(NULL, NULL, NULL,
                             NULL, NULL,
                             NULL, NULL, NULL);
}

ScenarioDefinition::ScenarioDefinition(
        const char * Re,
        const char * Ma,
        const char * Pr,
        const char * bulk_rho,
        const char * bulk_rhou,
        const char * alpha,
        const char * beta,
        const char * gamma)
    : IDefinition("Nondimensional scenario parameters"),
      Re(std::numeric_limits<real_t>::quiet_NaN()),
      Ma(std::numeric_limits<real_t>::quiet_NaN()),
      Pr(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rho(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rhou(std::numeric_limits<real_t>::quiet_NaN()),
      alpha(std::numeric_limits<real_t>::quiet_NaN()),
      beta(std::numeric_limits<real_t>::quiet_NaN()),
      gamma(std::numeric_limits<real_t>::quiet_NaN())
{
    this->initialize_options(Re, Ma, Pr,
                                bulk_rho, bulk_rhou,
                                alpha, beta, gamma);
}

void ScenarioDefinition::initialize_options(
        const char * default_Re,
        const char * default_Ma,
        const char * default_Pr,
        const char * default_bulk_rho,
        const char * default_bulk_rhou,
        const char * default_alpha,
        const char * default_beta,
        const char * default_gamma)
{
    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    std::auto_ptr<boost::program_options::typed_value<std::string> > p;

    // Re
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Re, "Re"));
    if (default_Re) p->default_value(default_Re);
    this->add_options()("Re", p.release(), "Reynolds number");

    // Ma
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Ma, "Ma"));
    if (default_Ma) p->default_value(default_Ma);
    this->add_options()("Ma", p.release(), "Mach number");

    // Pr
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &Pr, "Pr"));
    if (default_Pr) p->default_value(default_Pr);
    this->add_options()("Pr", p.release(), "Prandtl number");

    // bulk_rho
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &bulk_rho, "bulk_rho"));
    if (default_bulk_rho) p->default_value(default_bulk_rho);
    this->add_options()("bulk_rho", p.release(), "Bulk density target");

    // bulk_rhou
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &bulk_rhou, "bulk_rhou"));
    if (default_bulk_rhou) p->default_value(default_bulk_rhou);
    this->add_options()("bulk_rhou", p.release(), "Bulk momentum target");

    // alpha
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &alpha, "alpha"));
    if (default_alpha) p->default_value(default_alpha);
    this->add_options()("alpha", p.release(),
                        "Ratio of bulk to dynamic viscosity");

    // beta
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &beta, "beta"));
    if (default_beta) p->default_value(default_beta);
    this->add_options()("beta", p.release(), "Temperature power law exponent");

    // gamma
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_positive, _1, &gamma, "gamma"));
    if (default_gamma) p->default_value(default_gamma);
    this->add_options()("gamma", p.release(), "Ratio of specific heats");
}

} // namespace problem

} // namespace suzerain
