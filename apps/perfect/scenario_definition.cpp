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

#include "scenario_definition.hpp"

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling problem scenario parameters which are either
 * reference quantities or nondimensional parameters describing a particular
 * problem setup.
 */

namespace suzerain {

static void parse_positive(const std::string& s, real_t *t, const char *n)
{
    const real_t v = exprparse<real_t>(s, n);
    validation::ensure_positive(v, n);
    *t = v;
}

static void parse_nonnegative(const std::string& s, real_t *t, const char *n)
{
    const real_t v = exprparse<real_t>(s, n);
    validation::ensure_nonnegative(v, n);
    *t = v;
}

namespace perfect {

scenario_definition::scenario_definition()
    : definition_base("Nondimensional scenario parameters"),
      Re(std::numeric_limits<real_t>::quiet_NaN()),
      Ma(std::numeric_limits<real_t>::quiet_NaN()),
      Pr(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rho(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rho_u(std::numeric_limits<real_t>::quiet_NaN()),
      alpha(std::numeric_limits<real_t>::quiet_NaN()),
      beta(std::numeric_limits<real_t>::quiet_NaN()),
      gamma(std::numeric_limits<real_t>::quiet_NaN())
{
    this->initialize_options(NULL, NULL, NULL,
                             NULL, NULL,
                             NULL, NULL, NULL);
}

scenario_definition::scenario_definition(
        const char * Re,
        const char * Ma,
        const char * Pr,
        const char * bulk_rho,
        const char * bulk_rho_u,
        const char * alpha,
        const char * beta,
        const char * gamma)
    : definition_base("Nondimensional scenario parameters"),
      Re(std::numeric_limits<real_t>::quiet_NaN()),
      Ma(std::numeric_limits<real_t>::quiet_NaN()),
      Pr(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rho(std::numeric_limits<real_t>::quiet_NaN()),
      bulk_rho_u(std::numeric_limits<real_t>::quiet_NaN()),
      alpha(std::numeric_limits<real_t>::quiet_NaN()),
      beta(std::numeric_limits<real_t>::quiet_NaN()),
      gamma(std::numeric_limits<real_t>::quiet_NaN())
{
    this->initialize_options(Re, Ma, Pr,
                                bulk_rho, bulk_rho_u,
                                alpha, beta, gamma);
}

void scenario_definition::initialize_options(
        const char * default_Re,
        const char * default_Ma,
        const char * default_Pr,
        const char * default_bulk_rho,
        const char * default_bulk_rho_u,
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

    // bulk_rho_u
    p.reset(boost::program_options::value<std::string>(NULL));
    p->notifier(boost::bind(&parse_nonnegative, _1, &bulk_rho_u, "bulk_rho_u"));
    if (default_bulk_rho_u) p->default_value(default_bulk_rho_u);
    this->add_options()("bulk_rho_u", p.release(), "Bulk momentum target");

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

void save(const esio_handle h,
          const scenario_definition& scenario)
{
    DEBUG0("Storing scenario_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Re", &scenario.Re, 0,
            scenario.options().find("Re",false).description().c_str());

    esio_line_write(h, "Ma", &scenario.Ma, 0,
            scenario.options().find("Ma",false).description().c_str());

    esio_line_write(h, "Pr", &scenario.Pr, 0,
            scenario.options().find("Pr",false).description().c_str());

    esio_line_write(h, "bulk_rho", &scenario.bulk_rho, 0,
            scenario.options().find("bulk_rho",false).description().c_str());

    esio_line_write(h, "bulk_rho_u", &scenario.bulk_rho_u, 0,
            scenario.options().find("bulk_rho_u",false).description().c_str());

    esio_line_write(h, "alpha", &scenario.alpha, 0,
            scenario.options().find("alpha",false).description().c_str());

    esio_line_write(h, "beta", &scenario.beta, 0,
            scenario.options().find("beta",false).description().c_str());

    esio_line_write(h, "gamma", &scenario.gamma, 0,
            scenario.options().find("gamma",false).description().c_str());
}

void load(const esio_handle h,
          scenario_definition& scenario)
{
    DEBUG0("Loading scenario_definition parameters");

    esio_line_establish(h, 1, 0, 1); // All ranks load

    if (!(boost::math::isnan)(scenario.Re)) {
        INFO0("Overriding scenario using Re = " << scenario.Re);
    } else {
        esio_line_read(h, "Re", &scenario.Re, 0);
    }

    if (!(boost::math::isnan)(scenario.Ma)) {
        INFO0("Overriding scenario using Ma = " << scenario.Ma);
    } else {
        esio_line_read(h, "Ma", &scenario.Ma, 0);
    }

    if (!(boost::math::isnan)(scenario.Pr)) {
        INFO0("Overriding scenario using Pr = " << scenario.Pr);
    } else {
        esio_line_read(h, "Pr", &scenario.Pr, 0);
    }

    if (!(boost::math::isnan)(scenario.bulk_rho)) {
        INFO0("Overriding scenario using bulk_rho = " << scenario.bulk_rho);
    } else {
        esio_line_read(h, "bulk_rho", &scenario.bulk_rho, 0);
    }

    if (!(boost::math::isnan)(scenario.bulk_rho_u)) {
        INFO0("Overriding scenario using bulk_rho_u = " << scenario.bulk_rho_u);
    } else {
        esio_line_read(h, "bulk_rho_u", &scenario.bulk_rho_u, 0);
    }

    if (!(boost::math::isnan)(scenario.alpha)) {
        INFO0("Overriding scenario using alpha = " << scenario.alpha);
    } else {
        esio_line_read(h, "alpha", &scenario.alpha, 0);
    }

    if (!(boost::math::isnan)(scenario.beta)) {
        INFO0("Overriding scenario using beta = " << scenario.beta);
    } else {
        esio_line_read(h, "beta", &scenario.beta, 0);
    }

    if (!(boost::math::isnan)(scenario.gamma)) {
        INFO0("Overriding scenario using gamma = " << scenario.gamma);
    } else {
        esio_line_read(h, "gamma", &scenario.gamma, 0);
    }
}

} // namespace perfect

} // namespace suzerain
