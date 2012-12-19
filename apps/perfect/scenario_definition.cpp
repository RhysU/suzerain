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

/** @file
 * @copydoc scenario_definition.hpp
 */

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
    : Re        (std::numeric_limits<real_t>::quiet_NaN())
    , Ma        (std::numeric_limits<real_t>::quiet_NaN())
    , Pr        (std::numeric_limits<real_t>::quiet_NaN())
    , bulk_rho  (std::numeric_limits<real_t>::quiet_NaN())
    , bulk_rho_u(std::numeric_limits<real_t>::quiet_NaN())
    , alpha     (std::numeric_limits<real_t>::quiet_NaN())
    , beta      (std::numeric_limits<real_t>::quiet_NaN())
    , gamma     (std::numeric_limits<real_t>::quiet_NaN())
{
}

scenario_definition::scenario_definition(
        const real_t Re,
        const real_t Ma,
        const real_t Pr,
        const real_t bulk_rho,
        const real_t bulk_rho_u,
        const real_t alpha,
        const real_t beta,
        const real_t gamma)
    : Re        (Re        )
    , Ma        (Ma        )
    , Pr        (Pr        )
    , bulk_rho  (bulk_rho  )
    , bulk_rho_u(bulk_rho_u)
    , alpha     (alpha     )
    , beta      (beta      )
    , gamma     (gamma     )
{
}

// Descriptions used in options_description and possibly save/load.
static const char description_Re[]
        = "Reynolds number";

static const char description_Ma[]
        = "Mach number";

static const char description_Pr[]
        = "Prandtl number";

static const char description_bulk_rho[]
        = "Bulk density target";

static const char description_bulk_rho_u[]
        = "Bulk momentum target";

static const char description_alpha[]
        = "Ratio of bulk to dynamic viscosity";

static const char description_beta[]
        = "Temperature power law exponent";

static const char description_gamma[]
        = "Ratio of specific heats";


boost::program_options::options_description
scenario_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;

    options_description retval("Nondimensional scenario parameters");

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    auto_ptr<typed_value<string> > p;

    // Re
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Re, "Re"));
    if (!(boost::math::isnan)(Re)) {
        p->default_value(lexical_cast<string>(Re));
    }
    retval.add_options()("Re", p.release(), description_Re);

    // Ma
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Ma, "Ma"));
    if (!(boost::math::isnan)(Ma)) {
        p->default_value(lexical_cast<string>(Ma));
    }
    retval.add_options()("Ma", p.release(), description_Ma);

    // Pr
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Pr, "Pr"));
    if (!(boost::math::isnan)(Pr)) {
        p->default_value(lexical_cast<string>(Pr));
    }
    retval.add_options()("Pr", p.release(), description_Pr);

    // bulk_rho
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &bulk_rho, "bulk_rho"));
    if (!(boost::math::isnan)(bulk_rho)) {
        p->default_value(lexical_cast<string>(bulk_rho));
    }
    retval.add_options()("bulk_rho", p.release(), description_bulk_rho);

    // bulk_rho_u
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &bulk_rho_u, "bulk_rho_u"));
    if (!(boost::math::isnan)(bulk_rho_u)) {
        p->default_value(lexical_cast<string>(bulk_rho_u));
    }
    retval.add_options()("bulk_rho_u", p.release(), description_bulk_rho_u);

    // alpha
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &alpha, "alpha"));
    if (!(boost::math::isnan)(alpha)) {
        p->default_value(lexical_cast<string>(alpha));
    }
    retval.add_options()("alpha", p.release(), description_alpha);

    // beta
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &beta, "beta"));
    if (!(boost::math::isnan)(beta)) {
        p->default_value(lexical_cast<string>(beta));
    }
    retval.add_options()("beta", p.release(), description_beta);

    // gamma
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &gamma, "gamma"));
    if (!(boost::math::isnan)(gamma)) {
        p->default_value(lexical_cast<string>(gamma));
    }
    retval.add_options()("gamma", p.release(), description_gamma);

    return retval;
}

void save(const esio_handle h,
          const scenario_definition& scenario)
{
    DEBUG0("Storing scenario_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, "Re", &scenario.Re, 0, description_Re);
    esio_line_write(h, "Ma", &scenario.Ma, 0, description_Ma);
    esio_line_write(h, "Pr", &scenario.Pr, 0, description_Pr);
    esio_line_write(h, "bulk_rho",   &scenario.bulk_rho,
                    0, description_bulk_rho);
    esio_line_write(h, "bulk_rho_u", &scenario.bulk_rho_u,
                    0, description_bulk_rho_u);
    esio_line_write(h, "alpha", &scenario.alpha, 0, description_alpha);
    esio_line_write(h, "beta", &scenario.beta,   0, description_beta);
    esio_line_write(h, "gamma", &scenario.gamma, 0, description_gamma);
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
