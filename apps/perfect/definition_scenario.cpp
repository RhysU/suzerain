//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc definition_scenario.hpp
 */

#include "definition_scenario.hpp"

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

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
    , bulk_rho_E(std::numeric_limits<real_t>::quiet_NaN())
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
        const real_t bulk_rho_E,
        const real_t alpha,
        const real_t beta,
        const real_t gamma)
    : Re        (Re        )
    , Ma        (Ma        )
    , Pr        (Pr        )
    , bulk_rho  (bulk_rho  )
    , bulk_rho_u(bulk_rho_u)
    , bulk_rho_E(bulk_rho_E)
    , alpha     (alpha     )
    , beta      (beta      )
    , gamma     (gamma     )
{
}

scenario_definition::~scenario_definition()
{
}

// Strings used in options_description and populate/override/save/load.
static const char name_Re[]         = "Re";
static const char name_Ma[]         = "Ma";
static const char name_Pr[]         = "Pr";
static const char name_bulk_rho[]   = "bulk_rho";
static const char name_bulk_rho_u[] = "bulk_rho_u";
static const char name_bulk_rho_E[] = "bulk_rho_E";
static const char name_alpha[]      = "alpha";
static const char name_beta[]       = "beta";
static const char name_gamma[]      = "gamma";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_Re[]         = "Reynolds number";
static const char desc_Ma[]         = "Mach number";
static const char desc_Pr[]         = "Prandtl number";
static const char desc_bulk_rho[]   = "Bulk density target";
static const char desc_bulk_rho_u[] = "Bulk streamwise momentum target";
static const char desc_bulk_rho_E[] = "Bulk total energy target";
static const char desc_alpha[]      = "Ratio of bulk to dynamic viscosity";
static const char desc_beta[]       = "Temperature power law exponent";
static const char desc_gamma[]      = "Ratio of specific heats";

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
    p->notifier(bind(&parse_positive, _1, &Re, name_Re));
    if (!(boost::math::isnan)(Re)) {
        p->default_value(lexical_cast<string>(Re));
    }
    retval.add_options()(name_Re, p.release(), desc_Re);

    // Ma
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Ma, name_Ma));
    if (!(boost::math::isnan)(Ma)) {
        p->default_value(lexical_cast<string>(Ma));
    }
    retval.add_options()(name_Ma, p.release(), desc_Ma);

    // Pr
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &Pr, name_Pr));
    if (!(boost::math::isnan)(Pr)) {
        p->default_value(lexical_cast<string>(Pr));
    }
    retval.add_options()(name_Pr, p.release(), desc_Pr);

    // bulk_rho
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &bulk_rho, name_bulk_rho));
    if (!(boost::math::isnan)(bulk_rho)) {
        p->default_value(lexical_cast<string>(bulk_rho));
    }
    retval.add_options()(name_bulk_rho, p.release(), desc_bulk_rho);

    // bulk_rho_u
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &bulk_rho_u, name_bulk_rho_u));
    if (!(boost::math::isnan)(bulk_rho_u)) {
        p->default_value(lexical_cast<string>(bulk_rho_u));
    }
    retval.add_options()(name_bulk_rho_u, p.release(), desc_bulk_rho_u);

    // bulk_rho_E
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &bulk_rho_E, name_bulk_rho_E));
    if (!(boost::math::isnan)(bulk_rho_E)) {
        p->default_value(lexical_cast<string>(bulk_rho_E));
    }
    retval.add_options()(name_bulk_rho_E, p.release(), desc_bulk_rho_E);

    // alpha
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &alpha, name_alpha));
    if (!(boost::math::isnan)(alpha)) {
        p->default_value(lexical_cast<string>(alpha));
    }
    retval.add_options()(name_alpha, p.release(), desc_alpha);

    // beta
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &beta, name_beta));
    if (!(boost::math::isnan)(beta)) {
        p->default_value(lexical_cast<string>(beta));
    }
    retval.add_options()(name_beta, p.release(), desc_beta);

    // gamma
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &gamma, name_gamma));
    if (!(boost::math::isnan)(gamma)) {
        p->default_value(lexical_cast<string>(gamma));
    }
    retval.add_options()(name_gamma, p.release(), desc_gamma);

    return retval;
}

void
scenario_definition::populate(
        const scenario_definition& that,
        const bool verbose)
{
    using support::maybe_populate;
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(Re);
    CALL_MAYBE_POPULATE(Ma);
    CALL_MAYBE_POPULATE(Pr);
    CALL_MAYBE_POPULATE(bulk_rho);
    CALL_MAYBE_POPULATE(bulk_rho_u);
    CALL_MAYBE_POPULATE(bulk_rho_E);
    CALL_MAYBE_POPULATE(alpha);
    CALL_MAYBE_POPULATE(beta);
    CALL_MAYBE_POPULATE(gamma);
#undef CALL_MAYBE_POPULATE
}

void
scenario_definition::override(
        const scenario_definition& that,
        const bool verbose)
{
    using support::maybe_override;
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(Re);
    CALL_MAYBE_OVERRIDE(Ma);
    CALL_MAYBE_OVERRIDE(Pr);
    CALL_MAYBE_OVERRIDE(bulk_rho);
    CALL_MAYBE_OVERRIDE(bulk_rho_u);
    CALL_MAYBE_OVERRIDE(bulk_rho_E);
    CALL_MAYBE_OVERRIDE(alpha);
    CALL_MAYBE_OVERRIDE(beta);
    CALL_MAYBE_OVERRIDE(gamma);
#undef CALL_MAYBE_OVERRIDE
}

void
scenario_definition::save(
        const esio_handle h) const
{
    DEBUG0("Storing scenario_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, name_Re,         &this->Re,         0, desc_Re);
    esio_line_write(h, name_Ma,         &this->Ma,         0, desc_Ma);
    esio_line_write(h, name_Pr,         &this->Pr,         0, desc_Pr);
    esio_line_write(h, name_bulk_rho,   &this->bulk_rho,   0, desc_bulk_rho);
    esio_line_write(h, name_bulk_rho_u, &this->bulk_rho_u, 0, desc_bulk_rho_u);
    esio_line_write(h, name_bulk_rho_E, &this->bulk_rho_E, 0, desc_bulk_rho_E);
    esio_line_write(h, name_alpha,      &this->alpha,      0, desc_alpha);
    esio_line_write(h, name_beta,       &this->beta,       0, desc_beta);
    esio_line_write(h, name_gamma,      &this->gamma,      0, desc_gamma);
}

void
scenario_definition::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading scenario_definition parameters");

    // All ranks load
    esio_line_establish(h, 1, 0, 1);

    scenario_definition t;
    esio_line_read(h, name_Re,         &t.Re,         0);
    esio_line_read(h, name_Ma,         &t.Ma,         0);
    esio_line_read(h, name_Pr,         &t.Pr,         0);
    esio_line_read(h, name_bulk_rho,   &t.bulk_rho,   0);
    esio_line_read(h, name_bulk_rho_u, &t.bulk_rho_u, 0);
    if (ESIO_NOTFOUND != esio_line_size(h, name_bulk_rho_E, NULL)) {
        esio_line_read(h, name_bulk_rho_E, &t.bulk_rho_E, 0);
    } else {
        INFO0(desc_bulk_rho_E << " not found; defaulting to disabled");
        t.bulk_rho_E = std::numeric_limits<real_t>::quiet_NaN();
    }
    esio_line_read(h, name_alpha,      &t.alpha,      0);
    esio_line_read(h, name_beta,       &t.beta,       0);
    esio_line_read(h, name_gamma,      &t.gamma,      0);
    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace perfect

} // namespace suzerain
