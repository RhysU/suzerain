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
 * @copydoc definition_channel.hpp
 */

#include "definition_channel.hpp"

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

namespace reacting {

channel_definition::channel_definition()
    : bulk_rho  (std::numeric_limits<real_t>::quiet_NaN())
    , bulk_rho_u(std::numeric_limits<real_t>::quiet_NaN())
    , bulk_rho_E(std::numeric_limits<real_t>::quiet_NaN())
{
}

channel_definition::channel_definition(
        const real_t bulk_rho,
        const real_t bulk_rho_u,
        const real_t bulk_rho_E)
    : bulk_rho  (bulk_rho  )
    , bulk_rho_u(bulk_rho_u)
    , bulk_rho_E(bulk_rho_E)
{
}

// Strings used in options_description and populate/override/save/load.
static const char name_bulk_rho[]            = "bulk_rho";
static const char name_bulk_rho_u[]          = "bulk_rho_u";
static const char name_bulk_rho_E[]          = "bulk_rho_E";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_bulk_rho[]            = "Bulk density target";
static const char desc_bulk_rho_u[]          = "Bulk momentum target";
static const char desc_bulk_rho_E[]          = "Bulk total energy target";

boost::program_options::options_description
channel_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;

    options_description retval("Channel flow scenario parameters");

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    auto_ptr<typed_value<string> > p;

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

    return retval;
}

void
channel_definition::populate(
        const channel_definition& that,
        const bool verbose)
{
    using support::maybe_populate;
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(bulk_rho);
    CALL_MAYBE_POPULATE(bulk_rho_u);
    CALL_MAYBE_POPULATE(bulk_rho_E);
#undef CALL_MAYBE_POPULATE
}

void
channel_definition::override(
        const channel_definition& that,
        const bool verbose)
{
    using support::maybe_override;
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(bulk_rho);
    CALL_MAYBE_OVERRIDE(bulk_rho_u);
    CALL_MAYBE_OVERRIDE(bulk_rho_E);
#undef CALL_MAYBE_OVERRIDE
}

void
channel_definition::save(
        const esio_handle h) const
{
    DEBUG0("Storing channel_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    // scalars
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, name_bulk_rho,   &this->bulk_rho,   0, desc_bulk_rho);
    esio_line_write(h, name_bulk_rho_u, &this->bulk_rho_u, 0, desc_bulk_rho_u);
    esio_line_write(h, name_bulk_rho_E, &this->bulk_rho_E, 0, desc_bulk_rho_E);
}

void
channel_definition::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading channel_definition parameters");

    channel_definition t;

    // All ranks load

    // Scalars
    esio_line_establish(h, 1, 0, 1);
    esio_line_read(h, name_bulk_rho,   &t.bulk_rho,   0);
    esio_line_read(h, name_bulk_rho_u, &t.bulk_rho_u, 0);
    if (ESIO_NOTFOUND != esio_line_size(h, name_bulk_rho_E, NULL)) {
        esio_line_read(h, name_bulk_rho_E, &t.bulk_rho_E, 0);
    } else {
        INFO0(desc_bulk_rho_E << " not found; defaulting to disabled");
        t.bulk_rho_E = std::numeric_limits<real_t>::quiet_NaN();
    }

    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace reacting

} // namespace suzerain
