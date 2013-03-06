//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc channel_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "channel_definition.hpp"

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling channel flow scenario parameters 
 * which describe the particular channel flow being simulated.
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

namespace reacting {

channel_definition::channel_definition()
    : bulk_rho  (std::numeric_limits<real_t>::quiet_NaN())
    , bulk_rho_u(std::numeric_limits<real_t>::quiet_NaN())
{
}

channel_definition::channel_definition(
        const real_t bulk_rho,
        const real_t bulk_rho_u)
    : bulk_rho  (bulk_rho  )
    , bulk_rho_u(bulk_rho_u)
{
}

channel_definition::~channel_definition()
{
    // NOP
}

// Strings used in options_description and populate/override/save/load.
static const char name_bulk_rho[]   = "bulk_rho";
static const char name_bulk_rho_u[] = "bulk_rho_u";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_bulk_rho[]   = "Bulk density target";
static const char desc_bulk_rho_u[] = "Bulk momentum target";

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

    return retval;
}

void
channel_definition::populate(
        const channel_definition& that,
        const bool verbose)
{
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(bulk_rho);
    CALL_MAYBE_POPULATE(bulk_rho_u);
#undef CALL_MAYBE_POPULATE
}

void
channel_definition::override(
        const channel_definition& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(bulk_rho);
    CALL_MAYBE_OVERRIDE(bulk_rho_u);
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
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, name_bulk_rho,   &this->bulk_rho,   0, desc_bulk_rho);
    esio_line_write(h, name_bulk_rho_u, &this->bulk_rho_u, 0, desc_bulk_rho_u);
}

void
channel_definition::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading channel_definition parameters");

    // All ranks load
    esio_line_establish(h, 1, 0, 1);

    channel_definition t;
    esio_line_read(h, name_bulk_rho,   &t.bulk_rho,   0);
    esio_line_read(h, name_bulk_rho_u, &t.bulk_rho_u, 0);
    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace reacting

} // namespace suzerain
