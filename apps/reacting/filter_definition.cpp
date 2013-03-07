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
 * @copydoc filter_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "filter_definition.hpp"

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling problem defintion for filter
 * source term, e.g., strength coefficient.
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

filter_definition::filter_definition()
    : filter_phi(std::numeric_limits<real_t>::quiet_NaN())
{
}

filter_definition::filter_definition(
        const real_t filter_phi)
    : filter_phi(filter_phi)
{
}

filter_definition::~filter_definition()
{
    // NOP
}

// Strings used in options_description and populate/override/save/load.
static const char name_filter_phi[] = "filter_phi";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_filter_phi[] = "Filter source strength";

boost::program_options::options_description
filter_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;

    options_description retval("filter source parameters");

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    auto_ptr<typed_value<string> > p;

    // filter_phi
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &filter_phi, name_filter_phi));
    if (!(boost::math::isnan)(filter_phi)) {
        p->default_value(lexical_cast<string>(filter_phi));
    }
    retval.add_options()(name_filter_phi, p.release(), desc_filter_phi);

    return retval;
}

void
filter_definition::populate(
        const filter_definition& that,
        const bool verbose)
{
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(filter_phi);
#undef CALL_MAYBE_POPULATE
}

void
filter_definition::override(
        const filter_definition& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(filter_phi);
#undef CALL_MAYBE_OVERRIDE
}

void
filter_definition::save(
        const esio_handle h) const
{
    DEBUG0("Storing filter_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));

    esio_line_write(h, name_filter_phi, &this->filter_phi, 0, desc_filter_phi);
}

void
filter_definition::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading filter_definition parameters");

    // All ranks load
    esio_line_establish(h, 1, 0, 1);

    filter_definition t;
    esio_line_read(h, name_filter_phi, &t.filter_phi, 0);
    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace reacting

} // namespace suzerain
