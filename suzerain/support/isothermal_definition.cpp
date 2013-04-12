//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
 * @copydoc isothermal_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "isothermal_definition.hpp"

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

/** @file 
 * Provides classes handling problem defintion for isothermal
 * flow, e.g., bulk density and momentum.
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

isothermal_definition::isothermal_definition()
{
  // Call isothermal_specification constructor
  isothermal_specification();
}

// FIXME: Declare other constructors according to the 
//        parameteres passed


isothermal_definition::~isothermal_definition()
{
    // NOP
}

// Strings used in options_description and populate/override/save/load.
static const char name_lower_T[]   = "lower_T";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_lower_T[]   = "Input for temperature at lower boundary";

boost::program_options::options_description
isothermal_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;

    options_description retval("Isothermal boundary parameters");

    // Complicated add_options() calls done to allow changing the default value
    // displayed when the default is NaN.  NaN is used as a NOP value by client
    // code.  Validation routines used below all silently allow NaNs.

    auto_ptr<typed_value<string> > p;

    // lower_T
    p.reset(value<string>());
    p->notifier(bind(&parse_nonnegative, _1, &lower_T, name_lower_T));
    if (!(boost::math::isnan)(lower_T)) {
        p->default_value(lexical_cast<string>(lower_T));
    }
    retval.add_options()(name_lower_T, p.release(), desc_lower_T);

    return retval;
}

void
isothermal_definition::populate(
        const isothermal_definition& that,
        const bool verbose)
{
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(lower_T);
#undef CALL_MAYBE_POPULATE
}

void
isothermal_definition::override(
        const isothermal_definition& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(lower_T);
#undef CALL_MAYBE_OVERRIDE
}

void
isothermal_definition::save(
        const esio_handle h) const
{
    DEBUG0("Storing isothermal_definition parameters");

    // Only root writes data
    int procid;
    esio_handle_comm_rank(h, &procid);

    // scalars
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, name_lower_T,   &this->lower_T,   0, desc_lower_T);
}

void
isothermal_definition::load(
        const esio_handle h,
        const bool verbose)
{
    DEBUG0("Loading isothermal_definition parameters");

    isothermal_definition t;

    // All ranks load

    // Scalars
    esio_line_establish(h, 1, 0, 1);
    esio_line_read(h, name_lower_T,   &t.lower_T,   0);

    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace suzerain
