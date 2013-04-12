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

  
static void parse(const std::string& s, real_t *t, const char *n)
{
    const real_t v = exprparse<real_t>(s, n);
    *t = v;
}

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

// FIXME: Add a parse_bounded for species inputs,
//        maybe "ensure_zero_one"


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
static const char name_lower_T []  = "lower_T";
static const char name_lower_u []  = "lower_u";
static const char name_lower_v []  = "lower_v";
static const char name_lower_w []  = "lower_w";
static const char name_lower_cs[]  = "lower_cs";
static const char name_upper_T []  = "upper_T";
static const char name_upper_u []  = "upper_u";
static const char name_upper_v []  = "upper_v";
static const char name_upper_w []  = "upper_w";
static const char name_upper_cs[]  = "upper_cs";


// Descriptions used in options_description and populate/override/save/load.
static const char desc_lower_T []  = "Input for temperature at lower boundary";
static const char desc_lower_u []  = "Input for u-velocity at lower boundary";
static const char desc_lower_v []  = "Input for v-velocity at lower boundary";
static const char desc_lower_w []  = "Input for w-velocity at lower boundary";
static const char desc_lower_cs[]  = "Input for species mass fractions at \
                                      lower boundary";
static const char desc_upper_T []  = "Input for temperature at upper boundary";
static const char desc_upper_u []  = "Input for u-velocity at upper boundary";
static const char desc_upper_v []  = "Input for v-velocity at upper boundary";
static const char desc_upper_w []  = "Input for w-velocity at upper boundary";
static const char desc_upper_cs[]  = "Input for species mass fractions at \
                                      upper boundary";


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

    // FIXME: review which values will be accepted for each input

    // lower_T
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &lower_T, name_lower_T));
    if (!(boost::math::isnan)(lower_T)) {
        p->default_value(lexical_cast<string>(lower_T));
    }
    retval.add_options()(name_lower_T, p.release(), desc_lower_T);

    // lower_u
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &lower_u, name_lower_u));
    if (!(boost::math::isnan)(lower_u)) {
        p->default_value(lexical_cast<string>(lower_u));
    }
    retval.add_options()(name_lower_u, p.release(), desc_lower_u);

    // lower_v
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &lower_v, name_lower_v));
    if (!(boost::math::isnan)(lower_v)) {
        p->default_value(lexical_cast<string>(lower_v));
    }
    retval.add_options()(name_lower_v, p.release(), desc_lower_v);

    // lower_w
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &lower_w, name_lower_w));
    if (!(boost::math::isnan)(lower_w)) {
        p->default_value(lexical_cast<string>(lower_w));
    }
    retval.add_options()(name_lower_w, p.release(), desc_lower_w);

    // FIXME: Check validity of incoming values
    // Must be non-negative and sum to 1
    // lower_cs
    retval.add_options()(name_lower_cs, 
                         value(&lower_cs), 
                         desc_lower_cs);

    // upper_T
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &upper_T, name_upper_T));
    if (!(boost::math::isnan)(upper_T)) {
        p->default_value(lexical_cast<string>(upper_T));
    }
    retval.add_options()(name_upper_T, p.release(), desc_upper_T);

    // upper_u
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &upper_u, name_upper_u));
    if (!(boost::math::isnan)(upper_u)) {
        p->default_value(lexical_cast<string>(upper_u));
    }
    retval.add_options()(name_upper_u, p.release(), desc_upper_u);

    // upper_v
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &upper_v, name_upper_v));
    if (!(boost::math::isnan)(upper_v)) {
        p->default_value(lexical_cast<string>(upper_v));
    }
    retval.add_options()(name_upper_v, p.release(), desc_upper_v);

    // upper_w
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &upper_w, name_upper_w));
    if (!(boost::math::isnan)(upper_w)) {
        p->default_value(lexical_cast<string>(upper_w));
    }
    retval.add_options()(name_upper_w, p.release(), desc_upper_w);

    // FIXME: Check validity of incoming values
    // Must be non-negative and sum to 1
    // upper_cs
    retval.add_options()(name_upper_cs, 
                         value(&upper_cs), 
                         desc_upper_cs);

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
    CALL_MAYBE_POPULATE(lower_u);
    CALL_MAYBE_POPULATE(lower_v);
    CALL_MAYBE_POPULATE(lower_w);
    CALL_MAYBE_POPULATE(upper_T);
    CALL_MAYBE_POPULATE(upper_u);
    CALL_MAYBE_POPULATE(upper_v);
    CALL_MAYBE_POPULATE(upper_w);
#undef CALL_MAYBE_POPULATE
    if (this->lower_cs.size()!=0)
        this->lower_cs = that.lower_cs;
    if (this->upper_cs.size()!=0)
        this->upper_cs = that.upper_cs;
}

void
isothermal_definition::override(
        const isothermal_definition& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(lower_T);
    CALL_MAYBE_OVERRIDE(lower_u);
    CALL_MAYBE_OVERRIDE(lower_v);
    CALL_MAYBE_OVERRIDE(lower_w);
    CALL_MAYBE_OVERRIDE(upper_T);
    CALL_MAYBE_OVERRIDE(upper_u);
    CALL_MAYBE_OVERRIDE(upper_v);
    CALL_MAYBE_OVERRIDE(upper_w);
#undef CALL_MAYBE_OVERRIDE
    if (this->lower_cs.size()!=0)
        this->lower_cs = that.lower_cs;
    if (this->upper_cs.size()!=0)
        this->upper_cs = that.upper_cs;
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
    esio_line_write(h, name_lower_u,   &this->lower_u,   0, desc_lower_u);
    esio_line_write(h, name_lower_v,   &this->lower_v,   0, desc_lower_v);
    esio_line_write(h, name_lower_w,   &this->lower_w,   0, desc_lower_w);
    esio_line_write(h, name_upper_T,   &this->upper_T,   0, desc_upper_T);
    esio_line_write(h, name_upper_u,   &this->upper_u,   0, desc_upper_u);
    esio_line_write(h, name_upper_v,   &this->upper_v,   0, desc_upper_v);
    esio_line_write(h, name_upper_w,   &this->upper_w,   0, desc_upper_w);

    // mass fractions vector
    int Ns = this->lower_cs.size();
    esio_line_establish(h, Ns, 0, (procid == 0 ? Ns : 0));
    esio_line_write(h, name_lower_cs, 
                    this->lower_cs.data(), 1, 
                    desc_lower_cs);
    esio_line_write(h, name_upper_cs, 
                    this->upper_cs.data(), 1, 
                    desc_upper_cs);
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
    esio_line_read(h, name_lower_u,   &t.lower_u,   0);
    esio_line_read(h, name_lower_v,   &t.lower_v,   0);
    esio_line_read(h, name_lower_w,   &t.lower_w,   0);
    esio_line_read(h, name_upper_T,   &t.upper_T,   0);
    esio_line_read(h, name_upper_u,   &t.upper_u,   0);
    esio_line_read(h, name_upper_v,   &t.upper_v,   0);
    esio_line_read(h, name_upper_w,   &t.upper_w,   0);

    // Mass fractions vector
    int Ns;
    esio_line_size(h, name_lower_cs, &Ns);
    t.lower_cs.resize(Ns);

    esio_line_establish(h, Ns, 0, Ns);
    esio_line_read(h, name_lower_cs, 
                   t.lower_cs.data(), 1);

    esio_line_size(h, name_upper_cs, &Ns);
    t.upper_cs.resize(Ns);

    esio_line_establish(h, Ns, 0, Ns);
    esio_line_read(h, name_upper_cs, 
                   t.upper_cs.data(), 1);

    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace suzerain
