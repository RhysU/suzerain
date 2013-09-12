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

#include <suzerain/support/isothermal_definition.hpp>

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/common.hpp>
#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

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

static void parse_bounded(const std::string& s,
                          real_t *t,
                          real_t lower,
                          real_t upper,
                          const char *n)
{
    const real_t v = exprparse<real_t>(s, n);
    validation::ensure_bounded(v, lower, upper, true, true, n);
    *t = v;
}

isothermal_definition::isothermal_definition()
    : isothermal_specification()
{
}

isothermal_definition::isothermal_definition(
        real_t wall_T)
    : isothermal_specification(wall_T)
{
}

isothermal_definition::isothermal_definition(
        real_t wall_T,
        const std::vector<real_t>& wall_cs)
    : isothermal_specification(wall_T, wall_cs)
{
}

isothermal_definition::isothermal_definition(
        real_t wall_T,
        real_t inflow_velocity)
    : isothermal_specification(wall_T, inflow_velocity)
{
}

isothermal_definition::isothermal_definition(
        real_t wall_T,
        real_t inflow_velocity,
        const std::vector<real_t>& wall_cs)
    : isothermal_specification(wall_T, inflow_velocity, wall_cs)
{
}

isothermal_definition::isothermal_definition(
        real_t lower_T,
        real_t lower_v,
        real_t upper_T,
        real_t upper_v)
    : isothermal_specification(lower_T, lower_v, upper_T, upper_v)
{
}

isothermal_definition::isothermal_definition(
        real_t lower_T,
        real_t lower_v,
        const std::vector<real_t>& lower_cs,
        real_t upper_T,
        real_t upper_v,
        const std::vector<real_t>& upper_cs)
    : isothermal_specification(lower_T, lower_v, lower_cs,
                               upper_T, upper_v, upper_cs)
{
}

isothermal_definition::isothermal_definition(
        real_t lower_T,
        real_t lower_v,
        real_t lower_rho,
        const std::vector<real_t>& lower_cs,
        real_t upper_T,
        real_t upper_v,
        real_t upper_rho,
        const std::vector<real_t>& upper_cs)
    : isothermal_specification(lower_T, lower_v, lower_rho, lower_cs,
                               upper_T, upper_v, upper_rho, upper_cs)
{
}

// Strings used in options_description and populate/override/save/load
static const char name_lower_T  []  = "lower_T";
static const char name_lower_u  []  = "lower_u";
static const char name_lower_v  []  = "lower_v";
static const char name_lower_w  []  = "lower_w";
static const char name_lower_rho[]  = "lower_rho";
static const char name_lower_cs []  = "lower_cs";
static const char name_upper_T  []  = "upper_T";
static const char name_upper_u  []  = "upper_u";
static const char name_upper_v  []  = "upper_v";
static const char name_upper_w  []  = "upper_w";
static const char name_upper_rho[]  = "upper_rho";
static const char name_upper_cs []  = "upper_cs";

// Descriptions used in options_description and populate/override/save/load
static const char desc_lower_T  []  = "Temperature at lower boundary";
static const char desc_lower_u  []  = "Streamwise velocity at lower boundary";
static const char desc_lower_v  []  = "Wall-normal velocity at lower boundary";
static const char desc_lower_w  []  = "Spanwise velocity at lower boundary";
static const char desc_lower_rho[]  = "Density at lower boundary";
static const char desc_lower_cs []  = "Species mass fractions at lower boundary";
static const char desc_upper_T  []  = "Temperature at upper boundary";
static const char desc_upper_u  []  = "Streamwise velocity at upper boundary";
static const char desc_upper_v  []  = "Wall-normal velocity at upper boundary";
static const char desc_upper_w  []  = "Spanwise at upper boundary";
static const char desc_upper_rho[]  = "Density at upper boundary";
static const char desc_upper_cs []  = "Species mass fractions at upper boundary";

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

    // lower_rho
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &lower_rho, name_lower_rho));
    if (!(boost::math::isnan)(lower_rho)) {
        p->default_value(lexical_cast<string>(lower_rho));
    }
    retval.add_options()(name_lower_rho, p.release(), desc_lower_rho);

    // lower_cs
    // FIXME: suzerain::exprparse and check validity of incoming values
    // Must be in the range [0, 1] per parse_bounded.
    // Should sum to one, but floating point makes that tricky to check.
    (void) parse_bounded;
    retval.add_options()(name_lower_cs, value(&lower_cs), desc_lower_cs);

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

    // upper_rho
    p.reset(value<string>());
    p->notifier(bind(&parse, _1, &upper_rho, name_upper_rho));
    if (!(boost::math::isnan)(upper_rho)) {
        p->default_value(lexical_cast<string>(upper_rho));
    }
    retval.add_options()(name_upper_rho, p.release(), desc_upper_rho);

    // upper_cs
    // FIXME: suzerain::exprparse and check validity of incoming values
    // Must be in the range [0, 1] per parse_bounded.
    // Should sum to one, but floating point makes that tricky to check.
    (void) parse_bounded;
    retval.add_options()(name_upper_cs, value(&upper_cs), desc_upper_cs);

    return retval;
}

void
isothermal_definition::populate(
        const isothermal_specification& that,
        const bool verbose)
{
#define CALL_MAYBE_POPULATE(mem)                                             \
    maybe_populate(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_POPULATE(lower_T);
    CALL_MAYBE_POPULATE(lower_u);
    CALL_MAYBE_POPULATE(lower_v);
    CALL_MAYBE_POPULATE(lower_w);
    CALL_MAYBE_POPULATE(lower_rho);
    if (lower_cs.empty()) lower_cs = that.lower_cs;

    CALL_MAYBE_POPULATE(upper_T);
    CALL_MAYBE_POPULATE(upper_u);
    CALL_MAYBE_POPULATE(upper_v);
    CALL_MAYBE_POPULATE(upper_w);
    CALL_MAYBE_POPULATE(upper_rho);
    if (upper_cs.empty()) upper_cs = that.upper_cs;
#undef CALL_MAYBE_POPULATE
}

void
isothermal_definition::override(
        const isothermal_specification& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem)                                            \
    maybe_override(name_ ## mem, desc_ ## mem, this->mem, that.mem, verbose)
    CALL_MAYBE_OVERRIDE(lower_T);
    CALL_MAYBE_OVERRIDE(lower_u);
    CALL_MAYBE_OVERRIDE(lower_v);
    CALL_MAYBE_OVERRIDE(lower_w);
    CALL_MAYBE_OVERRIDE(lower_rho);
    if (!that.lower_cs.empty()) lower_cs = that.lower_cs;

    CALL_MAYBE_OVERRIDE(upper_T);
    CALL_MAYBE_OVERRIDE(upper_u);
    CALL_MAYBE_OVERRIDE(upper_v);
    CALL_MAYBE_OVERRIDE(upper_w);
    CALL_MAYBE_OVERRIDE(upper_rho);
    if (!that.upper_cs.empty()) upper_cs = that.upper_cs;
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
    esio_line_write(h, name_lower_T,   &lower_T,   0, desc_lower_T);
    esio_line_write(h, name_lower_u,   &lower_u,   0, desc_lower_u);
    esio_line_write(h, name_lower_v,   &lower_v,   0, desc_lower_v);
    esio_line_write(h, name_lower_w,   &lower_w,   0, desc_lower_w);
    esio_line_write(h, name_lower_rho, &lower_rho, 0, desc_lower_rho);
    esio_line_write(h, name_upper_T,   &upper_T,   0, desc_upper_T);
    esio_line_write(h, name_upper_u,   &upper_u,   0, desc_upper_u);
    esio_line_write(h, name_upper_v,   &upper_v,   0, desc_upper_v);
    esio_line_write(h, name_upper_w,   &upper_w,   0, desc_upper_w);
    esio_line_write(h, name_upper_rho, &upper_rho, 0, desc_upper_rho);

    // Lower mass fractions vector (only written when non-trivial)
    int Ns;
    if ((Ns = lower_cs.size())) {
        esio_line_establish(h, Ns, 0, (procid == 0 ? Ns : 0));
        esio_line_write(h, name_lower_cs, lower_cs.data(), 1, desc_lower_cs);
    } else {
        WARN0("No lower mass fractions saved because lower_cs.size() == 0");
    }

    // Upper mass fractions vector (ditto, could be of different length)
    if ((Ns = upper_cs.size())) {
        esio_line_establish(h, Ns, 0, (procid == 0 ? Ns : 0));
        esio_line_write(h, name_upper_cs, upper_cs.data(), 1, desc_upper_cs);
    } else {
        WARN0("No upper mass fractions saved because upper_cs.size() == 0");
    }
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

    // Backwards compatibility:
    // The old, unspecified default for apps/perfect was to use T_wall = 1.
    // The old, apps/reacting code saved lower_T = upper_T = T_wall in "/T_wall".
    // The loading logic honors those two older behaviors in a verbose manner.
    real_t default_T_wall = 1;
    if (ESIO_NOTFOUND != esio_line_size(h, "T_wall", NULL)) {
        esio_line_read(h, "T_wall", &default_T_wall, 0);
        INFO0("Found old-style wall temperature data = " << default_T_wall);
    }

#define LOAD_OR_DEFAULT(var, default_value)                              \
    do {                                                                 \
        if (ESIO_NOTFOUND == esio_line_size(h, name_##var, NULL)) {      \
            t.var = default_value;                                       \
            INFO0(desc_##var << " not found therefore using " << t.var); \
        } else {                                                         \
            esio_line_read(h, name_##var, &t.var, 0);                    \
        }                                                                \
    } while(0)

    LOAD_OR_DEFAULT(lower_T,   default_T_wall);
    LOAD_OR_DEFAULT(lower_u,                0);
    LOAD_OR_DEFAULT(lower_v,                0);
    LOAD_OR_DEFAULT(lower_w,                0);
    LOAD_OR_DEFAULT(lower_rho, std::numeric_limits<real_t>::quiet_NaN());
    LOAD_OR_DEFAULT(upper_T,   default_T_wall);
    LOAD_OR_DEFAULT(upper_u,                0);
    LOAD_OR_DEFAULT(upper_v,                0);
    LOAD_OR_DEFAULT(upper_w,                0);
    LOAD_OR_DEFAULT(upper_rho, 1);

#undef LOAD_OR_DEFAULT

    // Lower mass fractions vector (only present when non-trivial)
    int Ns;
    if (ESIO_NOTFOUND == esio_line_size(h, name_lower_cs, &Ns)) {
        INFO0(name_lower_cs << " not found therefore using 1 dilluter species");
        t.lower_cs.assign(1U, 1.0);
    } else {
        t.lower_cs.resize(Ns);
        esio_line_establish(h, Ns, 0, Ns);
        esio_line_read(h, name_lower_cs, t.lower_cs.data(), 0);
    }

    // Upper mass fractions vector (ditto, could be be of different length)
    if (ESIO_NOTFOUND == esio_line_size(h, name_upper_cs, &Ns)) {
        INFO0(name_upper_cs << " not found therefore using 1 dilluter species");
        t.upper_cs.assign(1U, 1.0);
    } else {
        t.upper_cs.resize(Ns);
        esio_line_establish(h, Ns, 0, Ns);
        esio_line_read(h, name_upper_cs, t.upper_cs.data(), 0);
    }

    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace support

} // namespace suzerain
