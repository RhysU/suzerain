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
 * Provides classes handling problem defintion for channel
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

namespace reacting {

channel_definition::channel_definition()
    : bulk_rho  (std::numeric_limits<real_t>::quiet_NaN())
    , bulk_rho_u(std::numeric_limits<real_t>::quiet_NaN())
    , T_wall    (std::numeric_limits<real_t>::quiet_NaN())
    , wall_mass_fractions()
{
}

channel_definition::channel_definition(
        const real_t bulk_rho,
        const real_t bulk_rho_u,
	const real_t T_wall,
        const std::vector<real_t> wall_mass_fractions)
    : bulk_rho  (bulk_rho  )
    , bulk_rho_u(bulk_rho_u)
    , T_wall    (T_wall)
    , wall_mass_fractions(wall_mass_fractions)
{
}

channel_definition::~channel_definition()
{
    // NOP
}

// Strings used in options_description and populate/override/save/load.
static const char name_bulk_rho[]   = "bulk_rho";
static const char name_bulk_rho_u[] = "bulk_rho_u";
static const char name_T_wall[]     = "T_wall";
static const char name_mass_fractions_wall[] = "wall_mass_frac";

// Descriptions used in options_description and populate/override/save/load.
static const char desc_bulk_rho[]   = "Bulk density target";
static const char desc_bulk_rho_u[] = "Bulk momentum target";
static const char desc_T_wall[]     = "Wall temperature";
static const char desc_mass_fractions_wall[] = "Wall mass fractions";

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

    // T_wall
    p.reset(value<string>());
    p->notifier(bind(&parse_positive, _1, &T_wall, name_T_wall));
    if (!(boost::math::isnan)(T_wall)) {
        p->default_value(lexical_cast<string>(T_wall));
    }
    retval.add_options()(name_T_wall, p.release(), desc_T_wall);

    // FIXME: Check validity of incoming values
    // Must be non-negative and sum to 1
    retval.add_options()(name_mass_fractions_wall, 
                         value(&wall_mass_fractions), 
                         desc_mass_fractions_wall);


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
    CALL_MAYBE_POPULATE(T_wall);
#undef CALL_MAYBE_POPULATE
    if (this->wall_mass_fractions.size()==0)
        this->wall_mass_fractions = that.wall_mass_fractions;
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
    CALL_MAYBE_OVERRIDE(T_wall);
#undef CALL_MAYBE_OVERRIDE
    if (this->wall_mass_fractions.size()!=0)
        this->wall_mass_fractions = that.wall_mass_fractions;
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
    esio_line_write(h, name_T_wall,     &this->T_wall,     0, desc_T_wall);

    // mass fractions vector
    esio_line_establish(h, this->wall_mass_fractions.size(), 
                        0, (procid == 0 ? 1 : 0));
    esio_line_write(h, name_mass_fractions_wall, 
                    this->wall_mass_fractions.data(), 1, 
                    desc_mass_fractions_wall);
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
    esio_line_read(h, name_T_wall,     &t.T_wall,     0);

    // Mass fractions vector
    int Ns;
    esio_line_size(h, name_mass_fractions_wall, &Ns);
    t.wall_mass_fractions.resize(Ns);

    esio_line_establish(h, Ns, 0, Ns);
    esio_line_read(h, name_mass_fractions_wall, 
                   t.wall_mass_fractions.data(), 1);

    this->populate(t, verbose);  // Prefer this to incoming
}

} // namespace reacting

} // namespace suzerain
