//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
 * @copydoc radialflow_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/radialflow_definition.hpp>

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

template<typename T>
static void parse_option(
        const std::string& s,
        T* value,
        const char* name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    *value = t;
}

template<typename T>
static void validate_option(
        const std::string& s,
        T* value,
        void (*validator)(T, const char*),
        const char* name)
{
#pragma warning(push,disable:2259)
    const T t = exprparse<real_t>(s, name);
#pragma warning(pop)
    validator(t, name);
    *value = t;
}

radialflow_definition::radialflow_definition(
            double deltae,
            double gam0,
            double Mae,
            double pexi,
            double Te)
    : radialflow_specification(deltae, gam0, Mae, pexi, Te)
{
}

// Strings used in options_description and populate/override/save/load.
static const char location   [] = "radialflow";
static const char name_deltae[] = "radialflow_deltae";
static const char name_gam0  [] = "radialflow_gam0";
static const char name_Mae   [] = "radialflow_Mae";
static const char name_pexi  [] = "radialflow_pexi";
static const char name_Te    [] = "radialflow_Te";
static const char * const attr_deltae = name_deltae + sizeof(location);
static const char * const attr_gam0   = name_gam0   + sizeof(location);
static const char * const attr_Mae    = name_Mae    + sizeof(location);
static const char * const attr_pexi   = name_pexi   + sizeof(location);
static const char * const attr_Te     = name_Te     + sizeof(location);

// Descriptions used in options_description and populate/override/save/load.
static const char desc_deltae[] = "Edge distance above the x-axis.";
static const char desc_gam0  [] = "Edge specific heat ratio."
                                  " Zero or NaN indicates an external value should be used.";
static const char desc_Mae   [] = "Edge Mach number."
                                  " Zero or NaN indicates an external value should be used.";
static const char desc_pexi  [] = "Edge pressure gradient parameter."
                                  " Zero or NaN disables the radial flow.";
static const char desc_Te    [] = "Edge temperature.";

boost::program_options::options_description
radialflow_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::options_description;
    using boost::program_options::typed_value;
    using boost::program_options::value;
    using std::auto_ptr;
    using std::string;
    using validation::ensure_nonnegative;
    using validation::ensure_positive;

    options_description retval(
            "Radial flow problem reference parameters"
            " for homogenized boundary layers with baseflow");

    auto_ptr<typed_value<string> > p;

    // deltae
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &deltae,
                     &ensure_positive<double>, name_deltae));
    if (!(boost::math::isnan)(deltae)) {
        p->default_value(lexical_cast<string>(deltae));
    }
    retval.add_options()(name_deltae, p.release(), desc_deltae);

    // gam0
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &gam0,
                     &ensure_positive<double>, name_gam0));
    if (!(boost::math::isnan)(gam0)) {
        p->default_value(lexical_cast<string>(gam0));
    }
    retval.add_options()(name_gam0, p.release(), desc_gam0);

    // Mae
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &Mae,
                     &ensure_positive<double>, name_Mae));
    if (!(boost::math::isnan)(Mae)) {
        p->default_value(lexical_cast<string>(Mae));
    }
    retval.add_options()(name_Mae, p.release(), desc_Mae);

    // pexi
    p.reset(value<string>());
    p->notifier(bind(&parse_option<double>, _1, &pexi, name_pexi));
    if (!(boost::math::isnan)(pexi)) {
        p->default_value(lexical_cast<string>(pexi));
    }
    retval.add_options()(name_pexi, p.release(), desc_pexi);

    // Te
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &Te,
                     &ensure_positive<double>, name_Te));
    if (!(boost::math::isnan)(Te)) {
        p->default_value(lexical_cast<string>(Te));
    }
    retval.add_options()(name_Te, p.release(), desc_Te);

    return retval;
}

void
radialflow_definition::populate(
        const radialflow_specification& that,
        const bool verbose)
{
    maybe_populate(name_deltae, desc_deltae, deltae, that.deltae, verbose);
    maybe_populate(name_gam0,   desc_gam0,   gam0,   that.gam0,   verbose);
    maybe_populate(name_Mae,    desc_Mae,    Mae,    that.Mae,    verbose);
    maybe_populate(name_pexi,   desc_pexi,   pexi,   that.pexi,   verbose);
    maybe_populate(name_Te,     desc_Te,     Te,     that.Te,     verbose);
}

void
radialflow_definition::override(
        const radialflow_specification& that,
        const bool verbose)
{
    maybe_override(name_deltae, desc_deltae, deltae, that.deltae, verbose);
    maybe_override(name_gam0,   desc_gam0,   gam0,   that.gam0,   verbose);
    maybe_override(name_Mae,    desc_Mae,    Mae,    that.Mae,    verbose);
    maybe_override(name_pexi,   desc_pexi,   pexi,   that.pexi,   verbose);
    maybe_override(name_Te,     desc_Te,     Te,     that.Te,     verbose);
}

void
radialflow_definition::save(
    const esio_handle h) const
{
    // Save nothing if there's nothing interesting to save
    if (this->trivial()) {
        return;
    }

    DEBUG0("Storing radialflow_definition parameters");

    // Write out the "container" holding all other settings
    const int one = 1;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &one, 0,
                    "Is a radial flow definition in use?");

    // Write out information as attributes within the container
    esio_attribute_write(h, location, attr_deltae, &deltae);
    esio_attribute_write(h, location, attr_gam0,   &gam0);
    esio_attribute_write(h, location, attr_Mae,    &Mae);
    esio_attribute_write(h, location, attr_pexi,   &pexi);
    esio_attribute_write(h, location, attr_Te,     &Te);
}

void
radialflow_definition::load(
    const esio_handle h,
    const bool verbose)
{
    radialflow_definition t;

    // Only proceed if a definition is active in the restart
    int in_use = 0;
    esio_line_establish(h, 1, 0, 1); // All ranks load any data
    if (ESIO_NOTFOUND != esio_line_size(h, location, NULL)) {
        esio_line_read(h, location, &in_use, 0);
    }
    if (!in_use) {
        return;
    }

    DEBUG0("Loading radialflow_definition parameters");

    // Read in information as attributes within the container
    esio_attribute_read(h, location, attr_deltae, &t.deltae);
    esio_attribute_read(h, location, attr_gam0,   &t.gam0);
    esio_attribute_read(h, location, attr_Mae,    &t.Mae);
    esio_attribute_read(h, location, attr_pexi,   &t.pexi);
    esio_attribute_read(h, location, attr_Te,     &t.Te);

    return this->populate(t, verbose);
}

} // end namespace support

} // end namespace suzerain
