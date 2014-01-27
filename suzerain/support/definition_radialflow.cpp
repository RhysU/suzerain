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
 * @copydoc definition_radialflow.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/definition_radialflow.hpp>

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

definition_radialflow::definition_radialflow(
            double deltae,
            double gamma,
            double Mae,
            double pexi)
    : specification_radialflow(deltae, gamma, Mae, pexi)
{
}

// Strings used in options_description and populate/override/save/load.
static const char location   [] = "radialflow";
static const char name_deltae[] = "radialflow_deltae";
static const char name_gamma [] = "radialflow_gamma";
static const char name_Mae   [] = "radialflow_Mae";
static const char name_pexi  [] = "radialflow_pexi";
static const char * const attr_deltae = name_deltae + sizeof(location);
static const char * const attr_gamma  = name_gamma  + sizeof(location);
static const char * const attr_Mae    = name_Mae    + sizeof(location);
static const char * const attr_pexi   = name_pexi   + sizeof(location);

// Descriptions used in options_description and populate/override/save/load.
static const char desc_deltae[] = "Edge distance above the x-axis.";
static const char desc_gamma [] = "Edge specific heat ratio."
                                  " Zero or NaN indicates an external value should be used.";
static const char desc_Mae   [] = "Edge Mach number."
                                  " Zero or NaN indicates an external value should be used.";
static const char desc_pexi  [] = "Edge pressure gradient parameter."
                                  " Zero or NaN disables the radial flow.";

boost::program_options::options_description
definition_radialflow::options_description()
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

    // gamma
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &gamma,
                     &ensure_positive<double>, name_gamma));
    if (!(boost::math::isnan)(gamma)) {
        p->default_value(lexical_cast<string>(gamma));
    }
    retval.add_options()(name_gamma, p.release(), desc_gamma);

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

    return retval;
}

void
definition_radialflow::populate(
        const specification_radialflow& that,
        const bool verbose)
{
    maybe_populate(name_deltae, desc_deltae, deltae, that.deltae, verbose);
    maybe_populate(name_gamma,  desc_gamma,  gamma,  that.gamma,  verbose);
    maybe_populate(name_Mae,    desc_Mae,    Mae,    that.Mae,    verbose);
    maybe_populate(name_pexi,   desc_pexi,   pexi,   that.pexi,   verbose);
}

void
definition_radialflow::override(
        const specification_radialflow& that,
        const bool verbose)
{
    maybe_override(name_deltae, desc_deltae, deltae, that.deltae, verbose);
    maybe_override(name_gamma,  desc_gamma,  gamma,  that.gamma,  verbose);
    maybe_override(name_Mae,    desc_Mae,    Mae,    that.Mae,    verbose);
    maybe_override(name_pexi,   desc_pexi,   pexi,   that.pexi,   verbose);
}

void
definition_radialflow::save(
    const esio_handle h) const
{
    // Save nothing if there's nothing interesting to save
    if (this->trivial()) {
        return;
    }

    DEBUG0("Storing definition_radialflow parameters");

    // Write out the "container" holding all other settings
    const int one = 1;
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(h, location, &one, 0,
                    "Is a radial flow definition in use?");

    // Write out information as attributes within the container
    esio_attribute_write(h, location, attr_deltae, &deltae);
    esio_attribute_write(h, location, attr_gamma,  &gamma);
    esio_attribute_write(h, location, attr_Mae,    &Mae);
    esio_attribute_write(h, location, attr_pexi,   &pexi);
}

void
definition_radialflow::load(
    const esio_handle h,
    const bool verbose)
{
    definition_radialflow t;

    // Only proceed if a definition is active in the restart
    int in_use = 0;
    esio_line_establish(h, 1, 0, 1); // All ranks load any data
    if (ESIO_NOTFOUND != esio_line_size(h, location, NULL)) {
        esio_line_read(h, location, &in_use, 0);
    }
    if (!in_use) {
        return;
    }

    DEBUG0("Loading definition_radialflow parameters");

    // Read in information as attributes within the container
    esio_attribute_read(h, location, attr_deltae, &t.deltae);
    esio_attribute_read(h, location, attr_gamma,  &t.gamma);
    esio_attribute_read(h, location, attr_Mae,    &t.Mae);
    esio_attribute_read(h, location, attr_pexi,   &t.pexi);

    return this->populate(t, verbose);
}

} // end namespace support

} // end namespace suzerain
