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
            double Ma0,
            double gam0,
            double rho1,
            double u1,
            double R1)
    : radialflow_specification(Ma0,
                                  gam0,
                                  rho1,
                                  u1,
                                  R1)
{
}

// Strings used in options_description and populate/override/save/load.
static const char location [] = "radialflow";
static const char name_Ma0 [] = "radialflow_Ma0";
static const char name_gam0[] = "radialflow_gam0";
static const char name_rho1[] = "radialflow_rho1";
static const char name_u1  [] = "radialflow_u1";
static const char name_R1  [] = "radialflow_R1";
static const char * const attr_Ma0  = name_Ma0  + sizeof(location);
static const char * const attr_gam0 = name_gam0 + sizeof(location);
static const char * const attr_rho1 = name_rho1 + sizeof(location);
static const char * const attr_u1   = name_u1   + sizeof(location);
static const char * const attr_R1   = name_R1   + sizeof(location);

// Descriptions used in options_description and populate/override/save/load.
static const char desc_Ma0 [] = "Reference Mach number, Ma_0";
static const char desc_gam0[] = "Reference specific heat ratio, gamma_0";
static const char desc_rho1[] = "Inner density, rho(R_1)";
static const char desc_u1  [] = "Inner radial velocity, u(R_1)";
static const char desc_R1  [] = "Inner radius of interest, R_1";

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
            "Radial nozzle problem parameters "
            " for homogenized boundary layers with baseflow");

    auto_ptr<typed_value<string> > p;

    // Ma0
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &Ma0,
                     &ensure_positive<double>, name_Ma0));
    if (!(boost::math::isnan)(Ma0)) {
        p->default_value(lexical_cast<string>(Ma0));
    }
    retval.add_options()(name_Ma0, p.release(), desc_Ma0);

    // gam0
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &gam0,
                     &ensure_positive<double>, name_gam0));
    if (!(boost::math::isnan)(gam0)) {
        p->default_value(lexical_cast<string>(gam0));
    }
    retval.add_options()(name_gam0, p.release(), desc_gam0);

    // rho1
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &rho1,
                     &ensure_positive<double>, name_rho1));
    if (!(boost::math::isnan)(rho1)) {
        p->default_value(lexical_cast<string>(rho1));
    }
    retval.add_options()(name_rho1, p.release(), desc_rho1);

    // u1
    p.reset(value<string>());
    p->notifier(bind(&parse_option<double>, _1, &u1, name_u1));
    if (!(boost::math::isnan)(u1)) {
        p->default_value(lexical_cast<string>(u1));
    }
    retval.add_options()(name_u1, p.release(), desc_u1);

    // R1
    p.reset(value<string>());
    p->notifier(bind(&validate_option<double>, _1, &R1,
                     &ensure_positive<double>, name_R1));
    if (!(boost::math::isnan)(R1)) {
        p->default_value(lexical_cast<string>(R1));
    }
    retval.add_options()(name_R1, p.release(), desc_R1);

    return retval;
}

void
radialflow_definition::populate(
        const radialflow_specification& that,
        const bool verbose)
{
    maybe_populate(name_Ma0,  desc_Ma0,  Ma0,  that.Ma0,  verbose);
    maybe_populate(name_gam0, desc_gam0, gam0, that.gam0, verbose);
    maybe_populate(name_rho1, desc_rho1, rho1, that.rho1, verbose);
    maybe_populate(name_u1,   desc_u1,   u1,   that.u1,   verbose);
    maybe_populate(name_R1,   desc_R1,   R1,   that.R1,   verbose);
}

void
radialflow_definition::override(
        const radialflow_specification& that,
        const bool verbose)
{
    maybe_override(name_Ma0,  desc_Ma0,  Ma0,  that.Ma0,  verbose);
    maybe_override(name_gam0, desc_gam0, gam0, that.gam0, verbose);
    maybe_override(name_rho1, desc_rho1, rho1, that.rho1, verbose);
    maybe_override(name_u1,   desc_u1,   u1,   that.u1,   verbose);
    maybe_override(name_R1,   desc_R1,   R1,   that.R1,   verbose);
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
                    "Is a radial nozzle definition in use?");

    // Write out information as attributes within the container
    esio_attribute_write(h, location, attr_Ma0,  &Ma0);
    esio_attribute_write(h, location, attr_gam0, &gam0);
    esio_attribute_write(h, location, attr_rho1, &rho1);
    esio_attribute_write(h, location, attr_u1,   &u1);
    esio_attribute_write(h, location, attr_R1,   &R1);
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
    esio_attribute_read(h, location, attr_Ma0,  &t.Ma0);
    esio_attribute_read(h, location, attr_gam0, &t.gam0);
    esio_attribute_read(h, location, attr_rho1, &t.rho1);
    esio_attribute_read(h, location, attr_u1,   &t.u1);
    esio_attribute_read(h, location, attr_R1,   &t.R1);

    return this->populate(t, verbose);
}

} // end namespace support

} // end namespace suzerain
