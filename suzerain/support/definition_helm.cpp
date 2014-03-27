//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
// Copyright (C) 2014 The PECOS Development Team
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
 * @copydoc definition_noise.hpp
 */

#include <suzerain/support/definition_helm.hpp>

#include <boost/version.hpp>
#include <esio/error.h>
#include <esio/esio.h>

#include <suzerain/exprparse.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

/** Helper used to parse string-based options */
template<typename T>
static void parse_option(const std::string& s,
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

definition_helm::definition_helm()
    : specification_helm(std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::quiet_NaN())
{
    // Fix setpoint within superclass storage
    this->r = std::numeric_limits<double>::quiet_NaN();
}

definition_helm::definition_helm(
        const double kp,
        const double Td,
        const double Tf,
        const double Ti,
        const double Tt,
        const double r)
    : specification_helm(kp, Td, Tf, Ti, Tt)
{
    // Preserve incoming setpoint within superclass storage
    this->r = r;
}

// Strings used in options_description and populate/override/save/load.
static const char location[] = "helm";
static const char desc_location []
    = "PID controller settings used to drive boundary layer thickness"
      " (measured using 99% of freestream velocity)"
      " by dynamically adjusting the Largo slow growth parameter."
      " Poorly chosen parameters will destabilize a simulation.";

static const char name_r [] = "helm_r";
static const char name_kp[] = "helm_kp";
static const char name_Td[] = "helm_Td";
static const char name_Tf[] = "helm_Tf";
static const char name_Ti[] = "helm_Ti";
static const char name_Tt[] = "helm_Tt";

static const char * const attr_r  = name_r  + sizeof(location);
static const char * const attr_kp = name_kp + sizeof(location);
static const char * const attr_Td = name_Td + sizeof(location);
static const char * const attr_Tf = name_Tf + sizeof(location);
static const char * const attr_Ti = name_Ti + sizeof(location);
static const char * const attr_Tt = name_Tt + sizeof(location);

static const char desc_r [] = "Reference value, often called the setpoint";
static const char desc_kp[] = "Proportional gain modifying P, I, and D terms.";
static const char desc_Td[] = "Time scale governing derivative action."
                              " Set to zero to disable derivative control.";
static const char desc_Tf[] = "Time scale filtering process observable for D."
                              " Set to infinity to disable filtering.";
static const char desc_Ti[] = "Time scale governing integral action."
                              " Set to infinity to disable integral control.";
static const char desc_Tt[] = "Time scale governing automatic reset."
                              " Set to infinity to disable automatic reset.";

boost::program_options::options_description
definition_helm::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::bool_switch;
    using boost::program_options::value;
    using std::string;
    using validation::ensure_nonnegative;
    using validation::ensure_positive;

    boost::program_options::options_description retval(desc_location);

    retval.add_options()
        (name_kp, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("GAIN")
#endif
         ->notifier(bind(&parse_option<double>, _1, &h.kp,
                         &ensure_nonnegative<double>, name_kp))
         ->default_value(lexical_cast<string>(h.kp)),
         desc_kp)
        (name_r, value<string>(NULL)
         ->notifier(bind(&parse_option<double>, _1, &this->r,
                         &ensure_nonnegative<double>, name_r))
         ->default_value(lexical_cast<string>(this->r)),
         desc_r)
        (name_Td, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option<double>, _1, &h.Td,
                         &ensure_nonnegative<double>, name_Td))
         ->default_value(lexical_cast<string>(h.Td)),
         desc_Td)
        (name_Tf, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option<double>, _1, &h.Tf,
                         &ensure_nonnegative<double>, name_Tf))
         ->default_value(lexical_cast<string>(h.Tf)),
         desc_Tf)
        (name_Ti, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option<double>, _1, &h.Ti,
                         &ensure_nonnegative<double>, name_Ti))
         ->default_value(lexical_cast<string>(h.Ti)),
         desc_Ti)
        (name_Tt, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option<double>, _1, &h.Tt,
                         &ensure_nonnegative<double>, name_Tt))
         ->default_value(lexical_cast<string>(h.Tt)),
         desc_Tt)
        ;

    return retval;
}

void
definition_helm::populate(
        const definition_helm& that,
        const bool verbose)
{
#define CALL_MAYBE_POPULATE(mem, loc)                                  \
    maybe_populate(name_ ## mem, desc_ ## mem, loc, that.loc, verbose)
    CALL_MAYBE_POPULATE(kp, h.kp);
    CALL_MAYBE_POPULATE(r , r   );
    CALL_MAYBE_POPULATE(Td, h.Td);
    CALL_MAYBE_POPULATE(Tf, h.Tf);
    CALL_MAYBE_POPULATE(Ti, h.Ti);
    CALL_MAYBE_POPULATE(Tt, h.Tt);
#undef CALL_MAYBE_POPULATE
}

void
definition_helm::override(
        const definition_helm& that,
        const bool verbose)
{
#define CALL_MAYBE_OVERRIDE(mem, loc)                                  \
    maybe_override(name_ ## mem, desc_ ## mem, loc, that.loc, verbose)
    CALL_MAYBE_OVERRIDE(kp, h.kp);
    CALL_MAYBE_OVERRIDE(r , r   );
    CALL_MAYBE_OVERRIDE(Td, h.Td);
    CALL_MAYBE_OVERRIDE(Tf, h.Tf);
    CALL_MAYBE_OVERRIDE(Ti, h.Ti);
    CALL_MAYBE_OVERRIDE(Tt, h.Tt);
#undef CALL_MAYBE_OVERRIDE
}

void
definition_helm::save(
        const esio_handle h) const
{
    DEBUG0("Storing definition_helm parameters");

    // Only root writes the containing location with
    // a 0 or 1 indicating whether or not helm was enabled.
    int procid;
    esio_handle_comm_rank(h, &procid);
    esio_line_establish(h, 1, 0, (procid == 0 ? 1 : 0));
    const int enabled_int = enabled;
    esio_line_write(h, location, &enabled_int, 0, desc_location);

    // Everyone writes the metadata
    esio_attribute_write(h, location, attr_kp, &this->h.kp);
    esio_attribute_write(h, location, attr_r,  &this->r   );
    esio_attribute_write(h, location, attr_Td, &this->h.Td);
    esio_attribute_write(h, location, attr_Tf, &this->h.Tf);
    esio_attribute_write(h, location, attr_Ti, &this->h.Ti);
    esio_attribute_write(h, location, attr_Tt, &this->h.Tt);
}

void
definition_helm::load(
    const esio_handle h,
    const bool verbose)
{
    definition_helm t;

    // When possible, load settings from restart.

    // All ranks load any data with short circuit when not present
    esio_line_establish(h, 1, 0, 1);
    if (ESIO_NOTFOUND == esio_line_size(h, location, NULL)) {
        return;
    }

    DEBUG0("Loading definition_helm parameters");
    esio_attribute_read(h, location, attr_kp, &t.h.kp);
    esio_attribute_read(h, location, attr_r,  &t.r   );
    esio_attribute_read(h, location, attr_Td, &t.h.Td);
    esio_attribute_read(h, location, attr_Tf, &t.h.Tf);
    esio_attribute_read(h, location, attr_Ti, &t.h.Ti);
    esio_attribute_read(h, location, attr_Tt, &t.h.Tt);

    // Prefer incoming to temporary
    this->populate(t, verbose);
}

} // end namespace support

} // end namespace suzerain
