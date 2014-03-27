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

#include <suzerain/exprparse.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

/** Helper for parsing and validating options. */
static void parse_option(const std::string& s,
                         double* value,
                         const char* name)
{

#pragma warning(push,disable:2259)
    const double t = exprparse<double>(s, name);
#pragma warning(pop)
    validation::ensure_nonnegative(t, name);
    *value = t;
}

definition_helm::definition_helm(
        const double kp,
        const double Td,
        const double Tf,
        const double Ti,
        const double Tt)
    : specification_helm(kp, Td, Tf, Ti, Tt)
{
}

static const char name_kp   [] = "helm_kp";
static const char name_Td   [] = "helm_Td";
static const char name_Tf   [] = "helm_Tf";
static const char name_Ti   [] = "helm_Ti";
static const char name_Tt   [] = "helm_Tt";

static const char desc_kp[] = "Proportional gain modifying P, I, and D terms.";
static const char desc_Td[] = "Time scale governing derivative action."
                              " Set to zero to disable derivative control.";
static const char desc_Tf[] = "Time scale filtering process observable for D."
                              " Set to infinity to disable filtering.";
static const char desc_Ti[] = "Time scale governing integral action."
                              " Set to infinity to disable integral control.";
static const char desc_Tt[] = "Time scale governing automatic reset."
                              " Set to infinity to disable automatic reset.";

static const char desc_location []
    = "PID controller settings used to drive boundary layer thickness"
      " (measured using 99% of freestream velocity)"
      " by dynamically adjusting the Largo slow growth parameter."
      " Poorly chosen parameters will destabilize a simulation.";

boost::program_options::options_description
definition_helm::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::bool_switch;
    using boost::program_options::value;
    using std::string;
    using validation::ensure_positive;

    boost::program_options::options_description retval(desc_location);

    retval.add_options()
        (name_kp, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("GAIN")
#endif
         ->notifier(bind(&parse_option, _1, &h.kp, name_kp))
         ->default_value(lexical_cast<string>(h.kp)),
         desc_kp)
        (name_Td, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option, _1, &h.Td, name_Td))
         ->default_value(lexical_cast<string>(h.Td)),
         desc_Td)
        (name_Tf, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option, _1, &h.Tf, name_Tf))
         ->default_value(lexical_cast<string>(h.Tf)),
         desc_Tf)
        (name_Ti, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option, _1, &h.Ti, name_Ti))
         ->default_value(lexical_cast<string>(h.Ti)),
         desc_Ti)
        (name_Tt, value<string>(NULL)
#if BOOST_VERSION >= 105000
         ->value_name("TIMESCALE")
#endif
         ->notifier(bind(&parse_option, _1, &h.Tt, name_Tt))
         ->default_value(lexical_cast<string>(h.Tt)),
         desc_Tt)
        ;

    return retval;
}

} // end namespace support

} // end namespace suzerain
