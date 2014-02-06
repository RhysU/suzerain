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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/definition_arsel.hpp>

#include <suzerain/exprparse.hpp>
#include <suzerain/validation.hpp>

namespace suzerain {

namespace support {

/** Helper used to parse size_t-based options */
static void parse_size_t(const std::string& s,
                         std::size_t* value,
                         const char* name)
{
    using std::floor;
#pragma warning(push,disable:2259)
    const real_t t = floor(exprparse<real_t>(s, name) + real_t(1)/2);
#pragma warning(pop)
    validation::ensure_nonnegative(t, name);
    *value = t;
}

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

definition_arsel::definition_arsel()
    : specification_arsel()
{
}

boost::program_options::options_description
definition_arsel::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::bool_switch;
    using boost::program_options::value;
    using std::string;
    using validation::ensure_positive;

    boost::program_options::options_description retval(
            "Automatic autocorrelation analysis by AR(p) models");

    retval.add_options()
        ("ar_minorder", value<string>(NULL)
         ->value_name("MIN")
         ->notifier(bind(&parse_size_t, _1, &minorder, "ar_minorder"))
         ->default_value(lexical_cast<string>(minorder)),
         "Consider only models of at least order AR(p=MIN)")
        ("ar_maxorder", value<string>(NULL)
         ->value_name("MAX")
         ->notifier(bind(&parse_size_t, _1, &maxorder, "ar_maxorder"))
         ->default_value(lexical_cast<string>(maxorder)),
         "Consider only models of at most order AR(p=MAX)")
        ("ar_criterion",
         value<std::string>()
         ->notifier(bind(&definition_arsel::criterion, this, _1))
         ->default_value(criterion()),
         "Employ the specified model selection criterion")
        ("ar_absolute_rho",
         bool_switch(&absrho)->default_value(absrho),
         "Integrate absolute autocorrelation when determining T0")
        ("ar_wlenT0", value<string>(NULL)
         ->value_name("WLEN")
         ->notifier(bind(&parse_option<real_t>, _1, &wlenT0,
                         &ensure_positive<real_t>, "ar_wlenT0"))
         ->default_value(lexical_cast<string>(wlenT0)),
         "Integrate for T0 until WLEN times the input length")
        ;

    return retval;
}

} // end namespace support

} // end namespace suzerain
