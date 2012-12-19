//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc statistics_definition.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/statistics_definition.hpp>

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

statistics_definition::statistics_definition(
        const std::string& destination,
        const std::size_t  retain,
        const real_t       dt,
        const std::size_t  nt)
    : destination(destination)
    , retain(retain)
    , dt(dt)
    , nt(nt)
{
}

boost::program_options::options_description
statistics_definition::options_description()
{
    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::bool_switch;
    using boost::program_options::options_description;
    using boost::program_options::value;
    using std::string;
    using validation::ensure_nonnegative;

    options_description retval("Statistics sampling parameters");

    retval.add_options()
    ("statistics_destination", value(&destination)
     ->default_value(destination),
     "Archiving destination to use when committing statistics files.  "
     "One or more #'s must be present and will be replaced by a sequence number.  "
     "Any trailing \"XXXXXX\" will be used to generate a unique template.")
    ("statistics_retain", value<string>(NULL)
     ->notifier(bind(&parse_size_t, _1, &retain, "statistics_retain"))
     ->default_value(lexical_cast<string>(retain)),
     "Maximum number of committed statistics files to retain")
    ("statistics_dt", value<string>(NULL)
     ->notifier(bind(&parse_option<real_t>, _1, &dt,
                     &ensure_nonnegative<real_t>, "statistics_dt"))
     ->default_value(lexical_cast<string>(dt)),
     "Maximum amount of simulation time between sampling statistics")
    ("statistics_nt", value<string>(NULL)
     ->notifier(bind(&parse_size_t, _1, &nt, "statistics_nt"))
     ->default_value(lexical_cast<string>(nt)),
     "Maximum number of time steps between sampling statistics")
    ;

    return retval;
}

} // namespace support

} // namespace suzerain
