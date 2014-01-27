//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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
 * @copydoc definition_restart.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/definition_restart.hpp>

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

definition_restart::definition_restart(
        const std::string& metadata,
        const std::string& uncommitted,
        const std::string& destination,
        const std::size_t  retain,
        const real_t       dt,
        const std::size_t  nt,
        const bool         final,
        const bool         physical)
    : metadata(metadata)
    , uncommitted(uncommitted)
    , destination(destination)
    , retain(retain)
    , dt(dt)
    , nt(nt)
    , final(final)
    , physical(physical)
{
}

boost::program_options::options_description
definition_restart::options_description()
{
    boost::program_options::options_description retval(
            "Restart-writing parameters");

    using boost::bind;
    using boost::lexical_cast;
    using boost::program_options::bool_switch;
    using boost::program_options::value;
    using std::string;
    using validation::ensure_nonnegative;

    retval.add_options()
    ("metadata", value(&metadata)
     ->default_value(metadata),
     "Path to use when saving common metadata for output files.  "
     "Any trailing \"XXXXXX\" will be used to generate a unique name.")
    ("uncommitted", value(&uncommitted)
     ->default_value(uncommitted),
     "Path to use when saving uncommitted output files.  "
     "Any trailing \"XXXXXX\" will be used to generate a unique name.")
    ("restart_destination", value(&destination)
     ->default_value(destination),
     "Archiving destination to use when committing restart files.  "
     "One or more #'s must be present and will be replaced by a sequence number.  "
     "Any trailing \"XXXXXX\" will be used to generate a unique template.")
    ("restart_retain", value<string>(NULL)
     ->notifier(bind(&parse_size_t, _1, &retain, "restart_retain"))
     ->default_value(lexical_cast<string>(retain)),
     "Maximum number of committed restart files to retain")
    ("restart_dt", value<string>(NULL)
     ->notifier(bind(&parse_option<real_t>, _1, &dt,
                     &ensure_nonnegative<real_t>, "restart_dt"))
     ->default_value(lexical_cast<string>(dt)),
     "Maximum amount of simulation time between restart files")
    ("restart_nt", value<string>(NULL)
     ->notifier(bind(&parse_size_t, _1, &nt, "restart_nt"))
     ->default_value(lexical_cast<string>(nt)),
     "Maximum number of time steps between restart files")
    ("restart_final", value<bool>(&final)
     ->default_value(final),
     "Should a final restart be written after advance successfully completes?")
    ("restart_physical", bool_switch(&physical),
     "Specify flag to save restart fields as primitive variables "
     "stored at collocation points in physical space")
    ;

    return retval;
}

} // end namespace support

} // end namespace suzerain
