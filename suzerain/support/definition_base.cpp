//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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
 * @copydoc definition_base.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/definition_base.hpp>

#include <suzerain/support/logging.hpp>

namespace suzerain {

namespace support {

// For maybe_XXX_impl to indicates a \c real_t is a default value
static bool default_value(const real_t& v)
{ return (boost::math::isnan)(v); }

// For maybe_XXX_impl to indicates an \c int is a default value
static bool default_value(const int& v)
{ return v == 0; }

// Compare and contrast maybe_override just below
// One translation-local template instantiated multiple times below
template<typename T>
static bool
maybe_populate_impl(
        const char* name,
        const char* description,
              T&    destination,
        const T&    source,
        const bool  verbose)
{
    if (default_value(destination)) {
        if (verbose) {
            // Populating a default isn't interesting, hence DEBUG0
            if (description) {
                DEBUG0("Populating " << name
                      << " (" << description << ") to be " << source);
            } else {
                DEBUG0("Populating " << name << " to be " << source);
            }
        }
        destination = source;
        return true;
    }

    if (verbose) {
        if (description) {
            // Retaining an existing setting is interesting, hence INFO0
            INFO0("Clutching onto " << name
                   << " (" << description << ") of " << destination);
        } else {
            INFO0("Clutching onto " << name << " of " << destination);
        }
    }
    return false;
}

// Compare and contrast maybe_populate just above
// One translation-local template instantiated multiple times below
template<typename T>
static bool
maybe_override_impl(
        const char* name,
        const char* description,
              T&    destination,
        const T&    source,
        const bool  verbose)
{
    if (!default_value(source)) {
#pragma warning(push,disable:1572)
        if (    verbose
             && source != destination
             && !default_value(destination)) {
#pragma warning(pop)
            // Overriding is value is interesting, hence INFO0
            if (description) {
                INFO0("Overriding " << name
                      << " (" << description << ") to be " << source);
            } else {
                INFO0("Overriding " << name << " to be " << source);
            }
        }
        destination = source;
        return true;
    }

    if (verbose) {
        // Retaining an existing setting isn't interesting, hence DEBUG0
        if (description) {
            DEBUG0("Retaining " << name
                   << " (" << description << ") as " << destination);
        } else {
            DEBUG0("Retaining " << name << " as " << destination);
        }
    }
    return false;
}

bool
definition_base::maybe_populate(
        const char*   name,
        const char*   description,
              real_t& destination,
        const real_t& source,
        const bool    verbose)
{
    return maybe_populate_impl(
            name, description, destination, source, verbose);
}

bool
definition_base::maybe_override(
        const char*   name,
        const char*   description,
              real_t& destination,
        const real_t& source,
        const bool    verbose)
{
    return maybe_override_impl(
            name, description, destination, source, verbose);
}

bool
definition_base::maybe_populate(
        const char* name,
        const char* description,
              int&  destination,
        const int&  source,
        const bool  verbose)
{
    return maybe_populate_impl(
            name, description, destination, source, verbose);
}

bool
definition_base::maybe_override(
        const char* name,
        const char* description,
              int&  destination,
        const int&  source,
        const bool  verbose)
{
    return maybe_override_impl(
            name, description, destination, source, verbose);
}

} // namespace support

} // namespace suzerain
