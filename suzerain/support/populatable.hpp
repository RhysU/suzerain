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

#ifndef SUZERAIN_SUPPORT_POPULATABLE_HPP
#define SUZERAIN_SUPPORT_POPULATABLE_HPP

/** @file
 * Provides \ref populatable and \ref maybe_populate.
 */

#include <suzerain/common.hpp>

#include <suzerain/support/logging.hpp>

namespace suzerain {

namespace support {

/**
 * Abstract interface indicating details may be populated from another instance.
 */
template <class Derived>
class populatable
{
public:

    /**
     * Populate any NaN members in \c this with values from \c that.
     *
     * Descendants should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when a value is retained?
     */
    virtual void populate(
            const Derived& that,
            const bool verbose = false) = 0;

    /** Virtual destructor to permit deleting through subclass. */
    virtual ~populatable() {}

};

/**@name Helpers for \ref populatable implementations */
/**@{*/

/**
 * If \c destination is NaN, populate it with the value from \c source.
 * When \c verbose, log an informative message using \c name and \c
 * description.
 *
 * @param name        Name to use for any logging.
 * @param description A short description of \c name to use for logging.
 *                    If \c NULL, no description will be logged.
 * @param destination Location to possibly assign from \c source.
 * @param source      Source data to possibly use.
 * @param verbose     Should a human-readable message be logged?
 *
 * @return \c True whenever \c destination was assigned from \c source.
 */
bool maybe_populate(const char*   name,
                    const char*   description,
                          real_t& destination,
                    const real_t& source,
                    const bool    verbose);

/**
 * If \c destination is zero, populate it with the value from \c source.
 * @copydetails maybe_populate(const char*,const char*,real_t&,const real_t&,const bool)
 */
bool maybe_populate(const char* name,
                    const char* description,
                          int&  destination,
                    const int&  source,
                    const bool  verbose);

/**@}*/

namespace internal {

// Compare and contrast maybe_override_impl in overridable.hpp
template<typename T>
static bool
maybe_populate_impl(
        const char* name,
        const char* description,
              T&    destination,
        const T&    source,
        const bool  verbose,
        bool (* const default_value)(const T& v))
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

} // namespace internal

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_POPULATABLE_HPP
