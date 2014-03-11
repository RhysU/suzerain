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

#ifndef SUZERAIN_SUPPORT_OVERRIDABLE_HPP
#define SUZERAIN_SUPPORT_OVERRIDABLE_HPP

/** @file
 * Provides \ref overridable and \ref maybe_override.
 */

#include <suzerain/common.hpp>

#include <suzerain/support/logging.hpp>

namespace suzerain {

namespace support {

/**
 * Abstract interface indicating data may be overridden from another instance.
 */
template <typename Derived>
class overridable
{
public:

    /**
     * Override members in \c this with non-NaN values from \c that.
     *
     * Descendants should override this method adding any desired functionality
     * either before or after invoking the superclass version.
     *
     * @param that    Instance from which information is taken.
     * @param verbose Should logging be emitted when an override occurs?
     */
    virtual void override(
            const Derived& that,
            const bool verbose = false) = 0;

    /** Virtual destructor to permit deleting through subclass. */
    virtual ~overridable() {}

};

/**@name Helpers for \ref overridable implementations */
/**@{*/

/**
 * If \c source is non-NaN, override \c destination with its value.
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
bool maybe_override(const char*   name,
                    const char*   description,
                          real_t& destination,
                    const real_t& source,
                    const bool    verbose);

/**
 * If \c source is non-zero, override \c destination with its value.
 * @copydetails maybe_override(const char*,const char*,real_t&,const real_t&,const bool)
 */
bool maybe_override(const char* name,
                    const char* description,
                          int&  destination,
                    const int&  source,
                    const bool  verbose);

/**
 * If \c source is not false, override \c destination with its value.
 * @copydetails maybe_override(const char*,const char*,real_t&,const real_t&,const bool)
 */
bool maybe_override(const char* name,
                    const char* description,
                          bool& destination,
                    const bool& source,
                    const bool  verbose);

/**@}*/

namespace internal {

// Compare and contrast maybe_populate_impl in populatable.hpp
template<typename T>
static bool
maybe_override_impl(
        const char* name,
        const char* description,
              T&    destination,
        const T&    source,
        const bool  verbose,
        bool (* const default_value)(const T& v))
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

} // namespace internal

} // namespace support

} // namespace suzerain

#endif // SUZERAIN_SUPPORT_OVERRIDABLE_HPP
