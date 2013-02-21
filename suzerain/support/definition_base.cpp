//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

// Compare and contrast maybe_override just below
bool
definition_base::maybe_populate(
        const char*   name,
        const char*   description,
              real_t& destination,
        const real_t& source,
        const bool    verbose)
{
    if ((boost::math::isnan)(destination)) {
        if (verbose) {
            if (description) {
                INFO0("Populating " << name
                      << " (" << description << ") to be " << source);
            } else {
                INFO0("Populating " << name << " to be " << source);
            }
        }
        destination = source;
        return true;
    }

    if (verbose) {
        if (description) {
            DEBUG0("Retaining " << name
                   << " (" << description << ") as " << destination);
        } else {
            DEBUG0("Retaining " << name << " as " << destination);
        }
    }
    return false;
}

// Compare and contrast maybe_populate just above
bool
definition_base::maybe_override(
        const char*   name,
        const char*   description,
              real_t& destination,
        const real_t& source,
        const bool    verbose)
{
    if (!(boost::math::isnan)(source)) {
#pragma warning(push,disable:1572)
        if (    verbose
             && source != destination
             && !(boost::math::isnan)(destination)) {
#pragma warning(pop)
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
        if (description) {
            DEBUG0("Retaining " << name
                   << " (" << description << ") as " << destination);
        } else {
            DEBUG0("Retaining " << name << " as " << destination);
        }
    }
    return false;
}

} // namespace support

} // namespace suzerain
