//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

/** @file
 * @copydoc version.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <string>
#include <sstream>

#include <suzerain/suzerain-revision.h>
#include <suzerain/version.hpp>

namespace suzerain {

std::string version(const std::string &application_name,
                    const std::string &application_version)
{
    std::ostringstream oss;

    oss << application_name;
    if (!application_version.empty()) {
        oss << ' ' << application_version;
    }
    if (   !application_name.empty()
        || (application_name.empty() && !application_version.empty())) {
        oss << " written using ";
    }
    oss << PACKAGE_NAME << ' ' << SUZERAIN_REVISION_STR
        << " ("
#if defined(__INTEL_COMPILER)
        << "Intel "
        << __INTEL_COMPILER << " " << __INTEL_COMPILER_BUILD_DATE
#elif defined(__GNUC__)
        << "GNU "
        << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__
#else
        << "unknown compiler"
#endif
        << ')';

    return oss.str();
}

} // namespace suzerain
