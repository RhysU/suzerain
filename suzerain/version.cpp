//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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
 * @copydoc version.hpp
 */

#include <suzerain/version.hpp>

#include <sstream>

#include <suzerain/suzerain-config.h>

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
    oss << SUZERAIN_PACKAGE_NAME << ' ' << REVISIONSTR
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
