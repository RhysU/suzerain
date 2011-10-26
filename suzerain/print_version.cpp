/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2011 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * print_version.cpp: dump a program version message to stdout
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <iostream>
#include <suzerain/suzerain-revision.h>
#include "print_version.hpp"

namespace suzerain {

void print_version(std::ostream &out,
                   const std::string &application_name,
                   const std::string &application_version)
{
    out << application_name;
    if (!application_version.empty()) {
        out << ' ' << application_version;
    }
    out << " written using "
        << PACKAGE_NAME << ' ' << SUZERAIN_REVISION_STR
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
        << ')'
        << std::endl;
}

} // namespace suzerain
