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
 * version.hpp: generate a nice program version message
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_VERSION_HPP
#define __SUZERAIN_VERSION_HPP

#include <iosfwd>
#include <string>

namespace suzerain {

/**
 * Generate a program version message.  Includes the application name (if
 * supplied), the application version (if supplied), the suzerain library
 * version, and the compiler used.
 *
 * @param application_name    Application name to be written.
 * @param application_version Application version to additionally output.
 *
 * @return A human-readable version string.
 */
std::string version(const std::string &application_name    = "",
                    const std::string &application_version = "");

} // namespace suzerain

#endif // __SUZERAIN_VERSION_HPP
