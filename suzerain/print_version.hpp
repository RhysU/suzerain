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
 * print_version.hpp: dump a program version message to stdout
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_PRINT_VERSION_HPP
#define __SUZERAIN_PRINT_VERSION_HPP

#include <iosfwd>
#include <string>

namespace suzerain {

/**
 * Print version information on the given stream.  Includes the application
 * name, the application version (if supplied), the suzerain library version,
 * and the compiler used.
 *
 * @param out Stream on which to write the information.
 * @param application_name Application name to be written.
 * @param application_version Application version string to output.
 */
void print_version(std::ostream &out,
                   const std::string &application_name,
                   const std::string &application_version = "");

} // namespace suzerain

#endif // __SUZERAIN_PRINT_VERSION_HPP
