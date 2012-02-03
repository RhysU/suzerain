//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// version.hpp: generate a nice program version message
// $Id$

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
