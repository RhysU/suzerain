//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
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
 * @copydoc ndx.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/ndx.hpp>

namespace suzerain {

namespace ndx {

const array<const char *, rho + 1u> identifier = {{
    "rho_E", "rho_u", "rho_v", "rho_w", "rho"
}};

const array<const char *, rho + 1u> description = {{
    "total energy",
    "streamwise momentum",
    "wall-normal momentum",
    "spanwise momentum",
    "density"
}};

} // namespace ndx

} // namespace suzerain
