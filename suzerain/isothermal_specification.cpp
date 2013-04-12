//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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
 * @copydoc isothermal_specification.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/isothermal_specification.hpp>

namespace suzerain {

isothermal_specification::isothermal_specification()
    : lower_u(std::numeric_limits<real_t>::quiet_NaN())
    , lower_v(std::numeric_limits<real_t>::quiet_NaN())
    , lower_w(std::numeric_limits<real_t>::quiet_NaN())
    , lower_T(std::numeric_limits<real_t>::quiet_NaN())
    , lower_mass_fractions(0)
    , upper_u(std::numeric_limits<real_t>::quiet_NaN())
    , upper_v(std::numeric_limits<real_t>::quiet_NaN())
    , upper_w(std::numeric_limits<real_t>::quiet_NaN())
    , upper_T(std::numeric_limits<real_t>::quiet_NaN())
    , upper_mass_fractions(0)
{
}

isothermal_specification::~isothermal_specification()
{
}

} // namespace suzerain
