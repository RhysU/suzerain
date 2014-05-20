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

/** @file
 * @copydoc operator_common_block.hpp
 */

#include "operator_common_block.hpp"

#include <suzerain/common.hpp>

namespace suzerain {

namespace perfect {

operator_common_block::operator_common_block()
    : linearization(static_cast<linearize::type>(0)) // Invalid is "false"
{
}

void
operator_common_block::set_zero(int Ny)
{
    ref.set_zero(Ny);
    sub.set_zero(Ny);
    imp.set_zero(Ny);
}

} // namespace perfect

} // namespace suzerain
