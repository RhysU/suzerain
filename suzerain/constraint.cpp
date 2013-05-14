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
 * @copydoc constraint.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/constraint.hpp>

#include <suzerain/bspline.hpp>

namespace suzerain {

constraint::constraint()
    : target(std::numeric_limits<real_t>::quiet_NaN())
{
    // NOP
}

constraint::constraint(const value_lower&, const real_t target, bspline &b)
    : target(target)
    , coeff(b.n())
{
    coeff.setZero();
    coeff.head<1>()[0] = 1;
}

constraint::constraint(const value_upper&, const real_t target, bspline &b)
    : target(target)
    , coeff(b.n())
{
    coeff.setZero();
    coeff.tail<1>()[0] = 1;
}

constraint::constraint(const value_bulk&, const real_t target, bspline &b)
    : target(target)
    , coeff(b.n())
{
    const real_t L = b.breakpoint(b.nbreak()-1) - b.breakpoint(0);
    b.integration_coefficients(0, coeff.data());
    coeff /= L;
}

} // namespace suzerain
