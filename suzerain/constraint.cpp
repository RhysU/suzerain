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

namespace constraint {

base::~base()
{
    // NOP
}

disabled::disabled(bspline &b)
    : base()
{
    SUZERAIN_UNUSED(b);
}

real_t
disabled::target() const
{
    return std::numeric_limits<real_t>::quiet_NaN();
}

bool
disabled::enabled() const
{
    return false;
}

lower::lower(bspline &b)
    : base()
{
    coeff.setZero(b.n());
    coeff.head<1>()[0] = 1;
}

upper::upper(bspline &b)
    : base()
{
    coeff.setZero(b.n());
    coeff.tail<1>()[0] = 1;
}

bulk::bulk(bspline &b)
    : base()
{
    coeff.resize(b.n());
    b.integration_coefficients(0, coeff.data());
    real_t L = b.breakpoint(b.nbreak()-1) - b.breakpoint(0);
    coeff /= L;
}

constant::constant(const real_t target)
    : base()
    , t(target)
{
    // NOP
}

real_t
constant::target() const
{
    return t;
}

bool
constant::enabled() const
{
    return !(boost::math::isnan)(t);
}

reference::reference(const real_t& target)
    : base()
    , t(&target)
{
    // NOP
}

real_t
reference::target() const
{
    return *t;
}

bool
reference::enabled() const
{
    return !(boost::math::isnan)(*t);
}

constant_lower::constant_lower(const real_t target, bspline &b)
    : constant(target)
    , lower(b)
{
    // NOP
}

constant_upper::constant_upper(const real_t target, bspline &b)
    : constant(target)
    , upper(b)
{
    // NOP
}

constant_bulk::constant_bulk(const real_t target, bspline &b)
    : constant(target)
    , bulk(b)
{
    // NOP
}

reference_lower::reference_lower(const real_t& target, bspline &b)
    : reference(target)
    , lower(b)
{
    // NOP
}

reference_upper::reference_upper(const real_t& target, bspline &b)
    : reference(target)
    , upper(b)
{
    // NOP
}

reference_bulk::reference_bulk(const real_t& target, bspline &b)
    : reference(target)
    , bulk(b)
{
    // NOP
}

} // namespace constraint

} // namespace suzerain
