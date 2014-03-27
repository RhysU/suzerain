//--------------------------------------------------------------------------
//
// Copyright (C) 2014 Rhys Ulerich
// Copyright (C) 2014 The PECOS Development Team
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
 * @copydoc specification_helm.hpp
 */

#include <suzerain/specification_helm.hpp>

#include <suzerain/common.hpp>

namespace suzerain {

specification_helm::specification_helm(
        const double kp,
        const double Td,
        const double Tf,
        const double Ti,
        const double Tt)
    : t(std::numeric_limits<double>::quiet_NaN())
    , v(std::numeric_limits<double>::quiet_NaN())
    , r(std::numeric_limits<double>::quiet_NaN())
{
    this->reset();
    h.kp = kp;
    h.Td = Td;
    h.Tf = Tf;
    h.Ti = Ti;
    h.Tt = Tt;
}

bool
specification_helm::enabled()
{
#pragma warning(push,disable:1572)
    return this->h.kp == 0;
#pragma warning(pop,disable:1572)
}

specification_helm *
specification_helm::reset()
{
    helm_reset(&this->h);
    return this;
}

specification_helm *
specification_helm::approach(
        const double t,
        const double v,
        const double r)
{
    this->t = t;
    this->v = v;
    this->r = r;
    helm_approach(&this->h);
    return this;
}

double
specification_helm::steady(
        const double t,
        const double u,
        const double y) const
{
    assert(!(boost::math::isnan)(this->t));
    assert(!(boost::math::isnan)(this->v));
    assert(!(boost::math::isnan)(this->r));

    // Invoke control logic
    const double dt = t - this->t;
    const double dv = helm_steady(&this->h, dt, r, u, v, y);

    // Update absolute state
    this->t  = t;
    this->v += dv;

    return dv;
}

} // namespace suzerain
