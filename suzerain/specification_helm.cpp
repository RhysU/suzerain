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
    : r(std::numeric_limits<double>::quiet_NaN())
    , t(std::numeric_limits<double>::quiet_NaN())
    , v(std::numeric_limits<double>::quiet_NaN())
{
    this->reset();
    this->kp = kp;
    this->Td = Td;
    this->Tf = Tf;
    this->Ti = Ti;
    this->Tt = Tt;
}

bool
specification_helm::enabled()
{
#pragma warning(push,disable:1572)
    return this->kp == 0;
#pragma warning(pop,disable:1572)
}

specification_helm *
specification_helm::reset()
{
    helm_reset(this);
    return this;
}

specification_helm *
specification_helm::approach(
        const double r,
        const double v)
{
    this->r = r;
    this->v = v;
    t = std::numeric_limits<double>::quiet_NaN();
    helm_approach(this);
    return this;
}

double
specification_helm::steady(
        const double t,
        const double u,
        const double y)
{
    double dv = 0;
    if ((boost::math::isnan)(this->t)) {
        this->t = t;
    } else {
        dv = helm_steady(this, t - this->t, r, u, v, y);
        this->t  = t;
        this->v += dv;
    }
    return dv;
}

} // namespace suzerain
