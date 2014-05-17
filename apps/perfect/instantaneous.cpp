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
 * @copydoc instantaneous.hpp
 */
#include "instantaneous.hpp"

#include "references.hpp"

namespace suzerain {

namespace perfect {

instantaneous::instantaneous()
{
}

instantaneous::~instantaneous()
{
}

void
instantaneous::set_zero(int Ny)
{
    super::setZero(Ny, q::count);
}

instantaneous&
instantaneous::operator=(const references& that)
{
    // Ensure adequate storage for the assignment (noticing transposed sizes)
    // FIXME Usage of q::count brittle wrt subclasses-- NoChange preferable
#ifndef NDEBUG
    setConstant(that.cols(), q::count, std::numeric_limits<Scalar>::quiet_NaN());
#else
    resize     (that.cols(), q::count);
#endif

    // Assign the data quantity-by-quantity
    // Notice target is stride-one while source is not
    this->rho  () = that.rho  ();
    this->p    () = that.p    ();
    this->p2   () = that.p2   ();
    this->u    () = that.u    ();
    this->v    () = that.v    ();
    this->w    () = that.w    ();
    this->uu   () = that.uu   ();
    this->uv   () = that.uv   ();
    this->uw   () = that.uw   ();
    this->vv   () = that.vv   ();
    this->vw   () = that.vw   ();
    this->ww   () = that.ww   ();
    this->rhou () = that.rhou ();
    this->rhov () = that.rhov ();
    this->rhow () = that.rhow ();
    this->rhoE () = that.rhoE ();
    this->rhouu() = that.rhouu();
    this->rhovv() = that.rhovv();
    this->rhoww() = that.rhoww();
    this->rhoEE() = that.rhoEE();

    return *this;
}

} // namespace perfect

} // namespace suzerain
