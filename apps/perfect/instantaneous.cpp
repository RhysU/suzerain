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

#include <suzerain/bspline.hpp>
#include <suzerain/samples.hpp>
#include <suzerain/summary.hpp>

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
instantaneous::set_zero(
        int Ny)
{
    super::setZero(Ny, q::count);
}

instantaneous&
instantaneous::operator=(
        const references& that)
{
    // Ensure adequate storage for the assignment (noticing transposed sizes)
    // FIXME Usage of q::count brittle wrt subclasses-- NoChange preferable
#ifndef NDEBUG
    setConstant(that.cols(), q::count, std::numeric_limits<Scalar>::quiet_NaN());
#else
    resize     (that.cols(), q::count);
#endif

    // Assign the data quantity-by-quantity
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
    this->rhouv() = that.rhouv();
    this->rhouw() = that.rhouw();
    this->rhovv() = that.rhovv();
    this->rhovw() = that.rhovw();
    this->rhoww() = that.rhoww();
    this->rhoEE() = that.rhoEE();

    return *this;
}

instantaneous&
instantaneous::operator=(
        const summary& that)
{
    // Ensure adequate storage for the assignment
    // FIXME Usage of q::count brittle wrt subclasses-- NoChange preferable
#ifndef NDEBUG
    setConstant(that.storage.rows(), q::count,
                std::numeric_limits<Scalar>::quiet_NaN());
#else
    resize     (that.storage.rows(), q::count);
#endif

    // Assign the data quantity-by-quantity
    this->rho()   = that.bar_rho();
    this->p()     = that.bar_p();
    this->p2()    = that.bar_p2();
    this->u()     = that.bar_u();
    this->v()     = that.bar_v();
    this->w()     = that.bar_w();
    this->uu()    = that.bar_u_u();
    this->uv()    = that.bar_u_v();
    this->uw()    = that.bar_u_w();
    this->vv()    = that.bar_v_v();
    this->vw()    = that.bar_v_w();
    this->ww()    = that.bar_w_w();
    this->rhou()  = that.bar_rho_u();
    this->rhov()  = that.bar_rho_v();
    this->rhow()  = that.bar_rho_w();
    this->rhoE()  = that.bar_rho_E();
    this->rhouu() = that.bar_rho_u_u();
    this->rhouv() = that.bar_rho_u_v();
    this->rhouw() = that.bar_rho_u_w();
    this->rhovv() = that.bar_rho_v_v();
    this->rhovw() = that.bar_rho_v_w();
    this->rhoww() = that.bar_rho_w_w();
    this->rhoEE() = that.bar_rho_E_E();

    return *this;
}

instantaneous&
instantaneous::copy_from(
        const bsplineop& cop,
        const samples& that)
{
    SUZERAIN_ENSURE(cop.get()->method == SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);

    // Ensure adequate storage for the assignment
    // FIXME Usage of q::count brittle wrt subclasses-- NoChange preferable
#ifndef NDEBUG
    setConstant(that.storage.rows(), q::count,
                std::numeric_limits<Scalar>::quiet_NaN());
#else
    resize     (that.storage.rows(), q::count);
#endif

    // Convert data, quantity-by-quantity, from coefficients to point values
    // Could perform multiple matvecs here in one shot , but why be error prone?
    // Additionally, in-place matvec operations require an auxiliary buffer,
    // so we wouldn't be saving much data movement in an amortizing operation.
    cop.accumulate(0, 1.0, that.rho()    .col(0).data(), 1, 0.0, this->rho()  .data(), 1);
    cop.accumulate(0, 1.0, that.p()      .col(0).data(), 1, 0.0, this->p()    .data(), 1);
    cop.accumulate(0, 1.0, that.p2()     .col(0).data(), 1, 0.0, this->p2()   .data(), 1);
    cop.accumulate(0, 1.0, that.u()      .col(0).data(), 1, 0.0, this->u()    .data(), 1);
    cop.accumulate(0, 1.0, that.u()      .col(1).data(), 1, 0.0, this->v()    .data(), 1);
    cop.accumulate(0, 1.0, that.u()      .col(2).data(), 1, 0.0, this->w()    .data(), 1);
    cop.accumulate(0, 1.0, that.u_u()    .col(0).data(), 1, 0.0, this->uu()   .data(), 1);
    cop.accumulate(0, 1.0, that.u_u()    .col(1).data(), 1, 0.0, this->uv()   .data(), 1);
    cop.accumulate(0, 1.0, that.u_u()    .col(2).data(), 1, 0.0, this->uw()   .data(), 1);
    cop.accumulate(0, 1.0, that.u_u()    .col(3).data(), 1, 0.0, this->vv()   .data(), 1);
    cop.accumulate(0, 1.0, that.u_u()    .col(4).data(), 1, 0.0, this->vw()   .data(), 1);
    cop.accumulate(0, 1.0, that.u_u()    .col(5).data(), 1, 0.0, this->ww()   .data(), 1);
    cop.accumulate(0, 1.0, that.rho_u()  .col(0).data(), 1, 0.0, this->rhou() .data(), 1);
    cop.accumulate(0, 1.0, that.rho_u()  .col(1).data(), 1, 0.0, this->rhov() .data(), 1);
    cop.accumulate(0, 1.0, that.rho_u()  .col(2).data(), 1, 0.0, this->rhow() .data(), 1);
    cop.accumulate(0, 1.0, that.rho_E()  .col(0).data(), 1, 0.0, this->rhoE() .data(), 1);
    cop.accumulate(0, 1.0, that.rho_u_u().col(0).data(), 1, 0.0, this->rhouu().data(), 1);
    cop.accumulate(0, 1.0, that.rho_u_u().col(1).data(), 1, 0.0, this->rhouv().data(), 1);
    cop.accumulate(0, 1.0, that.rho_u_u().col(2).data(), 1, 0.0, this->rhouw().data(), 1);
    cop.accumulate(0, 1.0, that.rho_u_u().col(3).data(), 1, 0.0, this->rhovv().data(), 1);
    cop.accumulate(0, 1.0, that.rho_u_u().col(4).data(), 1, 0.0, this->rhovw().data(), 1);
    cop.accumulate(0, 1.0, that.rho_u_u().col(5).data(), 1, 0.0, this->rhoww().data(), 1);
    cop.accumulate(0, 1.0, that.rho_E_E().col(0).data(), 1, 0.0, this->rhoEE().data(), 1);

    return *this;
}

} // namespace perfect

} // namespace suzerain
