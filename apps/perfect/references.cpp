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
 * @copydoc references.hpp
 */

#include "references.hpp"

#include <suzerain/rholut_imexop.h>

namespace suzerain {

namespace perfect {

references::references()
    : super(static_cast<super::Index>(q::count), 0)
{
}

references::~references()
{
}

void
references::set_zero(const int Ny)
{
    super::setZero(NoChange, Ny);
}

void
references::imexop_ref(
        suzerain_rholut_imexop_ref   &ref,
        suzerain_rholut_imexop_refld &ld)
{
    ref.ux         = u()         .data();
    ref.uy         = v()         .data();
    ref.uz         = w()         .data();
    ref.u2         = u2()        .data();
    ref.uxux       = uu()        .data();
    ref.uxuy       = uv()        .data();
    ref.uxuz       = uw()        .data();
    ref.uyuy       = vv()        .data();
    ref.uyuz       = vw()        .data();
    ref.uzuz       = ww()        .data();
    ref.nu         = nu()        .data();
    ref.nuux       = nu_u()      .data();
    ref.nuuy       = nu_v()      .data();
    ref.nuuz       = nu_w()      .data();
    ref.nuu2       = nu_u2()     .data();
    ref.nuuxux     = nu_uu()     .data();
    ref.nuuxuy     = nu_uv()     .data();
    ref.nuuxuz     = nu_uw()     .data();
    ref.nuuyuy     = nu_vv()     .data();
    ref.nuuyuz     = nu_vw()     .data();
    ref.nuuzuz     = nu_ww()     .data();
    ref.ex_gradrho = ex_gradrho().data();
    ref.ey_gradrho = ey_gradrho().data();
    ref.ez_gradrho = ez_gradrho().data();
    ref.e_divm     = e_divm()    .data();
    ref.e_deltarho = e_deltarho().data();

    const int inc = innerStride();
    ld.ux         = inc;
    ld.uy         = inc;
    ld.uz         = inc;
    ld.u2         = inc;
    ld.uxux       = inc;
    ld.uxuy       = inc;
    ld.uxuz       = inc;
    ld.uyuy       = inc;
    ld.uyuz       = inc;
    ld.uzuz       = inc;
    ld.nu         = inc;
    ld.nuux       = inc;
    ld.nuuy       = inc;
    ld.nuuz       = inc;
    ld.nuu2       = inc;
    ld.nuuxux     = inc;
    ld.nuuxuy     = inc;
    ld.nuuxuz     = inc;
    ld.nuuyuy     = inc;
    ld.nuuyuz     = inc;
    ld.nuuzuz     = inc;
    ld.ex_gradrho = inc;
    ld.ey_gradrho = inc;
    ld.ez_gradrho = inc;
    ld.e_divm     = inc;
    ld.e_deltarho = inc;
}

} // namespace perfect

} // namespace suzerain
