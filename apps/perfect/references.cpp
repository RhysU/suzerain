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

#include "references.hpp"

#include <suzerain/rholut_imexop.h>

namespace suzerain {

namespace perfect {

void
references::imexop_ref(
        suzerain_rholut_imexop_ref   &ref,
        suzerain_rholut_imexop_refld &ld)
{
    ref.ux         = ux().data();
    ref.uy         = uy().data();
    ref.uz         = uz().data();
    ref.u2         = u2().data();
    ref.uxux       = uxux().data();
    ref.uxuy       = uxuy().data();
    ref.uxuz       = uxuz().data();
    ref.uyuy       = uyuy().data();
    ref.uyuz       = uyuz().data();
    ref.uzuz       = uzuz().data();
    ref.nu         = nu().data();
    ref.nuux       = nu_ux().data();
    ref.nuuy       = nu_uy().data();
    ref.nuuz       = nu_uz().data();
    ref.nuu2       = nu_u2().data();
    ref.nuuxux     = nu_uxux().data();
    ref.nuuxuy     = nu_uxuy().data();
    ref.nuuxuz     = nu_uxuz().data();
    ref.nuuyuy     = nu_uyuy().data();
    ref.nuuyuz     = nu_uyuz().data();
    ref.nuuzuz     = nu_uzuz().data();
    ref.ex_gradrho = ex_gradrho().data();
    ref.ey_gradrho = ey_gradrho().data();
    ref.ez_gradrho = ez_gradrho().data();
    ref.e_divm     = e_divm().data();
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
