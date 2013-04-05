//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
 * @copydoc common_block.hpp
 */

#include "common_block.hpp"

#include <suzerain/common.hpp>

namespace suzerain {

namespace perfect {

operator_common_block::operator_common_block()
    : linearization(static_cast<linearize::type>(0))    // Invalid + "false"
    , slow_treatment(static_cast<slowgrowth::type>(0))  // Invalid + "false"
{
    // NOP
}

void
operator_common_block::imexop_ref(suzerain_rholut_imexop_ref   &ref,
                                  suzerain_rholut_imexop_refld &ld)
{
    ref.ux         = ref_ux().data();
    ref.uy         = ref_uy().data();
    ref.uz         = ref_uz().data();
    ref.u2         = ref_u2().data();
    ref.uxux       = ref_uxux().data();
    ref.uxuy       = ref_uxuy().data();
    ref.uxuz       = ref_uxuz().data();
    ref.uyuy       = ref_uyuy().data();
    ref.uyuz       = ref_uyuz().data();
    ref.uzuz       = ref_uzuz().data();
    ref.nu         = ref_nu().data();
    ref.nuux       = ref_nuux().data();
    ref.nuuy       = ref_nuuy().data();
    ref.nuuz       = ref_nuuz().data();
    ref.nuu2       = ref_nuu2().data();
    ref.nuuxux     = ref_nuuxux().data();
    ref.nuuxuy     = ref_nuuxuy().data();
    ref.nuuxuz     = ref_nuuxuz().data();
    ref.nuuyuy     = ref_nuuyuy().data();
    ref.nuuyuz     = ref_nuuyuz().data();
    ref.nuuzuz     = ref_nuuzuz().data();
    ref.ex_gradrho = ref_ex_gradrho().data();
    ref.ey_gradrho = ref_ey_gradrho().data();
    ref.ez_gradrho = ref_ez_gradrho().data();
    ref.e_divm     = ref_e_divm().data();
    ref.e_deltarho = ref_e_deltarho().data();

    const int inc = refs.colStride();
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
