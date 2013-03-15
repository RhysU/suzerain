//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_ZGBSV_BSMBSM_SOLVER_HPP
#define SUZERAIN_ZGBSV_BSMBSM_SOLVER_HPP

/** @file
 * Encapsulates solving BSMBSM problems per a \ref zgbsv_specification.
 * @see The definition of a BSMBSM problem in \ref bsmbsm.h.
 */

#include <suzerain/common.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/zgbsv_specification.hpp>

namespace suzerain {

// TODO Implement

struct zgbsv_bsmbsm_solver : public suzerain_bsmbsm
{

    zgbsv_bsmbsm_solver(const suzerain_bsmbsm& bsmbsm);

    typedef Map<ArrayXXc, Aligned, OuterStride<Dynamic> > A_type;

    typedef Map<ArrayXXc, Aligned, OuterStride<Dynamic> > LU_type;

    typedef Map<ArrayXc,  Aligned                       > x_type;

    A_type A;

    LU_type LU;

    ArrayXi ipiv;

    ArrayXc b;

    x_type x;

};

} // namespace suzerain

#endif // SUZERAIN_ZGBSV_BSMBSM_SOLVER_HPP
