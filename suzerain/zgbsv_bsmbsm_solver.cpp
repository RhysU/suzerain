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

/** @file
 * @copydoc zgbsv_bsmbsm_solver.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/zgbsv_bsmbsm_solver.hpp>

namespace suzerain {

// TODO Finish implementing

zgbsv_bsmbsm_solver::zgbsv_bsmbsm_solver(
        const zgbsv_specification& specification,
        const suzerain_bsmbsm&     bsmbsm)
    : zgbsv_specification(specification)
    , suzerain_bsmbsm(bsmbsm)
    , LU(KL + LD, N)
    , ipiv(N)
    , b(N)
    , A(KL + LU.data(), LD, N, KL + KU) // Aliases LU
    , x(b.data(), N)                    // Aliases b
{
    // NOP
}

zgbsv_bsmbsm_solver_zgbsv::zgbsv_bsmbsm_solver_zgbsv(
        const zgbsv_specification& specification,
        const suzerain_bsmbsm&    bsmbsm)
    : zgbsv_bsmbsm_solver(specification, bsmbsm)
{
    if (method() != zgbsv_specification::zgbsv) throw std::invalid_argument(
            "Invalid method() in zgbsv_bsmbsm_solver_zgbsv");
    assert(in_place() == true);
}

} // end namespace suzerain
