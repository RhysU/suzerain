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

zgbsv_bsmbsm_solver_zgbsvx::zgbsv_bsmbsm_solver_zgbsvx(
        const zgbsv_specification& specification,
        const suzerain_bsmbsm&    bsmbsm)
    : zgbsv_bsmbsm_solver(specification, bsmbsm)
    , r(N)       // Per zgbsvx requirements
    , c(N)       // Per zgbsvx requirements
    , work(2*N)  // Per zgbsvx requirements
    , rwork(N)   // Per zgbsvx requirements
    , A_(LD, N)  // Operator storage for out-of-place factorization
    , x_(N)      // Solution storage for out-of-place solution
{
    if (method() != zgbsv_specification::zgbsvx) throw std::invalid_argument(
            "Invalid method() in zgbsv_bsmbsm_solver_zgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(in_place() == false);
    new (&A) A_type(A_.data(), A_.rows(), A_.cols(), A_.colStride());
    new (&x) x_type(x_.data(), x_.rows());
}

zgbsv_bsmbsm_solver_zcgbsvx::zgbsv_bsmbsm_solver_zcgbsvx(
        const zgbsv_specification& specification,
        const suzerain_bsmbsm&    bsmbsm)
    : zgbsv_bsmbsm_solver(specification, bsmbsm)
    , work(2*N)  // Per zcgbsvx requirements
    , A_(LD, N)  // Operator storage for out-of-place factorization
    , x_(N)      // Solution storage for out-of-place solution
{
    if (method() != zgbsv_specification::zcgbsvx) throw std::invalid_argument(
            "Invalid method() in zgbsv_bsmbsm_solver_zcgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(in_place() == false);
    new (&A) A_type(A_.data(), A_.rows(), A_.cols(), A_.colStride());
    new (&x) x_type(x_.data(), x_.rows());
}

} // end namespace suzerain
