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
 * @copydoc bsmbsm_solver.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/bsmbsm_solver.hpp>

namespace suzerain {

// TODO Finish implementing

bsmbsm_solver::bsmbsm_solver(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& specification)
    : suzerain_bsmbsm(bsmbsm)
    , zgbsv_specification(specification)
    , LU(KL + LD, N)
    , ipiv(N)
    , b(N)
    , A(KL + LU.data(), LD, N, KL + KU) // Aliases LU
    , x(b.data(), N)                    // Aliases b
{
    // NOP
}

bsmbsm_solver_zgbsv::bsmbsm_solver_zgbsv(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& specification)
    : bsmbsm_solver(bsmbsm, specification)
{
    if (method() != zgbsv_specification::zgbsv) throw std::invalid_argument(
            "Invalid method() in bsmbsm_solver_zgbsv");
    assert(in_place() == true);
}

bsmbsm_solver_zgbsvx::bsmbsm_solver_zgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& specification)
    : bsmbsm_solver(bsmbsm, specification)
    , r(N)       // Per zgbsvx requirements
    , c(N)       // Per zgbsvx requirements
    , work(2*N)  // Per zgbsvx requirements
    , rwork(N)   // Per zgbsvx requirements
    , A_(LD, N)  // Operator storage for out-of-place factorization
    , x_(N)      // Solution storage for out-of-place solution
{
    if (method() != zgbsv_specification::zgbsvx) throw std::invalid_argument(
            "Invalid method() in bsmbsm_solver_zgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(in_place() == false);
    new (&A) A_type(A_.data(), A_.rows(), A_.cols(), A_.colStride());
    new (&x) x_type(x_.data(), x_.rows());
}

bsmbsm_solver_zcgbsvx::bsmbsm_solver_zcgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& specification)
    : bsmbsm_solver(bsmbsm, specification)
    , work(2*N)  // Per zcgbsvx requirements
    , A_(LD, N)  // Operator storage for out-of-place factorization
    , x_(N)      // Solution storage for out-of-place solution
{
    if (method() != zgbsv_specification::zcgbsvx) throw std::invalid_argument(
            "Invalid method() in bsmbsm_solver_zcgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(in_place() == false);
    new (&A) A_type(A_.data(), A_.rows(), A_.cols(), A_.colStride());
    new (&x) x_type(x_.data(), x_.rows());
}

} // end namespace suzerain
