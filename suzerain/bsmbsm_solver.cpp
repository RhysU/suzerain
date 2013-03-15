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
#include <suzerain/complex.hpp>

namespace suzerain {

// TODO Finish implementing

bsmbsm_solver::bsmbsm_solver(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec)
    : suzerain_bsmbsm(bsmbsm)
    , spec(spec)
    , LU(KL + LD, N)
    , Pb(N)
    , PAPT(KL + LU.data(), LD, N, KL + KU) // Aliases LU
    , Px(Pb.data(), N)                     // Aliases Pb
    , ipiv(N)
{
    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    LU  .setConstant(suzerain::complex::NaN<suzerain::complex_t>());
    Pb  .setConstant(suzerain::complex::NaN<suzerain::complex_t>());
    ipiv.setConstant(-12345);
#endif
}

bsmbsm_solver_zgbsv::bsmbsm_solver_zgbsv(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec)
    : bsmbsm_solver(bsmbsm, spec)
{
    if (spec.method() != zgbsv_specification::zgbsv)
        throw std::invalid_argument("Invalid spec in bsmbsm_solver_zgbsv");
    assert(spec.in_place() == true);
}

bsmbsm_solver_zgbsvx::bsmbsm_solver_zgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec)
    : bsmbsm_solver(bsmbsm, spec)
    , r(N)          // Per zgbsvx requirements
    , c(N)          // Per zgbsvx requirements
    , work(2*N)     // Per zgbsvx requirements
    , rwork(N)      // Per zgbsvx requirements
    , PAPT_(LD, N)  // Operator storage for out-of-place factorization
    , Px_(N)        // Solution storage for out-of-place solution
{
    if (spec.method() != zgbsv_specification::zgbsvx)
        throw std::invalid_argument("Invalid method in bsmbsm_solver_zgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(spec.in_place() == false);
    new (&PAPT) PAPT_type(PAPT_.data(), PAPT_.rows(),
                          PAPT_.cols(), PAPT_.colStride());
    new (&Px)   Px_type(Px_.data(), Px_.rows());

    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    r    .setConstant(std::numeric_limits<suzerain::real_t>::quiet_NaN());
    c    .setConstant(std::numeric_limits<suzerain::real_t>::quiet_NaN());
    work .setConstant(suzerain::complex::NaN<suzerain::complex_t>());
    rwork.setConstant(std::numeric_limits<suzerain::real_t>::quiet_NaN());
    PAPT_.setConstant(suzerain::complex::NaN<suzerain::complex_t>());
    Px_  .setConstant(suzerain::complex::NaN<suzerain::complex_t>());
#endif
}

bsmbsm_solver_zcgbsvx::bsmbsm_solver_zcgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec)
    : bsmbsm_solver(bsmbsm, spec)
    , work(2*N)     // Per zcgbsvx requirements
    , PAPT_(LD, N)  // Operator storage for out-of-place factorization
    , Px_(N)        // Solution storage for out-of-place solution
{
    if (spec.method() != zgbsv_specification::zcgbsvx)
        throw std::invalid_argument("Invalid method in bsmbsm_solver_zcgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(spec.in_place() == false);
    new (&PAPT) PAPT_type(PAPT_.data(), PAPT_.rows(),
                          PAPT_.cols(), PAPT_.colStride());
    new (&Px) Px_type(Px_.data(), Px_.rows());

    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    work .setConstant(suzerain::complex::NaN<suzerain::complex_t>());
    PAPT_.setConstant(suzerain::complex::NaN<suzerain::complex_t>());
    Px_  .setConstant(suzerain::complex::NaN<suzerain::complex_t>());
#endif
}

} // end namespace suzerain
