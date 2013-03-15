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
#include <suzerain/error.h>

namespace suzerain {

// TODO Finish implementing

bsmbsm_solver::bsmbsm_solver(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec,
        const int                  nrhs)
    : suzerain_bsmbsm(bsmbsm)
    , spec(spec)
    , LU(KL + LD, N)
    , PB(N, nrhs)
    , PAPT(KL + LU.data(), LD, N, KL + KU) // Aliases LU
    , PX(PB.data(), PB.rows(), PB.cols())  // Aliases PB
    , ipiv(N)
    , fact_(default_fact())
    , apprx_(0)
{
    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    LU  .setConstant(suzerain::complex::NaN<complex_double>());
    PB  .setConstant(suzerain::complex::NaN<complex_double>());
    ipiv.setConstant(-12345);
#endif
}

int
bsmbsm_solver::solve_internal(const char trans,
                              const int nrhs)
{
    if (SUZERAIN_UNLIKELY(nrhs < 1 || nrhs > PB.cols()))
        throw std::invalid_argument("Invalid nrhs supplied to solve()");
    const int info = solve_hook(trans, nrhs); // Invoke subclass-specific hook
    if (info == 0) return info;               // Eagerly return on success

    // Otherwise, loudly report any errors that occurred during the solve
    char buffer[128];
    if (info < 0) {
        snprintf(buffer, sizeof(buffer),
            "%s reported error in argument %d",
            spec.mname(), -info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    } else if (info <= N) {
        snprintf(buffer, sizeof(buffer),
            "%s reported singularity in PAP^T row %d"
            " corresponding to A row %d for state scalar %d",
            spec.mname(), info - 1, suzerain_bsmbsm_q(S, n, info - 1),
            suzerain_bsmbsm_q(S, n, info - 1) / n);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    } else {
        snprintf(buffer, sizeof(buffer),
            "%s reported unknown error %d", spec.mname(), info);
        SUZERAIN_ERROR(buffer, SUZERAIN_ESANITY);
    }
}

bsmbsm_solver_zgbsv::bsmbsm_solver_zgbsv(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec,
        const int                  nrhs)
    : bsmbsm_solver(bsmbsm, spec, nrhs)
{
    if (spec.method() != zgbsv_specification::zgbsv)
        throw std::invalid_argument("Invalid spec in bsmbsm_solver_zgbsv");
    assert(spec.in_place() == true);
}

int
bsmbsm_solver_zgbsv::solve_hook(
        const char trans,
        const int nrhs)
{
    SUZERAIN_UNUSED(trans);
    SUZERAIN_UNUSED(nrhs);

    return -1; // FIXME Implement
}

bsmbsm_solver_zgbsvx::bsmbsm_solver_zgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec,
        const int                  nrhs)
    : bsmbsm_solver(bsmbsm, spec, nrhs)
    , equed_('N')   // Default for non-factorized per zgbsvx
    , r_(N)         // Per zgbsvx requirements
    , c_(N)         // Per zgbsvx requirements
    , work_(N, 2)   // Per zgbsvx requirements
    , rwork_(N)     // Per zgbsvx requirements
    , rcond_(std::numeric_limits<double>::quiet_NaN())
    , ferr_(nrhs)   // Per zgbsvx requirements
    , berr_(nrhs)   // Per zgbsvx requirements
    , PAPT_(LD, N)  // Operator storage for out-of-place factorization
    , PX_(N, nrhs)  // Solution storage for out-of-place solution
{
    if (spec.method() != zgbsv_specification::zgbsvx)
        throw std::invalid_argument("Invalid method in bsmbsm_solver_zgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(spec.in_place() == false);
    new (&PAPT) PAPT_type(PAPT_.data(), PAPT_.rows(),
                          PAPT_.cols(), PAPT_.colStride());
    new (&PX)   PX_type(PX_.data(), PX_.rows(), PX_.cols());

    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    r_    .setConstant(std::numeric_limits<double>::quiet_NaN());
    c_    .setConstant(std::numeric_limits<double>::quiet_NaN());
    work_ .setConstant(suzerain::complex::NaN<complex_double>());
    rwork_.setConstant(std::numeric_limits<double>::quiet_NaN());
    PAPT_ .setConstant(suzerain::complex::NaN<complex_double>());
    PX_   .setConstant(suzerain::complex::NaN<complex_double>());
#endif
}

int
bsmbsm_solver_zgbsvx::solve_hook(
        const char trans,
        const int nrhs)
{
    SUZERAIN_UNUSED(trans);
    SUZERAIN_UNUSED(nrhs);

    return -1; // FIXME Implement
}

bsmbsm_solver_zcgbsvx::bsmbsm_solver_zcgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const zgbsv_specification& spec,
        const int                  nrhs)
    : bsmbsm_solver(bsmbsm, spec, nrhs)
    , afrob_(-1)    // Per zcgbsvx requirements
    , aiter_(nrhs)  // For nrhs distinct zcgbsvx invocations
    , siter_(nrhs)  // For nrhs distinct zcgbsvx invocations
    , diter_(nrhs)  // For nrhs distinct zcgbsvx invocations
    , tolsc_(nrhs)  // For nrhs distinct zcgbscx invocations
    , work_(N, 2)   // Per zcgbsvx requirements
    , res_(nrhs)    // For nrhs distinct zcgbscx invocations
    , PAPT_(LD, N)  // Operator storage for out-of-place factorization
    , PX_(N, nrhs)  // Solution storage for out-of-place solution
{
    if (spec.method() != zgbsv_specification::zcgbsvx)
        throw std::invalid_argument("Invalid method in bsmbsm_solver_zcgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(spec.in_place() == false);
    new (&PAPT) PAPT_type(PAPT_.data(), PAPT_.rows(),
                          PAPT_.cols(), PAPT_.colStride());
    new (&PX)   PX_type(PX_.data(), PX_.rows(), PX_.cols());

    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    work_.setConstant(suzerain::complex::NaN<complex_double>());
    PAPT_.setConstant(suzerain::complex::NaN<complex_double>());
    PX_  .setConstant(suzerain::complex::NaN<complex_double>());
#endif
}

int
bsmbsm_solver_zcgbsvx::solve_hook(
        const char trans,
        const int nrhs)
{
    SUZERAIN_UNUSED(trans);
    SUZERAIN_UNUSED(nrhs);

    return -1; // FIXME Implement
}

} // end namespace suzerain
