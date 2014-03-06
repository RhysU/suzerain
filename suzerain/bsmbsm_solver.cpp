//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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
 * @copydoc bsmbsm_solver.hpp
 */

#include <suzerain/bsmbsm_solver.hpp>

#include <suzerain/blas_et_al.h>
#include <suzerain/complex.hpp>
#include <suzerain/error.h>

namespace suzerain {

bsmbsm_solver*
bsmbsm_solver::build(
        const suzerain_bsmbsm&     bsmbsm,
        const specification_zgbsv& spec,
        const int                  nrhs)
{
    switch (spec.method()) {

    case specification_zgbsv::zgbsv:
        return new bsmbsm_solver_zgbsv(bsmbsm, spec, nrhs);

    case specification_zgbsv::zgbsvx:
        return new bsmbsm_solver_zgbsvx(bsmbsm, spec, nrhs);

    case specification_zgbsv::zcgbsvx:
        return new bsmbsm_solver_zcgbsvx(bsmbsm, spec, nrhs);

    }

    throw std::invalid_argument("bsmbsm_solver::build: unknown spec.method()");
}

bsmbsm_solver::bsmbsm_solver(
        const suzerain_bsmbsm&     bsmbsm,
        const specification_zgbsv& spec,
        const int                  nrhs)
    : suzerain_bsmbsm(bsmbsm)
    , spec(spec)
    , LU(KL + LD, N)
    , PB(N, nrhs)
    , PAPT(LU.data() + KL, LU.rows() - KL, LU.cols(), LU.colStride())  // Alias
    , PX(PB.data(), PB.rows(), PB.cols())                              // Alias
    , ipiv(N)
    , fact_(default_fact())
    , apprx_(false)
{
    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    LU  .setConstant(suzerain::complex::NaN<complex_double>());
    PB  .setConstant(suzerain::complex::NaN<complex_double>());
    ipiv.setConstant(-12345);
#endif
}

bsmbsm_solver&
bsmbsm_solver::supplied_PAPT()
{
    if (spec.reuse()) {
        apprx_ = fact_ != default_fact();
    } else {
        fact_ = default_fact();
    }
    return *this;
}

bool
bsmbsm_solver::apprx(const bool acceptable)
{
    const bool old = apprx_;
    if (spec.reuse()) {
        apprx_ = apprx_ && acceptable;
        if (!apprx_) {
            fact_ = default_fact();
        }
    }
    return old;
}

std::vector<std::string>
bsmbsm_solver::summarize_statistics() const
{
    std::vector<std::string> retval;
    std::ostringstream msg;
    msg << "Employed bsmbsm_solver per specification " << spec;
    retval.push_back(msg.str());
    return retval;
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
        const specification_zgbsv& spec,
        const int                  nrhs)
    : bsmbsm_solver(bsmbsm, spec, nrhs)
{
    if (spec.method() != specification_zgbsv::zgbsv)
        throw std::invalid_argument("Invalid spec in bsmbsm_solver_zgbsv");
    assert(spec.in_place() == true);
}

int
bsmbsm_solver_zgbsv::solve_hook(
        const char trans,
        const int nrhs)
{
    assert(apprx_ == false);
    // return suzerain_lapackext_zgbsv(&fact_, trans, N, KL, KU, nrhs,
    //                                 LU.data(), LU.colStride(), ipiv.data(),
    //                                 PB.data(), PB.colStride());

    // Logic below should be exactly equivalent to above commented out
    // call to suzerain_lapackext_zgbsv.  Inlined here to enable
    // separate timing of factorize and backsub steps.
    int info = 0;
    if (toupper(fact_) == 'N') {
        SUZERAIN_TIMER_SCOPED("bsmbsm_solver_zgbsv::suzerain_lapack_zgbtrf");
        info = suzerain_lapack_zgbtrf(N, N, KL, KU, 
                                      LU.data(), LU.colStride(), ipiv.data());
        fact_ = 'F';
    }
    if (!info) {
        SUZERAIN_TIMER_SCOPED("bsmbsm_solver_zgbsv::suzerain_lapack_zgbtrs");
        info = suzerain_lapack_zgbtrs(trans, N, KL, KU, nrhs, 
                                      LU.data(), LU.colStride(), ipiv.data(), 
                                      PB.data(), PB.colStride());
    }
    return info;
}

bsmbsm_solver_zgbsvx::bsmbsm_solver_zgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const specification_zgbsv& spec,
        const int                  nrhs)
    : bsmbsm_solver(bsmbsm, spec, nrhs)
    , equed_('N')   // Default for non-factorized per zgbsvx
    , rcwork_(N, 3) // Stores (r, c, rwork) x nrhs
    , rcond_(std::numeric_limits<double>::quiet_NaN())
    , err_(2, nrhs) // Stores (ferr, berr) x nrhs
    , work_(N, 2)   // Per zgbsvx requirements
    , PAPT_(LD, N)  // Operator storage for out-of-place factorization
    , PX_(N, nrhs)  // Solution storage for out-of-place solution
{
    if (spec.method() != specification_zgbsv::zgbsvx)
        throw std::invalid_argument("Invalid method in bsmbsm_solver_zgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(spec.in_place() == false);
    new (&PAPT) PAPT_type(PAPT_.data(), PAPT_.rows(),
                          PAPT_.cols(), PAPT_.colStride());
    new (&PX)   PX_type(PX_.data(), PX_.rows(), PX_.cols());

    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    rcwork_.setConstant(std::numeric_limits<double>::quiet_NaN());
    err_   .setConstant(std::numeric_limits<double>::quiet_NaN());
    work_  .setConstant(suzerain::complex::NaN<complex_double>());
    PAPT_  .setConstant(suzerain::complex::NaN<complex_double>());
    PX_    .setConstant(suzerain::complex::NaN<complex_double>());
#endif
}

const char * const bsmbsm_solver_zgbsvx::stats_names[
        bsmbsm_solver_zgbsvx::stats_type::static_size] = {
    "row equilibration (equed)",
    "column equilibration (equed)",
    "condition number (1/rcond)",
    "forward error bound (ferr)",
    "backward error bound (berr)"
};

std::vector<std::string>
bsmbsm_solver_zgbsvx::summarize_statistics() const
{
    std::vector<std::string> retval = bsmbsm_solver::summarize_statistics();

    if (stats.count() > 0) {  // Add operational details when solves occurred
        std::ostringstream msg;
        msg.precision(std::numeric_limits<double>::digits10/2 + 1);
        retval.reserve(retval.size() + stats_type::static_size);
        for (std::size_t i = spec.equil() ? 0 : 2; // Skip equed when disabled
             i < stats_type::static_size;
             ++i) {
            msg.str("");
            msg << "Min/avg/max/std of "
                << stats_names[i]
                << ": "
                << stats.min(i) << ", "
                << stats.avg(i) << ", "
                << stats.max(i) << ", "
                << stats.std(i);
            retval.push_back(msg.str());
        }
    }

    return retval;
}

int
bsmbsm_solver_zgbsvx::solve_hook(
        const char trans,
        const int nrhs)
{
#ifndef NDEBUG
    if (fact_ != 'F') {  // NaN r, c, and rwork on 'N' or 'E'...
        rcwork_ .setConstant(std::numeric_limits<double>::quiet_NaN());
    } else {             // ...but maintain r and c on 'F' per ZGBSVX
        rwork_().setConstant(std::numeric_limits<double>::quiet_NaN());
    }
    err_ .setConstant(std::numeric_limits<double>::quiet_NaN());
    work_.setConstant(suzerain::complex::NaN<complex_double>());
#endif

    // Perform the requested solve
    assert(apprx_ == false);
    const int info = suzerain_lapack_zgbsvx(fact_, trans, N, KL, KU, nrhs,
            PAPT.data(), PAPT.colStride(), LU.data(), LU.colStride(),
            ipiv.data(), &equed_, r_().data(), c_().data(), PB.data(),
            PB.colStride(), PX.data(), PX.colStride(), &rcond_,
            ferr_().data(), berr_().data(), work_.data(), rwork_().data());
    fact_ = 'F';

    // Track statistics for each right hand side per stats_type/stats_names
    double samples[stats_type::static_size] = {
        (equed_ == 'R' || equed_ == 'B'),
        (equed_ == 'C' || equed_ == 'B'),
        (1 / rcond_),
        /* ferr below */ std::numeric_limits<double>::quiet_NaN(),
        /* berr below */ std::numeric_limits<double>::quiet_NaN()
    };
    for (int j = 0; j < nrhs; ++j) {
        samples[3] = ferr_()[j];
        samples[4] = berr_()[j];
        stats(samples);
    }

    return info;
}

bsmbsm_solver_zcgbsvx::bsmbsm_solver_zcgbsvx(
        const suzerain_bsmbsm&     bsmbsm,
        const specification_zgbsv& spec,
        const int                  nrhs)
    : bsmbsm_solver(bsmbsm, spec, nrhs)
    , afrob_(-1)          // Per zcgbsvx requirements
    , iter_(2, nrhs)      // Stores (aiter, siter, diter) x nrhs
    , tolscres_(2, nrhs)  // Stores (tolsc, res) x nrhs
    , r_(N)               // Stores solution residual for each Pb.
    , PAPT_(LD, N)        // Operator storage for out-of-place factorization
    , PX_(N, nrhs)        // Solution storage for out-of-place solution
{
    if (spec.method() != specification_zgbsv::zcgbsvx)
        throw std::invalid_argument("Invalid method in bsmbsm_solver_zcgbsvx");

    // See Eigen "Changing the mapped array" documentation for details
    assert(spec.in_place() == false);
    new (&PAPT) PAPT_type(PAPT_.data(), PAPT_.rows(),
                          PAPT_.cols(), PAPT_.colStride());
    new (&PX)   PX_type(PX_.data(), PX_.rows(), PX_.cols());

    // Defensively set NaNs or NaN-like values on debug builds
#ifndef NDEBUG
    iter_    .setConstant(-12345);
    tolscres_.setConstant(std::numeric_limits<double>::quiet_NaN());
    r_       .setConstant(suzerain::complex::NaN<complex_double>());
    PAPT_    .setConstant(suzerain::complex::NaN<complex_double>());
    PX_      .setConstant(suzerain::complex::NaN<complex_double>());
#endif
}

bsmbsm_solver&
bsmbsm_solver_zcgbsvx::supplied_PAPT()
{
    afrob_ = -1;
    return bsmbsm_solver::supplied_PAPT();
}

const char * const bsmbsm_solver_zcgbsvx::stats_names[
        bsmbsm_solver_zcgbsvx::stats_type::static_size] = {
    "single precision LU (fact)",
    "double precision LU (fact)",
    "approx factorization (apprx)",
    "Frobenius norm (afrob)",
    "single precision refinement (siter)",
    "double precision refinement (diter)",
    "fraction of tolerance (tolsc)",
    "residual 2-norm (res)"
};

std::vector<std::string>
bsmbsm_solver_zcgbsvx::summarize_statistics() const
{
    std::vector<std::string> retval = bsmbsm_solver::summarize_statistics();

    if (stats.count() > 0) {  // Add operational details when solves occurred
        std::ostringstream msg;
        msg.precision(std::numeric_limits<double>::digits10/2 + 1);
        retval.reserve(retval.size() + stats_type::static_size);
        for (std::size_t i = 0; i < stats_type::static_size; ++i) {
            // Suppress by-definition uninteresting statistics
            if (i == 0 && (spec.siter() < 0 || spec.diter() < 0)) continue;
            if (i == 1 && (spec.siter() < 0 || spec.diter() < 0)) continue;
            if (i == 2 && !spec.reuse()    ) continue;
#pragma warning(push,disable:1572)
            if (i == 3 && spec.tolsc() == 0) continue;  // Implies afrob == -1
#pragma warning(pop)
            if (i == 4 && spec.siter() <  0) continue;
            if (i == 5 && spec.diter() <  0) continue;
            msg.str("");
            msg << "Min/avg/max/std of "
                << stats_names[i]
                << ": "
                << stats.min(i) << ", "
                << stats.avg(i) << ", "
                << stats.max(i) << ", "
                << stats.std(i);
            retval.push_back(msg.str());
        }
    }

    return retval;
}

int
bsmbsm_solver_zcgbsvx::solve_hook(
        const char trans,
        const int nrhs)
{
#ifndef NDEBUG
    iter_    .setConstant(-12345);
    tolscres_.setConstant(std::numeric_limits<double>::quiet_NaN());
    r_       .setConstant(suzerain::complex::NaN<complex_double>());
#endif

    // Perform the requested solve processing each right hand side in turn
    int info = 0, j = -1;
    while (!info && ++j < nrhs) {

        // Reset specification-related constants for each right hand side
        siter_()[j] = spec.siter();
        diter_()[j] = spec.diter();
        tolsc_()[j] = spec.tolsc();
        info = suzerain_lapackext_zcgbsvx(&fact_, &apprx_, spec.aiter(), trans,
                                          N, KL, KU, PAPT.data(), &afrob_,
                                          LU.data(), ipiv.data(),
                                          PB.data(), PX.data(),
                                          &siter_()[j], &diter_()[j],
                                          &tolsc_()[j], r_.data(), &res_()[j]);

        // Track statistics for each right hand side per stats_type/stats_names
        double samples[stats_type::static_size] = {
            (fact_ == 'S'), (fact_ == 'D'), apprx_, afrob_,
            siter_()[j], diter_()[j], tolsc_()[j], res_()[j]
        };
        stats(samples);
    }

    return info;
}

} // end namespace suzerain
