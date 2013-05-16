//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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
 * @copydoc constraint_treatment.hpp
 */

#include "constraint_treatment.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/constraint.hpp>
#include <suzerain/error.h>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

#include "common_block.hpp"

namespace suzerain {

namespace perfect {

constraint_treatment::constraint_treatment(
            const real_t& Ma,
            const pencil_grid& dgrid,
            operator_common_block& common)
    : implementation_defined()
    , Ma(Ma)
    , common(common)
    , rank_has_zero_zero_modes(dgrid.has_zero_zero_modes())
    , jacobiSvd(0, 0, Eigen::ComputeFullU | Eigen::ComputeFullV)
{
}

void constraint_treatment::apply_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const std::size_t substep_index) const
{
    return L->apply_mass_plus_scaled_operator(
            phi, state, substep_index);
}

void
constraint_treatment::accumulate_mass_plus_scaled_operator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        contiguous_state<4,complex_t> &output,
        const std::size_t substep_index) const
{
    return L->accumulate_mass_plus_scaled_operator(
            phi, input, beta, output, substep_index);
}

void
constraint_treatment::invert_mass_plus_scaled_operator(
        const complex_t& phi,
        multi_array::ref<complex_t,4>& state,
        const timestepper::lowstorage::method_interface<complex_t>& method,
        const real_t delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> *ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    SUZERAIN_TIMER_SCOPED("constraint_treatment");

    // Sidesteps assertions when local rank contains no wavespace information
    const int Ny = (int) state.shape()[1];
    if (SUZERAIN_UNLIKELY(0 == Ny)) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    // Any amount of incoming state is valid so long as there's enough there
    SUZERAIN_ENSURE(state.strides()[1] ==            1       );
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny      );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::e  );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::mx );
    SUZERAIN_ENSURE(state.shape()  [0] >  (unsigned) ndx::rho);

    // See channel_treatment writeup (redux) for information on the steps below.

    // Have a tantrum if caller expects us to compute any constraints
    SUZERAIN_ENSURE(ic0 == NULL);

    // On zero-zero rank, re-use ic0 to wrap cdata for BaseClass invocation
    if (rank_has_zero_zero_modes) {

        // Prepare RHS data for density, momentum, and total energy constraints
        //
        // Notice scaling by Mach^2 to cause mx-, my-, and mz-related forcing
        //     work to have nondimensional energy "units" because we will
        //     directly add the result to the total energy equation.
        const real_t Ma2 = Ma * Ma;
        cdata.setZero(state.shape()[0]*Ny, cdata.cols());
        cdata.col(ndx::e  ).segment(ndx::e   * Ny, Ny).setOnes();
        cdata.col(ndx::mx ).segment(ndx::e   * Ny, Ny).real() = Ma2*common.u();
        cdata.col(ndx::mx ).segment(ndx::mx  * Ny, Ny).setOnes();
        cdata.col(ndx::my ).segment(ndx::e   * Ny, Ny).real() = Ma2*common.v();
        cdata.col(ndx::my ).segment(ndx::my  * Ny, Ny).setOnes();
        cdata.col(ndx::mz ).segment(ndx::e   * Ny, Ny).real() = Ma2*common.w();
        cdata.col(ndx::mz ).segment(ndx::mz  * Ny, Ny).setOnes();
        cdata.col(ndx::rho).segment(ndx::rho * Ny, Ny).setOnes();

        // Wrap data into appropriately digestible format
        using std::size_t;
        const array<size_t,4> sizes = {{
                state.shape()[0], (size_t) Ny, (size_t) cdata.cols(), 1
        }};
        ic0 = new multi_array::ref<complex_t,4>(
                cdata.data(), sizes, storage::interleaved<4>());
    }

    // Delegate to wrapped operator for the solution procedure
    L->invert_mass_plus_scaled_operator(
            phi, state, method, delta_t, substep_index, ic0);

    // Clean up any integral constraint data wrapper we may have employed
    delete ic0;

    // Only the rank with zero-zero modes proceeds!
    if (!rank_has_zero_zero_modes) return;

    // Get an Eigen-friendly map of the zero-zero mode coefficients
    Map<VectorXc> mean(state.origin(), state.shape()[0]*Ny);

    // Solve the requested, possibly simultaneous constraint problem.  A fancy
    // decomposition is used for a simple 5x5 solve or to permit one or more
    // inactive constraints via least squares.  Least squares also adds
    // robustness whenever the constraints are somehow incompatible.
    //
    // 1) Assemble the matrix problem for simultaneous integral constraints
    Matrix5r cmat(Matrix5r::Zero());
    assert(static_cast<int>(ndx::e  ) < cmat.rows());
    assert(static_cast<int>(ndx::mx ) < cmat.rows());
    assert(static_cast<int>(ndx::my ) < cmat.rows());
    assert(static_cast<int>(ndx::mz ) < cmat.rows());
    assert(static_cast<int>(ndx::rho) < cmat.rows());
    for (int j = 0; j < cmat.cols(); ++j) {
        const constraint::base * const cj = operator[](j).get();
        if (cj) {
            for (int i = 0; i < cmat.rows(); ++i) {
                cmat(j,i) = cj->coeff.dot(cdata.col(i).segment(j*Ny,Ny).real());
            }
        }
    }
    // 2) Prepare initial values relative to desired targets when active.
    //    Otherwise, zero the associated row and column in matrix.
    Vector5r crhs(Vector5r::Zero());
    for (int j = 0; j < crhs.rows(); ++j) {
        const constraint::base * const cj = operator[](j).get();
        if (cj && cj->enabled()) {
            crhs[j] = cj->target()
                    - cj->coeff.dot(mean.segment(j*Ny,Ny).real());
        } else {
            cmat.col(j).setZero();
            cmat.row(j).setZero();
        }
    }
    // 3) Prepare matrix decomposition and solve using least squares
    Vector5r cphi = jacobiSvd.compute(cmat).solve(crhs);
    // 4) Add correctly scaled results to the mean state to satisfy constraints
    mean += (cdata * cphi.asDiagonal()).rowwise().sum();

    // The implicitly applied integral constraints, as coefficients, must be
    // averaged across each substep to permit accounting for their impact on
    // the Reynolds averaged equations using method_interface::iota as in
    //
    //    mean += iota * ((sample / delta_t) - mean).
    //
    // The delta_t accounts for step sizes already implicitly included in cphi.
    //
    // Notice mx-related forcing is NOT scaled by Mach^2 when tracked
    // because our post-processing routines will account for Mach^2 factor.
    cphi                 /= delta_t; // Henceforth includes 1/delta_t scaling!
    const real_t iota     = method.iota(substep_index);
    common.fx()          += iota * (
                                ArrayX1r::Constant(Ny, cphi[ndx::mx])
                              - common.fx()
                            );
    common.fy()          += iota * (
                                ArrayX1r::Constant(Ny, cphi[ndx::my])
                              - common.fy()
                            );
    common.fz()          += iota * (
                                ArrayX1r::Constant(Ny, cphi[ndx::mz])
                              - common.fz()
                            );
    common.f_dot_u()     += iota * (
                                cphi[ndx::mx] * common.u()
                              + cphi[ndx::my] * common.v()
                              + cphi[ndx::mz] * common.w()
                              - common.f_dot_u()
                            );
    common.qb()          += iota * (
                                ArrayX1r::Constant(Ny, cphi[ndx::e])
                              - common.qb()
                            );
    common.CrhoE()       += iota * (/* zero */ - common.CrhoE());
    common.Crhou()       += iota * (/* zero */ - common.Crhou());
    common.Crhov()       += iota * (/* zero */ - common.Crhov());
    common.Crhow()       += iota * (/* zero */ - common.Crhow());
    common.Crho()        += iota * (
                                ArrayX1r::Constant(Ny, cphi[ndx::rho])
                              - common.Crho()
                            );
    common.Crhou_dot_u() += iota * (/* zero */ - common.Crhou_dot_u());

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace perfect

} // namespace suzerain
