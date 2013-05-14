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
#include <suzerain/error.h>
#include <suzerain/grid_specification.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

#include "common_block.hpp"

namespace suzerain {

namespace perfect {

constraint_treatment::constraint_treatment(
            real_t& Ma,
            const grid_specification& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b,
            operator_common_block& common)
    : operator_base(grid, dgrid, cop, b)
    , Ma(Ma)
    , common(common)
    , rho(constraint())
    , mx(constraint())
    , e(constraint())
    , coeff_bulk(b.n())
    , coeff_rho(b.n())
    , coeff_mx(b.n())
    , coeff_e(b.n())
    , jacobiSvd(0, 0, Eigen::ComputeFullU | Eigen::ComputeFullV)
{
    // Precompute operator for finding bulk quantities from coefficients
    b.integration_coefficients(0, coeff_bulk.data());
    coeff_bulk /= grid.L.y();
}

constraint_treatment&
constraint_treatment::specify(const constraint& src,
                                    constraint& dst,
                                    VectorXr&   dstcoeff)
{
    dst = src;

    switch (src.what) {
    case constraint::nothing:
        dstcoeff.resize(0);
        break;
    case constraint::value_lower:
        dstcoeff.setZero(coeff_bulk.rows());
        dstcoeff.head<1>()[0] = 1;
        break;
    case constraint::value_upper:
        dstcoeff.setZero(coeff_bulk.rows());
        dstcoeff.tail<1>()[0] = 1;
        break;
    case constraint::value_bulk:
        dstcoeff = coeff_bulk;
        break;
    default:
        SUZERAIN_ERROR_VAL_UNIMPLEMENTED(*this);
    }

    return *this;
}

constraint_treatment&
constraint_treatment::specify_rho(const constraint& src)
{
    return this->specify(src, rho, coeff_rho);
}

constraint_treatment&
constraint_treatment::specify_rho_u(const constraint& src)
{
    return this->specify(src, mx, coeff_mx);
}

constraint_treatment&
constraint_treatment::specify_rho_E(const constraint& src)
{
    return this->specify(src, e, coeff_e);
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

    // Shorthand
    using std::size_t;
    const int Ny = this->dgrid.global_wave_extent.y();

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    // Any amount of incoming state is valid so long as there's enough there
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny      );
    SUZERAIN_ENSURE(state.strides()[1] ==             1      );
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
    if (dgrid.has_zero_zero_modes()) {

        // Prepare RHS data for density, momentum, and total energy constraints
        //
        // Notice scaling by Mach^2 to cause mx-related forcing
        //     work to have nondimensional energy "units" because we will
        //     directly add the result to the total energy equation.
        cdata.setZero(state.shape()[0]*Ny, cdata.cols());
        cdata.col(0).segment(ndx::rho * Ny, Ny).setOnes();
        cdata.col(1).segment(ndx::e   * Ny, Ny).real() = Ma * Ma * common.u();
        cdata.col(1).segment(ndx::mx  * Ny, Ny).setOnes();
        cdata.col(2).segment(ndx::e   * Ny, Ny).setOnes();

        // Wrap data into appropriately digestible format
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
    if (!dgrid.has_zero_zero_modes()) return;

    // Get an Eigen-friendly map of the zero-zero mode coefficients
    Map<VectorXc> mean(state.origin(), state.shape()[0]*Ny);

    // Solve the requested, possibly simultaneous constraint problem.  A fancy
    // decomposition is used for a simple 3x3 solve or to permit one or more
    // inactive constraints via least squares.  Least squares also adds
    // robustness if constraints incompatible.
    //
    // 1) Assemble the matrix problem for simultaneous integral constraints
    // Notice Matrix3r::operator<< depicts the matrix in row-major fashion
    Matrix3r cmat;
    cmat << coeff_rho.dot(cdata.col(0).segment(ndx::rho * Ny, Ny).real())
         ,  coeff_rho.dot(cdata.col(1).segment(ndx::rho * Ny, Ny).real())
         ,  coeff_rho.dot(cdata.col(2).segment(ndx::rho * Ny, Ny).real())
         ,  coeff_mx .dot(cdata.col(0).segment(ndx::mx  * Ny, Ny).real())
         ,  coeff_mx .dot(cdata.col(1).segment(ndx::mx  * Ny, Ny).real())
         ,  coeff_mx .dot(cdata.col(2).segment(ndx::mx  * Ny, Ny).real())
         ,  coeff_e  .dot(cdata.col(0).segment(ndx::e   * Ny, Ny).real())
         ,  coeff_e  .dot(cdata.col(1).segment(ndx::e   * Ny, Ny).real())
         ,  coeff_e  .dot(cdata.col(2).segment(ndx::e   * Ny, Ny).real())
         ;
    // 2) Prepare initial values relative to desired targets when active.
    // Otherwise, zero the associated row and column in matrix.
    Vector3r crhs;
    if (rho.enabled()) {
        crhs[0] = rho.target
                - coeff_rho.dot(mean.segment(ndx::rho * Ny, Ny).real());
    } else {
        crhs[0] = 0;
        cmat.col(0).setZero();
        cmat.row(0).setZero();
    }
    if (mx.enabled()) {
        crhs[1] = mx.target
                - coeff_mx. dot(mean.segment(ndx::mx  * Ny, Ny).real());
    } else {
        crhs[1] = 0;
        cmat.col(1).setZero();
        cmat.row(1).setZero();
    }
    if (e.enabled()) {
        crhs[2] = e.target
                - coeff_e.  dot(mean.segment(ndx::e   * Ny, Ny).real());
    } else {
        crhs[2] = 0;
        cmat.col(2).setZero();
        cmat.row(2).setZero();
    }
    // 3) Prepare matrix decomposition and solve using least squares
    Vector3r cphi = jacobiSvd.compute(cmat).solve(crhs);
    // 4) Add correctly scaled results to the mean state to satisfy constraints
    mean += cphi(0)*cdata.col(0) + cphi(1)*cdata.col(1) + cphi(2)*cdata.col(2);

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
                                ArrayX1r::Constant(Ny, cphi(1))
                              - common.fx()
                            );
    common.fy()          += iota * (/* zero */ - common.fy());
    common.fz()          += iota * (/* zero */ - common.fz());
    common.f_dot_u()     += iota * (
                                cphi(1) * common.u()
                              - common.f_dot_u()
                            );
    common.qb()          += iota * (
                                ArrayX1r::Constant(Ny, cphi(2))
                              - common.qb()
                            );
    common.CrhoE()       += iota * (/* zero */ - common.CrhoE());
    common.Crhou()       += iota * (/* zero */ - common.Crhou());
    common.Crhov()       += iota * (/* zero */ - common.Crhov());
    common.Crhow()       += iota * (/* zero */ - common.Crhow());
    common.Crho()        += iota * (
                                ArrayX1r::Constant(Ny, cphi(0))
                              - common.Crho()
                            );
    common.Crhou_dot_u() += iota * (/* zero */ - common.Crhou_dot_u());

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace perfect

} // namespace suzerain
