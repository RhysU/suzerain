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
#include <suzerain/grid_specification.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>

#include "common_block.hpp"
#include "scenario_definition.hpp"

namespace suzerain {

namespace perfect {

constraint_treatment::constraint_treatment(
            const scenario_definition& scenario,
            const grid_specification& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b,
            operator_common_block& common)
    : operator_base(grid, dgrid, cop, b)
    , scenario(scenario)
    , common(common)
    , jacobiSvd(2, 2, Eigen::ComputeFullU | Eigen::ComputeFullV)
{
    // Precomputed results only necessary on rank with zero-zero modes
    if (!dgrid.has_zero_zero_modes()) return;

    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid.L.y();
}

bool constraint_treatment::constrain_bulk_rho_u() const
{
    return !(boost::math::isnan)(scenario.bulk_rho_u);
}

bool constraint_treatment::constrain_bulk_rho() const
{
    return !(boost::math::isnan)(scenario.bulk_rho);
}

void constraint_treatment::apply_mass_plus_scaled_operator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const std::size_t substep_index) const
{
    return L->apply_mass_plus_scaled_operator(
            phi, state, substep_index);
}

void constraint_treatment::accumulate_mass_plus_scaled_operator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        contiguous_state<4,complex_t> &output,
        const std::size_t substep_index) const
{
    return L->accumulate_mass_plus_scaled_operator(
            phi, input, beta, output, substep_index);
}

void constraint_treatment::invert_mass_plus_scaled_operator(
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

        // Prepare data for bulk density and bulk momentum constraints
        //
        // Notice scaling by Mach^2 to cause bulk_rho_u-related forcing
        //     work to have nondimensional energy "units" because we will
        //     directly add the result to the total energy equation.
        cdata.setZero(state.shape()[0]*Ny, cdata.cols());
        cdata.col(0).segment(ndx::rho * Ny, Ny).setOnes();
        cdata.col(1).segment(ndx::e   * Ny, Ny).real() = scenario.Ma
                                                       * scenario.Ma
                                                       * common.u();
        cdata.col(1).segment(ndx::mx  * Ny, Ny).setOnes();

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

    // Assemble the matrix problem for simultaneous integral constraints
    // Notice Matrix2r::operator<< depicts the matrix in row-major fashion
    Vector2r crhs;
    crhs <<    scenario.bulk_rho
             - bulkcoeff.dot(mean.segment(ndx::rho * Ny, Ny).real())
         ,     scenario.bulk_rho_u
             - bulkcoeff.dot(mean.segment(ndx::mx  * Ny, Ny).real())
         ;
    Matrix2r cmat;
    cmat << bulkcoeff.dot(cdata.col(0).segment(ndx::rho * Ny, Ny).real())
         ,  bulkcoeff.dot(cdata.col(1).segment(ndx::rho * Ny, Ny).real())
         ,  bulkcoeff.dot(cdata.col(0).segment(ndx::mx  * Ny, Ny).real())
         ,  bulkcoeff.dot(cdata.col(1).segment(ndx::mx  * Ny, Ny).real())
         ;
    Vector2r cphi;

    // Solve the requested, possibly simultaneous constraint problem
    if (constrain_bulk_rho() && constrain_bulk_rho_u()) {
        jacobiSvd.compute(cmat);              // Fancy decomposition for
        cphi = jacobiSvd.solve(crhs);         // a robust 2x2 constraint solve
    } else if (constrain_bulk_rho()  ) {
        cphi(0) = crhs(0) / cmat(0,0);        // Solve 1x1 system
        cphi(1) = 0;                          // bulk_rho_u not applied
    } else if (constrain_bulk_rho_u()) {
        cphi(0) = 0;                          // bulk_rho not applied
        cphi(1) = crhs(1) / cmat(1,1);        // Solve 1x1 system
    } else {
        cphi.setZero();                       // Neither constraint applied
    }

    // Add scaled constraints to the mean state
    mean += cphi(0)*cdata.col(0) + cphi(1)*cdata.col(1);

    // The implicitly applied integral constraints, as coefficients, must be
    // averaged across each substep to permit accounting for their impact on
    // the Reynolds averaged equations using method_interface::iota similarly to
    //
    //    mean += iota * ((sample / delta_t) - mean).
    //
    // The delta_t accounts for step sizes already implicitly included in cphi.
    //
    // Notice bulk_rho_u-related forcing is NOT scaled by Mach^2 when tracked
    // because our post-processing routines will account for Mach^2 factor.
    cphi                 /= delta_t; // Henceforth includes 1/delta_t scaling!
    const real_t iota     = method.iota(substep_index);
    common.f()           += iota * (
                                ArrayX1r::Constant(Ny, cphi(1))
                              - common.f()
                            );
    common.f_dot_u()     += iota * (
                                cphi(1) * common.u()
                              - common.f_dot_u()
                            );
    common.qb()          += iota * (/* zero */ - common.qb());
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
