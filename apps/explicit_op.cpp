//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// explicit_op.hpp: Operators for channel simulation
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "nonlinear.hpp"
#include "explicit_op.hpp"

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>

#pragma warning(disable:383 1572)

namespace channel {

BsplineMassOperator::BsplineMassOperator(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop)
    : suzerain::OperatorBase<real_t>(scenario, grid, dgrid, b, bop),
      massluz(bop)
{
    SUZERAIN_UNUSED(scenario);
    SUZERAIN_UNUSED(grid);
    SUZERAIN_UNUSED(dgrid);
    SUZERAIN_UNUSED(b);

    massluz.factor_mass(bop);
}

void BsplineMassOperator::applyMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const component delta_t,
        const std::size_t substep_index) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Verify required assumptions
    SUZERAIN_ENSURE(state.strides()[1] == 1);
    SUZERAIN_ENSURE(state.shape()[1]   == static_cast<unsigned>(massluz.n()));
    SUZERAIN_ENSURE(suzerain::multi_array::is_contiguous(state));

    // Those assumptions holding, apply operator to each wall-normal pencil.
    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    bop.apply(0, nrhs, 1, state.data(), 1, state.shape()[1]);
}


void BsplineMassOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        suzerain::ContiguousState<4,complex_t> &output,
        const component delta_t,
        const std::size_t substep_index) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    SUZERAIN_ENSURE(output.isIsomorphic(input));

    const suzerain::multi_array::ref<complex_t,4> &x = input;  // Shorthand
    suzerain::ContiguousState<4,complex_t>        &y = output; // Shorthand
    const complex_t c_one = 1;

    // Sidesteps assertions triggered by dereferencing trivial input and output
    if (SUZERAIN_UNLIKELY(0U == x.shape()[1] * x.shape()[2])) return;

    // Loops go from slower to faster indices for ContiguousState<4,complex_t>
    typedef suzerain::ContiguousState<4,complex_t>::index index;
    for (index ix = x.index_bases()[0], iy = y.index_bases()[0];
        ix < static_cast<index>(x.index_bases()[0] + x.shape()[0]);
        ++ix, ++iy) {

        for (index lx = x.index_bases()[3], ly = y.index_bases()[3];
            lx < static_cast<index>(x.index_bases()[3] + x.shape()[3]);
            ++lx, ++ly) {

            bop.accumulate(0, x.shape()[2],
                    c_one,
                    &x[ix][x.index_bases()[1]][x.index_bases()[2]][lx],
                    x.strides()[1], x.strides()[2],
                    beta,
                    &y[iy][y.index_bases()[1]][y.index_bases()[2]][ly],
                    y.strides()[1], y.strides()[2]);
        }
    }
}

void BsplineMassOperator::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const component delta_t,
        const std::size_t substep_index,
        const real_t iota) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);
    SUZERAIN_UNUSED(iota);

    // Verify required assumptions
    SUZERAIN_ENSURE(state.strides()[1] == 1);
    SUZERAIN_ENSURE(state.shape()[1]   == static_cast<unsigned>(massluz.n()));
    SUZERAIN_ENSURE(suzerain::multi_array::is_contiguous(state));

    // Those assumptions holding, invert operator on each wall-normal pencil.
    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    massluz.solve(nrhs, state.data(), 1, state.shape()[1]);
}

BsplineMassOperatorIsothermal::BsplineMassOperatorIsothermal(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common)
    : BsplineMassOperator(scenario, grid, dgrid, b, bop),
      common(common)
{
    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= scenario.Ly;
}

// A helper class for implementing isothermal, no-slip boundary conditions
class IsothermalNoSlipFunctor
{
private:
    const ptrdiff_t field_stride;
    const real_t    inv_gamma_gamma1;

public:
    IsothermalNoSlipFunctor(ptrdiff_t field_stride, real_t gamma)
        : field_stride(field_stride),
          inv_gamma_gamma1(1 / (gamma * (gamma - 1)))
    {}

    void operator()(complex_t &rho) const
    {
        namespace ndx = channel::field::ndx;
        (&rho)[(ndx::rhou - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::rhov - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::rhow - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::rhoe - ndx::rho)*field_stride] = rho*inv_gamma_gamma1;
    }
};

void BsplineMassOperatorIsothermal::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const component delta_t,
        const std::size_t substep_index,
        const real_t iota) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    using Eigen::Map;
    using Eigen::ArrayXc;
    using Eigen::ArrayXr;
    namespace ndx = channel::field::ndx;
    const std::size_t Ny          = state.shape()[1];
    const std::size_t wall_lower  = 0;
    const std::size_t wall_upper  = Ny - 1;
    const bool has_zero_zero_modes = dgrid.has_zero_zero_modes();

    // See channel_treatment writeup for information on the steps below.
    // Steps appear out of order relative to the writeup (TODO fix writeup).

    // Means of the implicit momentum and energy forcing coefficients(!) are
    // maintained across each individual time step for sampling the statistics
    // /bar_f, /bar_f_dot_u, and /bar_qb using OperatorCommonBlock via
    // ILowStorageMethod::iota a la mean += iota*(sample - mean).  Note that
    // the sample must be divided by delta_t to account for step sizes.

    // channel_treatment step (1) done during nonlinear operator application
    // via shared OperatorCommonBlock storage space

    // channel_treatment step (8) sets no-slip conditions
    // on wall collocation points.
    //
    // channel_treatment step (9) sets isothermal conditions on wall
    // collocation points using e_wall = rho_wall / (gamma * (gamma - 1)).
    //
    // Possible pre-solve since L = 0.
    {
        // Prepare a state view of density locations at lower and upper walls
        using boost::multi_array_types::index_range;
        suzerain::multi_array::ref<complex_t,4>::array_view<3>::type view
                = state[boost::indices[ndx::rho]
                                      [index_range(wall_lower,
                                                   wall_upper + 1,
                                                   wall_upper - wall_lower)]
                                      [index_range()]
                                      [index_range()]];

        // Prepare functor setting pointwise BCs given density locations
        const IsothermalNoSlipFunctor bc_functor(
                state.strides()[0], scenario.gamma);

        // Apply the functor to all wall-only density locations
        suzerain::multi_array::for_each(view, bc_functor);
    }

    // Integral constraints enabled only when parameters are non-inf, non-NaN.
    // Allow disabling these to meet manufactured solution verification needs.
    const bool constrain_bulk_rhou
            = (boost::math::isnormal)(scenario.bulk_rhou);
    const bool constrain_bulk_rho
            = (boost::math::isnormal)(scenario.bulk_rho);

    // channel_treatment step (2) loads ones and mean streamwise velocity at
    // collocation points into the imaginary part of the constant (zero zero)
    // mode coefficients.  No forcing occurs at the lower or upper walls.
    if (constrain_bulk_rhou && has_zero_zero_modes) {

        Map<ArrayXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        mean_rhou.imag()[wall_lower] = 0;
        mean_rhou.imag().segment(1, Ny-2).setOnes();
        mean_rhou.imag()[wall_upper] = 0;

        Map<ArrayXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.imag()[wall_lower]      = 0;
        mean_rhoe.imag().segment(1, Ny-2) = common.u().segment(1, Ny-2);
        mean_rhoe.imag()[wall_upper]      = 0;
    }

    // B-spline numerics are not conservative when the mean density field has
    // any asymmetry.  Combat the bulk density's very, very small tendency to
    // drift by adding an integral constraint that the bulk density stay fixed
    // at a target value.  Approach follows that of the bulk momentum forcing.
    //
    // Forcing density at the walls as that no longer permits simultaneously
    // solving the zero-zero modes and the forcing constraint.  (at least not
    // without less-straightforward boundary treatments).  The complexity
    // arises from rhoe = rho / (gamma*(gamma-1)) not holding when rho at the
    // zero-zero wall modes has a non-zero imaginary part.
    if (constrain_bulk_rho && has_zero_zero_modes) {
        Map<ArrayXc> mean_rho(state[ndx::rho].origin(), Ny);
        mean_rho.imag()[wall_lower] = 0;
        mean_rho.imag().segment(1, Ny-2).setOnes();
        mean_rho.imag()[wall_upper] = 0;
    }

    // channel_treatment step (3) performs the usual operator solve
    base::invertMassPlusScaledOperator(
            phi, state, delta_t, substep_index, iota);

    if (constrain_bulk_rhou && has_zero_zero_modes) {

        // channel_treatment steps (4), (5), and (6) determine and
        // apply the appropriate bulk momentum forcing to achieve
        // a target value.
        Map<ArrayXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rhou.matrix());
        const real_t varphi = (scenario.bulk_rhou - bulk.real()) / bulk.imag();
        mean_rhou.real() += varphi * mean_rhou.imag();
        common.f() += iota*((varphi/delta_t)*mean_rhou.imag() - common.f());
        mean_rhou.imag().setZero();

        // channel_treatment step (7) accounts for the momentum forcing
        // within the total energy equation including the Mach squared
        // factor arising from the nondimensionalization choices.
        Map<ArrayXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.real() += (varphi*scenario.Ma*scenario.Ma) * mean_rhoe.imag();
        common.f_dot_u() += iota*(   (varphi/delta_t)*mean_rhoe.imag()
                                   - common.f_dot_u());
        mean_rhoe.imag().setZero();

        // channel_treatment steps (8) and (9) already performed above
    }

    // Complete the bulk density target forcing
    if (constrain_bulk_rho && has_zero_zero_modes) {

        Map<ArrayXc> mean_rho(state[ndx::rho].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rho.matrix());
        const real_t varphi = (scenario.bulk_rho - bulk.real()) / bulk.imag();
        mean_rho.real() += varphi * mean_rho.imag();
        mean_rho.imag().setZero();

        // Note that only the continuity equation is forced which neglects
        // contributions to the momentum and energy equations.  Neglecting
        // these terms is acceptable as the forcing is very small and vanishes
        // when the mean field is symmetric in the wall-normal direction.

        // Statistics for this quantity are not tracked.

    }

    // No volumetric energy forcing in performed current substep
    common.qb() += iota*(/* (Zero/delta_t) */ - common.qb());

    // State leaves method as coefficients in X, Y, and Z directions
}

std::vector<real_t> NonlinearOperator::applyOperator(
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{
    // Dispatch to implementation paying nothing for substep-related ifs
    if (substep_index == 0) {
        return channel::applyNonlinearOperator<true,  channel::linearize::none>
            (*this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    } else {
        return channel::applyNonlinearOperator<false, channel::linearize::none>
            (*this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    }
}

} // end namespace channel
