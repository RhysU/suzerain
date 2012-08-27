//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
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
// explicit_op.cpp: Fully explicit Navier--Stokes operators
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
    // TODO Speedup possible by not applying to dealiased modes
    // a la HybridIsothermalLinearOperator::applyMassPlusScaledOperator

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
    // TODO Speedup possible by not accumulating dealiased modes
    // a la HybridIsothermalLinearOperator::accumulateMassPlusScaledOperator

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

    // channel_treatment step (3) performs the usual operator solve
    base::invertMassPlusScaledOperator(
            phi, state, delta_t, substep_index, iota);

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
