//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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

#include "explicit_op.hpp"

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/state.hpp>

#include "nonlinear.hpp"

#pragma warning(disable:383 1572)

namespace suzerain { namespace perfect {

BsplineMassOperator::BsplineMassOperator(
        const problem::GridDefinition &grid,
        const pencil_grid &dgrid,
        bspline &b,
        const bsplineop &bop)
    : OperatorBase(grid, dgrid, b, bop),
      massluz(bop)
{
    SUZERAIN_UNUSED(grid);
    SUZERAIN_UNUSED(dgrid);
    SUZERAIN_UNUSED(b);

    massluz.factor_mass(bop);
}

void BsplineMassOperator::applyMassPlusScaledOperator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index) const
{
    // TODO Speedup possible by not applying to dealiased modes
    // a la HybridIsothermalLinearOperator::applyMassPlusScaledOperator

    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Verify required assumptions
    SUZERAIN_ENSURE(state.strides()[1] == 1);
    SUZERAIN_ENSURE(state.shape()[1]   == static_cast<unsigned>(massluz.n()));
    SUZERAIN_ENSURE(multi_array::is_contiguous(state));

    // Those assumptions holding, apply operator to each wall-normal pencil.
    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    bop.apply(0, nrhs, 1, state.data(), 1, state.shape()[1]);
}


void BsplineMassOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const multi_array::ref<complex_t,4> &input,
        const complex_t &beta,
        ContiguousState<4,complex_t> &output,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index) const
{
    // TODO Speedup possible by not accumulating dealiased modes
    // a la HybridIsothermalLinearOperator::accumulateMassPlusScaledOperator

    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    SUZERAIN_ENSURE(output.isIsomorphic(input));

    const multi_array::ref<complex_t,4> &x = input;  // Shorthand
    ContiguousState<4,complex_t>        &y = output; // Shorthand
    const complex_t c_one = 1;

    // Sidesteps assertions triggered by dereferencing trivial input and output
    if (SUZERAIN_UNLIKELY(0U == x.shape()[1] * x.shape()[2])) return;

    // Loops go from slower to faster indices for ContiguousState<4,complex_t>
    typedef ContiguousState<4,complex_t>::index index;
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
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> * const ic0) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(method);
    SUZERAIN_UNUSED(delta_t);
    SUZERAIN_UNUSED(substep_index);

    // Verify required assumptions
    SUZERAIN_ENSURE(state.strides()[1] == 1);
    SUZERAIN_ENSURE(state.shape()[1]   == static_cast<unsigned>(massluz.n()));
    SUZERAIN_ENSURE(multi_array::is_contiguous(state));

    // Those assumptions holding, invert operator on each wall-normal pencil
    massluz.solve(state.shape()[0]*state.shape()[2]*state.shape()[3],
                  state.data(), 1, state.shape()[1]);

    if (ic0) {

        // Likewise, verify required assumptions
        SUZERAIN_ENSURE(ic0->strides()[1] == 1);
        SUZERAIN_ENSURE(ic0->shape()[1]   == static_cast<unsigned>(massluz.n()));
        SUZERAIN_ENSURE(multi_array::is_contiguous(*ic0));

        // Likewise, invert operator on each wall-normal pencil
        massluz.solve(ic0->shape()[0]*ic0->shape()[2]*ic0->shape()[3],
                      ic0->data(), 1, ic0->shape()[1]);

    }
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
        namespace ndx = support::field::ndx;
        (&rho)[(ndx::mx - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::my - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::mz - ndx::rho)*field_stride] = 0;
        (&rho)[(ndx::e  - ndx::rho)*field_stride] = rho*inv_gamma_gamma1;
    }
};

void BsplineMassOperatorIsothermal::invertMassPlusScaledOperator(
        const complex_t &phi,
        multi_array::ref<complex_t,4> &state,
        const timestepper::lowstorage::IMethod<complex_t> &method,
        const component delta_t,
        const std::size_t substep_index,
        multi_array::ref<complex_t,4> * const ic0) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    namespace ndx = support::field::ndx;
    using boost::indices;
    typedef boost::multi_array_types::index_range range;
    typedef multi_array::ref<complex_t,4>::array_view<3>::type view_type;

    // Indexes only the first and last collocation point
    const std::size_t Ny         = state.shape()[1];
    const std::size_t wall_lower = 0;
    const std::size_t wall_upper = Ny - 1;
    range walls(wall_lower, wall_upper + 1, wall_upper - wall_lower);

    // channel_treatment step (8) sets no-slip conditions
    // on wall collocation points.
    //
    // channel_treatment step (9) sets isothermal conditions on wall
    // collocation points using e_wall = rho_wall / (gamma * (gamma - 1)).
    //
    // Possible pre-solve since L = 0.
    {
        // Prepare a state view of density locations at lower and upper walls
        view_type view = state[indices[ndx::rho][walls][range()][range()]];

        // Prepare functor setting pointwise BCs given density locations
        const IsothermalNoSlipFunctor bc_functor(
                state.strides()[0], scenario.gamma);

        // Apply the functor to all wall-only density locations
        multi_array::for_each(view, bc_functor);

        // Apply boundary conditions to any requested constraint problems
        if (ic0) {
            view = (*ic0)[indices[ndx::rho][walls][range()][range()]];
            SUZERAIN_ENSURE(state.strides()[0] == ic0->strides()[0]); // NB!
            multi_array::for_each(view, bc_functor);
        }
    }

    // channel_treatment step (3) performs the usual operator solve
    base::invertMassPlusScaledOperator(
            phi, state, method, delta_t, substep_index, ic0);

    // State leaves method as coefficients in X, Y, and Z directions
}

std::vector<real_t> NonlinearOperator::applyOperator(
            const real_t time,
            ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const
{
    // Dispatch to implementation paying nothing for substep-related ifs
    if (substep_index == 0) {
        return applyNonlinearOperator<true,  linearize::none>
            (this->scenario.alpha,
             this->scenario.beta,
             this->scenario.gamma,
             this->scenario.Ma,
             this->scenario.Pr,
             this->scenario.Re,
             *this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    } else {
        return applyNonlinearOperator<false, linearize::none>
            (this->scenario.alpha,
             this->scenario.beta,
             this->scenario.gamma,
             this->scenario.Ma,
             this->scenario.Pr,
             this->scenario.Re,
             *this, common, msoln, time, swave, evmaxmag_real, evmaxmag_imag);
    }
}

} /* namespace perfect */ } /* namespace suzerain */
