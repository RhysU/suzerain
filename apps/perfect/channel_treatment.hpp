//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
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
// channel_treatment.hpp: mixin providing channel problem forcing
// $Id$

#ifndef CHANNEL_TREATMENT_HPP
#define CHANNEL_TREATMENT_HPP

#include "nonlinear.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>

#include "perfect.hpp"

namespace suzerain { namespace perfect {

/**
 * A mixin providing channel problem treatment atop any
 * <tt>
 *    public suzerain::timestepper::INonlinearOperator<
 *          suzerain::ContiguousState<4,complex_t>
 *    >
 * </tt>.  During \ref invertMassPlusScaledOperator implicit momentum forcing
 * is applied following the section of <tt>writeups/channel_treatment.tex</tt>
 * titled "Enforcing a target bulk momentum via the linear operator" and using
 * information from OperatorCommonBlock::u() via instance provided at
 * construction time.
 *
 * Means of the implicit momentum and energy forcing coefficients are also
 * maintained across each individual time step for sampling the statistics
 * /bar_f, /bar_f_dot_u, and /bar_qb using OperatorCommonBlock.
 */
template< typename BaseClass >
class ChannelTreatment : public BaseClass
{
public:

    /**
     * Constructor delegating to BaseClass.
     *
     * BaseClass must make its constructor arguments available as member
     * variables under the same name as those found in this constructor.
     */
    ChannelTreatment(
            const ScenarioDefinition &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common);

    /**
     * Force the channel problem delegating to BaseClass when appropriate.
     * The BaseClass must perform the following steps.
     * <ul>
     * <li>
     *     channel_treatment step (3) performs the operator solve which for
     *     the implicit treatment must be combined with boundary conditions
     * </li><li>
     *     channel_treatment step (8) sets no-slip conditions
     *     on wall collocation points.
     * </li><li>
     * </li>
     *     channel_treatment step (9) sets isothermal conditions at walls
     *     using rho_wall = e_wall * gamma * (gamma - 1).
     * </ul>
     * The remainder of the channel treatment is done by this mixin.
     */
    virtual void invertMassPlusScaledOperator(
            const complex_t &phi,
            suzerain::multi_array::ref<complex_t,4> &state,
            const suzerain::timestepper::lowstorage::IMethod<complex_t> &method,
            const real_t delta_t,
            const std::size_t substep_index) const;

private:

    /** Should bulk streamwise momentum constraint be enforced? */
    bool constrain_bulk_rhou() const
    {
        // Yes when non-inf, non-NaN.
        return (boost::math::isnormal)(this->scenario.bulk_rhou);
    }

    /** Should bulk density constraint be enforced? */
    bool constrain_bulk_rho() const
    {
        // Yes when non-inf, non-NaN.
        return (boost::math::isnormal)(this->scenario.bulk_rho);
    }

    /** Precomputed, real-valued mass matrix factorization */
    boost::scoped_ptr<suzerain::bsplineop_lu> masslu;

    /** Precomputed integration coefficients */
    VectorXr bulkcoeff;

    /** Coefficients for function supported at non-wall points */
    ArrayXr interior;

    /** Dot product of bulkcoeff against interior */
    real_t bulkcoeff_dot_interior;
};

template< typename BaseClass >
ChannelTreatment<BaseClass>::ChannelTreatment(
            const ScenarioDefinition &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop,
            OperatorCommonBlock &common)
    : BaseClass(scenario, grid, dgrid, b, bop, common)
{
    // Precomputed results only necessary on rank with zero-zero modes
    if (!dgrid.has_zero_zero_modes()) return;

    masslu.reset(new suzerain::bsplineop_lu(bop));
    masslu->factor_mass(bop);

    // channel_treatment step (2) loads ones and local mean streamwise velocity
    // into non-wall collocation points.  We need an 'interior-only' forcing
    // profile.  The most pressure-gradient-like profile (as measured by
    // reproducing a constant \partial_y \tau_{x,y}) has zeros at the wall
    // collocation point and ones everywhere else.  Consequently, it has large
    // derivatives near the walls and will introduce small magnitude
    // oscillations near the channel center.  See Redmine ticket #2568 for
    // further background.
    interior.setOnes(b.n());                     // Set as collocation values
    interior.head<1>().setZero();                // Zero at lower wall
    interior.tail<1>().setZero();                // Zero at upper wall
    masslu->solve(1, interior.data(), 1, b.n()); // Convert to coefficients

    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid.L.y();
    bulkcoeff_dot_interior = bulkcoeff.dot(interior.matrix());

    // Release resources if not necessary beyond initialization
    if (!constrain_bulk_rhou()) masslu.reset();
}

template< typename BaseClass >
void ChannelTreatment<BaseClass>::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const suzerain::timestepper::lowstorage::IMethod<complex_t> &method,
        const real_t delta_t,
        const std::size_t substep_index) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    typedef Map<ArrayXr,0,Eigen::InnerStride<2> > ModesRealPart;
    using suzerain::inorder::wavenumber;
    namespace field = support::field;
    namespace ndx   = field::ndx;
    OperatorCommonBlock &common = this->common;
    const ScenarioDefinition &scenario = this->scenario;
    const suzerain::pencil_grid &dgrid = this->dgrid;
    const int Ny = dgrid.global_wave_extent.y();

    // Sidesteps assertions when local rank contains no wavespace information
    if (SUZERAIN_UNLIKELY(0U == state.shape()[1])) return;

    // Incoming state has wall-normal pencils of interleaved state scalars?
    SUZERAIN_ENSURE(state.shape()  [1] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.strides()[1] ==             1);
    SUZERAIN_ENSURE(state.strides()[0] == (unsigned) Ny);
    SUZERAIN_ENSURE(state.shape()  [0] ==  field::count);

    // Means of the implicit momentum and energy forcing coefficients(!) are
    // maintained across each individual time step for sampling the statistics
    // /bar_f, /bar_f_dot_u, and /bar_qb using OperatorCommonBlock via
    // IMethod::iota_alpha a la mean += iota_alpha*(sample - mean).  Note that
    // the sample must be divided by alpha*delta_t to account for step sizes.

    // See channel_treatment writeup for general information on the steps below.
    //
    // To handle either the fully explicit or the hybrid implicit/explicit
    // case, we must account for a non-trivial linear operator.  This implies
    // that the invertMassPlusScaledOperator() call below cannot be tricked
    // into solving our integral constraints while it solves the real-valued
    // zero-zero mode problem.  We must accommodate mass matrix solves.  As a
    // result, steps may appear out of order or differ compared to the writeup.
    //
    // B-spline numerics are not conservative when the mean density field has
    // any asymmetry.  Combat the bulk density's very, very small tendency to
    // drift by adding an integral constraint that the bulk density stay fixed
    // at a target value.  Approach follows that of the bulk momentum forcing.

    // channel_treatment step (1) done during nonlinear operator application
    // via shared OperatorCommonBlock storage space.

    // The BaseClass is responsible for all of the following steps:
    //
    // channel_treatment step (3) performs the operator solve which for the
    // implicit treatment must be combined with boundary conditions
    //
    // channel_treatment step (8) sets no-slip conditions
    // on wall collocation points.
    //
    // channel_treatment step (9) sets isothermal conditions at walls
    // using rho_wall = e_wall * gamma * (gamma - 1).
    BaseClass::invertMassPlusScaledOperator(
            phi, state, method, delta_t, substep_index);

    // Only the rank with zero-zero modes proceeds!
    if (!dgrid.has_zero_zero_modes()) return;

    // If necessary, constrain the bulk streamwise momentum.
    const real_t iota_alpha = method.iota_alpha(substep_index);
    const real_t alpha_dt   = method.alpha(substep_index)*delta_t;
    if (constrain_bulk_rhou()) {

        // channel_treatment step (3) was already performed for state.
        // Perform mass matrix solve for forcing work contribution.
        ArrayXr rhs = common.u();
        rhs[0]    = 0;
        rhs[Ny-1] = 0;
        masslu->solve(1, rhs.data(), 1, rhs.size());

        // channel_treatment steps (4), (5), and (6) determine and apply
        // the appropriate bulk momentum forcing to achieve a target value.
        ModesRealPart mean_rhou((real_t *)state[ndx::mx].origin(), Ny);
        const real_t observed = bulkcoeff.dot((mean_rhou.matrix()));
        const real_t varphi = (scenario.bulk_rhou - observed)
                            / bulkcoeff_dot_interior;
        mean_rhou += varphi*interior;

        // channel_treatment step (7) accounts for the momentum forcing
        // within the total energy equation including the Mach squared
        // factor arising from the nondimensionalization choices.
        ModesRealPart mean_rhoe((real_t *)state[ndx::e].origin(), Ny);
        mean_rhoe += (varphi*scenario.Ma*scenario.Ma)*rhs;

        // Track the forcing magnitude within statistics.
        common.f()       += iota_alpha
                          * ((varphi/alpha_dt)*interior - common.f()      );
        common.f_dot_u() += iota_alpha
                          * ((varphi/alpha_dt)*rhs      - common.f_dot_u());
    } else {

        // Track the lack of forcing within statistics.
        common.f()       += iota_alpha*(/* zero */ - common.f()      );
        common.f_dot_u() += iota_alpha*(/* zero */ - common.f_dot_u());

    }

    // If necessary, constraint the bulk density.
    if (constrain_bulk_rho()) {

        // We force only the continuity equation which neglects
        // associated contributions to the momentum and energy
        // equations.  The mass addition is of order the discrete
        // conservation error and vanishes when the mean density
        // field is symmetric in the wall-normal direction.
        ModesRealPart mean_rho((real_t *)state[ndx::rho].origin(), Ny);
        const real_t observed = bulkcoeff.dot((mean_rho.matrix()));
        const real_t varphi = (scenario.bulk_rho - observed)
                            / bulkcoeff_dot_interior;
        mean_rho += varphi*interior;

        // Statistics for the bulk density forcing are not tracked
    }

    // No volumetric energy forcing in performed current substep
    common.qb() += iota_alpha*(/* zero */ - common.qb());

    // State leaves method as coefficients in X, Y, and Z directions
}

} /* namespace perfect */ } /* namespace suzerain */

#endif /* CHANNEL_TREATMENT_HPP */
