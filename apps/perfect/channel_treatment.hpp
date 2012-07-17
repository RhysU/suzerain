//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
#include <suzerain/scenario_definition.hpp>

#include "../precision.hpp"

namespace channel {

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
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
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
            const real_t delta_t,
            const std::size_t substep_index,
            const real_t iota) const;

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
    Eigen::VectorXr bulkcoeff;

    // TODO Not strictly necessary beyond constructor.  Save some RAM.
    /** Coefficients for function supported at non-wall points */
    Eigen::ArrayXr interior;

    /** Dot product of bulkcoeff against interior */
    real_t bulkcoeff_dot_interior;
};

template< typename BaseClass >
ChannelTreatment<BaseClass>::ChannelTreatment(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
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

    // channel_treatment step (2) loads ones and local mean streamwise
    // velocity into non-wall collocation points.  Only non-wall points are
    // used to avoid upsetting Dirichlet boundary conditions.  More
    // thought will be necessary for Neumann or Robin conditions.
    interior.setZero(b.n());
    interior.segment(1, b.n()-2).setOnes();
    masslu->solve(1, interior.data(), 1, b.n());

    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= scenario.Ly;
    bulkcoeff_dot_interior = bulkcoeff.dot(interior.matrix());

    // Release resources if not necessary beyond initialization
    if (!constrain_bulk_rhou()) masslu.reset();
}

template< typename BaseClass >
void ChannelTreatment<BaseClass>::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::multi_array::ref<complex_t,4> &state,
        const real_t delta_t,
        const std::size_t substep_index,
        const real_t iota) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    typedef Eigen::Map<Eigen::ArrayXr,0,Eigen::InnerStride<2> > ModesRealPart;
    using suzerain::inorder::wavenumber;
    namespace field = channel::field;
    namespace ndx   = field::ndx;
    OperatorCommonBlock &common
            = this->common;
    const suzerain::problem::ScenarioDefinition<real_t> &scenario
            = this->scenario;
    const suzerain::pencil_grid &dgrid
            = this->dgrid;
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
    // ILowStorageMethod::iota a la mean += iota*(sample - mean).  Note that
    // the sample must be divided by delta_t to account for step sizes.

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
            phi, state, delta_t, substep_index, iota);

    // Only the rank with zero-zero modes proceeds!
    if (!dgrid.has_zero_zero_modes()) return;

    // If necessary, constrain the bulk streamwise momentum.
    if (constrain_bulk_rhou()) {

        // channel_treatment step (3) was already performed for state.
        // Perform mass matrix solve for forcing work contribution.
        Eigen::ArrayXr rhs = common.u();
        rhs[0]    = 0;
        rhs[Ny-1] = 0;
        masslu->solve(1, rhs.data(), 1, rhs.size());

        // channel_treatment steps (4), (5), and (6) determine and apply
        // the appropriate bulk momentum forcing to achieve a target value.
        ModesRealPart mean_rhou((real_t *)state[ndx::rhou].origin(), Ny);
        const real_t observed = bulkcoeff.dot((mean_rhou.matrix()));
        const real_t varphi = (scenario.bulk_rhou - observed)
                            / bulkcoeff_dot_interior;
        mean_rhou += varphi*interior;

        // channel_treatment step (7) accounts for the momentum forcing
        // within the total energy equation including the Mach squared
        // factor arising from the nondimensionalization choices.
        ModesRealPart mean_rhoe((real_t *)state[ndx::rhoe].origin(), Ny);
        mean_rhoe += (varphi*scenario.Ma*scenario.Ma)*rhs;

        // Track the forcing magnitude within statistics.
        common.f()       += iota*((varphi/delta_t)*interior - common.f());
        common.f_dot_u() += iota*((varphi/delta_t)*rhs - common.f_dot_u());
    } else {

        // Track the lack of forcing within statistics.
        common.f()       += iota*(/* (zero/delta_t) */ - common.f());
        common.f_dot_u() += iota*(/* (zero/delta_t) */ - common.f_dot_u());

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
    common.qb() += iota*(/* (zero/delta_t) */ - common.qb());

    // State leaves method as coefficients in X, Y, and Z directions
}

} // namespace channel

#endif /* CHANNEL_TREATMENT_HPP */
