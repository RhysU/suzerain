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
// explicit_op.hpp: Operators for channel_explicit
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/error.h>
#include <suzerain/math.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/state.hpp>
#include <suzerain/utility.hpp>

#include "explicit_op.hpp"

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
    if (constrain_bulk_rho && has_zero_zero_modes) {
        Map<ArrayXc> mean_rho(state[ndx::rho].origin(), Ny);
        mean_rho.imag().setConstant(1);
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

NonlinearOperator::NonlinearOperator(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        OperatorCommonBlock &common,
        const boost::shared_ptr<
            const channel::manufactured_solution>& msoln)
    : suzerain::OperatorBase<real_t>(scenario, grid, dgrid, b, bop),
      common(common),
      msoln(msoln)
{
    // NOP
}

// Notice that method is templated on being the zeroth substep to avoid runtime
// costs of "if (zeroth_substep) { ... }" logic and potentially aid the
// optimizer.  Penalty is code bloat from having two nearly identical versions
// of the method in the final binary and some extra typename keywords below.
template< bool zeroth_substep >
std::vector<real_t> NonlinearOperator::applyOperator(
    const real_t time,
    suzerain::ContiguousState<4,complex_t> &swave,
    const real_t evmaxmag_real,
    const real_t evmaxmag_imag) const
{
    namespace ndx = channel::field::ndx;
    using Eigen::Vector3r;
    using Eigen::Matrix3r;

    // State enters method as coefficients in X, Y, and Z directions

    // We need auxiliary scalar-field storage.  Prepare logical indices using a
    // struct for scoping (e.g. aux::rho_y).  Ordering will match usage below.
    struct aux { enum {
        rho_y, rho_yy, rho_x, rho_xx, rho_xz, rho_z, rho_zz, rho_xy, rho_yz,
        mx_y,  mx_yy,  mx_x,  mx_xx,  mx_xz,  mx_z,  mx_zz,  mx_xy,  mx_yz,
        my_y,  my_yy,  my_x,  my_xx,  my_xz,  my_z,  my_zz,  my_xy,  my_yz,
        mz_y,  mz_yy,  mz_x,  mz_xx,  mz_xz,  mz_z,  mz_zz,  mz_xy,  mz_yz,
        e_y, div_grad_e, e_x, e_z,
        count // Sentry
    }; };

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    boost::scoped_ptr<state_type> _auxw_ptr(
            allocate_padded_state<state_type>(aux::count, dgrid)); // RAII
    state_type &auxw = *_auxw_ptr;                                 // Shorthand

    // Sanity check incoming swave's and auxw's shape and contiguity
    assert(swave.shape()[0] == channel::field::count);
    assert(swave.shape()[1] == (unsigned) dgrid.local_wave_extent.y());
    assert(swave.shape()[2] == (unsigned) dgrid.local_wave_extent.x());
    assert(swave.shape()[3] == (unsigned) dgrid.local_wave_extent.z());
    assert((unsigned) swave.strides()[1] == 1u);
    assert((unsigned) swave.strides()[2] == swave.shape()[1]);
    assert((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    assert(std::equal(swave.shape() + 1, swave.shape() + 4,
                      auxw.shape() + 1));
    assert(std::equal(swave.strides() + 1, swave.strides() + 4,
                      auxw.strides() + 1));

    // Prepare common-block-like storage used to pass details from N to L.
    // Zeroing is done carefully as accumulated means must survive
    // from substep to substep while instantaneous profiles must not.
    common.storage.resize(/* Ny */ swave.shape()[1], Eigen::NoChange);
    common.u().setZero();

    // Maintain stable time step values to return to the caller
    boost::array<real_t, 2> delta_t_candidates = {{
            std::numeric_limits<real_t>::max(),
            std::numeric_limits<real_t>::max()
    }};
    real_t &convective_delta_t = delta_t_candidates[0];
    real_t &diffusive_delta_t  = delta_t_candidates[1];

    // Compute Y derivatives of density at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    diffwave_apply(0, 0, 1., swave, ndx::rho);
    bop_accumulate(1,    1., swave, ndx::rho, 0., auxw, aux::rho_y);
    bop_accumulate(2,    1., swave, ndx::rho, 0., auxw, aux::rho_yy);
    bop_apply     (0,    1., swave, ndx::rho);

    // Compute X- and Z- derivatives of density at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    diffwave_accumulate(1, 0, 1., swave, ndx::rho,   0., auxw, aux::rho_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rho,   0., auxw, aux::rho_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rho,   0., auxw, aux::rho_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rho,   0., auxw, aux::rho_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rho,   0., auxw, aux::rho_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::rho_y, 0., auxw, aux::rho_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::rho_y, 0., auxw, aux::rho_yz);

    // Compute Y derivatives of X momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    diffwave_apply(0, 0, 1., swave, ndx::rhou);
    bop_accumulate(1,    1., swave, ndx::rhou, 0., auxw, aux::mx_y);
    bop_accumulate(2,    1., swave, ndx::rhou, 0., auxw, aux::mx_yy);
    bop_apply     (0,    1., swave, ndx::rhou);

    // Compute X- and Z- derivatives of X momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    diffwave_accumulate(1, 0, 1., swave, ndx::rhou,  0., auxw, aux::mx_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhou,  0., auxw, aux::mx_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rhou,  0., auxw, aux::mx_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhou,  0., auxw, aux::mx_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhou,  0., auxw, aux::mx_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::mx_y,  0., auxw, aux::mx_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::mx_y,  0., auxw, aux::mx_yz);

    // Compute Y derivatives of Y momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    diffwave_apply(0, 0, 1., swave, ndx::rhov);
    bop_accumulate(1,    1., swave, ndx::rhov, 0., auxw, aux::my_y);
    bop_accumulate(2,    1., swave, ndx::rhov, 0., auxw, aux::my_yy);
    bop_apply     (0,    1., swave, ndx::rhov);

    // Compute X- and Z- derivatives of Y momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    diffwave_accumulate(1, 0, 1., swave, ndx::rhov,  0., auxw, aux::my_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhov,  0., auxw, aux::my_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rhov,  0., auxw, aux::my_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhov,  0., auxw, aux::my_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhov,  0., auxw, aux::my_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::my_y,  0., auxw, aux::my_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::my_y,  0., auxw, aux::my_yz);

    // Compute Y derivatives of Z momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    diffwave_apply(0, 0, 1., swave, ndx::rhow);
    bop_accumulate(1,    1., swave, ndx::rhow, 0., auxw, aux::mz_y);
    bop_accumulate(2,    1., swave, ndx::rhow, 0., auxw, aux::mz_yy);
    bop_apply     (0,    1., swave, ndx::rhow);

    // Compute X- and Z- derivatives of Z momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    diffwave_accumulate(1, 0, 1., swave, ndx::rhow,  0., auxw, aux::mz_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhow,  0., auxw, aux::mz_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rhow,  0., auxw, aux::mz_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhow,  0., auxw, aux::mz_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhow,  0., auxw, aux::mz_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::mz_y,  0., auxw, aux::mz_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::mz_y,  0., auxw, aux::mz_yz);

    // Compute Y derivatives of total energy at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    diffwave_apply(0, 0, 1., swave, ndx::rhoe);
    bop_accumulate(1,    1., swave, ndx::rhoe, 0., auxw, aux::e_y);
    bop_accumulate(2,    1., swave, ndx::rhoe, 0., auxw, aux::div_grad_e);
    bop_apply     (0,    1., swave, ndx::rhoe);

    // Compute X- and Z- derivatives of total energy at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    diffwave_accumulate(1, 0, 1., swave, ndx::rhoe, 0., auxw, aux::e_x       );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhoe, 1., auxw, aux::div_grad_e);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhoe, 0., auxw, aux::e_z       );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhoe, 1., auxw, aux::div_grad_e);

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    typename channel::physical_view<aux::count>::type auxp
        = channel::physical_view<aux::count>::create(dgrid, auxw);
    typename channel::physical_view<channel::field::count>::type sphys
        = channel::physical_view<channel::field::count>::create(dgrid, swave);
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (std::size_t i = 0; i < aux::count; ++i) {
        dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }

    // Retrieve constants and compute derived constants
    const real_t alpha            = scenario.alpha;
    const real_t beta             = scenario.beta;
    const real_t gamma            = scenario.gamma;
    const real_t Ma               = scenario.Ma;
    const real_t Pr               = scenario.Pr;
    const real_t Re               = scenario.Re;
    const real_t inv_Re           = 1 / Re;
    const real_t inv_Ma2          = 1 / (Ma * Ma);
    const real_t Ma2_over_Re      = (Ma * Ma) / Re;
    const real_t inv_Re_Pr_gamma1 = 1 / (Re * Pr * (gamma - 1));

    // Type of Boost.Accumulator to use for summation processes.
    // Kahan summation preferred when available as incremental cost is small
    // and we will add many small numbers to a large magnitude sum.
    typedef boost::accumulators::accumulator_set<
                real_t,
#if BOOST_VERSION >= 104700
                boost::accumulators::stats<boost::accumulators::tag::sum_kahan>
#else
                boost::accumulators::stats<boost::accumulators::tag::sum>
#endif
            > summing_accumulator_type;

    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary. Three traversals occur:
    //
    // (1) Computing mean velocity OR reference quantities and mean velocity
    //     (depending on which substep is being performed).
    // (2) Computing the nonlinear equation right hand sides.
    // (3) Computing any manufactured solution forcing (when enabled).
    //
    // The traversal pattern is best embodied by the third pass.  The first and
    // second passes use slightly more compact structure as only y(j) must be
    // known within them.  Study (3) and convince yourself that (1) and (2) are
    // equivalent but lack information on x(i) and z(k).

    // (1) Computing mean velocity OR reference quantities and mean velocity
    //     (depending on which substep is being performed).
    if (!zeroth_substep) { // MPI_Reduce mean velocity only

        // Sum streamwise velocities as a function of y(j) into common.u()
        size_t offset = 0;
        for (int j = dgrid.local_physical_start.y();
            j < dgrid.local_physical_end.y();
            ++j) {

            summing_accumulator_type ux;

            const size_t last_zxoffset = offset
                                       + dgrid.local_physical_extent.z()
                                       * dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {
                ux(sphys(ndx::rhou, offset)/sphys(ndx::rho, offset));
            } // end X // end Z

            // Store sum into common block in preparation for MPI Reduce
            common.u()[j] = boost::accumulators::sum(ux);

        } // end Y

        // Reduce and scale common.u() sums to obtain mean on zero-zero rank
        // Only zero-zero rank needs the information so Reduce is sufficient
        if (dgrid.has_zero_zero_modes()) {
            SUZERAIN_MPICHKR(MPI_Reduce(MPI_IN_PLACE, common.u().data(),
                    common.u().size(), suzerain::mpi::datatype<real_t>::value,
                    MPI_SUM, dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
            common.u() /= (   dgrid.global_physical_extent.x()
                            * dgrid.global_physical_extent.z());
        } else {
            Eigen::ArrayXr tmp;
            tmp.resizeLike(common.u());
            tmp.setZero();
            SUZERAIN_MPICHKR(MPI_Reduce(common.u().data(), tmp.data(),
                    common.u().size(), suzerain::mpi::datatype<real_t>::value,
                    MPI_SUM, dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
        }

    } else {               // MPI_Allreduce references and mean velocity

        // Sum reference quantities as a function of y(j) into common.ref_*
        size_t offset = 0;
        for (int j = dgrid.local_physical_start.y();
            j < dgrid.local_physical_end.y();
            ++j) {

            // See writeups/derivation.tex or rholut_imexop.h for definitions
            summing_accumulator_type ref_nu, ref_ux, ref_uy, ref_uz,
                                     ref_nuux, ref_nuuy, ref_nuuz,
                                     ref_m_gradrho, ref_ex_gradrho,
                                     ref_ey_gradrho, ref_ez_gradrho,
                                     ref_e_divm, ref_e_deltarho;

            const size_t last_zxoffset = offset
                                       + dgrid.local_physical_extent.z()
                                       * dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {

                // Unpack conserved state
                const real_t   rho(sphys(ndx::rho,  offset));
                const Vector3r m  (sphys(ndx::rhou, offset),
                                   sphys(ndx::rhov, offset),
                                   sphys(ndx::rhow, offset));
                const real_t   e  (sphys(ndx::rhoe, offset));

                // Compute quantities related to the equation of state
                real_t p, T, mu, lambda;
                suzerain::rholut::p_T_mu_lambda(
                    alpha, beta, gamma, Ma, rho, m, e, p, T, mu, lambda);

                // Accumulate reference quantities into running sums...
                const real_t nu = mu / rho;
                ref_nu(nu);

                // ...including simple velocity-related quantities...
                const Vector3r u = suzerain::rholut::u(rho, m);
                ref_ux(u.x());
                ref_uy(u.y());
                ref_uz(u.z());
                ref_nuux(nu*u.x());
                ref_nuuy(nu*u.y());
                ref_nuuz(nu*u.z());

                // ...and other expressions.
                ref_m_gradrho(u.squaredNorm());

                namespace rholut = suzerain::rholut;
                const Vector3r e_gradrho
                        = rholut::explicit_div_e_plus_p_u_refcoeff_grad_rho(
                                gamma, rho, m, e, p);
                ref_ex_gradrho(e_gradrho.x());
                ref_ey_gradrho(e_gradrho.y());
                ref_ez_gradrho(e_gradrho.z());

                ref_e_divm(
                        rholut::explicit_div_e_plus_p_u_refcoeff_div_m(
                            rho, e, p));

                ref_e_deltarho(
                        rholut::explicit_mu_div_grad_T_refcoeff_div_grad_rho(
                            gamma, mu, rho, e, p));

            } // end X // end Z

            // Store sum into common block in preparation for MPI Allreduce
            namespace accumulators = boost::accumulators;
            common.ref_nu        ()[j] = accumulators::sum(ref_nu        );
            common.ref_ux        ()[j] = accumulators::sum(ref_ux        );
            common.ref_uy        ()[j] = accumulators::sum(ref_uy        );
            common.ref_uz        ()[j] = accumulators::sum(ref_uz        );
            common.ref_nuux      ()[j] = accumulators::sum(ref_nuux      );
            common.ref_nuuy      ()[j] = accumulators::sum(ref_nuuy      );
            common.ref_nuuz      ()[j] = accumulators::sum(ref_nuuz      );
            common.ref_m_gradrho ()[j] = accumulators::sum(ref_m_gradrho );
            common.ref_ex_gradrho()[j] = accumulators::sum(ref_ex_gradrho);
            common.ref_ey_gradrho()[j] = accumulators::sum(ref_ey_gradrho);
            common.ref_ez_gradrho()[j] = accumulators::sum(ref_ez_gradrho);
            common.ref_e_divm    ()[j] = accumulators::sum(ref_e_divm    );
            common.ref_e_deltarho()[j] = accumulators::sum(ref_e_deltarho);

        } // end Y

        // Allreduce and scale common.refs() sums to obtain means on all ranks
        // Allreduce mandatory as all ranks need references for linearization
        SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, common.refs().data(),
                common.refs().size(), suzerain::mpi::datatype<real_t>::value,
                MPI_SUM, MPI_COMM_WORLD));
        common.refs() /= (   dgrid.global_physical_extent.x()
                           * dgrid.global_physical_extent.z());

        // Copy redundant mean streamwise velocity information into common.u()
        common.u() = common.ref_ux();
    }

    // (2) Computing the nonlinear equation right hand sides.
    size_t offset = 0;
    for (int j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Wall-normal grid spacing depends on wall-normal location
        const real_t one_over_delta_y_j = one_over_delta_y(j);

        const size_t last_zxoffset = offset
                                   + dgrid.local_physical_extent.z()
                                   * dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // Unpack density-related quantities
            const real_t   rho         ( sphys(ndx::rho,    offset));
            const Vector3r grad_rho    (  auxp(aux::rho_x,  offset),
                                          auxp(aux::rho_y,  offset),
                                          auxp(aux::rho_z,  offset));
            const real_t   div_grad_rho(  auxp(aux::rho_xx, offset)
                                        + auxp(aux::rho_yy, offset)
                                        + auxp(aux::rho_zz, offset));
            const Matrix3r grad_grad_rho;
            const_cast<Matrix3r&>(grad_grad_rho) <<
                                          auxp(aux::rho_xx, offset),
                                          auxp(aux::rho_xy, offset),
                                          auxp(aux::rho_xz, offset),
                                          auxp(aux::rho_xy, offset),
                                          auxp(aux::rho_yy, offset),
                                          auxp(aux::rho_yz, offset),
                                          auxp(aux::rho_xz, offset),
                                          auxp(aux::rho_yz, offset),
                                          auxp(aux::rho_zz, offset);

            // Unpack momentum-related quantities
            const Vector3r m    ( sphys(ndx::rhou, offset),
                                  sphys(ndx::rhov, offset),
                                  sphys(ndx::rhow, offset));
            const real_t   div_m(  auxp(aux::mx_x, offset)
                                 + auxp(aux::my_y, offset)
                                 + auxp(aux::mz_z, offset));
            const Matrix3r grad_m;
            const_cast<Matrix3r&>(grad_m) <<
                                        auxp(aux::mx_x,  offset),
                                        auxp(aux::mx_y,  offset),
                                        auxp(aux::mx_z,  offset),
                                        auxp(aux::my_x,  offset),
                                        auxp(aux::my_y,  offset),
                                        auxp(aux::my_z,  offset),
                                        auxp(aux::mz_x,  offset),
                                        auxp(aux::mz_y,  offset),
                                        auxp(aux::mz_z,  offset);
            const Vector3r div_grad_m(  auxp(aux::mx_xx, offset)
                                      + auxp(aux::mx_yy, offset)
                                      + auxp(aux::mx_zz, offset),
                                        auxp(aux::my_xx, offset)
                                      + auxp(aux::my_yy, offset)
                                      + auxp(aux::my_zz, offset),
                                        auxp(aux::mz_xx, offset)
                                      + auxp(aux::mz_yy, offset)
                                      + auxp(aux::mz_zz, offset));
            const Vector3r grad_div_m(  auxp(aux::mx_xx, offset)
                                      + auxp(aux::my_xy, offset)
                                      + auxp(aux::mz_xz, offset),
                                        auxp(aux::mx_xy, offset)
                                      + auxp(aux::my_yy, offset)
                                      + auxp(aux::mz_yz, offset),
                                        auxp(aux::mx_xz, offset)
                                      + auxp(aux::my_yz, offset)
                                      + auxp(aux::mz_zz, offset));

            // Unpack total energy-related quantities
            const real_t e        (sphys(ndx::rhoe,       offset));
            const Vector3r grad_e ( auxp(aux::e_x,        offset),
                                    auxp(aux::e_y,        offset),
                                    auxp(aux::e_z,        offset));
            const real_t div_grad_e(auxp(aux::div_grad_e, offset));

            // Compute velocity-related quantities
            const Vector3r u          = suzerain::rholut::u(
                                            rho, m);
            const real_t div_u        = suzerain::rholut::div_u(
                                            rho, grad_rho, m, div_m);
            const Matrix3r grad_u     = suzerain::rholut::grad_u(
                                            rho, grad_rho, m, grad_m);
            const Vector3r grad_div_u = suzerain::rholut::grad_div_u(
                                            rho, grad_rho, grad_grad_rho,
                                            m, div_m, grad_m, grad_div_m);
            const Vector3r div_grad_u = suzerain::rholut::div_grad_u(
                                            rho, grad_rho, div_grad_rho,
                                            m, grad_m, div_grad_m);

            // Compute quantities related to the equation of state
            real_t p, T, mu, lambda;
            Vector3r grad_p, grad_T, grad_mu, grad_lambda;
            suzerain::rholut::p_T_mu_lambda(
                alpha, beta, gamma, Ma,
                rho, grad_rho, m, grad_m, e, grad_e,
                p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
            const real_t div_grad_p = suzerain::rholut::div_grad_p(
                                        gamma, Ma,
                                        rho, grad_rho, div_grad_rho,
                                        m, grad_m, div_grad_m,
                                        e, grad_e, div_grad_e);
            const real_t div_grad_T = suzerain::rholut::div_grad_T(
                                        gamma,
                                        rho, grad_rho, div_grad_rho,
                                        p, grad_p, div_grad_p);

            // Compute quantities related to the viscous stress tensor
            const Matrix3r tau     = suzerain::rholut::tau(
                                        mu, lambda, div_u, grad_u);
            const Vector3r div_tau = suzerain::rholut::div_tau(
                                        mu, grad_mu, lambda, grad_lambda,
                                        div_u, grad_u, div_grad_u,
                                        grad_div_u);

            // Form continuity equation right hand side
            sphys(ndx::rho, offset) = - div_m
                ;

            // Form momentum equation right hand side
            const Vector3r momentum_rhs =
                - suzerain::rholut::div_u_outer_m(m, grad_m, u, div_u)
                - inv_Ma2 * grad_p
                + inv_Re * div_tau
                ;
            sphys(ndx::rhou, offset) = momentum_rhs.x();
            sphys(ndx::rhov, offset) = momentum_rhs.y();
            sphys(ndx::rhow, offset) = momentum_rhs.z();

            // Form energy equation right hand side
            sphys(ndx::rhoe, offset) =
                - suzerain::rholut::div_e_u(
                        e, grad_e, u, div_u
                    )
                - suzerain::rholut::div_p_u(
                        p, grad_p, u, div_u
                    )
                + inv_Re_Pr_gamma1 * suzerain::rholut::div_mu_grad_T(
                        grad_T, div_grad_T, mu, grad_mu
                    )
                + Ma2_over_Re * suzerain::rholut::div_tau_u<real_t>(
                        u, grad_u, tau, div_tau
                    )
                ;

            // Maintain the minimum observed stable time step, if necessary
            if (zeroth_substep) {
                namespace timestepper = suzerain::timestepper;
                // See convective_stability_criterion documentation for
                // why the magic number 4 modifies one_over_delta_y
                convective_delta_t = suzerain::math::minnan(
                        convective_delta_t,
                        timestepper::convective_stability_criterion(
                                u.x(), one_over_delta_x,
                                u.y(), one_over_delta_y_j / 4,
                                u.z(), one_over_delta_z,
                                evmaxmag_real,
                                std::sqrt(T) / Ma)); // a/u_0=sqrt(T*)/Ma
                const real_t nu = mu / rho;
                diffusive_delta_t = suzerain::math::minnan(
                        diffusive_delta_t,
                        timestepper::diffusive_stability_criterion(
                                one_over_delta_x,
                                one_over_delta_y_j,
                                one_over_delta_z,
                                Re, Pr, gamma, evmaxmag_imag,
                                nu, 0.0, alpha * nu, 0.0));
            }

        } // end X // end Z

    } // end Y

    // (3) Computing any manufactured solution forcing (when enabled).
    // Isolating this pass allows skipping the work when unnecessary
    if (msoln) {

        // Dereference the msoln smart pointer outside the compute loop
        const channel::manufactured_solution &ms = *msoln;

        offset = 0;
        for (int j = dgrid.local_physical_start.y();
             j < dgrid.local_physical_end.y();
             ++j) {

            const real_t y = this->y(j);

            for (int k = dgrid.local_physical_start.z();
                k < dgrid.local_physical_end.z();
                ++k) {

                const real_t z = this->z(k);

                for (int i = dgrid.local_physical_start.x();
                    i < dgrid.local_physical_end.x();
                    ++i, /* NB */ ++offset) {

                    const real_t x = this->x(i);

                    real_t Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoe;
                    ms.Q_conservative(x, y, z, time,
                                      Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoe);

                    sphys(ndx::rho,  offset) += Q_rho;
                    sphys(ndx::rhou, offset) += Q_rhou;
                    sphys(ndx::rhov, offset) += Q_rhov;
                    sphys(ndx::rhow, offset) += Q_rhow;
                    sphys(ndx::rhoe, offset) += Q_rhoe;

                } // end X

            } // end Z

        } // end Y

    } // end msoln

    // Collectively convert state to wave space using parallel FFTs
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        dgrid.transform_physical_to_wave(&sphys.coeffRef(i,0));
    }

    // Return the stable time step criteria separately on each rank.  The time
    // stepping logic must perform the Allreduce.  Delegating the Allreduce
    // responsibility allows reducing additional info with minimal overhead.
    return std::vector<real_t>(delta_t_candidates.begin(),
                               delta_t_candidates.end());

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

} // end namespace channel
