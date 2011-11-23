//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// explicit_op.hpp: Operators for channel_explicit
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

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

    massluz.form_mass(bop);
}

void BsplineMassOperator::applyMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::ContiguousState<4,complex_t> &state,
        const std::size_t substep_index) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(substep_index);

    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    assert(static_cast<unsigned>(massluz.n()) == state.shape()[1]);
    bop.apply(0, nrhs, 1, state.range().begin(), 1, state.shape()[1]);
}


void BsplineMassOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::ContiguousState<4,complex_t> &input,
        const complex_t &beta,
        suzerain::ContiguousState<4,complex_t> &output,
        const std::size_t substep_index) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(substep_index);
    const state_type &x   = input;  // Shorthand
    state_type &y         = output; // Shorthand
    const complex_t c_one = 1;
    assert(x.isIsomorphic(y));

    typedef state_type::index index;
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
        suzerain::ContiguousState<4,complex_t> &state,
        const std::size_t substep_index) const
{
    SUZERAIN_UNUSED(phi);
    SUZERAIN_UNUSED(substep_index);

    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    assert(static_cast<unsigned>(massluz.n()) == state.shape()[1]);
    massluz.solve(nrhs, state.range().begin(), 1, state.shape()[1]);
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

void BsplineMassOperatorIsothermal::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::ContiguousState<4,complex_t> &state,
        const std::size_t substep_index) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    using Eigen::Map;
    using Eigen::ArrayXc;
    using Eigen::ArrayXr;
    namespace ndx = channel::field::ndx;
    const std::size_t Ny         = state.shape()[1];
    const std::size_t wall_lower = 0;
    const std::size_t wall_upper = Ny - 1;

    // See channel_treatment writeup for information on the steps below.
    // Steps appear out of order relative to the writeup (TODO fix writeup).

    // Means of the implicit momentum and energy forcing coefficients(!) are
    // maintained across each individual time step for sampling the statistics
    // /bar_f, /bar_f_dot_u, and /bar_qb using OperatorCommonBlock.
    //
    // The accumulated means are updated in-place using expressions like
    //    mean = i/(i+1) * last + 1 / (i+1) * update
    // where is is the current substep index.  Values are reset on i = 0.
    const real_t prev_mean_coeff = real_t(substep_index) / (substep_index + 1);
    const real_t curr_mean_coeff = real_t(1)             / (substep_index + 1);
    common.storage.resize(Ny, Eigen::NoChange);        // Nondestructive on NOP

    // channel_treatment step (1) done during nonlinear operator application
    // via shared OperatorCommonBlock storage space

    // channel_treatment step (8) sets no-slip conditions
    // on wall collocation points.  Possible pre-solve since L = 0.
    assert(static_cast<int>(ndx::rhov) == static_cast<int>(ndx::rhou) + 1);
    assert(static_cast<int>(ndx::rhow) == static_cast<int>(ndx::rhov) + 1);
    for (std::size_t i = ndx::rhou; i <= ndx::rhow; ++i) {
        for (std::size_t k = 0; k < state.shape()[3]; ++k) {
            for (std::size_t j = 0; j < state.shape()[2]; ++j) {
                state[i][wall_lower][j][k] = 0;
                state[i][wall_upper][j][k] = 0;
            }
        }
    }

    // channel_treatment step (9) sets isothermal conditions on wall
    // collocation points using e_wall = rho_wall / (gamma * (gamma - 1)).
    // Possible pre-solve since L = 0.
    const real_t inv_gamma_gamma1
        = 1 / (scenario.gamma * (scenario.gamma - 1));
    for (std::size_t k = 0; k < state.shape()[3]; ++k) {
        for (std::size_t j = 0; j < state.shape()[2]; ++j) {
            state[ndx::rhoe][wall_lower][j][k]
                    = inv_gamma_gamma1 * state[ndx::rho][wall_lower][j][k];
            state[ndx::rhoe][wall_upper][j][k]
                    = inv_gamma_gamma1 * state[ndx::rho][wall_upper][j][k];
        }
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
    if (constrain_bulk_rhou && has_zero_zero_mode) {

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
    if (constrain_bulk_rho && has_zero_zero_mode) {
        Map<ArrayXc> mean_rho(state[ndx::rho].origin(), Ny);
        mean_rho.imag().setConstant(1);
    }

    // channel_treatment step (3) performs the usual operator solve
    base::invertMassPlusScaledOperator(phi, state, substep_index);

    if (constrain_bulk_rhou && has_zero_zero_mode) {

        // channel_treatment steps (4), (5), and (6) determine and
        // apply the appropriate bulk momentum forcing to achieve
        // a target value.
        Map<ArrayXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rhou.matrix());
        const real_t varphi = (scenario.bulk_rhou - bulk.real()) / bulk.imag();
        mean_rhou.real() += varphi * mean_rhou.imag();
        common.f()        = prev_mean_coeff * common.f()
                          + (curr_mean_coeff * varphi) * mean_rhou.imag();
        mean_rhou.imag().setZero();

        // channel_treatment step (7) accounts for the momentum forcing
        // within the total energy equation including the Mach squared
        // factor arising from the nondimensionalization choices.
        Map<ArrayXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        const real_t varphi_Ma_Ma = varphi * scenario.Ma * scenario.Ma;
        mean_rhoe.real() += varphi_Ma_Ma * mean_rhoe.imag();
        common.f_dot_u()  = prev_mean_coeff * common.f_dot_u()
                          + (curr_mean_coeff * varphi_Ma_Ma) * mean_rhoe.imag();
        mean_rhoe.imag().setZero();

        // channel_treatment steps (8) and (9) already performed above
    }

    // Complete the bulk density target forcing
    if (constrain_bulk_rho && has_zero_zero_mode) {

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
    common.qb() *= prev_mean_coeff;

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

std::vector<real_t> NonlinearOperator::applyOperator(
    const real_t time,
    suzerain::ContiguousState<4,complex_t> &swave,
    const real_t evmaxmag_real,
    const real_t evmaxmag_imag,
    const std::size_t substep_index) const
{
    namespace ndx = channel::field::ndx;

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

    // Prepare common-block-like storage used to pass details from N to L
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
    channel::physical_view<aux::count>::type auxp
        = channel::physical_view<aux::count>::create(dgrid, auxw);
    channel::physical_view<channel::field::count>::type sphys
        = channel::physical_view<channel::field::count>::create(dgrid, swave);
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        dgrid.transform_wave_to_physical(&sphys(i,0));
    }
    for (std::size_t i = 0; i < aux::count; ++i) {
        dgrid.transform_wave_to_physical(&auxp(i,0));
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

    // Working non-scalar storage used within following loop
    Eigen::Vector3r grad_rho;
    Eigen::Matrix3r grad_grad_rho;
    Eigen::Vector3r m;
    Eigen::Matrix3r grad_m;
    Eigen::Vector3r div_grad_m;
    Eigen::Vector3r grad_div_m;
    Eigen::Vector3r grad_e;
    Eigen::Vector3r u;
    Eigen::Matrix3r grad_u;
    Eigen::Vector3r grad_div_u, div_grad_u;
    Eigen::Vector3r grad_p, grad_T, grad_mu, grad_lambda;
    Eigen::Matrix3r tau;
    Eigen::Vector3r div_tau;
    Eigen::Vector3r momentum_rhs;

    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
    size_t offset = 0;
    for (int j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Used to accumulate mean quantities versus wall-normal position
        boost::accumulators::accumulator_set<
                real_t,
#if BOOST_VERSION >= 104700
                boost::accumulators::stats<boost::accumulators::tag::sum_kahan>
#else
                boost::accumulators::stats<boost::accumulators::tag::sum>
#endif
            > sum_u;

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Unpack density-related quantities
                const real_t rho          = sphys(ndx::rho, offset);
                grad_rho.x()              = auxp(aux::rho_x, offset);
                grad_rho.y()              = auxp(aux::rho_y, offset);
                grad_rho.z()              = auxp(aux::rho_z, offset);
                const real_t div_grad_rho = auxp(aux::rho_xx, offset)
                                          + auxp(aux::rho_yy, offset)
                                          + auxp(aux::rho_zz, offset);
                grad_grad_rho(0,0)        = auxp(aux::rho_xx, offset);
                grad_grad_rho(0,1)        = auxp(aux::rho_xy, offset);
                grad_grad_rho(0,2)        = auxp(aux::rho_xz, offset);
                grad_grad_rho(1,0)        = grad_grad_rho(0,1);
                grad_grad_rho(1,1)        = auxp(aux::rho_yy, offset);
                grad_grad_rho(1,2)        = auxp(aux::rho_yz, offset);
                grad_grad_rho(2,0)        = grad_grad_rho(0,2);
                grad_grad_rho(2,1)        = grad_grad_rho(1,2);
                grad_grad_rho(2,2)        = auxp(aux::rho_zz, offset);

                // Unpack momentum-related quantities
                m.x()              = sphys(ndx::rhou, offset);
                m.y()              = sphys(ndx::rhov, offset);
                m.z()              = sphys(ndx::rhow, offset);
                const real_t div_m = auxp(aux::mx_x, offset)
                                   + auxp(aux::my_y, offset)
                                   + auxp(aux::mz_z, offset);
                grad_m(0,0)        = auxp(aux::mx_x, offset);
                grad_m(0,1)        = auxp(aux::mx_y, offset);
                grad_m(0,2)        = auxp(aux::mx_z, offset);
                grad_m(1,0)        = auxp(aux::my_x, offset);
                grad_m(1,1)        = auxp(aux::my_y, offset);
                grad_m(1,2)        = auxp(aux::my_z, offset);
                grad_m(2,0)        = auxp(aux::mz_x, offset);
                grad_m(2,1)        = auxp(aux::mz_y, offset);
                grad_m(2,2)        = auxp(aux::mz_z, offset);
                div_grad_m.x()     = auxp(aux::mx_xx, offset)
                                   + auxp(aux::mx_yy, offset)
                                   + auxp(aux::mx_zz, offset);
                div_grad_m.y()     = auxp(aux::my_xx, offset)
                                   + auxp(aux::my_yy, offset)
                                   + auxp(aux::my_zz, offset);
                div_grad_m.z()     = auxp(aux::mz_xx, offset)
                                   + auxp(aux::mz_yy, offset)
                                   + auxp(aux::mz_zz, offset);
                grad_div_m.x()     = auxp(aux::mx_xx, offset)
                                   + auxp(aux::my_xy, offset)
                                   + auxp(aux::mz_xz, offset);
                grad_div_m.y()     = auxp(aux::mx_xy, offset)
                                   + auxp(aux::my_yy, offset)
                                   + auxp(aux::mz_yz, offset);
                grad_div_m.z()     = auxp(aux::mx_xz, offset)
                                   + auxp(aux::my_yz, offset)
                                   + auxp(aux::mz_zz, offset);

                // Unpack total energy-related quantities
                const real_t e          = sphys(ndx::rhoe, offset);
                grad_e.x()              = auxp(aux::e_x, offset);
                grad_e.y()              = auxp(aux::e_y, offset);
                grad_e.z()              = auxp(aux::e_z, offset);
                const real_t div_grad_e = auxp(aux::div_grad_e, offset);

                // Compute quantities based upon state.  Real-valued scalars
                // are declared inline.  Vector- and tensor-valued expressions
                // declared outside loop.
                u                  = suzerain::rholut::u(
                                        rho, m);
                sum_u(u.x());
                const real_t div_u = suzerain::rholut::div_u(
                                        rho, grad_rho, m, div_m);
                grad_u             = suzerain::rholut::grad_u(
                                        rho, grad_rho, m, grad_m);
                grad_div_u         = suzerain::rholut::grad_div_u(
                                        rho, grad_rho, grad_grad_rho,
                                        m, div_m, grad_m, grad_div_m);
                div_grad_u         = suzerain::rholut::div_grad_u(
                                        rho, grad_rho, div_grad_rho,
                                        m, grad_m, div_grad_m);
                real_t p, T, mu, lambda;
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
                tau     = suzerain::rholut::tau(
                            mu, lambda, div_u, grad_u);
                div_tau = suzerain::rholut::div_tau(
                            mu, grad_mu, lambda, grad_lambda,
                            div_u, grad_u, div_grad_u, grad_div_u);

                // Form continuity equation right hand side
                sphys(ndx::rho, offset) = - div_m
                    ;

                // Form momentum equation right hand side
                momentum_rhs =
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
                if (substep_index == 0) {
                    namespace timestepper = suzerain::timestepper;
                    // See convective_stability_criterion documentation for
                    // why the magic number 4 modifies one_over_delta_y
                    convective_delta_t = suzerain::math::minnan(
                            convective_delta_t,
                            timestepper::convective_stability_criterion(
                                    u.x(), one_over_delta_x,
                                    u.y(), one_over_delta_y(j) / 4,
                                    u.z(), one_over_delta_z,
                                    evmaxmag_real,
                                    std::sqrt(T) / Ma)); // a/u_0=sqrt(T*)/Ma
                    const real_t nu = mu / rho;
                    diffusive_delta_t = suzerain::math::minnan(
                            diffusive_delta_t,
                            timestepper::diffusive_stability_criterion(
                                    one_over_delta_x,
                                    one_over_delta_y(j),
                                    one_over_delta_z,
                                    Re, Pr, gamma, evmaxmag_imag,
                                    nu, 0.0, alpha * nu, 0.0));
                }

            } // end X

        } // end Z

        // Store sum(s) into common block in preparation for MPI reduction
        common.u()[j] = boost::accumulators::sum(sum_u);

    } // end Y

    // Reduce and scale common.u() sums to obtain mean quantities on rank zero
    if (has_zero_zero_mode) {
        assert(suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0);
        SUZERAIN_MPICHKR(MPI_Reduce(MPI_IN_PLACE, common.u().data(),
                    common.u().size(), suzerain::mpi::datatype<real_t>::value,
                    MPI_SUM, 0, MPI_COMM_WORLD));
        common.u() /= (   dgrid.global_physical_extent.x()
                        * dgrid.global_physical_extent.z());
    } else {
        Eigen::ArrayXr tmp;
        tmp.resizeLike(common.u());
        SUZERAIN_MPICHKR(MPI_Reduce(common.u().data(), tmp.data(),
                    common.u().size(), suzerain::mpi::datatype<real_t>::value,
                    MPI_SUM, 0, MPI_COMM_WORLD));
    }

    // If active, add manufactured solution forcing in a second pass.
    // Isolating this pass allows us to quickly skip the associated work when
    // not using a manufactured solution.  Same traversal idiom as above.
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
        dgrid.transform_physical_to_wave(&sphys(i,0));
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
