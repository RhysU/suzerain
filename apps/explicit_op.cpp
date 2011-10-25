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
        suzerain::ContiguousState<4,complex_t> &state) const
{
    SUZERAIN_UNUSED(phi);

    assert(static_cast<unsigned>(
                state.strides()[1]) == 1);
    assert(static_cast<unsigned>(
                state.strides()[2]) == state.shape()[1] * state.strides()[1]);
    assert(static_cast<unsigned>(
                state.strides()[3]) == state.shape()[2] * state.strides()[2]);
    assert(static_cast<unsigned>(
                state.strides()[0]) == state.shape()[3] * state.strides()[3]);

    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    assert(static_cast<unsigned>(massluz.n()) == state.shape()[1]);
    bop.apply(0, nrhs, 1, state.range().begin(), 1, state.shape()[1]);
}


void BsplineMassOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::ContiguousState<4,complex_t> &input,
        const complex_t &beta,
        suzerain::ContiguousState<4,complex_t> &output) const
{
    SUZERAIN_UNUSED(phi);
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
        suzerain::ContiguousState<4,complex_t> &state) const
{
    SUZERAIN_UNUSED(phi);

    assert(static_cast<unsigned>(
                state.strides()[1]) == 1);
    assert(static_cast<unsigned>(
                state.strides()[2]) == state.shape()[1] * state.strides()[1]);
    assert(static_cast<unsigned>(
                state.strides()[3]) == state.shape()[2] * state.strides()[2]);
    assert(static_cast<unsigned>(
                state.strides()[0]) == state.shape()[3] * state.strides()[3]);

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

void BsplineMassOperatorIsothermal::applyMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::ContiguousState<4,complex_t> &state) const
{
    // State enters method as coefficients in X, Y, and Z directions

    save_mean_state_at_collocation_points(state);

    return base::applyMassPlusScaledOperator(phi, state);

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

void BsplineMassOperatorIsothermal::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::ContiguousState<4,complex_t> &input,
        const complex_t &beta,
        suzerain::ContiguousState<4,complex_t> &output) const
{
    // State enters method as coefficients in X, Y, and Z directions

    save_mean_state_at_collocation_points(input);

    return base::accumulateMassPlusScaledOperator(phi, input, beta, output);

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

void BsplineMassOperatorIsothermal::save_mean_state_at_collocation_points(
        const state_type &state) const
{
    if (!has_zero_zero_mode) return;

    using Eigen::InnerStride;
    using Eigen::Map;
    using Eigen::VectorXr;

    // Copy mean coefficients.  Ugly syntax works around an Intel 11.1 20100806
    // (l_cproc_p_11.1.073) O3 segfault which is likely an optimizer bug.
    saved_mean_rho = Map<const VectorXr,0,InnerStride<2> >(
            reinterpret_cast<const real_t *>(
                state[channel::field::ndx::rho].origin()), state.shape()[1]);
    saved_mean_rhou = Map<const VectorXr,0,InnerStride<2> >(
            reinterpret_cast<const real_t *>(
                state[channel::field::ndx::rhou].origin()), state.shape()[1]);

    // Convert mean coefficients to mean collocation point values
    bop.apply(0, 1, 1.0, saved_mean_rho.data(),  1, state.shape()[1]);
    bop.apply(0, 1, 1.0, saved_mean_rhou.data(), 1, state.shape()[1]);
}

void BsplineMassOperatorIsothermal::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::ContiguousState<4,complex_t> &state) const
{
    // State enters method as coefficients in X and Z directions
    // State enters method as collocation point values in Y direction

    // Shorthand
    using Eigen::Map;
    using Eigen::VectorXc;
    using Eigen::VectorXr;
    namespace ndx = channel::field::ndx;
    const std::size_t Ny         = state.shape()[1];
    const std::size_t wall_lower = 0;
    const std::size_t wall_upper = Ny - 1;

    // See channel_treatment writeup for information on the steps below.
    // Steps appear out of order relative to the writeup (TODO fix writeup).

    // channel_treatment step (1) done during operator apply/accumulate
    // within save_mean_state_at_collocation_points();

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

    // channel_treatment step (2) loads mean state at collocation points
    // into the imaginary part of the constant (zero zero) mode coefficients
    if (constrain_bulk_rhou && has_zero_zero_mode) {

        saved_mean_rho[wall_lower] = 0;  // No forcing at lower wall
        saved_mean_rho[wall_upper] = 0;  // No forcing at upper wall
        Map<VectorXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        mean_rhou.imag() = saved_mean_rho;

        saved_mean_rhou[wall_lower] = 0;  // No forcing at lower wall
        saved_mean_rhou[wall_upper] = 0;  // No forcing at upper wall
        Map<VectorXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.imag() = saved_mean_rhou;
    }

    // B-spline numerics are not conservative when the mean density field has
    // any asymmetry.  Combat the bulk density's very, very small tendency to
    // drift by adding an integral constraint that the bulk density stay fixed
    // at a target value.  Approach follows that of the bulk momentum forcing.
    if (constrain_bulk_rho && has_zero_zero_mode) {
        Map<VectorXc> mean_rho(state[ndx::rho].origin(), Ny);
        mean_rho.imag().setConstant(1);
    }

    // channel_treatment step (3) performs the usual operator solve
    base::invertMassPlusScaledOperator(phi, state);

    if (constrain_bulk_rhou && has_zero_zero_mode) {

        // channel_treatment steps (4), (5), and (6) determine and
        // apply the appropriate bulk momentum forcing to achieve
        // a target value.
        Map<VectorXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rhou);
        const real_t varphi = (scenario.bulk_rhou - bulk.real()) / bulk.imag();
        mean_rhou.real() += varphi * mean_rhou.imag();
        mean_rhou.imag() = VectorXr::Zero(Ny);

        // channel_treatment step (7) accounts for the momentum forcing
        // within the total energy equation including the Mach squared
        // factor arising from the nondimensionalization choices.
        Map<VectorXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.real() += (varphi*scenario.Ma*scenario.Ma)*mean_rhoe.imag();
        mean_rhoe.imag() = VectorXr::Zero(Ny);

        // channel_treatment steps (8) and (9) already performed above
    }

    // Complete the bulk density target forcing
    if (constrain_bulk_rho && has_zero_zero_mode) {

        Map<VectorXc> mean_rho(state[ndx::rho].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rho);
        const real_t varphi = (scenario.bulk_rho - bulk.real()) / bulk.imag();
        mean_rho.real() += varphi * mean_rho.imag();
        mean_rho.imag() = VectorXr::Zero(Ny);

        // Note that only the continuity equation is forced which neglects
        // contributions to the momentum and energy equations.  Neglecting
        // these terms is acceptable as the forcing is very small and vanishes
        // when the mean field is symmetric in the wall-normal direction.

    }

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
    const bool delta_t_requested) const
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
                if (delta_t_requested) {
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

    } // end Y


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

NonlinearOperatorIsothermal::NonlinearOperatorIsothermal(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        OperatorCommonBlock &common,
        const boost::shared_ptr<
            const channel::manufactured_solution>& msoln)
    : NonlinearOperator(scenario, grid, dgrid, b, bop, common, msoln)
{
    // Precompute operator for finding bulk quantities from coefficients
    bulkcoeff.resize(b.n());
    b.integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= scenario.Ly;
}

std::vector<real_t> NonlinearOperatorIsothermal::applyOperator(
    const real_t time,
    suzerain::ContiguousState<4,complex_t> &swave,
    const real_t evmaxmag_real,
    const real_t evmaxmag_imag,
    const bool delta_t_requested) const
{
    using Eigen::Map;
    using Eigen::VectorXc;
    namespace ndx = channel::field::ndx;

    // Shorthand for the wall-normal size and wall indices
    const std::size_t Ny         = swave.shape()[1];
    const std::size_t wall_lower = 0;
    const std::size_t wall_upper = Ny - 1;

    // Compute and store quantities used later to implement forcing
    real_t bulk_density = std::numeric_limits<real_t>::quiet_NaN();
    if (has_zero_zero_mode) {

        // Save mean density values at collocation points
        rho_fm = Map<VectorXc>(swave[ndx::rho].origin(), Ny).real();
        bop.apply(0, 1, 1.0, rho_fm.data(), 1, Ny);

        // Save bulk density
        bulk_density = bulkcoeff.dot(rho_fm.real());

        // Save mean X momentum values at collocation points
        fm_dot_m = Map<VectorXc>(swave[ndx::rhou].origin(), Ny).real();
        bop.apply(0, 1, 1.0, fm_dot_m.data(), 1, Ny);
    }

    // Apply an operator that cares nothing about the boundaries.
    // Operator application turns coefficients in X, Y, and Z into
    // coefficients in X and Z but COLLOCATION POINT VALUES IN Y.
    const std::vector<real_t> delta_t_candidates = base::applyOperator(
            time, swave, evmaxmag_real, evmaxmag_imag, delta_t_requested);

    // Set no-slip condition for momentum on walls per writeup step (3)
    // Condition achieved by removing time evolution at walls
    assert(static_cast<int>(ndx::rhov) == static_cast<int>(ndx::rhou) + 1);
    assert(static_cast<int>(ndx::rhow) == static_cast<int>(ndx::rhov) + 1);
    for (std::size_t i = ndx::rhou; i <= ndx::rhow; ++i) {
        for (std::size_t k = 0; k < swave.shape()[3]; ++k) {
            for (std::size_t j = 0; j < swave.shape()[2]; ++j) {
                swave[i][wall_lower][j][k] = 0;
                swave[i][wall_upper][j][k] = 0;
            }
        }
    }

    // Set isothermal condition on walls per writeup step (4)
    // Condition achieved by removing time evolution at walls
    // independent of changes due to local density evolution
    const real_t inv_gamma_gamma1
        = 1 / (scenario.gamma * (scenario.gamma - 1));
    for (std::size_t k = 0; k < swave.shape()[3]; ++k) {
        for (std::size_t j = 0; j < swave.shape()[2]; ++j) {
            swave[ndx::rhoe][wall_lower][j][k]
                = inv_gamma_gamma1 * swave[ndx::rho][wall_lower][j][k];
            swave[ndx::rhoe][wall_upper][j][k]
                = inv_gamma_gamma1 * swave[ndx::rho][wall_upper][j][k];
        }
    }

    // Apply f_{m_x} to mean x-momentum, mean energy at non-wall locations
    if (has_zero_zero_mode) {

        // Compute temporary per writeup implementation step (5)
        Map<VectorXc> mean_rhou(swave[ndx::rhou].origin(), Ny);
        const real_t alpha = bulkcoeff.dot(mean_rhou.real()) / bulk_density;

        // Apply to non-wall mean x-momentum right hand side per step (6)
        rho_fm.head<1>()[0] = 0;
        rho_fm.tail<1>()[0] = 0;
        mean_rhou.real() -= alpha * rho_fm;

        // Apply to non-wall mean energy right hand side per step (7)
        fm_dot_m.head<1>()[0] = 0;
        fm_dot_m.tail<1>()[0] = 0;
        Eigen::Map<Eigen::VectorXc> mean_rhoe(swave[ndx::rhoe].origin(), Ny);
        mean_rhoe.real() -= alpha * fm_dot_m;
    }

    // Return the time step found by the BC-agnostic operator
    return delta_t_candidates;
}

} // end namespace channel
