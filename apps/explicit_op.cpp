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
#include <suzerain/orthonormal.hpp>
#include <suzerain/utility.hpp>
#include <suzerain/state.hpp>

#include "explicit_op.hpp"

#pragma warning(disable:383 1572)

namespace channel {

OperatorBase::OperatorBase(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop)
    : scenario(scenario),
      grid(grid),
      dgrid(dgrid),
      bop(bop),
      has_zero_zero_mode(    dgrid.local_wave_start.x() == 0
                          && dgrid.local_wave_start.z() == 0),
      one_over_delta_x(scenario.Lx / grid.N.x() /* !dN.x() */),
      one_over_delta_z(scenario.Lz / grid.N.z() /* !dN.z() */),
      y_(boost::extents[boost::multi_array_types::extent_range(
              dgrid.local_physical_start.y(), dgrid.local_physical_end.y())]),
      one_over_delta_y_(boost::extents[boost::multi_array_types::extent_range(
              dgrid.local_physical_start.y(), dgrid.local_physical_end.y())])
{
    // Compute y collocation point locations and spacing local to this rank
    for (int j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        y_[j] = b.collocation_point(j);
        const int jm = (j == 0        ) ? 1         : j - 1;
        const int jp = (j == b.n() - 1) ? b.n() - 2 : j + 1;
        const real_t delta_y = std::min(
                std::abs(b.collocation_point(jm) - y_[j]),
                std::abs(b.collocation_point(jp) - y_[j]));
        one_over_delta_y_[j] = 1.0 / delta_y;
    }
}

BsplineMassOperator::BsplineMassOperator(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop)
    : OperatorBase(scenario, grid, dgrid, b, bop),
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
        suzerain::NoninterleavedState<4,complex_t> &state) const
{
    SUZERAIN_UNUSED(phi);

    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    assert(1 == state.strides()[1]);
    assert(static_cast<unsigned>(massluz.n()) == state.shape()[1]);
    bop.apply(0, nrhs, 1, state.memory_begin(), 1, state.strides()[2]);
}


void BsplineMassOperator::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::NoninterleavedState<4,complex_t> &input,
        const complex_t &beta,
        suzerain::NoninterleavedState<4,complex_t> &output) const
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
        suzerain::NoninterleavedState<4,complex_t> &state) const
{
    SUZERAIN_UNUSED(phi);

    const int nrhs = state.shape()[0]*state.shape()[2]*state.shape()[3];
    assert(1 == state.strides()[1]);
    assert(static_cast<unsigned>(massluz.n()) == state.shape()[1]);
    massluz.solve(nrhs, state.memory_begin(), 1, state.strides()[2]);
}

BsplineMassOperatorIsothermal::BsplineMassOperatorIsothermal(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop)
    : BsplineMassOperator(scenario, grid, dgrid, b, bop)
{
    if (has_zero_zero_mode) {

        // Precompute operator for finding bulk quantities from coefficients
        bulkcoeff.resize(b.n());
        b.integration_coefficients(0, bulkcoeff.data());
        bulkcoeff /= scenario.Ly;

        // Precompute M^{-1}*(e_{lower wall}, e_{upper_wall} for applying
        // collocation point conditions in coefficient space
        massinv_elower = Eigen::VectorXr::Unit(b.n(), 0);
        massinv_eupper = Eigen::VectorXr::Unit(b.n(), b.n() - 1);
        suzerain::bsplineop_lu masslu(bop);
        masslu.form_mass(bop);
        masslu.solve(1, massinv_elower.data(), 1, b.n());
        masslu.solve(1, massinv_eupper.data(), 1, b.n());
    }
}

void BsplineMassOperatorIsothermal::applyMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::NoninterleavedState<4,complex_t> &state) const
{
    // State enters method as coefficients in X, Y, and Z directions

    save_mean_state_at_collocation_points(state);

    return base::applyMassPlusScaledOperator(phi, state);

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

void BsplineMassOperatorIsothermal::accumulateMassPlusScaledOperator(
        const complex_t &phi,
        const suzerain::NoninterleavedState<4,complex_t> &input,
        const complex_t &beta,
        suzerain::NoninterleavedState<4,complex_t> &output) const
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

    // Copy mean coefficients
    saved_mean_rho = Eigen::Map<const Eigen::VectorXc>(
        state[channel::field::ndx::rho].origin(), state.shape()[1]).real();
    saved_mean_rhou = Eigen::Map<const Eigen::VectorXc>(
        state[channel::field::ndx::rhou].origin(), state.shape()[1]).real();

    // Convert mean coefficients to mean collocation point values
    bop.apply(0, 1, 1.0, saved_mean_rho.data(),  1, state.shape()[1]);
    bop.apply(0, 1, 1.0, saved_mean_rhou.data(), 1, state.shape()[1]);
}

void BsplineMassOperatorIsothermal::invertMassPlusScaledOperator(
        const complex_t &phi,
        suzerain::NoninterleavedState<4,complex_t> &state) const
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

    // See channel_treatment document for information on the steps below

    // Channel treatment step (1) done during operator apply/accumulate
    // within save_mean_state_at_collocation_points();

    // Channel treatment step (2) loads mean state at collocation points
    // into the imaginary part of the constant mode coefficients
    if (has_zero_zero_mode) {
        saved_mean_rho[wall_lower] = 0;  // No forcing at lower wall
        saved_mean_rho[wall_upper] = 0;  // No forcing at upper wall
        Map<VectorXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        mean_rhou.imag() = saved_mean_rho;

        saved_mean_rhou[wall_lower] = 0;  // No forcing at lower wall
        saved_mean_rhou[wall_upper] = 0;  // No forcing at upper wall
        Map<VectorXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.imag() = saved_mean_rhou;
    }

    // Channel treatment step (3) performs the usual operator solve
    base::invertMassPlusScaledOperator(phi, state);

    if (has_zero_zero_mode) {

        // Channel treatment steps (4), (5), and (6) determine and
        // apply the appropriate bulk momentum forcing using the
        // Mach number as the target value
        Map<VectorXc> mean_rhou(state[ndx::rhou].origin(), Ny);
        const complex_t bulk = bulkcoeff.cast<complex_t>().dot(mean_rhou);
        const real_t phi = (scenario.Ma - bulk.real()) / bulk.imag();
        mean_rhou.real() += phi * mean_rhou.imag();
        mean_rhou.imag() = VectorXr::Zero(Ny);

        // Channel treatment step (7) accounts for the momentum forcing
        // within the total energy equation
        Map<VectorXc> mean_rhoe(state[ndx::rhoe].origin(), Ny);
        mean_rhoe.real() += phi * mean_rhoe.imag();
        mean_rhoe.imag() = VectorXr::Zero(Ny);

        // Channel treatment step (8) enforces no-slip at the walls
        // which is complicated by being in coefficient space.
        // Zero the collocation point value at both walls.
        assert(static_cast<int>(ndx::rhov) == static_cast<int>(ndx::rhou) + 1);
        assert(static_cast<int>(ndx::rhow) == static_cast<int>(ndx::rhov) + 1);
        for (std::size_t i = ndx::rhou; i <= ndx::rhow; ++i) {
            for (std::size_t k = 0; k < state.shape()[3]; ++k) {
                for (std::size_t j = 0; j < state.shape()[2]; ++j) {

                    const complex_t m_lower = state[i][wall_lower][j][k];
                    const complex_t m_upper = state[i][wall_upper][j][k];

                    Map<VectorXc> m_jk(&state[i][0][j][k], Ny);
                    m_jk -= (   massinv_elower*m_lower
                              + massinv_eupper*m_upper).cast<complex_t>();
                }
            }
        }

        // Channel treatment step (9) enforces isothermal wall by
        // setting e_wall = rho_wall / (gamma * (gamma - 1)) again
        // dealing with imposing collocation point restrictions in
        // coefficient space.  Uses that only the wall coefficients
        // contribute to rho_wall's collocation point value.
        const real_t inv_gamma_gamma1
            = 1 / (scenario.gamma * (scenario.gamma - 1));
        for (std::size_t k = 0; k < state.shape()[3]; ++k) {
            for (std::size_t j = 0; j < state.shape()[2]; ++j) {

                    Map<VectorXc> rhoe_jk(&state[ndx::rhoe][0][j][k], Ny);

                    const complex_t factor_lower
                        = inv_gamma_gamma1*state[ndx::rho][wall_lower][j][k]
                        -                  rhoe_jk[wall_lower];
                    const complex_t factor_upper
                        = inv_gamma_gamma1*state[ndx::rho][wall_upper][j][k]
                        -                  rhoe_jk[wall_upper];
                    rhoe_jk += (   massinv_elower*factor_lower
                                 + massinv_eupper*factor_upper).cast<complex_t>();
            }
        }

    }

    // State leaves method as coefficients in X, Y, and Z directions
}

NonlinearOperator::NonlinearOperator(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop)
    : OperatorBase(scenario, grid, dgrid, b, bop),
      auxw(suzerain::to_yxz(static_cast<std::size_t>(aux::count),
                            dgrid.local_wave_extent),
              suzerain::prepend(dgrid.local_wave_storage(),
                  suzerain::strides_cm(
                      suzerain::to_yxz(dgrid.local_wave_extent))))
{
    // NOP
}

real_t NonlinearOperator::applyOperator(
    suzerain::NoninterleavedState<4,complex_t> &swave,
    const real_t evmaxmag_real,
    const real_t evmaxmag_imag,
    const bool delta_t_requested) const
{
    namespace ndx = channel::field::ndx;

    // State enters method as coefficients in X, Y, and Z directions

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

    real_t delta_t_candidates[2] = { std::numeric_limits<real_t>::max(),
                                     std::numeric_limits<real_t>::max()  };
    real_t &convective_delta_t = delta_t_candidates[0];
    real_t &diffusive_delta_t  = delta_t_candidates[1];

    // Compute Y derivatives of density at collocation points
    bop_accumulate(1, 1., swave, ndx::rho, 0., auxw, aux::rho_y);
    bop_accumulate(2, 1., swave, ndx::rho, 0., auxw, aux::rho_yy);
    bop_apply     (0, 1., swave, ndx::rho);

    // Compute X- and Z- derivatives of density at collocation points
    diffwave_accumulate(1, 0, 1., swave, ndx::rho,   0., auxw, aux::rho_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rho,   0., auxw, aux::rho_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rho,   0., auxw, aux::rho_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rho,   0., auxw, aux::rho_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rho,   0., auxw, aux::rho_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::rho_y, 0., auxw, aux::rho_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::rho_y, 0., auxw, aux::rho_yz);

    // Compute Y derivatives of X momentum at collocation points
    bop_accumulate(1, 1., swave, ndx::rhou, 0., auxw, aux::mx_y);
    bop_accumulate(2, 1., swave, ndx::rhou, 0., auxw, aux::mx_yy);
    bop_apply     (0, 1., swave, ndx::rhou);

    // Compute X- and Z- derivatives of X momentum at collocation points
    diffwave_accumulate(1, 0, 1., swave, ndx::rhou,  0., auxw, aux::mx_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhou,  0., auxw, aux::mx_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rhou,  0., auxw, aux::mx_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhou,  0., auxw, aux::mx_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhou,  0., auxw, aux::mx_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::mx_y,  0., auxw, aux::mx_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::mx_y,  0., auxw, aux::mx_yz);

    // Compute Y derivatives of Y momentum at collocation points
    bop_accumulate(1, 1., swave, ndx::rhov, 0., auxw, aux::my_y);
    bop_accumulate(2, 1., swave, ndx::rhov, 0., auxw, aux::my_yy);
    bop_apply     (0, 1., swave, ndx::rhov);

    // Compute X- and Z- derivatives of Y momentum at collocation points
    diffwave_accumulate(1, 0, 1., swave, ndx::rhov,  0., auxw, aux::my_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhov,  0., auxw, aux::my_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rhov,  0., auxw, aux::my_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhov,  0., auxw, aux::my_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhov,  0., auxw, aux::my_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::my_y,  0., auxw, aux::my_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::my_y,  0., auxw, aux::my_yz);

    // Compute Y derivatives of Z momentum at collocation points
    bop_accumulate(1, 1., swave, ndx::rhow, 0., auxw, aux::mz_y);
    bop_accumulate(2, 1., swave, ndx::rhow, 0., auxw, aux::mz_yy);
    bop_apply     (0, 1., swave, ndx::rhow);

    // Compute X- and Z- derivatives of Z momentum at collocation points
    diffwave_accumulate(1, 0, 1., swave, ndx::rhow,  0., auxw, aux::mz_x );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhow,  0., auxw, aux::mz_xx);
    diffwave_accumulate(1, 1, 1., swave, ndx::rhow,  0., auxw, aux::mz_xz);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhow,  0., auxw, aux::mz_z );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhow,  0., auxw, aux::mz_zz);
    diffwave_accumulate(1, 0, 1., auxw,  aux::mz_y,  0., auxw, aux::mz_xy);
    diffwave_accumulate(0, 1, 1., auxw,  aux::mz_y,  0., auxw, aux::mz_yz);

    // Compute Y derivatives of total energy at collocation points
    bop_accumulate(1, 1., swave, ndx::rhoe, 0., auxw, aux::e_y);
    bop_accumulate(2, 1., swave, ndx::rhoe, 0., auxw, aux::div_grad_e);
    bop_apply     (0, 1., swave, ndx::rhoe);

    // Compute X- and Z- derivatives of total energy at collocation points
    diffwave_accumulate(1, 0, 1., swave, ndx::rhoe, 0., auxw, aux::e_x       );
    diffwave_accumulate(2, 0, 1., swave, ndx::rhoe, 1., auxw, aux::div_grad_e);
    diffwave_accumulate(0, 1, 1., swave, ndx::rhoe, 0., auxw, aux::e_z       );
    diffwave_accumulate(0, 2, 1., swave, ndx::rhoe, 1., auxw, aux::div_grad_e);

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    Eigen::Map<
            Eigen::Array<real_t, aux::count,
                         Eigen::Dynamic, Eigen::RowMajor>,
            Eigen::Unaligned, // Defensive
            Eigen::OuterStride<Eigen::Dynamic>
        > auxp(reinterpret_cast<real_t *>(auxw.origin()),
               aux::count, dgrid.local_physical_extent.prod(),
               Eigen::OuterStride<>(
                   auxw.strides()[0] * sizeof(complex_t)/sizeof(real_t)));
    Eigen::Map<
            Eigen::Array<real_t, channel::field::count,
                         Eigen::Dynamic, Eigen::RowMajor>,
            Eigen::Unaligned, // Defensive
            Eigen::OuterStride<Eigen::Dynamic>
        > sphys(reinterpret_cast<real_t *>(swave.origin()),
               channel::field::count, dgrid.local_physical_extent.prod(),
               Eigen::OuterStride<>(
                   swave.strides()[0] * sizeof(complex_t)/sizeof(real_t)));
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        dgrid.transform_wave_to_physical(&sphys(i,0));
    }
    for (std::size_t i = 0; i < aux::count; ++i) {
        dgrid.transform_wave_to_physical(&auxp(i,0));
    }

    // Retrieve constants and compute derived constants
    const real_t beta             = scenario.beta;
    const real_t gamma            = scenario.gamma;
    const real_t Pr               = scenario.Pr;
    const real_t Re               = scenario.Re;
    const real_t inv_Re           = 1 / Re;
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
                                   + auxp(aux::my_z, offset);
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

                // Shorten computational kernel names
                namespace orthonormal = suzerain::orthonormal;

                // Compute quantities based upon state.  Real-valued scalars
                // are declared inline.  Vector- and tensor-valued expressions
                // declared outside loop.
                u                  = orthonormal::rhome::u(
                                        rho, m);
                const real_t div_u = orthonormal::rhome::div_u(
                                        rho, grad_rho, m, div_m);
                grad_u             = orthonormal::rhome::grad_u(
                                        rho, grad_rho, m, grad_m);
                grad_div_u         = orthonormal::rhome::grad_div_u(
                                        rho, grad_rho, grad_grad_rho,
                                        m, div_m, grad_m, grad_div_m);
                div_grad_u         = orthonormal::rhome::div_grad_u(
                                        rho, grad_rho, div_grad_rho,
                                        m, grad_m, div_grad_m);
                real_t p, T, mu, lambda;
                orthonormal::rhome::p_T_mu_lambda(
                    beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
                    p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
                const real_t div_grad_p = orthonormal::rhome::div_grad_p(
                                            gamma,
                                            rho, grad_rho, div_grad_rho,
                                            m, grad_m, div_grad_m,
                                            e, grad_e, div_grad_e);
                const real_t div_grad_T = orthonormal::rhome::div_grad_T(
                                            gamma,
                                            rho, grad_rho, div_grad_rho,
                                            p, grad_p, div_grad_p);
                tau     = orthonormal::tau(
                            mu, lambda, div_u, grad_u);
                div_tau = orthonormal::div_tau(
                            mu, grad_mu, lambda, grad_lambda,
                            div_u, grad_u, div_grad_u, grad_div_u);

                // Form continuity equation right hand side
                sphys(ndx::rho, offset) = - div_m
                    ;

                // Form momentum equation right hand side
                momentum_rhs =
                    - orthonormal::div_u_outer_m(m, grad_m, u, div_u)
                    - grad_p
                    + inv_Re * div_tau
                    ;
                sphys(ndx::rhou, offset) = momentum_rhs.x();
                sphys(ndx::rhov, offset) = momentum_rhs.y();
                sphys(ndx::rhow, offset) = momentum_rhs.z();

                // Form energy equation right hand side
                sphys(ndx::rhoe, offset) =
                    - orthonormal::div_e_u(
                            e, grad_e, u, div_u
                        )
                    - orthonormal::div_p_u(
                            p, grad_p, u, div_u
                        )
                    + inv_Re_Pr_gamma1 * orthonormal::div_mu_grad_T(
                            grad_T, div_grad_T, mu, grad_mu
                        )
                    + inv_Re * orthonormal::div_tau_u<real_t>(
                            u, grad_u, tau, div_tau
                        )
                    ;

                // Maintain the minimum observed stable time step, if necessary
                if (delta_t_requested) {
                    namespace timestepper = suzerain::timestepper;
                    convective_delta_t = suzerain::math::minnan(
                            convective_delta_t,
                            timestepper::convective_stability_criterion(
                                    u.x(), one_over_delta_x,
                                    u.y(), one_over_delta_y(j),
                                    u.z(), one_over_delta_z,
                                    evmaxmag_real,
                                    std::sqrt(T)) /* nondimen a == sqrt(T) */);
                    diffusive_delta_t = suzerain::math::minnan(
                            diffusive_delta_t,
                            timestepper::diffusive_stability_criterion(
                                    one_over_delta_x,
                                    one_over_delta_y(j),
                                    one_over_delta_z,
                                    Re, Pr, gamma, evmaxmag_imag, mu / rho));
                }

            } // end X

        } // end Z

    } // end Y


    // Collectively convert state to wave space using parallel FFTs
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        dgrid.transform_physical_to_wave(&sphys(i,0));
    }

    // Perform Allreduce on stable time step sizes when necessary
    // Note delta_t_candidates aliases {convective,diffusive}_delta_t
    if (delta_t_requested) {
        SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, delta_t_candidates,
                    sizeof(delta_t_candidates)/sizeof(delta_t_candidates[0]),
                    suzerain::mpi::datatype<real_t>::value,
                    MPI_MIN, MPI_COMM_WORLD));
    }

    // Return minimum of either time step criterion, accounting for NaNs
    return suzerain::math::minnan(convective_delta_t, diffusive_delta_t);

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

NonlinearOperatorIsothermal::NonlinearOperatorIsothermal(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop)
    : NonlinearOperator(scenario, grid, dgrid, b, bop)
{
    // Precompute operator for finding bulk quantities from coefficients
    if (has_zero_zero_mode) {
        bulkcoeff.resize(b.n());
        b.integration_coefficients(0, bulkcoeff.data());
        bulkcoeff /= scenario.Ly;
    }
}

real_t NonlinearOperatorIsothermal::applyOperator(
    suzerain::NoninterleavedState<4,complex_t> &swave,
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
    const real_t delta_t = base::applyOperator(
            swave, evmaxmag_real, evmaxmag_imag, delta_t_requested);

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
    return delta_t;
}

} // end namespace channel
