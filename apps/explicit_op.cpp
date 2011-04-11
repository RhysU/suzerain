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
// explicit_op.hpp: Nonlinear operators for channel_explicit
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
#include <suzerain/multi_array.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/utility.hpp>

#include "explicit_op.hpp"

#pragma warning(disable:383 1572)

// Helper function to turn two values into a full Collection
template<typename T>
static boost::array<T,2> make_collection_of_two(const T& x, const T& y)
{
    boost::array<T,2> retval = {{ x, y }};
    return retval;
}

namespace channel {

NonlinearOperatorBase::NonlinearOperatorBase(
            const suzerain::problem::ScenarioDefinition<real_t> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop)
    : scenario(scenario),
      grid(grid),
      dgrid(dgrid),
      bop(bop),
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

NonlinearOperator::NonlinearOperator(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        const suzerain::bsplineop_luz &massluz)
    : NonlinearOperatorBase(scenario, grid, dgrid, b, bop),
      massluz(massluz),
      auxw(suzerain::to_yxz(static_cast<std::size_t>(aux::count),
                            dgrid.local_wave_extent),
              suzerain::prepend(dgrid.local_wave_storage(),
                  suzerain::strides_cm(
                      suzerain::to_yxz(dgrid.local_wave_extent)))),
      auxp(reinterpret_cast<real_t *>(auxw.origin()),
           make_collection_of_two<std::size_t>(
               aux::count, dgrid.local_physical_extent.prod()),
           make_collection_of_two<std::size_t>(
               auxw.strides()[0] * sizeof(complex_t)/sizeof(real_t), 0),
           suzerain::storage::noninterleaved<2>())
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

    // All state enters routine as coefficients in X, Y, and Z directions

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


    // Collectively convert swave to physical space using parallel FFTs
    // In physical space, we'll employ a view to reshape the four-dimensional
    // data stored (F, X, Z, Y) a 2D (F, X*Z*Y) layout.  Reducing the
    // dimensionality encourages linear access and eases indexing overhead.
    suzerain::multi_array::ref<real_t,2> sphys(
            reinterpret_cast<real_t *>(swave.origin()),
            make_collection_of_two<std::size_t>(
                swave.shape()[0],
                dgrid.local_physical_extent.prod()),
            make_collection_of_two<std::size_t>(
                swave.strides()[0] * sizeof(complex_t)/sizeof(real_t), 0),
            suzerain::storage::noninterleaved<2>());
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        dgrid.transform_wave_to_physical(sphys[i].origin());
    }
    for (std::size_t i = 0; i < aux::count; ++i) {
        dgrid.transform_wave_to_physical(auxp[i].origin());
    }

    // Retrieve constants and compute derived constants
    const real_t beta             = scenario.beta;
    const real_t gamma            = scenario.gamma;
    const real_t Pr               = scenario.Pr;
    const real_t Re               = scenario.Re;
    const real_t inv_Re           = 1 / Re;
    const real_t inv_Re_Pr_gamma1 = 1 / (Re * Pr * (gamma - 1));

    // Temporary storage used within following loop
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

    // Used to track local \frac{1}{\Delta{}y} for time criterion.
    // In physical space storage is X Z Y with Y direction slowest
    // Maintain where we are relative to local rank's starting y offset
    // FIXME Post #1492 work to simplify
    const size_t stride_y = dgrid.local_physical_extent.x()
                          * dgrid.local_physical_extent.z();
    size_t ndx_y = stride_y * dgrid.local_physical_start.y();

//  for (std::size_t j = dgrid.local_physical_start.y();
//       j < dgrid.local_physical.end.y(); ++j) {
//      const real_t yj = STARTHERE
//      SUZERAIN_UNUSED(yj);

//      for (std::size_t k = dgrid.local_physical_start.z();
//          k < dgrid.local_physical.end.z(); ++k) {
//          const real_t zk = k * scenario.Lz / grid.dN.z
//                          - scenario.Lz / 2;
//          SUZERAIN_UNUSED(zk);

//          for (std::size_t i = dgrid.local_physical_start.x();
//              i < dgrid.local_physical.end.x(); ++i) {
//              const real_t xi = i * scenario.Lx / grid.dN.x
//                              - scenario.Lx / 2;
//              SUZERAIN_UNUSED(xi);

//          } // end X

//      } // end Z

//  } // end Y

    // Walk physical space swave storage linearly
    for (// Loop initialization
         real_t *p_rho        = sphys[ndx::rho].origin(),
                *p_rho_x      = auxp[aux::rho_x].origin(),
                *p_rho_y      = auxp[aux::rho_y].origin(),
                *p_rho_z      = auxp[aux::rho_z].origin(),
                *p_rho_xx     = auxp[aux::rho_xx].origin(),
                *p_rho_xy     = auxp[aux::rho_xy].origin(),
                *p_rho_xz     = auxp[aux::rho_xz].origin(),
                *p_rho_yy     = auxp[aux::rho_yy].origin(),
                *p_rho_yz     = auxp[aux::rho_yz].origin(),
                *p_rho_zz     = auxp[aux::rho_zz].origin(),
                *p_mx         = sphys[ndx::rhou].origin(),
                *p_mx_x       = auxp[aux::mx_x].origin(),
                *p_mx_y       = auxp[aux::mx_y].origin(),
                *p_mx_z       = auxp[aux::mx_z].origin(),
                *p_mx_xx      = auxp[aux::mx_xx].origin(),
                *p_mx_xy      = auxp[aux::mx_xy].origin(),
                *p_mx_xz      = auxp[aux::mx_xz].origin(),
                *p_mx_yy      = auxp[aux::mx_yy].origin(),
                *p_mx_yz      = auxp[aux::mx_yz].origin(),
                *p_mx_zz      = auxp[aux::mx_zz].origin(),
                *p_my         = sphys[ndx::rhov].origin(),
                *p_my_x       = auxp[aux::my_x].origin(),
                *p_my_y       = auxp[aux::my_y].origin(),
                *p_my_z       = auxp[aux::my_z].origin(),
                *p_my_xx      = auxp[aux::my_xx].origin(),
                *p_my_xy      = auxp[aux::my_xy].origin(),
                *p_my_xz      = auxp[aux::my_xz].origin(),
                *p_my_yy      = auxp[aux::my_yy].origin(),
                *p_my_yz      = auxp[aux::my_yz].origin(),
                *p_my_zz      = auxp[aux::my_zz].origin(),
                *p_mz         = sphys[ndx::rhow].origin(),
                *p_mz_x       = auxp[aux::mz_x].origin(),
                *p_mz_y       = auxp[aux::mz_y].origin(),
                *p_mz_z       = auxp[aux::mz_z].origin(),
                *p_mz_xx      = auxp[aux::mz_xx].origin(),
                *p_mz_xy      = auxp[aux::mz_xy].origin(),
                *p_mz_xz      = auxp[aux::mz_xz].origin(),
                *p_mz_yy      = auxp[aux::mz_yy].origin(),
                *p_mz_yz      = auxp[aux::mz_yz].origin(),
                *p_mz_zz      = auxp[aux::mz_zz].origin(),
                *p_e          = sphys[ndx::rhoe].origin(),
                *p_e_x        = auxp[aux::e_x].origin(),
                *p_e_y        = auxp[aux::e_y].origin(),
                *p_e_z        = auxp[aux::e_z].origin(),
                *p_div_grad_e = auxp[aux::div_grad_e].origin();
         // Loop completion test
         p_rho != sphys[ndx::rho].origin() + dgrid.local_physical_extent.prod();
         // Loop increment
         ++p_rho,
         ++p_rho_x,
         ++p_rho_y,
         ++p_rho_z,
         ++p_rho_xx,
         ++p_rho_xy,
         ++p_rho_xz,
         ++p_rho_yy,
         ++p_rho_yz,
         ++p_rho_zz,
         ++p_mx,
         ++p_mx_x,
         ++p_mx_y,
         ++p_mx_z,
         ++p_mx_xx,
         ++p_mx_xy,
         ++p_mx_xz,
         ++p_mx_yy,
         ++p_mx_yz,
         ++p_mx_zz,
         ++p_my,
         ++p_my_x,
         ++p_my_y,
         ++p_my_z,
         ++p_my_xx,
         ++p_my_xy,
         ++p_my_xz,
         ++p_my_yy,
         ++p_my_yz,
         ++p_my_zz,
         ++p_mz,
         ++p_mz_x,
         ++p_mz_y,
         ++p_mz_z,
         ++p_mz_xx,
         ++p_mz_xy,
         ++p_mz_xz,
         ++p_mz_yy,
         ++p_mz_yz,
         ++p_mz_zz,
         ++p_e,
         ++p_e_x,
         ++p_e_y,
         ++p_e_z,
         ++p_div_grad_e,
         ++ndx_y) {

        // Prepare local density-related quantities
        const real_t rho          = *p_rho;
        grad_rho[0]               = *p_rho_x;
        grad_rho[1]               = *p_rho_y;
        grad_rho[2]               = *p_rho_z;
        const real_t div_grad_rho = *p_rho_xx + *p_rho_yy + *p_rho_zz;
        grad_grad_rho(0,0)        = *p_rho_xx;
        grad_grad_rho(0,1)        = *p_rho_xy;
        grad_grad_rho(0,2)        = *p_rho_xz;
        grad_grad_rho(1,0)        = grad_grad_rho(0,1);
        grad_grad_rho(1,1)        = *p_rho_yy;
        grad_grad_rho(1,2)        = *p_rho_yz;
        grad_grad_rho(2,0)        = grad_grad_rho(0,2);
        grad_grad_rho(2,1)        = grad_grad_rho(1,2);
        grad_grad_rho(2,2)        = *p_rho_zz;

        // Prepare local momentum-related quantities
        m[0]               = *p_mx;
        m[1]               = *p_my;
        m[2]               = *p_mz;
        const real_t div_m = *p_mx_x + *p_my_y + *p_my_z;
        grad_m(0,0)        = *p_mx_x;
        grad_m(0,1)        = *p_mx_y;
        grad_m(0,2)        = *p_mx_z;
        grad_m(1,0)        = *p_my_x;
        grad_m(1,1)        = *p_my_y;
        grad_m(1,2)        = *p_my_z;
        grad_m(2,0)        = *p_mz_x;
        grad_m(2,1)        = *p_mz_y;
        grad_m(2,2)        = *p_mz_z;
        div_grad_m[0]      = *p_mx_xx + *p_mx_yy + *p_mx_zz;
        div_grad_m[1]      = *p_my_xx + *p_my_yy + *p_my_zz;
        div_grad_m[2]      = *p_mz_xx + *p_mz_yy + *p_mz_zz;
        grad_div_m[0]      = *p_mx_xx + *p_my_xy + *p_mz_xz;
        grad_div_m[1]      = *p_mx_xy + *p_my_yy + *p_mz_yz;
        grad_div_m[2]      = *p_mx_xz + *p_my_yz + *p_mz_zz;

        // Prepare local total energy-related quantities
        const real_t e          = *p_e;
        grad_e[0]               = *p_e_x;
        grad_e[1]               = *p_e_y;
        grad_e[2]               = *p_e_z;

        // Prepare quantities derived from local state and its derivatives
        u                  = suzerain::orthonormal::rhome::u(rho, m);
        const real_t div_u = suzerain::orthonormal::rhome::div_u(
                                rho, grad_rho, m, div_m);
        grad_u             = suzerain::orthonormal::rhome::grad_u(
                                rho, grad_rho, m, grad_m);
        grad_div_u         = suzerain::orthonormal::rhome::grad_div_u(
                                rho, grad_rho, grad_grad_rho,
                                m, div_m, grad_m, grad_div_m);
        div_grad_u         = suzerain::orthonormal::rhome::div_grad_u(
                                rho, grad_rho, div_grad_rho,
                                m, grad_m, div_grad_m);
        real_t p, T, mu, lambda;
        suzerain::orthonormal::rhome::p_T_mu_lambda(
            beta, gamma, rho, grad_rho, m, grad_m, e, grad_e,
            p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
        const real_t div_grad_p = suzerain::orthonormal::rhome::div_grad_p(
                                    gamma,
                                    rho, grad_rho, div_grad_rho,
                                    m, grad_m, div_grad_m,
                                    e, grad_e, *p_div_grad_e);
        const real_t div_grad_T = suzerain::orthonormal::rhome::div_grad_T(
                                    gamma,
                                    rho, grad_rho, div_grad_rho,
                                    p, grad_p, div_grad_p);
        tau     = suzerain::orthonormal::tau(
                    mu, lambda, div_u, grad_u);
        div_tau = suzerain::orthonormal::div_tau(
                    mu, grad_mu, lambda, grad_lambda,
                    div_u, grad_u, div_grad_u, grad_div_u);

        // Continuity equation right hand side
        *p_rho = - div_m
            ;

        // Momentum equation right hand side
        Eigen::Vector3r momentum =
            - suzerain::orthonormal::div_u_outer_m(m, grad_m, u, div_u)
            - grad_p
            + inv_Re * div_tau
            ;
        *p_mx = momentum[0];
        *p_my = momentum[1];
        *p_mz = momentum[2];

        // Energy equation right hand side
        *p_e = - suzerain::orthonormal::div_e_u(
                    e, grad_e, u, div_u
                 )
               - suzerain::orthonormal::div_p_u(
                    p, grad_p, u, div_u
                 )
               + inv_Re_Pr_gamma1 * suzerain::orthonormal::div_mu_grad_T(
                    grad_T, div_grad_T, mu, grad_mu
                 )
               + inv_Re * suzerain::orthonormal::div_tau_u<real_t>(
                    u, grad_u, tau, div_tau
                 )
               ;

        // Maintain the minimum observed stable time step
        convective_delta_t = suzerain::math::minnan(
                convective_delta_t,
                suzerain::timestepper::convective_stability_criterion(
                        u.x(), one_over_delta_x,
                        u.y(), one_over_delta_y(ndx_y / stride_y),
                        u.z(), one_over_delta_z,
                        evmaxmag_real,
                        std::sqrt(T)) /* nondimensional a = sqrt(T) */);
        diffusive_delta_t = suzerain::math::minnan(
                diffusive_delta_t,
                suzerain::timestepper::diffusive_stability_criterion(
                        one_over_delta_x,
                        one_over_delta_y(ndx_y / stride_y),
                        one_over_delta_z,
                        Re, Pr, gamma, evmaxmag_imag, mu / rho));
    }

    // Collectively convert state to wave space using parallel FFTs
    for (std::size_t i = 0; i < channel::field::count; ++i) {
        dgrid.transform_physical_to_wave(sphys[i].origin());
    }

    // Convert collocation point values to Bspline coefficients
    massluz.solve(swave.shape()[2]*swave.shape()[3],
                  swave[ndx::rho].origin(), 1, swave.shape()[1]);
    massluz.solve(swave.shape()[2]*swave.shape()[3],
                  swave[ndx::rhou].origin(), 1, swave.shape()[1]);
    massluz.solve(swave.shape()[2]*swave.shape()[3],
                  swave[ndx::rhov].origin(), 1, swave.shape()[1]);
    massluz.solve(swave.shape()[2]*swave.shape()[3],
                  swave[ndx::rhow].origin(), 1, swave.shape()[1]);
    massluz.solve(swave.shape()[2]*swave.shape()[3],
                  swave[ndx::rhoe].origin(), 1, swave.shape()[1]);

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

    // All state leaves routine as coefficients in X, Y, and Z directions
}

NonlinearOperatorWithBoundaryConditions::NonlinearOperatorWithBoundaryConditions(
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::pencil_grid &dgrid,
        suzerain::bspline &b,
        const suzerain::bsplineop &bop,
        const suzerain::bsplineop_luz &massluz)
    : NonlinearOperator(scenario, grid, dgrid, b, bop, massluz)
{
    // Obtain integration coefficients
    bintcoeff.resize(b.n(), 1);
    b.integration_coefficients(0, bintcoeff.data());
}

real_t NonlinearOperatorWithBoundaryConditions::applyOperator(
    suzerain::NoninterleavedState<4,complex_t> &swave,
    const real_t evmaxmag_real,
    const real_t evmaxmag_imag,
    const bool delta_t_requested) const
{
    namespace ndx = channel::field::ndx;

    // Special handling occurs only on rank holding the "zero-zero" mode
    const bool zero_zero_rank =    (dgrid.local_wave_start.x() == 0)
                                && (dgrid.local_wave_start.z() == 0);

    // Precompute and store quantities necessary for BC implementation
    Eigen::ArrayXr original_state_mx;
    real_t bulk_density = std::numeric_limits<real_t>::quiet_NaN();
    if (zero_zero_rank) {
        // Save a copy of the constant x-momentum modes
        original_state_mx.resize(swave.shape()[1]);
        for (std::size_t i = 0; i < swave.shape()[1]; ++i) {
            original_state_mx[i]
                = suzerain::complex::real(swave[ndx::rhou][i][0][0]);
        }

        // Compute the bulk density so we can hold it constant in time
        bulk_density = suzerain::blas::dot(
                bintcoeff.size(), bintcoeff.data(), 1,
                reinterpret_cast<real_t *>(&swave[ndx::rho][0][0][0]),
                sizeof(complex_t)/sizeof(real_t));
    }

    // Apply an operator that cares nothing about the boundaries
    const real_t delta_t = base::applyOperator(
            swave, evmaxmag_real, evmaxmag_imag, delta_t_requested);

    // Indices that will be useful as shorthand
    const std::size_t lower_wall = 0;                    // index of wall
    const std::size_t upper_wall = swave.shape()[1] - 1; // index of wall

    // Add f_rho to mean density per writeup step (2)
    if (zero_zero_rank) {
        const real_t f_rho = suzerain::blas::dot(
                bintcoeff.size(), bintcoeff.data(), 1,
                reinterpret_cast<real_t *>(&swave[ndx::rho][0][0][0]),
                sizeof(complex_t)/sizeof(real_t)) / scenario.Ly;
        for (std::size_t i = 0; i < swave.shape()[1]; ++i) {
            swave[ndx::rho][i][0][0] -= f_rho;
        }
    }

    // Set no-slip condition for momentum on walls per writeup step (3)
    assert(static_cast<int>(ndx::rhov) == static_cast<int>(ndx::rhou) + 1);
    assert(static_cast<int>(ndx::rhow) == static_cast<int>(ndx::rhov) + 1);
    for (std::size_t i = ndx::rhou; i <= ndx::rhow; ++i) {
        for (std::size_t k = 0; k < swave.shape()[3]; ++k) {
            for (std::size_t j = 0; j < swave.shape()[2]; ++j) {
                swave[i][lower_wall][j][k] = 0;
                swave[i][upper_wall][j][k] = 0;
            }
        }
    }

    // Set isothermal condition on walls per writeup step (4)
    const real_t inv_gamma_gamma1
        = 1 / (scenario.gamma * (scenario.gamma - 1));
    for (std::size_t k = 0; k < swave.shape()[3]; ++k) {
        for (std::size_t j = 0; j < swave.shape()[2]; ++j) {
            swave[ndx::rhoe][lower_wall][j][k]
                = inv_gamma_gamma1 * swave[ndx::rho][lower_wall][j][k];
            swave[ndx::rhoe][upper_wall][j][k]
                = inv_gamma_gamma1 * swave[ndx::rho][upper_wall][j][k];
        }
    }

    // Apply f_{m_x} to mean x-momentum, mean energy
    if (zero_zero_rank) {

        // Compute temporary per writeup implementation step (5)
        real_t alpha = suzerain::blas::dot(
                bintcoeff.size(), bintcoeff.data(), 1,
                reinterpret_cast<real_t *>(&swave[ndx::rhou][0][0][0]),
                sizeof(complex_t)/sizeof(real_t)) / scenario.Ly;

        // Apply to non-wall mean x-momentum right hand side per step (6)
        for (std::size_t i = lower_wall + 1; i < upper_wall; ++i) {
            swave[ndx::rhou][i][0][0] -= alpha;
        }

        // Apply to non-wall mean energy right hand side per step (7)
        alpha /= bulk_density;
        for (std::size_t i = lower_wall + 1; i < upper_wall; ++i) {
            swave[ndx::rhoe][i][0][0] -= alpha * original_state_mx[i];
        }
    }

    // Return the time step found by the BC-agnostic operator
    return delta_t;
}

} // end namespace channel
