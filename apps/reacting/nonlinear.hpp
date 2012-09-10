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
// nonlinear.hpp: Building blocks for nonlinear Navier--Stokes operators
// $Id$

#ifndef NONLINEAR_HPP
#define NONLINEAR_HPP

#include "nonlinear_fwd.hpp"

#include <suzerain/error.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/state.hpp>

#include "../support.hpp"

#pragma warning(disable:383 1572)

namespace channel {

template<bool ZerothSubstep,
         linearize::type Linearize,
         class ManufacturedSolution>
std::vector<real_t> applyNonlinearOperator(
            const suzerain::OperatorBase<real_t> &o,
            OperatorCommonBlock &common,
            const boost::shared_ptr<const ManufacturedSolution>& msoln,
            const real_t time,
            suzerain::ContiguousState<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag)
{
    GRVY_TIMER_BEGIN("applyNonlinearOperator");

    // Shorthand
    typedef suzerain::ContiguousState<4,complex_t> state_type;
    namespace ndx = channel::field::ndx;
    using Eigen::Vector3r;
    using Eigen::Matrix3r;
    using std::size_t;

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
    typename boost::scoped_ptr<state_type> _auxw_ptr(
            allocate_padded_state<state_type>(aux::count, o.dgrid)); // RAII
    state_type &auxw = *_auxw_ptr;                                   // Brevity

    // Sanity check incoming swave's and auxw's shape and contiguity
    assert(swave.shape()[0] == channel::field::count);
    assert(swave.shape()[1] == (unsigned) o.dgrid.local_wave_extent.y());
    assert(swave.shape()[2] == (unsigned) o.dgrid.local_wave_extent.x());
    assert(swave.shape()[3] == (unsigned) o.dgrid.local_wave_extent.z());
    assert((unsigned) swave.strides()[1] == 1u);
    assert((unsigned) swave.strides()[2] == swave.shape()[1]);
    assert((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    assert(std::equal(swave.shape() + 1, swave.shape() + 4,
                      auxw.shape() + 1));
    assert(std::equal(swave.strides() + 1, swave.strides() + 4,
                      auxw.strides() + 1));

    // Prepare common-block-like storage used to pass details from N to L.
    // Zeroing is done carefully as accumulated means and reference quantities
    // must survive from nonzero substep to substep while instant profiles do not.
    if (ZerothSubstep) common.setZero(/* Ny */ swave.shape()[1]);
    common.u().setZero();

    // Maintain stable time step values to return to the caller
    boost::array<real_t, 2> delta_t_candidates = {{
            std::numeric_limits<real_t>::max(),
            std::numeric_limits<real_t>::max()
    }};
    real_t &convective_delta_t = delta_t_candidates[0];
    real_t &diffusive_delta_t  = delta_t_candidates[1];

    // GRVY_TIMER_{BEGIN,END} pairs for differentiation done in OperatorBase

    // Compute Y derivatives of density at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.diffwave_apply(0, 0, 1, swave, ndx::rho);
    o.bop_accumulate(1,    1, swave, ndx::rho, 0, auxw, aux::rho_y);
    o.bop_accumulate(2,    1, swave, ndx::rho, 0, auxw, aux::rho_yy);
    o.bop_apply     (0,    1, swave, ndx::rho);

    // Compute X- and Z- derivatives of density at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::rho,   0, auxw, aux::rho_x );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::rho,   0, auxw, aux::rho_xx);
    o.diffwave_accumulate(1, 1, 1, swave, ndx::rho,   0, auxw, aux::rho_xz);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::rho,   0, auxw, aux::rho_z );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::rho,   0, auxw, aux::rho_zz);
    o.diffwave_accumulate(1, 0, 1, auxw,  aux::rho_y, 0, auxw, aux::rho_xy);
    o.diffwave_accumulate(0, 1, 1, auxw,  aux::rho_y, 0, auxw, aux::rho_yz);

    // Compute Y derivatives of X momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.diffwave_apply(0, 0, 1, swave, ndx::rhou);
    o.bop_accumulate(1,    1, swave, ndx::rhou, 0, auxw, aux::mx_y);
    o.bop_accumulate(2,    1, swave, ndx::rhou, 0, auxw, aux::mx_yy);
    o.bop_apply     (0,    1, swave, ndx::rhou);

    // Compute X- and Z- derivatives of X momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::rhou,  0, auxw, aux::mx_x );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::rhou,  0, auxw, aux::mx_xx);
    o.diffwave_accumulate(1, 1, 1, swave, ndx::rhou,  0, auxw, aux::mx_xz);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::rhou,  0, auxw, aux::mx_z );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::rhou,  0, auxw, aux::mx_zz);
    o.diffwave_accumulate(1, 0, 1, auxw,  aux::mx_y,  0, auxw, aux::mx_xy);
    o.diffwave_accumulate(0, 1, 1, auxw,  aux::mx_y,  0, auxw, aux::mx_yz);

    // Compute Y derivatives of Y momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.diffwave_apply(0, 0, 1, swave, ndx::rhov);
    o.bop_accumulate(1,    1, swave, ndx::rhov, 0, auxw, aux::my_y);
    o.bop_accumulate(2,    1, swave, ndx::rhov, 0, auxw, aux::my_yy);
    o.bop_apply     (0,    1, swave, ndx::rhov);

    // Compute X- and Z- derivatives of Y momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::rhov,  0, auxw, aux::my_x );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::rhov,  0, auxw, aux::my_xx);
    o.diffwave_accumulate(1, 1, 1, swave, ndx::rhov,  0, auxw, aux::my_xz);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::rhov,  0, auxw, aux::my_z );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::rhov,  0, auxw, aux::my_zz);
    o.diffwave_accumulate(1, 0, 1, auxw,  aux::my_y,  0, auxw, aux::my_xy);
    o.diffwave_accumulate(0, 1, 1, auxw,  aux::my_y,  0, auxw, aux::my_yz);

    // Compute Y derivatives of Z momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.diffwave_apply(0, 0, 1, swave, ndx::rhow);
    o.bop_accumulate(1,    1, swave, ndx::rhow, 0, auxw, aux::mz_y);
    o.bop_accumulate(2,    1, swave, ndx::rhow, 0, auxw, aux::mz_yy);
    o.bop_apply     (0,    1, swave, ndx::rhow);

    // Compute X- and Z- derivatives of Z momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::rhow,  0, auxw, aux::mz_x );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::rhow,  0, auxw, aux::mz_xx);
    o.diffwave_accumulate(1, 1, 1, swave, ndx::rhow,  0, auxw, aux::mz_xz);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::rhow,  0, auxw, aux::mz_z );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::rhow,  0, auxw, aux::mz_zz);
    o.diffwave_accumulate(1, 0, 1, auxw,  aux::mz_y,  0, auxw, aux::mz_xy);
    o.diffwave_accumulate(0, 1, 1, auxw,  aux::mz_y,  0, auxw, aux::mz_yz);

    // Compute Y derivatives of total energy at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.diffwave_apply(0, 0, 1, swave, ndx::rhoe);
    o.bop_accumulate(1,    1, swave, ndx::rhoe, 0, auxw, aux::e_y);
    o.bop_accumulate(2,    1, swave, ndx::rhoe, 0, auxw, aux::div_grad_e);
    o.bop_apply     (0,    1, swave, ndx::rhoe);

    // Compute X- and Z- derivatives of total energy at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::rhoe, 0, auxw, aux::e_x       );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::rhoe, 1, auxw, aux::div_grad_e);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::rhoe, 0, auxw, aux::e_z       );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::rhoe, 1, auxw, aux::div_grad_e);

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    typename channel::physical_view<aux::count>::type auxp
        = channel::physical_view<aux::count>::create(o.dgrid, auxw);
    typename channel::physical_view<channel::field::count>::type sphys
        = channel::physical_view<channel::field::count>::create(o.dgrid, swave);
    for (size_t i = 0; i < channel::field::count; ++i) {
        GRVY_TIMER_BEGIN("transform_wave_to_physical");
        o.dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
        GRVY_TIMER_END("transform_wave_to_physical");
    }
    for (size_t i = 0; i < aux::count; ++i) {
        GRVY_TIMER_BEGIN("transform_wave_to_physical");
        o.dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
        GRVY_TIMER_END("transform_wave_to_physical");
    }

    // Retrieve constants and compute derived constants before inner loops
    const real_t alpha            = o.scenario.alpha;
    const real_t alpha13          = alpha + real_t(1)/real_t(3);
    const real_t beta             = o.scenario.beta;
    const real_t gamma            = o.scenario.gamma;
    const real_t Ma               = o.scenario.Ma;
    const real_t Pr               = o.scenario.Pr;
    const real_t Re               = o.scenario.Re;
    const real_t inv_Re           = 1 / Re;
    const real_t inv_Ma2          = 1 / (Ma * Ma);
    const real_t Ma2_over_Re      = (Ma * Ma) / Re;
    const real_t inv_Re_Pr_gamma1 = 1 / (Re * Pr * (gamma - 1));
    const real_t lambda1_x        = o.lambda1_x;
    const real_t lambda1_z        = o.lambda1_z;
    const real_t lambda2_x        = o.lambda2_x;
    const real_t lambda2_z        = o.lambda2_z;

    // Type of Boost.Accumulator to use for summation processes.
    // Kahan summation preferred when available as incremental cost is small
    // and we will add many small numbers to a large magnitude sum.
    // During debugging, also make the number of samples available.
    typedef boost::accumulators::accumulator_set<
                real_t,
                boost::accumulators::stats<
#if BOOST_VERSION >= 104700
                    boost::accumulators::tag::sum_kahan
#else
                    boost::accumulators::tag::sum
#endif
#ifndef NDEBUG
                    , boost::accumulators::tag::count
#endif
                >
            > summing_accumulator_type;

    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
    //
    // Three traversals occur:
    // (1) Computing reference quantities and mean velocity OR mean quantities
    //     (depending on linearization and which substep is being performed).
    // (2) Computing the nonlinear equation right hand sides.
    // (3) Computing any manufactured solution forcing (when enabled).
    //
    // The traversal pattern is best embodied by the third pass.  The first and
    // second passes use slightly more compact structure as only y(j) must be
    // known within them.  Study (3) and convince yourself that (1) and (2) are
    // equivalent but lack information on x(i) and z(k).

    // Traversal:
    // (1) Computing reference quantities and mean velocity OR mean velocity
    //     (depending on linearization and which substep is being performed).
    //
    // Always gather all reference quantities when using implicit solves.
    // Profiling indicates the overhead is tiny and it keeps the code readable.
    if (    ZerothSubstep
         && Linearize != linearize::none) {  // References and mean velocity

        GRVY_TIMER_BEGIN("reference quantities");

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t leftNotOnRank = o.dgrid.local_physical_start.y();
        if (leftNotOnRank) common.refs.leftCols(leftNotOnRank).setZero();

        // Sum reference quantities as a function of y(j) into common.ref_*
        // See writeups/derivation.tex or rholut_imexop.h for definitions
        size_t offset = 0;
        for (int j = o.dgrid.local_physical_start.y();
            j < o.dgrid.local_physical_end.y();
            ++j) {

            // Prepare logical indices using a struct for scoping (e.g. ref::ux).
            struct ref { enum { ux, uy, uz, u2,
                                uxux, uxuy, uxuz, uyuy, uyuz, uzuz,
                                nu, nuux, nuuy, nuuz, nuu2,
                                nuuxux, nuuxuy, nuuxuz, nuuyuy, nuuyuz, nuuzuz,
                                ex_gradrho, ey_gradrho, ez_gradrho,
                                e_divm, e_deltarho,
                                count // Sentry
            }; };

            // An array of summing_accumulator_type holds all running sums.
            // This gives nicer construction and allows looping over results.
            summing_accumulator_type acc[ref::count];

            const size_t last_zxoffset = offset
                                       + o.dgrid.local_physical_extent.z()
                                       * o.dgrid.local_physical_extent.x();
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

                // ...including simple velocity-related quantities...
                const Vector3r u = suzerain::rholut::u(rho, m);
                acc[ref::ux](u.x());
                acc[ref::uy](u.y());
                acc[ref::uz](u.z());
                acc[ref::u2](u.squaredNorm());
                acc[ref::uxux](u.x()*u.x());
                acc[ref::uxuy](u.x()*u.y());
                acc[ref::uxuz](u.x()*u.z());
                acc[ref::uyuy](u.y()*u.y());
                acc[ref::uyuz](u.y()*u.z());
                acc[ref::uzuz](u.z()*u.z());

                // ...including simple viscosity-related quantities...
                const real_t nu = mu / rho;
                acc[ref::nu](nu);
                acc[ref::nuux](nu*u.x());
                acc[ref::nuuy](nu*u.y());
                acc[ref::nuuz](nu*u.z());
                acc[ref::nuu2](nu*u.squaredNorm());
                acc[ref::nuuxux](nu*u.x()*u.x());
                acc[ref::nuuxuy](nu*u.x()*u.y());
                acc[ref::nuuxuz](nu*u.x()*u.z());
                acc[ref::nuuyuy](nu*u.y()*u.y());
                acc[ref::nuuyuz](nu*u.y()*u.z());
                acc[ref::nuuzuz](nu*u.z()*u.z());

                // ...and other, more complicated expressions.
                namespace rholut = suzerain::rholut;
                const Vector3r e_gradrho
                        = rholut::explicit_div_e_plus_p_u_refcoeff_grad_rho(
                                gamma, rho, m, e, p);
                acc[ref::ex_gradrho](e_gradrho.x());
                acc[ref::ey_gradrho](e_gradrho.y());
                acc[ref::ez_gradrho](e_gradrho.z());

                acc[ref::e_divm](
                        rholut::explicit_div_e_plus_p_u_refcoeff_div_m(
                            rho, e, p));

                acc[ref::e_deltarho](
                        rholut::explicit_mu_div_grad_T_refcoeff_div_grad_rho(
                            gamma, mu, rho, e, p));

            } // end X // end Z

#ifndef NDEBUG
            // Ensure that all accumulators saw a consistent number of samples
            const size_t expected = boost::accumulators::count(acc[0]);
            for (size_t k = 1; k < sizeof(acc)/sizeof(acc[0]); ++k) {
                const size_t observed = boost::accumulators::count(acc[k]);
                assert(expected == observed);
            }
#endif

            // Store sums into common block in preparation for MPI Allreduce
            using boost::accumulators::sum;
            common.ref_ux        ()[j] = sum(acc[ref::ux        ]);
            common.ref_uy        ()[j] = sum(acc[ref::uy        ]);
            common.ref_uz        ()[j] = sum(acc[ref::uz        ]);
            common.ref_u2        ()[j] = sum(acc[ref::u2        ]);
            common.ref_uxux      ()[j] = sum(acc[ref::uxux      ]);
            common.ref_uxuy      ()[j] = sum(acc[ref::uxuy      ]);
            common.ref_uxuz      ()[j] = sum(acc[ref::uxuz      ]);
            common.ref_uyuy      ()[j] = sum(acc[ref::uyuy      ]);
            common.ref_uyuz      ()[j] = sum(acc[ref::uyuz      ]);
            common.ref_uzuz      ()[j] = sum(acc[ref::uzuz      ]);
            common.ref_nu        ()[j] = sum(acc[ref::nu        ]);
            common.ref_nuux      ()[j] = sum(acc[ref::nuux      ]);
            common.ref_nuuy      ()[j] = sum(acc[ref::nuuy      ]);
            common.ref_nuuz      ()[j] = sum(acc[ref::nuuz      ]);
            common.ref_nuu2      ()[j] = sum(acc[ref::nuu2      ]);
            common.ref_nuuxux    ()[j] = sum(acc[ref::nuuxux    ]);
            common.ref_nuuxuy    ()[j] = sum(acc[ref::nuuxuy    ]);
            common.ref_nuuxuz    ()[j] = sum(acc[ref::nuuxuz    ]);
            common.ref_nuuyuy    ()[j] = sum(acc[ref::nuuyuy    ]);
            common.ref_nuuyuz    ()[j] = sum(acc[ref::nuuyuz    ]);
            common.ref_nuuzuz    ()[j] = sum(acc[ref::nuuzuz    ]);
            common.ref_ex_gradrho()[j] = sum(acc[ref::ex_gradrho]);
            common.ref_ey_gradrho()[j] = sum(acc[ref::ey_gradrho]);
            common.ref_ez_gradrho()[j] = sum(acc[ref::ez_gradrho]);
            common.ref_e_divm    ()[j] = sum(acc[ref::e_divm    ]);
            common.ref_e_deltarho()[j] = sum(acc[ref::e_deltarho]);

        } // end Y

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t rightNotOnRank = common.refs.cols()
                                    - o.dgrid.local_physical_end.y();
        if (rightNotOnRank) common.refs.rightCols(rightNotOnRank).setZero();

        // Allreduce and scale common.refs sums to obtain means on all ranks
        // Allreduce mandatory as all ranks need references for linearization
        SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, common.refs.data(),
                common.refs.size(), suzerain::mpi::datatype<real_t>::value,
                MPI_SUM, MPI_COMM_WORLD));
        common.refs /= (   o.dgrid.global_physical_extent.x()
                         * o.dgrid.global_physical_extent.z());

        // Copy mean streamwise velocity information into common.u()
        common.u() = common.ref_ux();

        GRVY_TIMER_END("reference quantities");

    } else {                                 // Mean velocity profile only

        GRVY_TIMER_BEGIN("mean velocity profile");

        // Zero all reference quantities on fully-explicit zeroth substep
        if (ZerothSubstep && Linearize == linearize::none) {
            common.refs.setZero();
        }

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t topNotOnRank = o.dgrid.local_physical_start.y();
        if (topNotOnRank) common.u().topRows(topNotOnRank).setZero();

        // Sum streamwise velocities as a function of y(j) into common.u()
        size_t offset = 0;
        for (int j = o.dgrid.local_physical_start.y();
            j < o.dgrid.local_physical_end.y();
            ++j) {

            summing_accumulator_type ux;

            const size_t last_zxoffset = offset
                                       + o.dgrid.local_physical_extent.z()
                                       * o.dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {
                ux(sphys(ndx::rhou, offset)/sphys(ndx::rho, offset));
            } // end X // end Z

            // Store sum into common block in preparation for MPI Reduce
            common.u()[j] = boost::accumulators::sum(ux);

        } // end Y

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t bottomNotOnRank = common.u().rows()
                                     - o.dgrid.local_physical_end.y();
        if (bottomNotOnRank) common.u().bottomRows(bottomNotOnRank).setZero();

        // Reduce and scale common.u() sums to obtain mean on zero-zero rank
        // Only zero-zero rank needs the information so Reduce is sufficient
        if (o.dgrid.has_zero_zero_modes()) {
            SUZERAIN_MPICHKR(MPI_Reduce(MPI_IN_PLACE, common.u().data(),
                    common.u().size(), suzerain::mpi::datatype<real_t>::value,
                    MPI_SUM, o.dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
            common.u() /= (   o.dgrid.global_physical_extent.x()
                            * o.dgrid.global_physical_extent.z());
        } else {
            Eigen::ArrayXr tmp;
            tmp.resizeLike(common.u());
            tmp.setZero();
            SUZERAIN_MPICHKR(MPI_Reduce(common.u().data(), tmp.data(),
                    common.u().size(), suzerain::mpi::datatype<real_t>::value,
                    MPI_SUM, o.dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
        }

        GRVY_TIMER_END("mean velocity profile");

    }

    // Traversal:
    // (2) Computing the nonlinear equation right hand sides.
    GRVY_TIMER_BEGIN("nonlinear right hand sides");
    size_t offset = 0;
    for (int j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {

        // Wall-normal operator eigenvalue estimates depend on location
        const real_t lambda1_y = o.lambda1_y(j);
        const real_t lambda2_y = o.lambda2_y(j);

        // Unpack appropriate wall-normal reference quantities
        const Vector3r ref_u              (common.ref_ux        ()[j],
                                           common.ref_uy        ()[j],
                                           common.ref_uz        ()[j]);
        const real_t   ref_u2             (common.ref_u2        ()[j]);
        const Matrix3r ref_uu;
        const_cast<Matrix3r&>(ref_uu) <<   common.ref_uxux      ()[j],
                                           common.ref_uxuy      ()[j],
                                           common.ref_uxuz      ()[j],
                                           common.ref_uxuy      ()[j],
                                           common.ref_uyuy      ()[j],
                                           common.ref_uyuz      ()[j],
                                           common.ref_uxuz      ()[j],
                                           common.ref_uyuz      ()[j],
                                           common.ref_uzuz      ()[j];
        const real_t   ref_nu             (common.ref_nu        ()[j]);
        const Vector3r ref_nuu            (common.ref_nuux      ()[j],
                                           common.ref_nuuy      ()[j],
                                           common.ref_nuuz      ()[j]);
        const real_t   ref_nuu2           (common.ref_nuu2      ()[j]);
        const Matrix3r ref_nuuu;
        const_cast<Matrix3r&>(ref_nuuu) << common.ref_nuuxux    ()[j],
                                           common.ref_nuuxuy    ()[j],
                                           common.ref_nuuxuz    ()[j],
                                           common.ref_nuuxuy    ()[j],
                                           common.ref_nuuyuy    ()[j],
                                           common.ref_nuuyuz    ()[j],
                                           common.ref_nuuxuz    ()[j],
                                           common.ref_nuuyuz    ()[j],
                                           common.ref_nuuzuz    ()[j];
        const Vector3r ref_e_gradrho      (common.ref_ex_gradrho()[j],
                                           common.ref_ey_gradrho()[j],
                                           common.ref_ez_gradrho()[j]);
        const real_t   ref_e_divm         (common.ref_e_divm    ()[j]);
        const real_t   ref_e_deltarho     (common.ref_e_deltarho()[j]);

        // Iterate across the j-th ZX plane
        const size_t last_zxoffset = offset
                                   + o.dgrid.local_physical_extent.z()
                                   * o.dgrid.local_physical_extent.x();
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

            // FORM CONTINUITY EQUATION RIGHT HAND SIDE
            //
            // Implicit continuity equation handling requires zeroing RHS in
            // anticipation of possible manufactured solution forcing.  See
            // subsequent transform_physical_to_wave if you monkey around here.
            switch (Linearize) {
                case linearize::none:
                    sphys(ndx::rho, offset) = - div_m; // Explicit convection
                    break;
                case linearize::rhome:
                    sphys(ndx::rho, offset) = 0;       // Implicit convection
                    break;
            }

            // FORM MOMENTUM EQUATION RIGHT HAND SIDE
            Vector3r momentum_rhs =
                // Explicit viscous term
                  inv_Re * div_tau
                ;
            switch (Linearize) {
                case linearize::none:
                    momentum_rhs +=
                        // Explicit convective term
                        - suzerain::rholut::div_u_outer_m(m, grad_m, u, div_u)
                        // Explicit pressure term
                        - inv_Ma2 * grad_p
                        ;
                    break;
                case linearize::rhome:
                    momentum_rhs +=
                        // Explicit convective term less implicit portion
                        - suzerain::rholut::explicit_div_rho_inverse_m_outer_m(
                                grad_rho, div_m, grad_m, u, ref_u, ref_uu)
                        // Explicit pressure less implicit pressure terms
                        - inv_Ma2 * suzerain::rholut::explicit_grad_p(
                                gamma, Ma, rho, grad_rho, m, grad_m,
                                ref_u2, ref_u)
                        // Subtract implicit portions of viscous terms per
                        // suzerain::rholut::explicit_mu_div_grad_u and
                        // suzerain::rholut::explicit_mu_plus_lambda_grad_div_u
                        - inv_Re * (
                            ref_nu*(div_grad_m + alpha13*grad_div_m)
                          - ref_nuu*div_grad_rho
                          - alpha13*grad_grad_rho*ref_nuu
                        )
                        ;
                    break;
            }
            sphys(ndx::rhou, offset) = momentum_rhs.x();
            sphys(ndx::rhov, offset) = momentum_rhs.y();
            sphys(ndx::rhow, offset) = momentum_rhs.z();

            // FORM ENERGY EQUATION RIGHT HAND SIDE
            sphys(ndx::rhoe, offset) =
                // Explicit viscous work term
                + Ma2_over_Re * suzerain::rholut::div_tau_u<real_t>(
                        u, grad_u, tau, div_tau
                    )
                ;
            switch (Linearize) {
                case linearize::none:
                    sphys(ndx::rhoe, offset) +=
                        // Explicit convective and acoustic terms
                        - suzerain::rholut::div_e_u(
                                e, grad_e, u, div_u
                            )
                        - suzerain::rholut::div_p_u(
                                p, grad_p, u, div_u
                            )
                        // Explicit energy diffusion terms
                        + inv_Re_Pr_gamma1 * suzerain::rholut::div_mu_grad_T(
                                grad_T, div_grad_T, mu, grad_mu
                            )
                        ;
                        // No need to adjust explicit viscous work term
                    break;
                case linearize::rhome:
                    sphys(ndx::rhoe, offset) +=
                        // Explicit convective/acoustic less implicit portion
                        - suzerain::rholut::explicit_div_e_plus_p_u(
                                gamma, Ma, rho, grad_rho,
                                m, div_m, grad_m, e, grad_e, p,
                                ref_e_divm, ref_e_gradrho, ref_u)
                        // Explicit portion of energy diffusion terms
                        + inv_Re_Pr_gamma1 * (
                              grad_mu.dot(grad_T)
                            + suzerain::rholut::explicit_mu_div_grad_T(
                                 gamma, Ma, mu, rho, grad_rho, div_grad_rho, m,
                                 grad_m, div_grad_m, e, div_grad_e, p, grad_p,
                                 ref_nu, ref_nuu, ref_e_deltarho)
                        )
                        // Subtract implicit portions of viscous work terms per
                        // rholut::explicit_u_dot_mu_div_grad_u and
                        // rholut::explicit_u_dot_mu_plus_lambda_grad_div_u
                        - Ma2_over_Re * (
                              ref_nuu.dot(div_grad_m)
                            - ref_nuu2*div_grad_rho
                            + alpha13*(
                                ref_nuu.dot(grad_div_m)
                              - grad_grad_rho.cwiseProduct(ref_nuuu).sum()
                            )
                        )
                        ;
                    break;
            }

            // Determine the minimum observed stable time step when necessary
            if (ZerothSubstep) {

                // Implicit acoustic handling sets the effective sound speed to
                // zero when computing the convective_stability_criterion.
                real_t a;
                switch (Linearize) {
                    case linearize::none:
                        a = std::sqrt(T) / Ma; // a/u_0 = sqrt(T*)/Ma
                        break;
                    case linearize::rhome:
                        a = 0;
                        break;
                }

                // Strictly speaking, linearize::rhome's implicit treatment of
                // the convective terms in all equations "entitles" us to
                // compute convective stability using the difference between u
                // and ref_u.  We computed those terms implicitly to buy an
                // incremental bump in the stability safety factor (which
                // should have already incorporated into evmaxmag_imag).  While
                // the resulting time steps are in general too large for
                // reasonable accuracy, we expect to be limited by diffusive
                // stability and so we do take (u - ref_u) here.
                convective_delta_t = suzerain::math::minnan(
                        convective_delta_t,
                        suzerain::timestepper::convective_stability_criterion(
                                u.x(), a, lambda1_x,
                                u.y(), a, lambda1_y,
                                u.z(), a, lambda1_z,
                                evmaxmag_imag));

                // The diffusive stability uses ref_nu which has been
                // previously set as appropriate for Linearize.  This form
                // assumes "isotropic linearization" across X, Y, and Z.
                const real_t nu = mu / rho;
                diffusive_delta_t = suzerain::math::minnan(
                        diffusive_delta_t,
                        suzerain::timestepper::diffusive_stability_criterion(
                                lambda2_x, lambda2_y, lambda2_z,
                                Re, Pr, gamma, evmaxmag_real,
                                nu, ref_nu, alpha*nu, alpha*ref_nu));
            }

        } // end X // end Z

    } // end Y
    GRVY_TIMER_END("nonlinear right hand sides");

    // Traversal:
    // (3) Computing any manufactured solution forcing (when enabled).
    // Isolating this pass allows skipping the work when unnecessary
    if (msoln) {
        GRVY_TIMER_BEGIN("manufactured forcing");

        // Dereference the msoln smart pointer outside the compute loop
        const channel::manufactured_solution &ms = *msoln;

        offset = 0;
        for (int j = o.dgrid.local_physical_start.y();
             j < o.dgrid.local_physical_end.y();
             ++j) {

            const real_t y = o.y(j);

            for (int k = o.dgrid.local_physical_start.z();
                k < o.dgrid.local_physical_end.z();
                ++k) {

                const real_t z = o.z(k);

                for (int i = o.dgrid.local_physical_start.x();
                    i < o.dgrid.local_physical_end.x();
                    ++i, /* NB */ ++offset) {

                    const real_t x = o.x(i);

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

        GRVY_TIMER_END("manufactured forcing");
    } // end msoln

    // Collectively convert state to wave space using parallel FFTs
    for (size_t i = 0; i < channel::field::count; ++i) {

        if (Linearize == linearize::rhome && i == ndx::rho && !msoln) {

            // When density equation is handled fully implicitly AND no
            // manufactured solution is employed, save some communications by
            // using that the density right hand side is identically zero.
            assert(suzerain::multi_array::is_contiguous(swave[i]));
            std::memset(swave[i].origin(), 0,
                    sizeof(complex_t)*o.dgrid.local_wave_extent.prod());

        } else {

            // Otherwise, bring the field back to wave space...
            GRVY_TIMER_BEGIN("transform_physical_to_wave");
            o.dgrid.transform_physical_to_wave(&sphys.coeffRef(i,0));
            GRVY_TIMER_END("transform_physical_to_wave");
            // ...and zero wavenumbers present only for dealiasing to
            // prevent "leakage" of dealiasing modes to other routines.
            o.zero_dealiasing_modes(swave, i);

        }

    }

    GRVY_TIMER_END("applyNonlinearOperator");

    // Return the stable time step criteria separately on each rank.  The time
    // stepping logic must perform the Allreduce.  Delegating the Allreduce
    // responsibility allows reducing additional info with minimal overhead.
    return std::vector<real_t>(delta_t_candidates.begin(),
                               delta_t_candidates.end());

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

} // namespace channel

#endif  /* NONLINEAR_HPP */
