//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_PERFECT_NAVIER_STOKES_HPP
#define SUZERAIN_PERFECT_NAVIER_STOKES_HPP

/** @file
 * Implementation of nonlinear Navier--Stokes spatial operators.
 */

#include <largo/largo.h>

#include <suzerain/error.h>
#include <suzerain/lowstorage.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/largo_definition.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/timers.h>

#include "common_block.hpp"
#include "largo_state.hpp"
#include "linearize_type.hpp"

#pragma warning(push, disable:280 383 1572)

namespace suzerain {

namespace perfect {

/**
 * A complete Navier&ndash;Stokes \c apply_operator implementation.  The
 * implementation is provided as a common building block for
 * <tt>lowstorage::nonlinear_operator< contiguous_state<4,complex_t> ></tt>
 * subclasses allowing varying numbers of passive scalars or varying hybrid
 * implicit/explicit treatment.  Such subclasses feature an overwhelming amount
 * of redundancy and are error prone to create.  This implementation allows
 * writing the "futsy" bits once and then sharing the logic repeatedly.
 * Templating allows for compile-time branching amongst different
 * implementation choices rather than paying runtime cost for such flexibility.
 *
 * Computation follows the "Numerical Considerations" section of
 * writeups/derivation.tex.  No boundary conditions are applied.  The
 * instantaneous wall-normal velocity is averaged across the streamwise and
 * spanwise directions and stored into <tt>common.u()</tt>.
 *
 * \param alpha Parameter controlling bulk viscosity
 *        per \f$\alpha\f$ in namespace \ref rholut.
 * \param beta Temperature power law exponent
 *        per \f$\beta\f$ in namespace \ref rholut.
 * \param gamma Constant ratio of specific heats
 *        per \f$\gamma\f$ in namespace \ref rholut.
 * \param Ma Mach number
 *        per \f$\textrm{Ma}\f$ in namespace \ref rholut.
 * \param Pr Prandtl number
 *        per \f$\textrm{Pr}\f$ in namespace \ref rholut.
 * \param Re Reynolds number
 *        per \f$\textrm{Re}\f$ in namespace \ref rholut.
 * \param o Provides access to discretization and parallel decomposition
 *        operational details.
 * \param common Shared storage for interaction with an linear_operator
 *        implementation providing forcing and boundary conditions.
 * \param sg Slow growth forcing to be applied via the Largo library.
 * \param msoln If \c msoln evaluates to \c true in a boolean context,
 *        then it will be used to provide manufactured forcing terms.
 * \param time Simulation time at which the operator should be applied.
 *        This allows time-dependent forcing (e.g. from \c msoln).
 * \param swave State to which the operator should be applied.  On
 *        entry, it must be coefficients in the X, Y, and Z directions.
 *        on exit, it must be coefficients in the X and Z directions but
 *        collocation point values in the Y direction.
 * \param method Low-storage timestepping scheme used to compute a stable
 *        time step when <tt>ZerothSubstep == true</tt>.
 *
 * \tparam ZerothSubstep Should one-time activities taking place at the
 *         beginning of a Runge-Kutta step be performed?  Examples include
 *         computing a stable time step size and also computing reference
 *         quantities for linearization.
 * \tparam Linearize What type of hybrid implicit/explicit linearization
 *         is employed?
 * \tparam ManufacturedSolution What manufactured solution should be used to
 *         provide additional forcing (when enabled)?
 *
 * @return A vector of stable timestep sizes according to different criteria
 *         per lowstorage::nonlinear_operator::apply_operator.
 *
 * @see lowstorage::nonlinear_operator for the (slightly different)
 *      interface that an actual operator would provide.
 */
template <bool ZerothSubstep,
          linearize::type Linearize,
          class ManufacturedSolution>
std::vector<real_t> apply_navier_stokes_spatial_operator(
            const real_t alpha,
            const real_t beta,
            const real_t gamma,
            const real_t Ma,
            const real_t Pr,
            const real_t Re,
            const operator_base &o,
            operator_common_block &common,
            support::largo_definition& sg,
            const shared_ptr<const ManufacturedSolution>& msoln,
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method)
{
    // State enters method as coefficients in X, Y, and Z directions
    SUZERAIN_TIMER_SCOPED("apply_navier_stokes_spatial_operator");

    // Shorthand
    typedef contiguous_state<4,complex_t> state_type;
    using std::abs;
    using std::equal;
    using std::max;
    using std::min;
    using std::size_t;
    using std::sqrt;

    // Compile-time template parameters are used to reduce jumps at runtime.
    // Ensure that someone didn't hand in a mismatched common block.
    SUZERAIN_ENSURE(common.linearization  == Linearize);

    SUZERAIN_ENSURE(!sg.formulation.enabled()); // FIXME Redmine #2495

    // FIXME Ticket #2477 retrieve linearization-dependent CFL information
    // Afterwards, change the stable time step computation accordingly
    const real_t evmaxmag_real = method.evmaxmag_real();
    const real_t evmaxmag_imag = method.evmaxmag_imag();

    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    enum { swave_count = 5 };
    assert(static_cast<int>(ndx::e  ) < swave_count);
    assert(static_cast<int>(ndx::mx ) < swave_count);
    assert(static_cast<int>(ndx::my ) < swave_count);
    assert(static_cast<int>(ndx::mz ) < swave_count);
    assert(static_cast<int>(ndx::rho) < swave_count);

    // We need auxiliary scalar-field storage.  Prepare logical indices using a
    // struct for scoping (e.g. aux::rho_y).  Ordering will match usage below.
    // TODO Only linearize::rhome_y needs e_yy.  Avoid overhead in other cases.
    struct aux { enum {
        e_y,   e_yy,   div_grad_e, e_x, e_z,
        mx_y,  mx_yy,  mx_x,  mx_xx,  mx_xz,  mx_z,  mx_zz,  mx_xy,  mx_yz,
        my_y,  my_yy,  my_x,  my_xx,  my_xz,  my_z,  my_zz,  my_xy,  my_yz,
        mz_y,  mz_yy,  mz_x,  mz_xx,  mz_xz,  mz_z,  mz_zz,  mz_xy,  mz_yz,
        rho_y, rho_yy, rho_x, rho_xx, rho_xz, rho_z, rho_zz, rho_xy, rho_yz,
        count // Sentry
    }; };

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    scoped_ptr<state_type> _auxw_ptr(
            support::allocate_padded_state<state_type>(
                aux::count, o.dgrid)); // RAII
    state_type &auxw = *_auxw_ptr;                                   // Brevity

    // Sanity check incoming swave's and auxw's shape and contiguity
    SUZERAIN_ENSURE(swave.shape()[0] == swave_count);
    SUZERAIN_ENSURE(swave.shape()[1] == (unsigned) o.dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(swave.shape()[2] == (unsigned) o.dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(swave.shape()[3] == (unsigned) o.dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE((unsigned) swave.strides()[1] == 1u);
    SUZERAIN_ENSURE((unsigned) swave.strides()[2] == swave.shape()[1]);
    SUZERAIN_ENSURE((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    SUZERAIN_ENSURE(equal(swave.shape()   + 1, swave.shape()   + 4, auxw.shape()   + 1));
    SUZERAIN_ENSURE(equal(swave.strides() + 1, swave.strides() + 4, auxw.strides() + 1));

    // Prepare common-block-like storage used to pass details from N to L.
    // Zeroing is done carefully as accumulated means and reference quantities
    // must survive across nonzero substeps while instant profiles must not.
    // We do not modify common.implicits as other logic populates its content.
    if (ZerothSubstep) {
        common.set_zero(/* Ny */ swave.shape()[1]);
    } else {
        common.means.setZero();
    }

    // Maintain stable time step values to return to the caller:
    //
    //   convtotal_* is convective stability using total velocity.
    //               It may or may not have acoustics included as necessary.
    //   convfluct_* is convective stability using only velocity fluctuations.
    //               It permits measuring any advantage in treating mean
    //               convection implicitly.
    //   diffusive_* is diffusive stability linearized however is appropriate.
    //   *_xyz_*     accounts for restrictions in all three directions.
    //               It must be more restrictive than any single direction.
    //   *_x_*       accounts for restrictions in streamwise direction only.
    //   *_y_*       accounts for restrictions in wall-normal direction only
    //   *_z_*       accounts for restrictions in spanwise direction only
    //
    // Results reported by the code are reported in this order.
    array<real_t, 12> delta_t_candidates;
    std::fill(delta_t_candidates.begin(), delta_t_candidates.end(),
              std::numeric_limits<real_t>::max());
    real_t &convtotal_xyz_delta_t = delta_t_candidates[ 0];
    real_t &convfluct_xyz_delta_t = delta_t_candidates[ 1];
    real_t &diffusive_xyz_delta_t = delta_t_candidates[ 2];
    real_t &convtotal_x_delta_t   = delta_t_candidates[ 3];
    real_t &convfluct_x_delta_t   = delta_t_candidates[ 4];
    real_t &diffusive_x_delta_t   = delta_t_candidates[ 5];
    real_t &convtotal_y_delta_t   = delta_t_candidates[ 6];
    real_t &convfluct_y_delta_t   = delta_t_candidates[ 7];
    real_t &diffusive_y_delta_t   = delta_t_candidates[ 8];
    real_t &convtotal_z_delta_t   = delta_t_candidates[ 9];
    real_t &convfluct_z_delta_t   = delta_t_candidates[10];
    real_t &diffusive_z_delta_t   = delta_t_candidates[11];

    // Compute Y derivatives of total energy at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.zero_dealiasing_modes(  swave, ndx::e);
    o.bop_accumulate(1,    1, swave, ndx::e, 0, auxw, aux::e_y);
    o.bop_accumulate(2,    1, swave, ndx::e, 0, auxw, aux::e_yy);
    std::memcpy(auxw[aux::div_grad_e].origin(), auxw[aux::e_yy].origin(),
                static_cast<size_t>(   auxw[aux::div_grad_e + 1].origin()
                                     - auxw[aux::div_grad_e    ].origin())
                                   * sizeof(complex_t));
    o.bop_apply     (0,    1, swave, ndx::e);

    // Compute X- and Z- derivatives of total energy at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::e, 0, auxw, aux::e_x       );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::e, 1, auxw, aux::div_grad_e);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::e, 0, auxw, aux::e_z       );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::e, 1, auxw, aux::div_grad_e);

    // Compute Y derivatives of X momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.zero_dealiasing_modes(  swave, ndx::mx);
    o.bop_accumulate(1,    1, swave, ndx::mx, 0, auxw, aux::mx_y);
    o.bop_accumulate(2,    1, swave, ndx::mx, 0, auxw, aux::mx_yy);
    o.bop_apply     (0,    1, swave, ndx::mx);

    // Compute X- and Z- derivatives of X momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::mx,   0, auxw, aux::mx_x );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::mx,   0, auxw, aux::mx_xx);
    o.diffwave_accumulate(1, 1, 1, swave, ndx::mx,   0, auxw, aux::mx_xz);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::mx,   0, auxw, aux::mx_z );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::mx,   0, auxw, aux::mx_zz);
    o.diffwave_accumulate(1, 0, 1, auxw,  aux::mx_y, 0, auxw, aux::mx_xy);
    o.diffwave_accumulate(0, 1, 1, auxw,  aux::mx_y, 0, auxw, aux::mx_yz);

    // Compute Y derivatives of Y momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.zero_dealiasing_modes(  swave, ndx::my);
    o.bop_accumulate(1,    1, swave, ndx::my, 0, auxw, aux::my_y);
    o.bop_accumulate(2,    1, swave, ndx::my, 0, auxw, aux::my_yy);
    o.bop_apply     (0,    1, swave, ndx::my);

    // Compute X- and Z- derivatives of Y momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::my,   0, auxw, aux::my_x );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::my,   0, auxw, aux::my_xx);
    o.diffwave_accumulate(1, 1, 1, swave, ndx::my,   0, auxw, aux::my_xz);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::my,   0, auxw, aux::my_z );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::my,   0, auxw, aux::my_zz);
    o.diffwave_accumulate(1, 0, 1, auxw,  aux::my_y, 0, auxw, aux::my_xy);
    o.diffwave_accumulate(0, 1, 1, auxw,  aux::my_y, 0, auxw, aux::my_yz);

    // Compute Y derivatives of Z momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.zero_dealiasing_modes(  swave, ndx::mz);
    o.bop_accumulate(1,    1, swave, ndx::mz, 0, auxw, aux::mz_y);
    o.bop_accumulate(2,    1, swave, ndx::mz, 0, auxw, aux::mz_yy);
    o.bop_apply     (0,    1, swave, ndx::mz);

    // Compute X- and Z- derivatives of Z momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, swave, ndx::mz,   0, auxw, aux::mz_x );
    o.diffwave_accumulate(2, 0, 1, swave, ndx::mz,   0, auxw, aux::mz_xx);
    o.diffwave_accumulate(1, 1, 1, swave, ndx::mz,   0, auxw, aux::mz_xz);
    o.diffwave_accumulate(0, 1, 1, swave, ndx::mz,   0, auxw, aux::mz_z );
    o.diffwave_accumulate(0, 2, 1, swave, ndx::mz,   0, auxw, aux::mz_zz);
    o.diffwave_accumulate(1, 0, 1, auxw,  aux::mz_y, 0, auxw, aux::mz_xy);
    o.diffwave_accumulate(0, 1, 1, auxw,  aux::mz_y, 0, auxw, aux::mz_yz);

    // Compute Y derivatives of density at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    o.zero_dealiasing_modes(  swave, ndx::rho);
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

    // Now that conserved state is Fourier in X and Z but collocation in Y,
    // if slow growth forcing is active...
    std::vector<field_L2xz> rms  (0);
    std::vector<field_L2xz> rms_y(0);
    if (sg.formulation.enabled()) {
        SUZERAIN_TIMER_SCOPED("root-mean-square of state");
        // ...collectively compute L^2_{xz} of state at each collocation point
        rms = compute_field_L2xz(swave, o.grid, o.dgrid);

        // ...allocate space to additionally house wall-normal derivatives
        rms_y.resize(rms.size());

        // ...rescale results to convert to root-mean-square (RMS) values,
        // and use the B-spline basis to compute wall-normal derivatives.
        const real_t rms_adjust = 1 / std::sqrt(o.grid.L.x() * o.grid.L.z());
        ArrayXr tmp;
        for (size_t i = 0; i < rms.size(); ++i) {
            rms[i].mean *= rms_adjust;
            tmp = rms[i].mean;
            o.masslu()->solve(1, tmp.data(), 1, tmp.size());
            rms_y[i].mean.resizeLike(tmp);
            o.cop.accumulate(1, 1, tmp.data(),           1,
                                0, rms_y[i].mean.data(), 1);

            rms[i].fluctuating *= rms_adjust;
            tmp = rms[i].fluctuating;
            o.masslu()->solve(1, tmp.data(), 1, tmp.size());
            rms_y[i].fluctuating.resizeLike(tmp);
            o.cop.accumulate(1, 1, tmp.data(),                  1,
                                0, rms_y[i].fluctuating.data(), 1);
        }
    }

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    physical_view<aux::count>  auxp(o.dgrid, auxw);
    physical_view<swave_count> sphys(o.dgrid, swave);
    for (size_t i = 0; i < swave_count; ++i) {
        o.dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (size_t i = 0; i < aux::count; ++i) {
        o.dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }

    // Compute derived constants before inner loops
    const real_t alpha13          = alpha + real_t(1)/real_t(3);
    const real_t alpha43          = alpha + real_t(4)/real_t(3);
    const real_t inv_Re           = 1 / Re;
    const real_t inv_Ma2          = 1 / (Ma * Ma);
    const real_t Ma2_over_Re      = (Ma * Ma) / Re;
    const real_t gamma1           = gamma - 1;
    const real_t inv_Re_Pr_gamma1 = 1 / (Re * Pr * gamma1);
    const real_t gamma_over_Pr    = gamma / Pr;
    const real_t lambda1_x        = o.lambda1_x;
    const real_t lambda1_z        = o.lambda1_z;
    const real_t maxdiffconst     = inv_Re*max(gamma/Pr, max(real_t(1), alpha));
    const real_t md_lambda2_x     = maxdiffconst * o.lambda2_x;
    const real_t md_lambda2_z     = maxdiffconst * o.lambda2_z;

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
    // (1) Computing reference quantities and miscellaneous moments OR
    //     just reduced moments depending on the substep being performed.
    // (2) Computing the nonlinear equation right hand sides.
    // (3) Computing any manufactured solution forcing (when enabled).
    //
    // The traversal pattern is best embodied by the third pass.  The first and
    // second passes use slightly more compact structure as only y(j) must be
    // known within them.  Study (3) and convince yourself that (1) and (2) are
    // equivalent but lack information on x(i) and z(k).

    // Traversal:
    // (1) Computing reference quantities and miscellaneous moments OR
    //     just reduced moments depending on the substep being performed.
    if (ZerothSubstep) {

        SUZERAIN_TIMER_SCOPED("reference quantities");

        // To avoid accumulating garbage, must zero y(j) not present on rank.
        // Clearing everything is expected to be a bit more performant.
        common.refs.setZero();

        // Sum reference quantities as a function of y(j) into common.ref_*
        // See writeups/derivation.tex or rholut_imexop.h for definitions
        for (int offset = 0, j = o.dgrid.local_physical_start.y();
            j < o.dgrid.local_physical_end.y();
            ++j) {

            // Prepare logical indices using struct for scoping (e.g. ref::ux).
            struct ref { enum { rho, p, p2, T, a,
                                ux, uy, uz, u2,
                                uxux, uxuy, uxuz, uyuy, uyuz, uzuz,
                                nu, nuux, nuuy, nuuz, nuu2,
                                nuuxux, nuuxuy, nuuxuz, nuuyuy, nuuyuz, nuuzuz,
                                ex_gradrho, ey_gradrho, ez_gradrho,
                                e_divm, e_deltarho,
                                rhouxux, rhouyuy, rhouzuz, rhoEE,
                                count // Sentry
            }; };

            // An array of summing_accumulator_type holds all running sums.
            // This gives nicer construction and allows looping over results.
            summing_accumulator_type acc[ref::count];

            const int last_zxoffset = offset
                                    + o.dgrid.local_physical_extent.z()
                                    * o.dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {

                // Unpack conserved state
                const real_t   e  (sphys(ndx::e,   offset));
                const Vector3r m  (sphys(ndx::mx,  offset),
                                   sphys(ndx::my,  offset),
                                   sphys(ndx::mz,  offset));
                const real_t   rho(sphys(ndx::rho, offset));

                // Compute quantities related to the equation of state
                real_t p, T, mu, lambda;
                rholut::p_T_mu_lambda(
                    alpha, beta, gamma, Ma, rho, m, e, p, T, mu, lambda);

                // Accumulate reference quantities into running sums...
                acc[ref::rho](rho);
                acc[ref::p  ](p);
                acc[ref::p2 ](p*p);
                acc[ref::T  ](T);
                acc[ref::a  ](sqrt(T));

                // ...including simple velocity-related quantities...
                const Vector3r u = rholut::u(rho, m);
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

                // ...other, more complicated expressions...
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

                // ...and, lastly, details needed for slow growth forcing.
                acc[ref::rhouxux](m.x()* m.x() / rho);
                acc[ref::rhouyuy](m.y()* m.y() / rho);
                acc[ref::rhouzuz](m.z()* m.z() / rho);
                acc[ref::rhoEE  ](e    * e     / rho);

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
            common.ref_rho       ()[j] = sum(acc[ref::rho       ]);
            common.ref_p         ()[j] = sum(acc[ref::p         ]);
            common.ref_p2        ()[j] = sum(acc[ref::p2        ]);
            common.ref_T         ()[j] = sum(acc[ref::T         ]);
            common.ref_a         ()[j] = sum(acc[ref::a         ]);
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
            common.ref_rhouxux   ()[j] = sum(acc[ref::rhouxux   ]);
            common.ref_rhouyuy   ()[j] = sum(acc[ref::rhouyuy   ]);
            common.ref_rhouzuz   ()[j] = sum(acc[ref::rhouzuz   ]);
            common.ref_rhoEE     ()[j] = sum(acc[ref::rhoEE     ]);

        } // end Y

        // Allreduce and scale common.refs sums to obtain means on all ranks
        SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, common.refs.data(),
                common.refs.size(), mpi::datatype<real_t>::value,
                MPI_SUM, MPI_COMM_WORLD));
        common.refs *= o.dgrid.chi();

        // Copy mean velocity information into common.{u, v, w}()
        common.u() = common.ref_ux();
        common.v() = common.ref_uy();
        common.w() = common.ref_uz();

        // Copy mean Reynolds stresses into common.{uu, uv, uw, vv, vw, ww}()
        common.uu() = common.ref_uxux();
        common.uv() = common.ref_uxuy();
        common.uw() = common.ref_uxuz();
        common.vv() = common.ref_uyuy();
        common.vw() = common.ref_uyuz();
        common.ww() = common.ref_uzuz();

        // Copy mean information additionally required for slow growth forcing
        common.rho  () = common.ref_rho    ();
        common.rhouu() = common.ref_rhouxux();
        common.rhovv() = common.ref_rhouyuy();
        common.rhoww() = common.ref_rhouzuz();
        common.rhoEE() = common.ref_rhoEE  ();
        common.p    () = common.ref_p      ();
        common.p2   () = common.ref_p2     ();

    } else {

        SUZERAIN_TIMER_SCOPED("instantaneous moments");

        // To avoid accumulating garbage, must zero y(j) not present on rank.
        // Clearing everything is expected to be a bit more performant.
        common.means.setZero();

        // Accumulate velocity moments into common storage as function of y(j)
        for (int offset = 0, j = o.dgrid.local_physical_start.y();
            j < o.dgrid.local_physical_end.y();
            ++j) {

            // Prepare logical indices using struct for scoping (e.g. q::u).
            struct q { enum { u,  v,  w, uu, uv, uw, vv, vw, ww,
                              rho, rhouu, rhovv, rhoww, rhoEE, p, p2,
                              count // Sentry
            }; };

            // An array of summing_accumulator_type holds all running sums.
            // This gives nicer construction and allows looping over results.
            summing_accumulator_type acc[q::count];

            const int last_zxoffset = offset
                                    + o.dgrid.local_physical_extent.z()
                                    * o.dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {

                // Unpack conserved state
                const real_t   e  (sphys(ndx::e,   offset));
                const Vector3r m  (sphys(ndx::mx,  offset),
                                   sphys(ndx::my,  offset),
                                   sphys(ndx::mz,  offset));
                const real_t   rho(sphys(ndx::rho, offset));

                // Compute derived quantities
                const Vector3r u  (m / rho);
                real_t p;
                rholut::p(alpha, beta, gamma, Ma, rho, m, e, p);

                // Accumulate pointwise information
                acc[q::u    ](u.x());
                acc[q::v    ](u.y());
                acc[q::w    ](u.z());
                acc[q::uu   ](u.x() * u.x());
                acc[q::uv   ](u.x() * u.y());
                acc[q::uw   ](u.x() * u.z());
                acc[q::vv   ](u.y() * u.y());
                acc[q::vw   ](u.y() * u.z());
                acc[q::ww   ](u.z() * u.z());
                acc[q::rho  ](rho);
                acc[q::rhouu](u.x() * u.x() * rho);
                acc[q::rhovv](u.y() * u.y() * rho);
                acc[q::rhoww](u.z() * u.z() * rho);
                acc[q::rhoEE](e     * e     / rho);
                acc[q::p    ](p);
                acc[q::p2   ](p * p);

            } // end X // end Z

#ifndef NDEBUG
            // Ensure that all accumulators saw a consistent number of samples
            const size_t expected = boost::accumulators::count(acc[0]);
            for (size_t k = 1; k < sizeof(acc)/sizeof(acc[0]); ++k) {
                const size_t observed = boost::accumulators::count(acc[k]);
                assert(expected == observed);
            }
#endif

            // Store sum into common block in preparation for MPI Allreduce
            using boost::accumulators::sum;
            common.u    ()[j] = sum(acc[q::u    ]);
            common.v    ()[j] = sum(acc[q::v    ]);
            common.w    ()[j] = sum(acc[q::w    ]);
            common.uu   ()[j] = sum(acc[q::uu   ]);
            common.uv   ()[j] = sum(acc[q::uv   ]);
            common.uw   ()[j] = sum(acc[q::uw   ]);
            common.vv   ()[j] = sum(acc[q::vv   ]);
            common.vw   ()[j] = sum(acc[q::vw   ]);
            common.ww   ()[j] = sum(acc[q::ww   ]);
            common.rho  ()[j] = sum(acc[q::rho  ]);
            common.rhouu()[j] = sum(acc[q::rhouu]);
            common.rhovv()[j] = sum(acc[q::rhovv]);
            common.rhoww()[j] = sum(acc[q::rhoww]);
            common.rhoEE()[j] = sum(acc[q::rhoEE]);
            common.p    ()[j] = sum(acc[q::p    ]);
            common.p2   ()[j] = sum(acc[q::p2   ]);

        } // end Y

        // Allreduce and scale sums to produce the mean on all ranks
        SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE,
                    common.means.data(), common.means.size(),
                    mpi::datatype<real_t>::value, MPI_SUM, MPI_COMM_WORLD));
        common.means *= o.dgrid.chi();

    }

    // If necessary, perform globally-relevant initialization calls to Largo.
    largo_state basewall;
    if (sg.formulation.enabled()) {
        largo_state dy, dx;
        sg.get_baseflow(0.0, basewall.as_is(), dy.as_is(), dx.as_is());

        largo_state grDA, grDArms;
        if (!basewall.trivial()) {
            grDA.mx = - basewall.u() * dx.mx / (0.0 - basewall.u());
        }

        largo_init(sg.workspace, sg.grdelta,
                   grDA.rescale(inv_Ma2), grDArms.rescale(inv_Ma2));
    }

    // Traversal:
    // (2) Computing the nonlinear equation right hand sides.
    SUZERAIN_TIMER_BEGIN("nonlinear right hand sides");
    for (int offset = 0, j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {

        // Wall-normal operator eigenvalue estimates depend on location
        const real_t lambda1_y    = o.lambda1_y(j);
        const real_t md_lambda2_y = maxdiffconst * o.lambda2_y(j);

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

        // If necessary, perform Largo base flow and Y-dependent invocations
        if (sg.formulation.enabled()) {
            largo_state base, dy, dx;
            sg.get_baseflow(o.y(j), base.as_is(), dy.as_is(), dx.as_is());

            // When a nontrivial inviscid base flow is present, compute
            // Euler residual so that it might be eradicated.  The base
            // flow is intended to be stationary so any residual it is
            // a numerical artifact and not part of the intended problem.
            largo_state dt, src;
            if (!base.trivial()) {

                // Compute time derivative of base flow using wall information
                dt.rho = basewall.u() * dx.rho;
                dt.mx  = basewall.u() * dx.mx;
                dt.my  = basewall.u() * dx.my;
                dt.mz  = basewall.u() * dx.mz;
                dt.e   = basewall.u() * dx.e;

                // Compute pressure-related quantities
                real_t P, dyP, dxP;
                sg.get_baseflow_pressure(o.y(j), P, dyP, dxP); // as_is()

                // Compute inviscid base flow residual from Euler equations
                // Residual permits only non-trivial wall-normal derivatives
                const double u = base.u();
                const double v = base.v();
                const double w = base.w();
                const double H = (base.e + P) / base.rho;
                src.rho = dt.rho + dy.my;
                src.mx  = dt.mx  + v*(dy.mx - u*dy.rho) + u*dy.my;
                src.my  = dt.my  + v*(dy.my - v*dy.rho) + v*dy.my + inv_Ma2*dyP;
                src.mz  = dt.mz  + v*(dy.mz - w*dy.rho) + w*dy.my;
                src.e   = dt.e   + H*(dy.my - v*dy.rho) + v*(dy.e + dyP);
            }

            largo_prestep_baseflow(sg.workspace, base.rescale(inv_Ma2),
                                   dy.rescale(inv_Ma2), dt.rescale(inv_Ma2),
                                   dx.rescale(inv_Ma2), src.rescale(inv_Ma2));

            // FIXME #2495 innery
        }

        // Iterate across the j-th ZX plane
        const int last_zxoffset = offset
                                + o.dgrid.local_physical_extent.z()
                                * o.dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // Unpack total energy-related quantities
            const real_t e        (sphys(ndx::e,          offset));
            const Vector3r grad_e ( auxp(aux::e_x,        offset),
                                    auxp(aux::e_y,        offset),
                                    auxp(aux::e_z,        offset));
            const real_t e_yy      (auxp(aux::e_yy,       offset));
            const real_t div_grad_e(auxp(aux::div_grad_e, offset));

            // Unpack momentum-related quantities
            const Vector3r m    ( sphys(ndx::mx,   offset),
                                  sphys(ndx::my,   offset),
                                  sphys(ndx::mz,   offset));
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
            const real_t   mx_yy     (  auxp(aux::mx_yy, offset));
            const real_t   my_yy     (  auxp(aux::my_yy, offset));
            const real_t   mz_yy     (  auxp(aux::mz_yy, offset));
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

            // Unpack density-related quantities
            const real_t   rho         ( sphys(ndx::rho,    offset));
            const Vector3r grad_rho    (  auxp(aux::rho_x,  offset),
                                          auxp(aux::rho_y,  offset),
                                          auxp(aux::rho_z,  offset));
            const real_t   rho_yy      (  auxp(aux::rho_yy, offset));
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

            // Compute velocity-related quantities
            const Vector3r u          = rholut::u(
                                            rho, m);
            const real_t div_u        = rholut::div_u(
                                            rho, grad_rho, m, div_m);
            const Matrix3r grad_u     = rholut::grad_u(
                                            rho, grad_rho, m, grad_m);
            const Vector3r grad_div_u = rholut::grad_div_u(
                                            rho, grad_rho, grad_grad_rho,
                                            m, div_m, grad_m, grad_div_m);
            const Vector3r div_grad_u = rholut::div_grad_u(
                                            rho, grad_rho, div_grad_rho,
                                            m, grad_m, div_grad_m);

            // Compute quantities related to the equation of state
            real_t p, T, mu, lambda;
            Vector3r grad_p, grad_T, grad_mu, grad_lambda;
            rholut::p_T_mu_lambda(alpha, beta, gamma, Ma,
                                  rho, grad_rho, m, grad_m, e, grad_e,
                                  p, grad_p, T, grad_T,
                                  mu, grad_mu, lambda, grad_lambda);
            const real_t nu         = mu / rho;
            const real_t div_grad_p = rholut::div_grad_p(
                                        gamma, Ma,
                                        rho, grad_rho, div_grad_rho,
                                        m, grad_m, div_grad_m,
                                        e, grad_e, div_grad_e);
            const real_t div_grad_T = rholut::div_grad_T(
                                        gamma,
                                        rho, grad_rho, div_grad_rho,
                                        p, grad_p, div_grad_p);

            // Compute quantities related to the viscous stress tensor
            const Matrix3r tau     = rholut::tau(
                                        mu, lambda, div_u, grad_u);
            const Vector3r div_tau = rholut::div_tau(
                                        mu, grad_mu, lambda, grad_lambda,
                                        div_u, grad_u, div_grad_u,
                                        grad_div_u);

            // FORM ENERGY EQUATION RIGHT HAND SIDE
            sphys(ndx::e, offset) =
                // Explicit viscous work term
                + Ma2_over_Re * rholut::div_tau_u<real_t>(
                        u, grad_u, tau, div_tau
                    )
                ;
            switch (Linearize) {
            case linearize::rhome_xyz:
                sphys(ndx::e, offset) +=
                    // Explicit convective/acoustic less implicit portion
                    - rholut::explicit_div_e_plus_p_u(
                            gamma, Ma, rho, grad_rho,
                            m, div_m, grad_m, e, grad_e, p,
                            ref_e_divm, ref_e_gradrho, ref_u)
                    // Explicit portion of energy diffusion terms
                    + inv_Re_Pr_gamma1 * (
                          grad_mu.dot(grad_T)
                        + rholut::explicit_mu_div_grad_T(
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

            case linearize::rhome_y:
            {
                // Build up colored terms from perfect gas writeup figure 2
                const real_t term_rho_y
                    = - ref_e_gradrho.y();
                const real_t term_rho_yy
                    = gamma * inv_Re_Pr_gamma1 * ref_e_deltarho
                    - Ma2_over_Re * (ref_nuu2 + alpha13 * ref_nuuu(1,1));
                const real_t term_mx_yy
                    = Ma2_over_Re * (1 - gamma_over_Pr) * ref_nuu.x();
                const real_t term_my_y
                    = - ref_e_divm;
                const real_t term_my_yy
                    = Ma2_over_Re * (alpha43 - gamma_over_Pr) * ref_nuu.y();
                const real_t term_mz_yy
                    = Ma2_over_Re * (1 - gamma_over_Pr) * ref_nuu.z();
                const real_t term_e_y
                    = - gamma * ref_u.y();
                const real_t term_e_yy
                    = inv_Re * gamma_over_Pr * ref_nu;

                // Subtract terms scaled by appropriate state derivatives
                sphys(ndx::e, offset) -=
                      term_rho_y  * grad_rho.y()
                    + term_rho_yy * rho_yy
                    + term_mx_yy  * mx_yy
                    + term_my_y   * grad_m(1,1)
                    + term_my_yy  * my_yy
                    + term_mz_yy  * mz_yy
                    + term_e_y    * grad_e.y()
                    + term_e_yy   * e_yy
                    ;

                // ...
                // Fall through!
                // ...
            }

            case linearize::none:
                sphys(ndx::e, offset) +=
                    // Explicit convective and acoustic terms
                    - rholut::div_e_u(
                            e, grad_e, u, div_u
                        )
                    - rholut::div_p_u(
                            p, grad_p, u, div_u
                        )
                    // Explicit energy diffusion terms
                    + inv_Re_Pr_gamma1 * rholut::div_mu_grad_T(
                            grad_T, div_grad_T, mu, grad_mu
                        )
                    ;
                    // No need to adjust explicit viscous work term
                break;

            default:
                SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();
            }

            // FORM MOMENTUM EQUATION RIGHT HAND SIDE
            Vector3r momentum_rhs =
                // Explicit viscous term
                  inv_Re * div_tau
                ;
            switch (Linearize) {
            case linearize::rhome_xyz:
                momentum_rhs +=
                    // Explicit convective term less implicit portion
                    - rholut::explicit_div_rho_inverse_m_outer_m(
                            grad_rho, div_m, grad_m, u, ref_u, ref_uu)
                    // Explicit pressure less implicit pressure terms
                    - inv_Ma2 * rholut::explicit_grad_p(
                            gamma, Ma, rho, grad_rho, m, grad_m,
                            ref_u2, ref_u)
                    // Subtract implicit portions of viscous terms per
                    // rholut::explicit_mu_div_grad_u and
                    // rholut::explicit_mu_plus_lambda_grad_div_u
                    - inv_Re * (
                        ref_nu*(div_grad_m + alpha13*grad_div_m)
                      - ref_nuu*div_grad_rho
                      - alpha13*grad_grad_rho*ref_nuu
                    )
                    ;
                break;

            case linearize::rhome_y:
                // Subtract colored terms from perfect gas writeup figure 2
                momentum_rhs.x() -=
                               ref_uu(0,1) * grad_rho.y()
                    - inv_Re * ref_nuu.x() * rho_yy
                    -          ref_u.y()   * grad_m(0,1)
                    + inv_Re * ref_nu      * mx_yy
                    -          ref_u.x()   * grad_m(1,1)
                    ;
                momentum_rhs.y() -=
                      (ref_uu(1,1) - gamma1 / 2 * ref_u2) * grad_rho.y()
                    - alpha43 * inv_Re * ref_nuu.y()      * rho_yy
                    + gamma1 * ref_u.x()                  * grad_m(0,1)
                    + (gamma1 - 2) * ref_u.y()            * grad_m(1,1)
                    + alpha43 * inv_Re * ref_nu           * my_yy
                    + gamma1 * ref_u.z()                  * grad_m(2,1)
                    - gamma1 * inv_Ma2                    * grad_e.y()
                    ;
                momentum_rhs.z() -=
                               ref_uu(1,2) * grad_rho.y()
                    - inv_Re * ref_nuu.z() * rho_yy
                    -          ref_u.z()   * grad_m(1,1)
                    -          ref_u.y()   * grad_m(2,1)
                    + inv_Re * ref_nu      * mz_yy
                    ;

                // ...
                // Fall through!
                // ...

            case linearize::none:
                momentum_rhs +=
                    // Explicit convective term
                    - rholut::div_u_outer_m(m, grad_m, u, div_u)
                    // Explicit pressure term
                    - inv_Ma2 * grad_p
                    ;
                break;

            default:
                SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();
            }
            sphys(ndx::mx, offset) = momentum_rhs.x();
            sphys(ndx::my, offset) = momentum_rhs.y();
            sphys(ndx::mz, offset) = momentum_rhs.z();

            // FORM CONTINUITY EQUATION RIGHT HAND SIDE
            //
            // Implicit continuity equation handling requires zeroing RHS in
            // anticipation of possible manufactured solution forcing.  See
            // subsequent transform_physical_to_wave if you monkey around here.
            switch (Linearize) {
            case linearize::rhome_xyz:    // Fully implicit convection
                sphys(ndx::rho, offset) = 0;
                break;

            case linearize::rhome_y:      // Implicit Y but Explicit X, Z
                sphys(ndx::rho, offset) = - grad_m(0,0) - grad_m(2,2);
                break;

            case linearize::none:         // Full explicit convection
                sphys(ndx::rho, offset) = - div_m;
                break;

            default:
                SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();
            }

            // Determine the minimum observed stable time step when necessary
            // This logic used to call some canned routines, but additional
            // monitoring requirements forced inlining many of these details.
            // See delta_t_candidates declaration (above) for descriptions.
            //
            // details on these computations in the three-directional case.
            if (ZerothSubstep) {

                // minnan(...) calls only required on *_xyz_delta_t as the
                // three-dimensional criterion should collect any NaNs
                // occurring within the direction-dependent criteria.  min(...)
                // is presumably a touch faster and is preferred when possible.
                //
                // maxnan(...) truncate non-positive values but prefer NaNs.
                using math::minnan;
                using math::maxnan;

                // See lowstorage::convective_stability_criterion
                const real_t a = sqrt(T) / Ma;  // Because a/u_0 = sqrt(T*)/Ma
                real_t       ua_l1_x,       ua_l1_y,       ua_l1_z;
                real_t fluct_ua_l1_x, fluct_ua_l1_y, fluct_ua_l1_z;
                switch (Linearize) {
                // Implicit acoustics sets the effective sound speed
                // to zero within the convective_stability_criterion.
                // Fluctuating velocity is taken relative to references.
                case linearize::rhome_xyz:
                    ua_l1_x       = abs(u.x()            ) * lambda1_x;
                    ua_l1_y       = abs(u.y()            ) * lambda1_y;
                    ua_l1_z       = abs(u.z()            ) * lambda1_z;
                    fluct_ua_l1_x = abs(u.x() - ref_u.x()) * lambda1_x;
                    fluct_ua_l1_y = abs(u.y() - ref_u.y()) * lambda1_y;
                    fluct_ua_l1_z = abs(u.z() - ref_u.z()) * lambda1_z;
                    break;

                // Explicit treatment forces including acoustics
                // in stability and has a zero reference velocity.
                case linearize::none:
                    ua_l1_x       = (abs(u.x()) + a) * lambda1_x;
                    ua_l1_y       = (abs(u.y()) + a) * lambda1_y;
                    ua_l1_z       = (abs(u.z()) + a) * lambda1_z;
                    fluct_ua_l1_x = ua_l1_x;
                    fluct_ua_l1_y = ua_l1_y;
                    fluct_ua_l1_z = ua_l1_z;
                    break;

                // Wall-normal-only implicit acoustics and convection
                // is nothing but a hybrid of the above two cases.
                case linearize::rhome_y:
                    ua_l1_x       = (abs(u.x()) + a) * lambda1_x;
                    ua_l1_y       = (abs(u.y())    ) * lambda1_y;
                    ua_l1_z       = (abs(u.z()) + a) * lambda1_z;
                    fluct_ua_l1_x = ua_l1_x;
                    fluct_ua_l1_y = abs(u.y() - ref_u.y()) * lambda1_y;
                    fluct_ua_l1_z = ua_l1_z;
                    break;

                default:
                    SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();
                }
                convtotal_xyz_delta_t = minnan(convtotal_xyz_delta_t,
                        evmaxmag_imag / (ua_l1_x + ua_l1_y + ua_l1_z));
                convtotal_x_delta_t   = min   (convtotal_x_delta_t,
                        evmaxmag_imag / ua_l1_x);
                convtotal_y_delta_t   = min   (convtotal_y_delta_t,
                        evmaxmag_imag / ua_l1_y);
                convtotal_z_delta_t   = min   (convtotal_z_delta_t,
                        evmaxmag_imag / ua_l1_z);
                convfluct_xyz_delta_t = minnan(convfluct_xyz_delta_t,
                        evmaxmag_imag / (  fluct_ua_l1_x
                                         + fluct_ua_l1_y
                                         + fluct_ua_l1_z));
                convfluct_x_delta_t   = min   (convfluct_x_delta_t,
                        evmaxmag_imag / fluct_ua_l1_x);
                convfluct_y_delta_t   = min   (convfluct_y_delta_t,
                        evmaxmag_imag / fluct_ua_l1_y);
                convfluct_z_delta_t   = min   (convfluct_z_delta_t,
                        evmaxmag_imag / fluct_ua_l1_z);

                // See lowstorage::diffusive_stability_criterion
                // Antidiffusive locations might be ignored when linearized.
                // Hence we compute criteria within the switch statement.
                //
                // We already absorbed maxdiffconst within md_lambda2_{x,y,z}
                // to account for viscous, bulk viscous, and thermal effects.
                real_t diffusivity = nu;
                switch (Linearize) {
                // Implicit diffusion permits removing a reference value.
                // Antidiffusive (nu - ref_nu) is fine and not computed.
                case linearize::rhome_xyz:
                    diffusivity -= ref_nu;        // Compute sign wrt ref.
                    if (diffusivity <= 0) break;  // NaN => false, proceed
                    // ...
                    // Fall through!
                    // ...

                // Explicit treatment forces a zero reference diffusivity
                case linearize::none:
                    diffusive_xyz_delta_t = minnan(diffusive_xyz_delta_t,
                              evmaxmag_real / (   diffusivity
                                                * (   md_lambda2_x
                                                    + md_lambda2_y
                                                    + md_lambda2_z)));
                    diffusive_x_delta_t   = min   (diffusive_x_delta_t,
                            evmaxmag_real / (diffusivity * md_lambda2_x));
                    diffusive_y_delta_t   = min   (diffusive_y_delta_t,
                            evmaxmag_real / (diffusivity * md_lambda2_y));
                    diffusive_z_delta_t   = min   (diffusive_z_delta_t,
                            evmaxmag_real / (diffusivity * md_lambda2_z));
                    break;

                // Wall-normal implicit diffusion permits removing a
                // reference value in only the Y direction.  Notice
                // antidiffusive (nu - ref_nu) is fine but requires
                // chomping to zero to avoid circumventing the XZ checks,
                // however maxnan is required in case ref_nu is NaN.
                // Notice also that we /want/ to avoid NaN's arising from
                // diffusivity_y == 0 in diffusive_y_delta_t computation.
                case linearize::rhome_y:
                {
                    const real_t diffusivity_y
                            = maxnan(diffusivity - ref_nu, real_t(0));
                    diffusive_xyz_delta_t = minnan(diffusive_xyz_delta_t,
                              evmaxmag_real
                            / (   diffusivity   * md_lambda2_x
                                + diffusivity_y * md_lambda2_y
                                + diffusivity   * md_lambda2_z));
                    diffusive_x_delta_t   = min   (diffusive_x_delta_t,
                            evmaxmag_real / (diffusivity   * md_lambda2_x));
                    diffusive_y_delta_t   = min   (diffusive_y_delta_t,
                            evmaxmag_real / (diffusivity_y * md_lambda2_y));
                    diffusive_z_delta_t   = min   (diffusive_z_delta_t,
                            evmaxmag_real / (diffusivity   * md_lambda2_z));
                    break;
                }

                default:
                    SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();
                }
            }

        } // end X // end Z

    } // end Y
    SUZERAIN_TIMER_END("nonlinear right hand sides");

    // Traversal:
    // (3) Computing any manufactured solution forcing (when enabled).
    // Isolating this pass allows skipping the work when unnecessary
    if (msoln) {
        SUZERAIN_TIMER_SCOPED("manufactured forcing");

        // Dereference the msoln smart pointer outside the compute loop
        const ManufacturedSolution &ms = *msoln;

        for (int offset = 0, j = o.dgrid.local_physical_start.y();
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

                    real_t Q_rho, Q_rho_u, Q_rho_v, Q_rho_w, Q_rho_E;
                    ms.Q_conservative(x, y, z, time, Q_rho,
                                      Q_rho_u, Q_rho_v, Q_rho_w, Q_rho_E);

                    sphys(ndx::e,   offset) += Q_rho_E;
                    sphys(ndx::mx,  offset) += Q_rho_u;
                    sphys(ndx::my,  offset) += Q_rho_v;
                    sphys(ndx::mz,  offset) += Q_rho_w;
                    sphys(ndx::rho, offset) += Q_rho;

                } // end X

            } // end Z

        } // end Y

    } // end msoln

    // Collectively convert state to wave space using parallel FFTs
    for (size_t i = 0; i < swave_count; ++i) {

        if (Linearize == linearize::rhome_xyz && i == ndx::rho && !msoln) {

            // When density equation is handled fully implicitly AND no
            // manufactured solution is employed, save some communications by
            // using that the density right hand side is identically zero.
            assert(multi_array::is_contiguous(swave[i]));
            std::memset(swave[i].origin(), 0,
                    sizeof(complex_t)*o.dgrid.local_wave_extent.prod());

        } else {

            // Otherwise, bring the field back to wave space...
            o.dgrid.transform_physical_to_wave(&sphys.coeffRef(i,0));
            // ...and zero wavenumbers present only for dealiasing to
            // prevent "leakage" of dealiasing modes to other routines.
            o.zero_dealiasing_modes(swave, i);

        }

    }

    // Return the stable time step criteria separately on each rank.  The time
    // stepping logic must perform the Allreduce.  Delegating the Allreduce
    // responsibility allows reducing additional info with minimal overhead.
    return std::vector<real_t>(delta_t_candidates.begin(),
                               delta_t_candidates.end());

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

} // namespace perfect

} // namespace suzerain

#pragma warning(pop)

#endif  /* SUZERAIN_PERFECT_NAVIER_STOKES_HPP */
