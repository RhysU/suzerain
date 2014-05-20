//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#include <suzerain/error.h>
#include <suzerain/largo_state.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/math.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/slowgrowth.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/timers.h>
#include <suzerain/utility.hpp>

#include "definition_scenario.hpp"
#include "operator_common_block.hpp"
#include "perfect.hpp"

#pragma warning(push, disable:280 383 1572)

namespace suzerain {

namespace perfect {

/**
 * A complete Navier&ndash;Stokes \c apply_operator implementation.  The
 * implementation is provided as a common building block for
 * <tt>lowstorage::operator_nonlinear< contiguous_state<4,complex_t> ></tt>
 * subclasses allowing varying numbers of passive scalars or varying hybrid
 * implicit/explicit treatment.  Such subclasses feature an overwhelming amount
 * of redundancy and are error prone to create.  This implementation allows
 * writing the "futsy" bits once and then sharing the logic repeatedly.
 * Templating allows for compile-time branching amongst different
 * implementation choices rather than paying runtime cost for such flexibility.
 *
 * Computation follows the "Numerical Considerations" section of
 * writeups/perfect_gas.tex.  No boundary conditions are applied.  Gathers
 * information required per \ref operator_common_block documentation.
 *
 * \param scenario Nondimensional scenario parameters.
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
 * \param substep_index The index of the current substep being taken.
 *
 * \tparam ZerothSubstep Should one-time activities taking place at the
 *         beginning of a Runge-Kutta step be performed?  Examples include
 *         computing a stable time step size and also computing reference
 *         quantities for linearization.
 * \tparam Linearize What type of hybrid implicit/explicit linearization
 *         is employed?
 * \tparam SlowGrowth Is slow growth forcing being applied?
 * \tparam ManufacturedSolution What manufactured solution should be used to
 *         provide additional forcing (when enabled)?
 *
 * @return A vector of stable timestep sizes according to different criteria
 *         per lowstorage::operator_nonlinear::apply_operator.
 *
 * @see lowstorage::operator_nonlinear for the (slightly different)
 *      interface that an actual operator would provide.
 */
template <bool ZerothSubstep,
          linearize::type Linearize,
          bool SlowGrowth,
          class ManufacturedSolution>
std::vector<real_t> apply_navier_stokes_spatial_operator(
            const definition_scenario &scenario,
            const operator_base &o,
            operator_common_block &common,
            const specification_largo &sg,
            const shared_ptr<const ManufacturedSolution> &msoln,
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index)
{
    // State enters method as coefficients in X, Y, and Z directions
    SUZERAIN_TIMER_SCOPED("apply_navier_stokes_spatial_operator");

    // Shorthand
    const real_t alpha = scenario.alpha;
    const real_t beta  = scenario.beta;
    const real_t gamma = scenario.gamma;
    const real_t Ma    = scenario.Ma;
    const real_t Pr    = scenario.Pr;
    const real_t Re    = scenario.Re;
    typedef contiguous_state<4,complex_t> state_type;
    using std::abs;
    using std::equal;
    using std::max;
    using std::min;
    using std::numeric_limits;
    using std::size_t;
    using std::sqrt;

    // Compile-time template parameters are used to reduce jumps at runtime.
    // Ensure someone didn't hand in a mismatched substep, common block or sg.
    SUZERAIN_ENSURE(ZerothSubstep == (substep_index == 0));
    SUZERAIN_ENSURE(Linearize == common.linearization);
    SUZERAIN_ENSURE(SlowGrowth == sg.formulation.enabled());

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
    SUZERAIN_ENSURE((int) swave.shape()[0] == swave_count);
    SUZERAIN_ENSURE((int) swave.shape()[1] == o.dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE((int) swave.shape()[2] == o.dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE((int) swave.shape()[3] == o.dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE(swave.strides()[1] == (int) 1);
    SUZERAIN_ENSURE(swave.strides()[2] == (int) swave.shape()[1]);
    SUZERAIN_ENSURE(swave.strides()[3] == (int) ( swave.shape()[1]
                                                 *swave.shape()[2]));
    SUZERAIN_ENSURE(equal(swave.shape()   + 1, swave.shape()   + 4,
                          auxw.shape()    + 1));
    SUZERAIN_ENSURE(equal(swave.strides() + 1, swave.strides() + 4,
                          auxw.strides()  + 1));
    const size_t Ny = swave.shape()[1];

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
              numeric_limits<real_t>::max());
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

    // Initialize slow growth treatment if necessary
    slowgrowth sg_treater(sg, Ma);
    if (SlowGrowth) sg_treater.initialize(substep_index);

    // With conserved state Fourier in X and Z but collocation in Y,
    // gather any RMS-like information necessary for slow growth forcing.
    if (SlowGrowth) sg_treater.gather_wavexz(o, swave);

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
        collect_references(scenario, o.grid, o.dgrid, sphys, common.ref);
        common.sub = common.ref;
    } else {
        collect_instantaneous(scenario, o.grid, o.dgrid, sphys, common.sub);
    }

    // Slow growth requires mean conserved state at collocation points,
    // mean pressure and pressure fluctuation profiles, "rhoqq" values for
    // q \in {u, v, w, E}, and wall-normal derivatives of everything.
    if (SlowGrowth) sg_treater.gather_physical_cons(o, common.sub);
    if (SlowGrowth) sg_treater.gather_physical_rqq (o, common.sub);

    // Compute inviscid-only constants before any Y, X, or Z inner loops
    const real_t inv_Ma2   = 1 / (Ma * Ma);
    const real_t gamma1    = gamma - 1;
    const real_t lambda1_x = o.lambda1_x;
    const real_t lambda1_z = o.lambda1_z;

    // Traversal:
    // (2) Computing the nonlinear equation right hand sides.
    SUZERAIN_TIMER_BEGIN("nonlinear right hand sides");
    for (int offset = 0, j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {

        // Prepare any y-dependent slow growth computation
        if (SlowGrowth) sg_treater.inner_y(j, o.y(j));

        // Precompute subtractive slow growth velocity for wall-normal
        // convective stability criterion (Euler_Eigensystem_3D_Temporal.nb)
        // assuming not one bit of slow growth is handled implicitly.
        const real_t y_grdelta = SlowGrowth ? o.y(j)*sg.grdelta : 0;

        // Unpack appropriate wall-normal reference quantities
        const Vector3r ref_u              (common.ref.u         ()[j],
                                           common.ref.v         ()[j],
                                           common.ref.w         ()[j]);
        const real_t   ref_u2             (common.ref.u2        ()[j]);
        const Matrix3r ref_uu;
        const_cast<Matrix3r&>(ref_uu) <<   common.ref.uu        ()[j],
                                           common.ref.uv        ()[j],
                                           common.ref.uw        ()[j],
                                           common.ref.uv        ()[j],
                                           common.ref.vv        ()[j],
                                           common.ref.vw        ()[j],
                                           common.ref.uw        ()[j],
                                           common.ref.vw        ()[j],
                                           common.ref.ww        ()[j];
        const real_t   ref_nu             (common.ref.nu        ()[j]);
        const Vector3r ref_nuu            (common.ref.nu_u      ()[j],
                                           common.ref.nu_v      ()[j],
                                           common.ref.nu_w      ()[j]);
        const real_t   ref_nuu2           (common.ref.nu_u2     ()[j]);
        const Matrix3r ref_nuuu;
        const_cast<Matrix3r&>(ref_nuuu) << common.ref.nu_uu     ()[j],
                                           common.ref.nu_uv     ()[j],
                                           common.ref.nu_uw     ()[j],
                                           common.ref.nu_uv     ()[j],
                                           common.ref.nu_vv     ()[j],
                                           common.ref.nu_vw     ()[j],
                                           common.ref.nu_uw     ()[j],
                                           common.ref.nu_vw     ()[j],
                                           common.ref.nu_ww     ()[j];
        const Vector3r ref_e_gradrho      (common.ref.ex_gradrho()[j],
                                           common.ref.ey_gradrho()[j],
                                           common.ref.ez_gradrho()[j]);
        const real_t   ref_e_divm         (common.ref.e_divm    ()[j]);
        const real_t   ref_e_deltarho     (common.ref.e_deltarho()[j]);

        // Redmine #2983 disables mu and lambda on the freestream boundary
        // to attempt adhering to Poinsot and Lele subsonic NRBC conditions
        const bool locallyviscous = !(o.grid.one_sided() && j+1U == Ny);

        // Compute viscous-related constants before X or Z inner loops
        // (as this permits easily making a single XZ plan inviscid)
        const real_t alpha13          = alpha + real_t(1)/real_t(3);
        const real_t alpha43          = alpha + real_t(4)/real_t(3);
        const real_t inv_Re           = static_cast<int>(locallyviscous) / Re;
        const real_t Ma2_over_Re      = (Ma * Ma) * inv_Re;
        const real_t inv_Re_Pr_gamma1 = 1 / (Pr * gamma1) * inv_Re;
        const real_t gamma_over_Pr    = gamma / Pr;
        const real_t maxdiffconst     = inv_Re
                                      * max(gamma/Pr, max(real_t(1), alpha));
        const real_t md_lambda2_x     = maxdiffconst * o.lambda2_x;
        const real_t md_lambda2_z     = maxdiffconst * o.lambda2_z;

        // Wall-normal operator eigenvalue estimates depend on location
        const real_t lambda1_y    = o.lambda1_y(j);
        const real_t md_lambda2_y = maxdiffconst * o.lambda2_y(j);

        // Prepare to iterate across the j-th ZX plane
        const int last_zxoffset = offset
                                + o.dgrid.local_physical_extent.z()
                                * o.dgrid.local_physical_extent.x();

        // Iterate across the j-th ZX plane
        for (; offset < last_zxoffset; ++offset) {

            // Unpack local conserved state
            const real_t   e      (sphys(ndx::e,   offset));
            const Vector3r m      (sphys(ndx::mx,  offset),
                                   sphys(ndx::my,  offset),
                                   sphys(ndx::mz,  offset));
            const real_t   rho    (sphys(ndx::rho, offset));

            // Unpack total energy-related derivatives
            const Vector3r grad_e  (auxp(aux::e_x,        offset),
                                    auxp(aux::e_y,        offset),
                                    auxp(aux::e_z,        offset));
            const real_t e_yy      (auxp(aux::e_yy,       offset));
            const real_t div_grad_e(auxp(aux::div_grad_e, offset));

            // Unpack momentum-related derivatives
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

            // Unpack density-related derivatives
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

            // Compute velocity-related quantities as well as
            // quantities related to the viscous stress tensor
            const Vector3r u          = rholut::u(
                                            rho, m);
            const Matrix3r grad_u     = rholut::grad_u(
                                            rho, grad_rho, m, grad_m);
            const real_t   div_u      = grad_u.trace();
            const Matrix3r tau        = rholut::tau(
                                            mu, lambda, div_u, grad_u);
            const real_t div_m        = grad_m.trace();
            const Vector3r grad_div_u = rholut::grad_div_u(
                                            rho, grad_rho, grad_grad_rho,
                                            m, div_m, grad_m, grad_div_m);
            const Vector3r div_grad_u = rholut::div_grad_u(
                                            rho, grad_rho, div_grad_rho,
                                            m, grad_m, div_grad_m);
            const Vector3r div_tau    = rholut::div_tau(
                                            mu, grad_mu, lambda, grad_lambda,
                                            div_u, grad_u, div_grad_u,
                                            grad_div_u);

            // COMPUTE ANY SLOW GROWTH FORCING APPLIED TO THE FLOW
            //
            // Accumulation into each scalar right hand side is delayed
            // to improve striding when accessing sphys buffers below.
            largo_state slowgrowth_f;
            if (SlowGrowth) {
                largo_state state(e, m.x(), m.y(), m.z(), rho, p);
                sg_treater.inner_xz(state, slowgrowth_f);
            }

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
            sphys(ndx::e, offset) += slowgrowth_f.e;

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
            sphys(ndx::mx, offset) = momentum_rhs.x() + slowgrowth_f.mx;
            sphys(ndx::my, offset) = momentum_rhs.y() + slowgrowth_f.my;
            sphys(ndx::mz, offset) = momentum_rhs.z() + slowgrowth_f.mz;

            // FORM CONTINUITY EQUATION RIGHT HAND SIDE
            //
            // Implicit continuity equation handling requires zeroing RHS in
            // anticipation of possible manufactured solution forcing and/or
            // slow growth forcing.  See subsequent transform_physical_to_wave
            // if you monkey around here.
            switch (Linearize) {
            case linearize::rhome_xyz:    // Fully implicit convection
                sphys(ndx::rho, offset) = 0;
                break;

            case linearize::rhome_y:      // Implicit Y but Explicit X, Z
                sphys(ndx::rho, offset) = - grad_m(0,0) - grad_m(2,2);
                break;

            case linearize::none:         // Fully explicit convection
                sphys(ndx::rho, offset) = - div_m;
                break;

            default:
                SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();
            }
            sphys(ndx::rho, offset) += slowgrowth_f.rho;

            // Determine the minimum observed stable time step when necessary
            // This logic used to call some canned routines, but additional
            // monitoring requirements forced inlining many of these details.
            // See delta_t_candidates declaration (above) for descriptions.
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
                const real_t a   = sqrt(T) / Ma;      // As a/u_0 = sqrt(T*)/Ma
                const real_t u_y = u.y() - y_grdelta; // Slow growth adjustment
                real_t       ua_l1_x,       ua_l1_y,       ua_l1_z;
                real_t fluct_ua_l1_x, fluct_ua_l1_y, fluct_ua_l1_z;
                switch (Linearize) {
                // Implicit acoustics sets the effective sound speed
                // to zero within the convective_stability_criterion.
                // Fluctuating velocity is taken relative to references.
                case linearize::rhome_xyz:
                    ua_l1_x       = abs(u.x()            ) * lambda1_x;
                    ua_l1_y       = abs(u_y              ) * lambda1_y;
                    ua_l1_z       = abs(u.z()            ) * lambda1_z;
                    fluct_ua_l1_x = abs(u.x() - ref_u.x()) * lambda1_x;
                    fluct_ua_l1_y = abs(u_y   - ref_u.y()) * lambda1_y;
                    fluct_ua_l1_z = abs(u.z() - ref_u.z()) * lambda1_z;
                    break;

                // Explicit treatment forces including acoustics
                // in stability and has a zero reference velocity.
                case linearize::none:
                    ua_l1_x       = (abs(u.x()) + a) * lambda1_x;
                    ua_l1_y       = (abs(u_y  ) + a) * lambda1_y;
                    ua_l1_z       = (abs(u.z()) + a) * lambda1_z;
                    fluct_ua_l1_x = ua_l1_x;
                    fluct_ua_l1_y = ua_l1_y;
                    fluct_ua_l1_z = ua_l1_z;
                    break;

                // Wall-normal-only implicit acoustics and convection
                // is nothing but a hybrid of the above two cases.
                case linearize::rhome_y:
                    ua_l1_x       = (abs(u.x()) + a) * lambda1_x;
                    ua_l1_y       = (abs(u_y  )    ) * lambda1_y;
                    ua_l1_z       = (abs(u.z()) + a) * lambda1_z;
                    fluct_ua_l1_x = ua_l1_x;
                    fluct_ua_l1_y = abs(u_y   - ref_u.y()) * lambda1_y;
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

        if (   Linearize == linearize::rhome_xyz  // Fully implicit density
            && i == ndx::rho                      // Looking at density RHS
            && !SlowGrowth                        // No slow growth forcing
            && !msoln) {                          // No manufactured solution

            // When the above conditions are all met, save some communications
            // by using that the density right hand side is identically zero.
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
