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

#ifndef SUZERAIN_REACTING_NONLINEAR_OPERATOR_HPP
#define SUZERAIN_REACTING_NONLINEAR_OPERATOR_HPP

/** @file
 * Implementation of Nonlinear Navier--Stokes spatial operators
 * declared within nonlinear_operator_fwd.hpp.
 */

#include "nonlinear_operator_fwd.hpp"

#include <suzerain/error.h>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/support/logging.hpp>

#include "reacting.hpp"

#pragma warning(disable:280 383 1572)

namespace suzerain { namespace reacting {

  // Indices for state fields
  struct ndx { enum {
      e,
      mx,
      my,
      mz,
      rho,
      species
    }; };

  // Indices for directions
  struct dir { enum {
      y,
      x,
      z,
      count
    }; };

  // Indices for auxiliary storage.
  //
  // NOTE: Could eliminate 1 field here by accumulating T derivatives
  // back into T storage.  Makes code more complex for minimal gain,
  // so I didn't do it.
  //
  struct aux { enum {
      kap     = 0,
      T       = 1,
      gT      = 2,
      e       = 2 + 1*dir::count,
      mx      = 2 + 2*dir::count,
      my      = 2 + 3*dir::count,
      mz      = 2 + 4*dir::count,
      rho     = 2 + 5*dir::count,
      species = 2 + 6*dir::count
    }; };


template <bool ZerothSubstep,
          linearize::type Linearize,
          class ManufacturedSolution,
	  class ConstitutiveModels>
std::vector<real_t> apply_navier_stokes_spatial_operator(
            const operator_base &o,
            operator_common_block &common,
            const shared_ptr<const ManufacturedSolution>& msoln,
            const ConstitutiveModels& cmods,
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag)
{
    SUZERAIN_TIMER_SCOPED("applyNonlinearOperator");

    const size_t Ns = cmods.Ns();

    assert(Ns>0);

    // Shorthand
    typedef contiguous_state<4,complex_t> state_type;
    using std::abs;
    using std::equal;
    using std::max;
    using std::min;
    using std::size_t;
    using std::sqrt;

    // State enters method as coefficients in X, Y, and Z directions


   // Get total number of conserved state fields
    //const int state_count = Ns+4; // = (Ns-1) species + 5 (rho, ru, rv, rw, rE)
    const size_t state_count = Ns+4; // = (Ns-1) species + 5 (rho, ru, rv, rw, rE)

    // Get total number of fields in aux storage
    const size_t aux_count = (aux::species      +   // everything *except* species derivatives
                              dir::count*(Ns-1) );  // spatial derivatives of species

    // FIXME: Will have to change for reacting!!
    // We are only prepared to handle rho_E, rho_u, rho_v, rho_w, rho!
    enum { swave_count = 5 };
    assert(static_cast<int>(ndx::e  ) < swave_count);
    assert(static_cast<int>(ndx::mx ) < swave_count);
    assert(static_cast<int>(ndx::my ) < swave_count);
    assert(static_cast<int>(ndx::mz ) < swave_count);
    assert(static_cast<int>(ndx::rho) < swave_count);

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    scoped_ptr<state_type> _auxw_ptr(
	    support::allocate_padded_state<state_type>(
		aux_count, o.dgrid));                               // RAII
    state_type &auxw = *_auxw_ptr;                                  // Brevity

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
    if (ZerothSubstep) common.set_zero(/* Ny */ swave.shape()[1]);
    common.u().setZero();

    // Maintain stable time step values to return to the caller:
    //
    //   convtotal_* is convective stability using total velocity.
    //               It may or may not have acoustics included as necessary.
    //   convfluct_* is convective stability using only velocity fluctuations.
    //               It permits measuring any advantage in treating mean
    //               convection implicitly.
    //   diffusive_* is diffusive stability linearized however is appropriate.
    //   *_xyz_*     accounts for restrictions in all three directions.
    //               It must be more retrictive than any single direction.
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

    //********************************************************************
    //
    // From coefficients, compute values at collocation points for
    // all state and derivatives of state *except* derivatives of
    // total energy.
    //
    // This is the first wave to physical pass.  There will be another
    // to get the temperature gradient.
    //

    // NOTE: The indexing here *assumes* that the total energy is
    // stored first.

    // Energy (no derivatives)
    o.zero_dealiasing_modes(swave, ndx::e);
    o.bop_apply   (0,    1, swave, ndx::e);

    // Everything else (all spatial derivatives)
    for (size_t var = ndx::mx; var<state_count; ++var) {

      // Indexing note: aux::mx is beginning of derivatives of
      // conserved state in aux ordering.  Things appearing ahead of
      // aux::mx are not used until later.

      // Compute Y derivatives of variable var at collocation points
      // Zero wavenumbers present only for dealiasing along the way
      o.zero_dealiasing_modes(swave, var);
      o.bop_accumulate(1,    1, swave, var, 0, auxw, aux::e + dir::count*var + dir::y);
      o.bop_apply     (0,    1, swave, var);

      // Compute X- and Z- derivatives of variable var at collocation points
      // Zeros wavenumbers present only for dealiasing in the target storage
      o.diffwave_accumulate(1, 0, 1, swave, var,  0, auxw, aux::e + dir::count*var + dir::x );
      o.diffwave_accumulate(0, 1, 1, swave, var,  0, auxw, aux::e + dir::count*var + dir::z );
    }

    physical_view<> auxp (o.dgrid, auxw );
    physical_view<> sphys(o.dgrid, swave);
    for (size_t i = 0; i < state_count; ++i) {
      o.dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (size_t i = aux::mx; i < aux_count; ++i) {
      o.dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }

    //
    // Done with first wave to physical pass.
    //
    //********************************************************************

    // // Compute derived constants before inner loops
    // const real_t alpha13          = alpha + real_t(1)/real_t(3);
    // const real_t inv_Re           = 1 / Re;
    // const real_t inv_Ma2          = 1 / (Ma * Ma);
    // const real_t Ma2_over_Re      = (Ma * Ma) / Re;
    // const real_t inv_Re_Pr_gamma1 = 1 / (Re * Pr * (gamma - 1));
    const real_t lambda1_x        = o.lambda1_x;
    const real_t lambda1_z        = o.lambda1_z;
    const real_t lambda2_x        = o.lambda2_x;
    const real_t lambda2_z        = o.lambda2_z;
    //const real_t maxdiffconst     = inv_Re*max(gamma/Pr, max(real_t(1), alpha));

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

    //********************************************************************
    //
    // Physical space is traversed linearly using a single offset 'offset'.
    // The three loop structure is present to provide the global absolute
    // positions x(i), y(j), and z(k) where necessary.
    //
    // Four traversals occur:
    // (1) Computing reference quantities and mean velocity OR mean quantities
    //     (depending on linearization and which substep is being performed).
    //
    // (2) Computing as much of the sources and fluxes as possible
    //     given *only* conserved state and conserved state gradients.
    //     In theory this is everything.  In practice, we don't do the
    //     heat flux in this pass b/c the temperature gradient isn't
    //     available.
    //
    // (3) (After going to wave space with T and back down to physical
    //     space with grad(T)), Compute heat flux contribution.
    //
    // (4) Computing any manufactured solution forcing (when enabled).
    //
    // The traversal pattern is best embodied by the last pass.  The
    // first passes use slightly more compact structure as only y(j)
    // must be known within them.  Study (4) and convince yourself
    // that the rest are equivalent but lack information on x(i)
    // and z(k).
    //
    // TODO: Combine some of these traversals (e.g., 1 and 2)

    // Traversal:
    // (1) Computing reference quantities and mean velocity OR mean velocity
    //     (depending on linearization and which substep is being performed).
    //
    // Always gather all reference quantities when using implicit solves.
    // Profiling indicates the overhead is tiny and it keeps the code readable.
    if (    ZerothSubstep
         && Linearize != linearize::none) {  // References and mean velocity

        SUZERAIN_TIMER_SCOPED("reference quantities");
	
	// FIXME: Have to add this code back (w/ appropriate mods) to
	// support implicit reacting
	FATAL0("suzerain::reacting::apply_navier_stokes_spatial_operator only supports linearize::none");

    } else {                                 // Mean velocity profile only

        SUZERAIN_TIMER_SCOPED("mean velocity profile");

	// NOTE: I'm leaving this block alone w/out really looking at it.
	// TODO: Figure out if/how it must change for reacting.

        // Zero all reference quantities on fully-explicit zeroth substep
        if (ZerothSubstep && Linearize == linearize::none) {
            common.refs.setZero();
        }

        // Logic below here absolutely requires common.{u,v,w}() be housed
        // in the leftmost three columns of common.means.  Be sure.
        assert(common.means.data() == &common.u()[0]);
        assert(&common.u()[0] + common.means.rows() == &common.v()[0]);
        assert(&common.v()[0] + common.means.rows() == &common.w()[0]);

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t topNotOnRank = o.dgrid.local_physical_start.y();
        if (topNotOnRank) {
            common.means.leftCols<3>().topRows(topNotOnRank).setZero();
        }

        // Sum velocities as a function of y(j) into common.{u,v,w}()
        for (int offset = 0, j = o.dgrid.local_physical_start.y();
            j < o.dgrid.local_physical_end.y();
            ++j) {

            summing_accumulator_type ux;
            summing_accumulator_type uy;
            summing_accumulator_type uz;

            const int last_zxoffset = offset
                                    + o.dgrid.local_physical_extent.z()
                                    * o.dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {
                const real_t inv_rho = 1 / sphys(ndx::rho, offset);
                ux(inv_rho * sphys(ndx::mx, offset));
                uy(inv_rho * sphys(ndx::my, offset));
                uz(inv_rho * sphys(ndx::mz, offset));
            } // end X // end Z

            // Store sum into common block in preparation for MPI Reduce
            common.u()[j] = boost::accumulators::sum(ux);
            common.v()[j] = boost::accumulators::sum(uy);
            common.w()[j] = boost::accumulators::sum(uz);

        } // end Y

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t bottomNotOnRank = common.means.leftCols<3>().rows()
                                     - o.dgrid.local_physical_end.y();
        if (bottomNotOnRank) {
            common.means.leftCols<3>().bottomRows(bottomNotOnRank).setZero();
        }

        // Reduce, scale common.{u,v,w}() sums to obtain mean on zero-zero rank
        // Only zero-zero rank needs the information so Reduce is sufficient
        if (o.dgrid.has_zero_zero_modes()) {
            SUZERAIN_MPICHKR(MPI_Reduce(MPI_IN_PLACE,
                        common.means.leftCols<3>().data(),
                        common.means.leftCols<3>().size(),
                        mpi::datatype<real_t>::value, MPI_SUM,
                        o.dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
            common.means.leftCols<3>()
                    /= (   o.dgrid.global_physical_extent.x()
                         * o.dgrid.global_physical_extent.z());
        } else {
            ArrayXXr tmp;
            tmp.resizeLike(common.means.leftCols<3>());
            tmp.setZero();
            SUZERAIN_MPICHKR(MPI_Reduce(common.means.leftCols<3>().data(),
                        tmp.data(), common.means.leftCols<3>().size(),
                        mpi::datatype<real_t>::value, MPI_SUM,
                        o.dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
        }

    }

    // Traversal:
    // (2) Computing the nonlinear equation right hand sides.
    SUZERAIN_TIMER_BEGIN("nonlinear right hand sides (pass 1)");
    for (int offset = 0, j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {

        // FIXME: Put stable time step calculation back in.
        // // Wall-normal operator eigenvalue estimates depend on location
        // const real_t lambda1_y = o.lambda1_y(j);
        // const real_t lambda2_y = o.lambda2_y(j);

        // Iterate across the j-th ZX plane
        const int last_zxoffset = offset
                                + o.dgrid.local_physical_extent.z()
                                * o.dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // TODO: Add physics calculations!!!

	    // FIXME: Put stable time step calculation back in.

	    //
            // // Determine the minimum observed stable time step when necessary
            // // This logic used to call some canned routines, but additional
            // // monitoring requirements forced inlining many of these details.
            // // See delta_t_candidates declaration (above) for descriptions.
            // //
            // // details on these computations in the three-directional case.
            // if (ZerothSubstep) {

            //     // minnan(...) calls only required on *_xyz_delta_t as the
            //     // three-dimensional criterion should collect any NaNs
            //     // occurring within the direction-dependent criteria.  min(...)
            //     // is presumably a touch faster and is preferred when possible.
            //     using math::minnan;

            //     // See timestepper::convective_stability_criterion
            //     const real_t a = sqrt(T) / Ma;  // Because a/u_0 = sqrt(T*)/Ma
            //     real_t       ua_l1_x,       ua_l1_y,       ua_l1_z;
            //     real_t fluct_ua_l1_x, fluct_ua_l1_y, fluct_ua_l1_z;
            //     switch (Linearize) {
            //         default:
            //             SUZERAIN_ERROR_REPORT("Unimplemented!",
            //                                   SUZERAIN_ESANITY);

            //         // Explicit treatment forces including acoustics
            //         // in stability and has a zero reference velocity.
            //         case linearize::none:
            //             ua_l1_x       = (abs(u.x()) + a) * lambda1_x;
            //             ua_l1_y       = (abs(u.y()) + a) * lambda1_y;
            //             ua_l1_z       = (abs(u.z()) + a) * lambda1_z;
            //             fluct_ua_l1_x = ua_l1_x;
            //             fluct_ua_l1_y = ua_l1_y;
            //             fluct_ua_l1_z = ua_l1_z;
            //             break;

            //         // Implicit acoustics sets the effective sound speed
            //         // to zero within the convective_stability_criterion.
            //         // Fluctuating velocity is taken relative to references.
            //         case linearize::rhome:
            //             ua_l1_x       = abs(u.x()            ) * lambda1_x;
            //             ua_l1_y       = abs(u.y()            ) * lambda1_y;
            //             ua_l1_z       = abs(u.z()            ) * lambda1_z;
            //             fluct_ua_l1_x = abs(u.x() - ref_u.x()) * lambda1_x;
            //             fluct_ua_l1_y = abs(u.y() - ref_u.y()) * lambda1_y;
            //             fluct_ua_l1_z = abs(u.z() - ref_u.z()) * lambda1_z;
            //             break;
            //     }
            //     convtotal_xyz_delta_t = minnan(convtotal_xyz_delta_t,
            //             evmaxmag_imag / (ua_l1_x + ua_l1_y + ua_l1_z));
            //     convtotal_x_delta_t   = min   (convtotal_x_delta_t,
            //             evmaxmag_imag / ua_l1_x);
            //     convtotal_y_delta_t   = min   (convtotal_y_delta_t,
            //             evmaxmag_imag / ua_l1_y);
            //     convtotal_z_delta_t   = min   (convtotal_z_delta_t,
            //             evmaxmag_imag / ua_l1_z);
            //     convfluct_xyz_delta_t = minnan(convfluct_xyz_delta_t,
            //             evmaxmag_imag / (  fluct_ua_l1_x
            //                              + fluct_ua_l1_y
            //                              + fluct_ua_l1_z));
            //     convfluct_x_delta_t   = min   (convfluct_x_delta_t,
            //             evmaxmag_imag / fluct_ua_l1_x);
            //     convfluct_y_delta_t   = min   (convfluct_y_delta_t,
            //             evmaxmag_imag / fluct_ua_l1_y);
            //     convfluct_z_delta_t   = min   (convfluct_z_delta_t,
            //             evmaxmag_imag / fluct_ua_l1_z);

            //     // See timestepper::diffusive_stability_criterion
            //     // Antidiffusive locations might be ignored when linearized.
            //     // Hence we compute criteria within the switch statment.
            //     const real_t nu = mu / rho;
            //     real_t diffusivity;
            //     switch (Linearize) {
            //         default:
            //             SUZERAIN_ERROR_REPORT("Unimplemented!",
            //                                   SUZERAIN_ESANITY);

            //         // Explicit treatment forces a zero reference diffusivity
            //         case linearize::none:
            //             diffusivity = maxdiffconst * nu;
            //             diffusive_xyz_delta_t = minnan(diffusive_xyz_delta_t,
            //                       evmaxmag_real
            //                     / diffusivity
            //                     / (lambda2_x + lambda2_y + lambda2_z));
            //             diffusive_x_delta_t   = min   (diffusive_x_delta_t,
            //                     evmaxmag_real / diffusivity / lambda2_x);
            //             diffusive_y_delta_t   = min   (diffusive_y_delta_t,
            //                     evmaxmag_real / diffusivity / lambda2_y);
            //             diffusive_z_delta_t   = min   (diffusive_z_delta_t,
            //                     evmaxmag_real / diffusivity / lambda2_z);
            //             break;

            //         // Implicit diffusion permits removing a reference value.
            //         // Antidiffusive (nu - ref_nu) is fine and not computed.
            //         case linearize::rhome:
            //             diffusivity = nu - ref_nu;    // Compute sign wrt ref.
            //             if (diffusivity <= 0) break;  // NaN => false, proceed
            //             diffusivity *= maxdiffconst;  // Rescale as necessary.
            //             diffusive_xyz_delta_t = minnan(diffusive_xyz_delta_t,
            //                       evmaxmag_real
            //                     / diffusivity
            //                     / (lambda2_x + lambda2_y + lambda2_z));
            //             diffusive_x_delta_t   = min   (diffusive_x_delta_t,
            //                     evmaxmag_real / diffusivity / lambda2_x);
            //             diffusive_y_delta_t   = min   (diffusive_y_delta_t,
            //                     evmaxmag_real / diffusivity / lambda2_y);
            //             diffusive_z_delta_t   = min   (diffusive_z_delta_t,
            //                     evmaxmag_real / diffusivity / lambda2_z);
            //             break;
            //     } // end switch(Linearize)
            // } // end if(ZerothSubstep}

        } // end X // end Z

    } // end Y
    SUZERAIN_TIMER_END("nonlinear right hand sides (pass 1)");

    //------------------------------------------------------------------
    //
    // At this point, we have done everything except the heat flux.
    // To do this, we will
    //
    // (1) Take temperature from physical to wave space.
    // (2) Differentiate temperature in wave space and bring grad(T)
    //     back to physical space.
    //
    // This is the second transform stage.


    // (1) Temperature from physical to wave space
    o.dgrid.transform_physical_to_wave(&auxp.coeffRef(aux::T,0));


    // TODO: Apply inverse mass matrix to get to pure coefficient space

    // ...and zero wavenumbers present only for dealiasing
    o.zero_dealiasing_modes(auxw, aux::T);


    // (2) Get derivatives.
    //
    // Compute Y derivatives of variable var at collocation points
    o.bop_accumulate(1,    1, auxw, aux::T, 0, auxw, aux::gT + dir::y);

    // Compute X- and Z- derivatives of variable var at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    o.diffwave_accumulate(1, 0, 1, auxw, aux::T,  0, auxw, aux::gT + dir::x );
    o.diffwave_accumulate(0, 1, 1, auxw, aux::T,  0, auxw, aux::gT + dir::z );

    // FFTs to get to physical space
    o.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::y,0));
    o.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::x,0));
    o.dgrid.transform_wave_to_physical(&auxp.coeffRef(aux::gT+dir::z,0));

    // Now have temperature gradient at collocation points.
    //-----------------------------------------------------------------

    // Traversal:
    // (3) (After going back to wave space with T and back down the
    //      grad(T)), Compute heat flux contribution.
    SUZERAIN_TIMER_BEGIN("nonlinear right hand sides (pass 2)");
    for (int offset = 0, j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {

      // Iterate across the j-th ZX plane
      const int last_zxoffset = offset
                              + o.dgrid.local_physical_extent.z()
                              * o.dgrid.local_physical_extent.x();
      for (; offset < last_zxoffset; ++offset) {

        // Thermal conductivity
        const real_t   kap       ( auxp(aux::kap     ,   offset));

        // Extract grad(T)
        const Vector3r grad_T    ( auxp(aux::gT+dir::x,  offset),
                                   auxp(aux::gT+dir::y,  offset),
                                   auxp(aux::gT+dir::z,  offset));

        // Accumulate into energy fluxes
        //
        // NOTE: Sign correct for fluxes appearing on the LHS---i.e.,
        // U_t + div(F) = S(U).
        auxp(aux::e +dir::x, offset) -= kap*grad_T.x();
        auxp(aux::e +dir::y, offset) -= kap*grad_T.y();
        auxp(aux::e +dir::z, offset) -= kap*grad_T.z();


      } // end x,z
    } // end y
    SUZERAIN_TIMER_END("nonlinear right hand sides (pass 2)");

    // Traversal:
    // (4) Computing any manufactured solution forcing (when enabled).
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
    		    // FIXME(?): nsctpl.hpp does not provide this
    		    // method (where nsctpl_rholut did). Could add it.
    		    // For now however, I'm choosing to refactor here
    		    // to call individual src fcns to avoid drift from
    		    // masa src.
    		    //
    		    ms.Q_conservative(x, y, z, time, Q_rho,
    				      Q_rho_u, Q_rho_v, Q_rho_w, Q_rho_E);

    		    // Q_rho   = ms.Q_rho (x, y, z, time);
    		    // Q_rho_u = ms.Q_rhou(x, y, z, time);
    		    // Q_rho_v = ms.Q_rhov(x, y, z, time);
    		    // Q_rho_w = ms.Q_rhow(x, y, z, time);
    		    // Q_rho_E = ms.Q_rhoe(x, y, z, time);

                    sphys(ndx::e,   offset) += Q_rho_E;
                    sphys(ndx::mx,  offset) += Q_rho_u;
                    sphys(ndx::my,  offset) += Q_rho_v;
                    sphys(ndx::mz,  offset) += Q_rho_w;
                    sphys(ndx::rho, offset) += Q_rho;

                } // end X

            } // end Z

        } // end Y

    } // end msoln


    // TODO: Anything related to BCs that must be done in physical
    // space.

    // At this point we have accumulated sources (into sphys/swave)
    // and fluxes (into auxp/auxw) in physical space.  Now we need to
    // get back to wave space to finish accumulating the RHS.


    // Collectively convert source to wave space using parallel FFTs
    for (size_t i = 0; i < state_count; ++i) {
      o.dgrid.transform_physical_to_wave(&sphys.coeffRef(i,0));
    }

    // Collectively convert fluxes to wave space using parallel FFTs
    for (size_t i = aux::e; i < aux_count; ++i) {
      o.dgrid.transform_physical_to_wave(&auxp.coeffRef(i,0));
    }


    // Now we have sources and fluxes in wave space.  More
    // specifically, we have X,Z coefficients and Y collocation pts.
    // The last piece is to apply the divergence to the fluxes and
    // accumulate into the sources.

    // Note that the divergence is never formed explicitly.  We
    // accumulate the Y, X, and Z derivatives into the source
    // separately.
    {

      suzerain::bsplineop_luz massluz(o.cop);
      //const complex_t scale_factor = grid.dN.x() * grid.dN.z();
      const complex_t scale_factor = 1.0; //grid.dN.x() * grid.dN.z();
      massluz.opform(1, &scale_factor, o.cop);
      massluz.factor();

      //suzerain::bsplineop_lu boplu(o.cop);

      // FIXME: Need form mass matrix call
      // FIXME: Move these ops outside this function

      for (size_t i = 0; i < state_count; ++i) {

	// FIXME: This doesn't compile... input types incorrect
        // Apply inverse mass matrix to Y flux to get to pure
        // coefficient representation
        //o.bop_solve(boplu, auxw, aux::e + dir::count*i + dir::y);
	o.bop_solve(massluz, auxw, aux::e + dir::count*i + dir::y);

        //o.diffwave_apply(0, 0, 1, auxw, aux::e + dir::count*i + dir::y);

        // Accumulate Y derivative of Y flux (at collocation pts in Y
        // and coefficients in X,Z) into source.  Note that alpha = -1
        // b/c we need to subtract divergence from source
        o.bop_accumulate(1,   -1, auxw , aux::e + dir::count*i + dir::y,
                               1, swave, i );

        // Accumulate X and Z derivatives of X and Z fluxes into source

        // alpha = -1 b/c we need to subtract divergence from source
        o.diffwave_accumulate(1, 0, -1, auxw , aux::e + dir::count*i + dir::x,
                                     1, swave, i );

        // alpha = -1 b/c we need to subtract divergence from source
        o.diffwave_accumulate(0, 1, -1, auxw , aux::e + dir::count*i + dir::z,
                                     1, swave, i );

        // and zero wavenumbers present only for dealiasing to
        // prevent "leakage" of dealiasing modes to other routines.
        o.zero_dealiasing_modes(swave, i);

      } // end for

    } // end accumulate

    // Return the stable time step criteria separately on each rank.  The time
    // stepping logic must perform the Allreduce.  Delegating the Allreduce
    // responsibility allows reducing additional info with minimal overhead.
    return std::vector<real_t>(delta_t_candidates.begin(),
                               delta_t_candidates.end());

    // State leaves method as coefficients in X and Z directions
    // State leaves method as collocation point values in Y direction
}

} /* namespace reacting */ } /* namespace suzerain */

#endif  /* SUZERAIN_REACTING_NONLINEAR_OPERATOR_HPP */
