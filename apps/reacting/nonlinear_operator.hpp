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
#include "filter_definition.hpp"
#include "reacting_ndx.hpp"

#pragma warning(disable:280 383 1572)

namespace suzerain { namespace reacting {

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
          filter::type Filter,
          class ManufacturedSolution,
          class ConstitutiveModels>
std::vector<real_t> apply_navier_stokes_spatial_operator(
            const operator_base &o,
            operator_common_block &common,
            const filter_definition &fsdef,
            const shared_ptr<const ManufacturedSolution>& msoln,
            const ConstitutiveModels& cmods,
            const suzerain::bsplineop_luz& massluz,
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
    const size_t state_count = Ns+4; // = (Ns-1) species + 5 (rho, ru, rv, rw, rE)

    // Get total number of fields in aux storage
    const size_t aux_count = (aux::species      +   // everything *except* species derivatives
                              dir::count*(Ns-1) );  // spatial derivatives of species

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    scoped_ptr<state_type> _auxw_ptr(
            support::allocate_padded_state<state_type>(
                aux_count, o.dgrid));                               // RAII
    state_type &auxw = *_auxw_ptr;                                  // Brevity

    // Auxiliary storage for filter source:
    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    // Number of fields for filter source is equal to number of variables
    scoped_ptr<state_type> _fsrcw_ptr(
            support::allocate_padded_state<state_type>(
                state_count, o.dgrid));                               // RAII
    state_type &fsrcw = *_fsrcw_ptr;                                  // Brevity

    // Sanity check incoming swave's and auxw's shape and contiguity
    SUZERAIN_ENSURE(swave.shape()[0] == state_count);
    SUZERAIN_ENSURE(swave.shape()[1] == (unsigned) o.dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(swave.shape()[2] == (unsigned) o.dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(swave.shape()[3] == (unsigned) o.dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE((unsigned) swave.strides()[1] == 1u);
    SUZERAIN_ENSURE((unsigned) swave.strides()[2] == swave.shape()[1]);
    SUZERAIN_ENSURE((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    SUZERAIN_ENSURE(equal(swave.shape()   + 1, swave.shape()   + 4, auxw.shape()   + 1));
    SUZERAIN_ENSURE(equal(swave.strides() + 1, swave.strides() + 4, auxw.strides() + 1));
    SUZERAIN_ENSURE(equal(swave.shape()   + 1, swave.shape()   + 4, fsrcw.shape()   + 1));
    SUZERAIN_ENSURE(equal(swave.strides() + 1, swave.strides() + 4, fsrcw.strides() + 1));

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

    if (Filter == filter::viscous) {
        //**************************************************************
        // Begin computing viscous filter.  There are five steps.  The
        // first four must done here, while the state (swave) is still
        // entirely coefficients.  The last step (scaling by reference
        // quantities) must be done after the reference quantities are
        // computed below (of course).
        //
        // Each step is performed independently for each state variable
        for (size_t var = ndx::e; var<state_count; ++var) {

            // 1.) Accumulate D_1 swave in to fsrcw
            o.bop_accumulate(1, 1, swave, var,
                             0, fsrcw, var );

            // 2.) fsrcw <- M \ fsrcw
            o.bop_solve(massluz, fsrcw, var);

            // 3.) fsrcw <- D_1 fsrcw + 0 * fsrcw
            o.bop_apply(1, 1, fsrcw, var);

            // 4.) fsrcw <- D_2 swave - fsrcw
            o.bop_accumulate(2,  1, swave, var,
                             -1, fsrcw, var );

        }
        // Now, fsrcw contains D_2 * swave - D_1 * M \ (D_1 * swave).
        //
        // To complete the filter source, we need to scale it by the
        // appropriate reference.
        //**************************************************************
    }


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


    if (Filter == filter::cook) {
        // FIXME: to start testing, use directly a user defined value;
        //        to do is to figure out proper coefficients for each variable;
        //        see if it works to set the coefficient for rho = 0;

        const complex_t f_phi = fsdef.filter_phi;

        // Filter source: compute for energy
        fsdef.source_accumulate(o.grid, o.dgrid, f_phi, swave, ndx::e,
                                                 0.,    fsrcw, ndx::e);
    }

    // Everything else (all spatial derivatives)
    for (size_t var = ndx::mx; var<state_count; ++var) {

      // Indexing note: aux::mx is beginning of derivatives of
      // conserved state in aux ordering.  Things appearing ahead of
      // aux::mx are not used until later.

      // Compute Y derivatives of variable var at collocation points
      // Zero wavenumbers present only for dealiasing along the way
      o.zero_dealiasing_modes(swave, var);
      o.bop_accumulate(1,    1, swave, var,
                             0, auxw , aux::e + dir::count*var + dir::y);
      o.bop_apply     (0,    1, swave, var);

      // Compute X- and Z- derivatives of variable var at collocation points
      // Zeros wavenumbers present only for dealiasing in the target storage
      o.diffwave_accumulate(1, 0, 1, swave, var,
                                  0, auxw , aux::e + dir::count*var + dir::x );
      o.diffwave_accumulate(0, 1, 1, swave, var,
                                  0, auxw , aux::e + dir::count*var + dir::z );

      if (Filter == filter::cook) {
          const complex_t f_phi = fsdef.filter_phi;
          // Filter source: compute for variable var
          fsdef.source_accumulate(o.grid, o.dgrid, f_phi, swave, var,
                                                   0.,    fsrcw, var);
      }
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

    // Compute derived constants before inner loops
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

    // Before any traversals, instantiate Eigen variable length arrays to store species
    // specific info
    VectorXr species(Ns); // species densities
    VectorXr Ds     (Ns); // mass diffusivities
    VectorXr hs     (Ns); // species enthalpies
    VectorXr om     (Ns); // reaction source terms
    VectorXr cs     (Ns); // species mass fractions

    Matrix3Xr grad_species (3,Ns); // spatial derivatives of species densities
    Matrix3Xr grad_cs      (3,Ns); // spatial derivatives of species mass fractions
    Matrix3Xr sdiff        (3,Ns); // diffusive fluxes for species equations


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
    // (1) Computing reference quantities and mean velocity OR mean
    //     velocity (depending on linearization, filtering, and which
    //     substep is being performed).
    //
    // Always gather all reference quantities when using implicit solves.
    // Profiling indicates the overhead is tiny and it keeps the code readable.
    //
    // Also need reference quantities if we are using the viscous filter option.
    if (    ZerothSubstep
        && (Linearize != linearize::none || Filter == filter::viscous) ) {

        // Gather reference profiles and mean velocity

        SUZERAIN_TIMER_SCOPED("reference quantities");

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t leftNotOnRank = o.dgrid.local_physical_start.y();
        if (leftNotOnRank) common.refs.leftCols(leftNotOnRank).setZero();

        // Sum reference quantities as a function of y(j) into common.ref_*
        // See writeups/derivation.tex or rholut_imexop.h for definitions
        for (int offset = 0, j = o.dgrid.local_physical_start.y();
            j < o.dgrid.local_physical_end.y();
            ++j) {

            // Prepare logical indices using a struct for scoping (e.g. ref::ux).
            struct ref { enum { ux, uy, uz, uxuy, uzuy,
                                p_ru, p_rw, p_rE,
                                vp_ru, vp_rw, vp_rE,
                                Cmy_rho, Ce_rho, Ce_rv,
                                nu, korCv, Ds,
                                count // Sentry
            }; };

            // An array of summing_accumulator_type holds all running sums.
            // This gives nicer construction and allows looping over results.
            summing_accumulator_type acc[ref::count];

            const int last_zxoffset = offset
                                    + o.dgrid.local_physical_extent.z()
                                    * o.dgrid.local_physical_extent.x();
            for (; offset < last_zxoffset; ++offset) {

                // First, unpack conserved state...

                // ... mixture density
                const real_t   rho         ( sphys(ndx::rho,    offset));
                const real_t  irho = 1.0/rho;

                // ... momentum
                const Vector3r m    ( sphys(ndx::mx, offset),
                                      sphys(ndx::my, offset),
                                      sphys(ndx::mz, offset));

                // ... total energy
                const real_t e        (sphys(ndx::e,       offset));

                // ... species densities
                // NOTE: In species vector, idx 0 is the dilluter (the species
                // that is not explicitly part of the state vector)
                species(0) = rho;
                for (unsigned int s=1; s<Ns; ++s) {
                    species(s) = sphys(ndx::species + s - 1, offset);

                    // dilluter density = rho_0 = rho - sum_{s=1}^{Ns-1} rho_s
                    species(0)        -= species(s);
                }

                // Second, compute required derived quantities...

                // ... mass fractions
                for (unsigned int s=0; s<Ns; ++s) {
                    cs(s) = irho * species(s);
                }

                // ... velocity
                const Vector3r u     = suzerain::rholut::u(rho, m);

                // ... pressure-related quantities
                real_t p, p_rho, p_rsum, p_e, mu, korCv;
                Vector3r p_m;
                cmods.evaluate_pressure_derivs_and_trans(
                    e, m, rho, species, cs,
                    p, p_rho, p_rsum, p_m, p_e, mu, korCv, Ds);

                real_t H = irho*(e + p);

                // Finally, accumulate reference quantities into running sums...

                // ...including simple velocity-related quantities...
                acc[ref::ux     ](u.x());
                acc[ref::uy     ](u.y());
                acc[ref::uz     ](u.z());
                acc[ref::uxuy   ](u.x()*u.y());
                acc[ref::uzuy   ](u.z()*u.y());

                // ...and pressure/equation of state related quantities...
                acc[ref::p_ru](p_m.x());
                acc[ref::p_rw](p_m.z());
                acc[ref::p_rE](p_e    );

                acc[ref::vp_ru](u.y()*p_m.x());
                acc[ref::vp_rw](u.y()*p_m.z());
                acc[ref::vp_rE](u.y()*p_e    );

                acc[ref::Cmy_rho](u.y()*u.y() - p_rho - p_rsum);
                acc[ref::Ce_rho ](u.y()*(H - p_rho));
                acc[ref::Ce_rv  ](-H - u.y()*p_m.y());

                // ...and viscous-term-related quantities
                acc[ref::nu   ](mu*irho);
                acc[ref::korCv](korCv);
                acc[ref::Ds   ](Ds[0]); // Yes, b/c constant Lewis number!

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
            common.ref_uxuy      ()[j] = sum(acc[ref::uxuy      ]);
            common.ref_uzuy      ()[j] = sum(acc[ref::uzuy      ]);
            common.ref_p_ru      ()[j] = sum(acc[ref::p_ru      ]);
            common.ref_p_rw      ()[j] = sum(acc[ref::p_rw      ]);
            common.ref_p_rE      ()[j] = sum(acc[ref::p_rE      ]);
            common.ref_vp_ru     ()[j] = sum(acc[ref::vp_ru     ]);
            common.ref_vp_rw     ()[j] = sum(acc[ref::vp_rw     ]);
            common.ref_vp_rE     ()[j] = sum(acc[ref::vp_rE     ]);
            common.ref_Cmy_rho   ()[j] = sum(acc[ref::Cmy_rho   ]);
            common.ref_Ce_rho    ()[j] = sum(acc[ref::Ce_rho    ]);
            common.ref_Ce_rv     ()[j] = sum(acc[ref::Ce_rv     ]);
            common.ref_nu        ()[j] = sum(acc[ref::nu        ]);
            common.ref_korCv     ()[j] = sum(acc[ref::korCv     ]);
            common.ref_Ds        ()[j] = sum(acc[ref::Ds        ]);

        } // end Y

        // Zero y(j) not present on this rank to avoid accumulating garbage
        const size_t rightNotOnRank = common.refs.cols()
                                    - o.dgrid.local_physical_end.y();
        if (rightNotOnRank) common.refs.rightCols(rightNotOnRank).setZero();

        // Allreduce and scale common.refs sums to obtain means on all ranks
        // Allreduce mandatory as all ranks need references for linearization
        SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, common.refs.data(),
                common.refs.size(), mpi::datatype<real_t>::value,
                MPI_SUM, MPI_COMM_WORLD));
        common.refs *= o.dgrid.chi();

        // Copy mean velocity information into common.{u, v, w}()
        common.u() = common.ref_ux();
        common.v() = common.ref_uy();
        common.w() = common.ref_uz();

    } else {                                 // Mean velocity profile only

        SUZERAIN_TIMER_SCOPED("mean velocity profile");

        // NOTE: I'm leaving this block alone w/out really looking at it.
        // TODO: Figure out if/how it must change for reacting.

        // Zero all reference quantities on fully-explicit zeroth
        // substep when not using viscous filter
        if (   ZerothSubstep
            && Linearize == linearize::none
            && Filter    != filter::viscous ) {
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
            common.means.leftCols<3>() *= o.dgrid.chi();
        } else {
            ArrayXXr tmp;
            tmp.resizeLike(common.means.leftCols<3>());
            tmp.setZero();
            SUZERAIN_MPICHKR(MPI_Reduce(common.means.leftCols<3>().data(),
                        tmp.data(), common.means.leftCols<3>().size(),
                        mpi::datatype<real_t>::value, MPI_SUM,
                        o.dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
        }

    } // end traversal (1)

    if (Filter == filter::viscous) {
        // After first traversal, have gathered reference profiles.  So,
        // we can complete the viscous filter source calculation by
        // multiplying by the appropriate reference quantities.

        const std::size_t Ny = fsrcw.shape()[1];
        const std::size_t Nplane = fsrcw.shape()[2]*fsrcw.shape()[3];

        // Energy
        {
            Map<MatrixXXc> F(fsrcw[ndx::e].origin(), Ny, Nplane);
            const VectorXr& D(common.ref_korCv());
            MatrixXXc tmp(Ny, Nplane);
            tmp = D.asDiagonal()*F;
            F = tmp;
        }

        // x-momentum
        {
            Map<MatrixXXc> F(fsrcw[ndx::mx].origin(), Ny, Nplane);
            const VectorXr& D(common.ref_nu());
            MatrixXXc tmp(Ny, Nplane);
            tmp = D.asDiagonal()*F;
            F = tmp;
        }

        // y-momentum
        {
            Map<MatrixXXc> F(fsrcw[ndx::my].origin(), Ny, Nplane);
            const VectorXr& D(common.ref_nu());
            MatrixXXc tmp(Ny, Nplane);
            tmp = (cmods.alpha + 4.0/3.0) * D.asDiagonal()*F;
            F = tmp;
        }

        // z-momentum
        {
            Map<MatrixXXc> F(fsrcw[ndx::mz].origin(), Ny, Nplane);
            const VectorXr& D(common.ref_nu());
            MatrixXXc tmp(Ny, Nplane);
            tmp = D.asDiagonal()*F;
            F = tmp;
        }

        // mass
        {
            Map<MatrixXXc> F(fsrcw[ndx::rho].origin(), Ny, Nplane);
            F *= 0;
        }

        // species
        {
            for (unsigned int s=1; s<Ns; ++s) {
                Map<MatrixXXc> F(fsrcw[ndx::rho+s].origin(), Ny, Nplane);
                const VectorXr& D(common.ref_Ds());
                MatrixXXc tmp(Ny, Nplane);
                tmp = D.asDiagonal()*F;
                F = tmp;
            }
        }

        // Done with viscous filter source
    }

    // Traversal:
    // (2) Computing the nonlinear equation right hand sides.
    SUZERAIN_TIMER_BEGIN("nonlinear right hand sides (pass 1)");
    for (int offset = 0, j = o.dgrid.local_physical_start.y();
         j < o.dgrid.local_physical_end.y();
         ++j) {

        // Wall-normal operator eigenvalue estimates depend on location
        const real_t lambda1_y = o.lambda1_y(j);
        const real_t lambda2_y = o.lambda2_y(j);

        // Iterate across the j-th ZX plane
        const int last_zxoffset = offset
                                + o.dgrid.local_physical_extent.z()
                                * o.dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // Unpack reference quantities (used in dt calc)
            const Vector3r ref_u              (common.ref_ux        ()[j],
                                               common.ref_uy        ()[j],
                                               common.ref_uz        ()[j]);

            const real_t ref_nu   (common.ref_nu()   [j]);
            const real_t ref_korCv(common.ref_korCv()[j]);
            const real_t ref_Ds   ( (Ds.size()>1) ? common.ref_Ds()[j] : 0.0);


            // Unpack density-related quantities
            const real_t   rho         ( sphys(ndx::rho,    offset));

            const real_t  irho = 1.0/rho;

            const Vector3r grad_rho    (  auxp(aux::rho+dir::x,  offset),
                                          auxp(aux::rho+dir::y,  offset),
                                          auxp(aux::rho+dir::z,  offset));

            // Unpack momentum-related quantities
            const Vector3r m    ( sphys(ndx::mx, offset),
                                  sphys(ndx::my, offset),
                                  sphys(ndx::mz, offset));
            const real_t   div_m(  auxp(aux::mx+dir::x, offset)
                                 + auxp(aux::my+dir::y, offset)
                                 + auxp(aux::mz+dir::z, offset));
            const Matrix3r grad_m;
            const_cast<Matrix3r&>(grad_m) <<
                auxp(aux::mx+dir::x,  offset),
                auxp(aux::mx+dir::y,  offset),
                auxp(aux::mx+dir::z,  offset),
                auxp(aux::my+dir::x,  offset),
                auxp(aux::my+dir::y,  offset),
                auxp(aux::my+dir::z,  offset),
                auxp(aux::mz+dir::x,  offset),
                auxp(aux::mz+dir::y,  offset),
                auxp(aux::mz+dir::z,  offset);

            // Unpack total energy-related quantities
            const real_t e        (sphys(ndx::e,       offset));

            // Unpack species variables
            // NOTE: In species vector, idx 0 is the dilluter (the species
            // that is not explicitly part of the state vector)
            species(0) = rho;

            grad_species(0,0) = grad_rho(0);
            grad_species(1,0) = grad_rho(1);
            grad_species(2,0) = grad_rho(2);

            for (unsigned int s=1; s<Ns; ++s) {
                species(s) = sphys(ndx::species + s - 1, offset);

                grad_species(0,s) = auxp(aux::species + s - 1 + dir::x, offset);
                grad_species(1,s) = auxp(aux::species + s - 1 + dir::y, offset);
                grad_species(2,s) = auxp(aux::species + s - 1 + dir::z, offset);

                // dilluter density = rho_0 = rho - sum_{s=1}^{Ns-1} rho_s
                species(0)        -= species(s);

                grad_species(0,0) -= grad_species(0,s);
                grad_species(1,0) -= grad_species(1,s);
                grad_species(2,0) -= grad_species(2,s);
            }

            // Compute mass fractions and mass fraction gradients
            for (unsigned int s=0; s<Ns; ++s) {

                cs(s) = irho * species(s);

                grad_cs(0,s) = irho * (grad_species(0,s) - cs(s) * grad_rho(0));
                grad_cs(1,s) = irho * (grad_species(1,s) - cs(s) * grad_rho(1));
                grad_cs(2,s) = irho * (grad_species(2,s) - cs(s) * grad_rho(2));

            }


            // Compute velocity-related quantities
            const Vector3r u     = suzerain::rholut::u(rho, m);

            const real_t div_u   =
                suzerain::rholut::div_u(rho, grad_rho, m, div_m);

            const Matrix3r grad_u=
                suzerain::rholut::grad_u(rho, grad_rho, m, grad_m);


            // Compute temperature, pressure, mass diffusivities,
            // viscosity, thermal conductivity, species enthalpies, and
            // reaction source terms
            real_t T, p, mu, kap, a, Cp;
            cmods.evaluate(e, m, rho, species, cs,
                           T, p, Ds, mu, kap, hs, om, a, Cp);

            const real_t lam = (cmods.alpha - 2.0/3.0)*mu;

            // Compute quantities related to the viscous stress tensor
            const Matrix3r tau = suzerain::rholut::tau(mu, lam, div_u, grad_u);

            // Place to sum species fluxes from Fick's model
            Vector3r sdifftot (0.0, 0.0, 0.0);

            // Compute Fick's model contribution to diffusive fluxes and sum
            for (unsigned int s=0; s<Ns; ++s) {

                sdiff(0,s) = rho * Ds(s) * grad_cs(0,s);
                sdiff(1,s) = rho * Ds(s) * grad_cs(1,s);
                sdiff(2,s) = rho * Ds(s) * grad_cs(2,s);

                sdifftot(0) += sdiff(0,s);
                sdifftot(1) += sdiff(1,s);
                sdifftot(2) += sdiff(2,s);

            }

            // Subtract off cs*sdifftot to get SCEBD fluxes
            for (unsigned int s=0; s<Ns; ++s) {

                sdiff(0,s) -= cs(s) * sdifftot(0);
                sdiff(1,s) -= cs(s) * sdifftot(1);
                sdiff(2,s) -= cs(s) * sdifftot(2);

            }


            // Source terms get accumulated into state storage
            //
            // NOTE: Sign correct for source appearing on the RHS---i.e.,
            // U_t + div(F) = S.
            //
            // TODO: Add slow growth terms.
            sphys(ndx::e  , offset) = 0.0;
            sphys(ndx::mx , offset) = 0.0;
            sphys(ndx::my , offset) = 0.0;
            sphys(ndx::mz , offset) = 0.0;
            sphys(ndx::rho, offset) = 0.0;

            for (unsigned int s=1; s<Ns; ++s) {
                sphys(ndx::rho+s, offset) = om(s);
            }


            // Fluxes get accumulated into auxp
            //
            // NOTE: Sign correct for fluxes appearing on the LHS---i.e.,
            // U_t + div(F) = S(U).
            //
            // TODO: "Eigenify" these calcs where appropriate

            Vector3r vwork = tau*u;

            //----------------------------------------------------------------
            // ENERGY                    = u     * (rho*H) - viscous work ...
            auxp(aux::e +dir::x, offset) = u.x() * (e + p) - vwork.x() ;
            auxp(aux::e +dir::y, offset) = u.y() * (e + p) - vwork.y() ;
            auxp(aux::e +dir::z, offset) = u.z() * (e + p) - vwork.z() ;

            // ... - species enthalpy term
            if (Ns>1) {
                // NOTE: If Ns=1, we should have sdiff(*,0) = 0.0.  Thus,
                // this loop would be entered but shouldn't do anything.
                for (unsigned int s=0; s<Ns; ++s) {
                    auxp(aux::e +dir::x, offset) -= sdiff(0,s) * hs(s);
                    auxp(aux::e +dir::y, offset) -= sdiff(1,s) * hs(s);
                    auxp(aux::e +dir::z, offset) -= sdiff(2,s) * hs(s);
                }
            }
            // NOTE: the heat flux is accumulated below (in traversal 3)
            //----------------------------------------------------------------

            //----------------------------------------------------------------
            // MOMENTUM                  = convection     - viscous    +pressure
            auxp(aux::mx+dir::x, offset) = u.x() * m.x()  -  tau(0,0)  +  p ;
            auxp(aux::mx+dir::y, offset) = u.x() * m.y()  -  tau(1,0)       ;
            auxp(aux::mx+dir::z, offset) = u.x() * m.z()  -  tau(2,0)       ;

            auxp(aux::my+dir::x, offset) = u.y() * m.x()  -  tau(0,1)       ;
            auxp(aux::my+dir::y, offset) = u.y() * m.y()  -  tau(1,1)  +  p ;
            auxp(aux::my+dir::z, offset) = u.y() * m.z()  -  tau(2,1)       ;

            auxp(aux::mz+dir::x, offset) = u.z() * m.x()  -  tau(0,2)       ;
            auxp(aux::mz+dir::y, offset) = u.z() * m.y()  -  tau(1,2)       ;
            auxp(aux::mz+dir::z, offset) = u.z() * m.z()  -  tau(2,2)  +  p ;
            //----------------------------------------------------------------

            //----------------------------------------------------------------
            // mass                       = mass flux
            auxp(aux::rho+dir::x, offset) = m.x();
            auxp(aux::rho+dir::y, offset) = m.y();
            auxp(aux::rho+dir::z, offset) = m.z();
            //----------------------------------------------------------------

            //----------------------------------------------------------------
            // species
            for (unsigned int s=0; s<Ns-1; ++s) {
                // NOTE: species(0) is the diluter!

                //                                = convection    - diffusion
                auxp(aux::species+dir::x, offset) = cs(s+1)*m.x() - sdiff(0,s+1);
                auxp(aux::species+dir::y, offset) = cs(s+1)*m.y() - sdiff(1,s+1);
                auxp(aux::species+dir::z, offset) = cs(s+1)*m.z() - sdiff(2,s+1);

            }
            //----------------------------------------------------------------------

            // Finally, put temperature and thermal conductivity data
            // into auxp for later use
            auxp(aux::kap, offset) = kap;
            auxp(aux::T  , offset) = T;

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
                using math::minnan;

                // See timestepper::convective_stability_criterion
                //
                // NOTE: Speed of sound a computed by cmods.evaluate above
                real_t       ua_l1_x,       ua_l1_y,       ua_l1_z;
                real_t fluct_ua_l1_x, fluct_ua_l1_y, fluct_ua_l1_z;
                switch (Linearize) {
                default:
                    SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
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

                // See timestepper::diffusive_stability_criterion
                // Antidiffusive locations might be ignored when linearized.
                // Hence we compute criteria within the switch statment.
                const real_t nu  = mu  /  rho;
                const real_t kd  = kap / (rho*Cp); // thermal diffusivity
                const real_t Ds0 = ( Ds.size()>1 ? Ds[0] : 0.0 );

                real_t diffusivity;
                switch (Linearize) {
                default:
                    SUZERAIN_ERROR_VAL_UNIMPLEMENTED(std::vector<real_t>());
                    break;

                // Explicit treatment forces a zero reference diffusivity
                case linearize::none:
                    // NB: species diffusivities ok b/c constant Le
                    diffusivity = max(Ds0,
                                  max(kd,
                                  max(real_t(1), cmods.alpha)*nu));

                    diffusive_xyz_delta_t = minnan(diffusive_xyz_delta_t,
                                                   evmaxmag_real
                                                   / diffusivity
                                                   / (lambda2_x + lambda2_y + lambda2_z));
                    diffusive_x_delta_t   = min   (diffusive_x_delta_t,
                                                   evmaxmag_real / diffusivity / lambda2_x);
                    diffusive_y_delta_t   = min   (diffusive_y_delta_t,
                                                   evmaxmag_real / diffusivity / lambda2_y);
                    diffusive_z_delta_t   = min   (diffusive_z_delta_t,
                                                   evmaxmag_real / diffusivity / lambda2_z);
                    break;

                // Implicit diffusion permits removing a reference value.
                // Antidiffusive (nu - ref_nu) is fine and not computed.
                case linearize::rhome_y:
                    // Compute sign wrt ref.
                    diffusivity = max(Ds0-ref_Ds,
                                  max(kd-ref_korCv,
                                  max(real_t(1), cmods.alpha)*(nu-ref_nu)));

                    if (diffusivity <= 0) break;  // NaN => false, proceed
                    //diffusivity *= maxdiffconst;  // Rescale as necessary.
                    diffusive_xyz_delta_t = minnan(diffusive_xyz_delta_t,
                                                   evmaxmag_real
                                                   / diffusivity
                                                   / (lambda2_x + lambda2_y + lambda2_z));
                    diffusive_x_delta_t   = min   (diffusive_x_delta_t,
                                                   evmaxmag_real / diffusivity / lambda2_x);
                    diffusive_y_delta_t   = min   (diffusive_y_delta_t,
                                                   evmaxmag_real / diffusivity / lambda2_y);
                    diffusive_z_delta_t   = min   (diffusive_z_delta_t,
                                                   evmaxmag_real / diffusivity / lambda2_z);
                    break;


                } // end switch(Linearize)
            } // end if(ZerothSubstep}

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


    // (1a) Temperature from physical to wave space (still collocation in y)
    o.dgrid.transform_physical_to_wave(&auxp.coeffRef(aux::T,0));

    // (1b) Apply inverse mass matrix to get to pure coefficient space
    o.bop_solve(massluz, auxw, aux::T);

    // ... and zero wavenumbers present only for dealiasing
    // while simultaneously normalizing FFT
    o.diffwave_apply(0, 0, o.dgrid.chi(), auxw, aux::T);

    // (2) Get derivatives.
    //
    // Compute Y derivatives of variable var at collocation points
    o.bop_accumulate(1,    1, auxw, aux::T, 0, auxw, aux::gT + dir::y);

    // Get to collocation points in preparation for computing x and z
    // derivatives
    o.bop_apply     (0,    1, auxw, aux::T);

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
            const real_t   kap       ( auxp(aux::kap,   offset));

            // Extract grad(T)
            const Vector3r grad_T    ( auxp(aux::gT+dir::x, offset),
                                       auxp(aux::gT+dir::y, offset),
                                       auxp(aux::gT+dir::z, offset));

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

                    // Would prefer
                    //   ms.Q_conservative(x, y, z, time, Q_rho,
                    //                     Q_rho_u, Q_rho_v, Q_rho_w, Q_rho_E);
                    // but see MASA ticket #2794

                    Q_rho   = ms.Q_rho (x, y, z, time);
                    Q_rho_u = ms.Q_rhou(x, y, z, time);
                    Q_rho_v = ms.Q_rhov(x, y, z, time);
                    Q_rho_w = ms.Q_rhow(x, y, z, time);
                    Q_rho_E = ms.Q_rhoe(x, y, z, time);

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

        for (size_t i = 0; i < state_count; ++i) {

            // Apply inverse mass matrix to Y flux to get to pure
            // coefficient representation
            o.bop_solve(massluz, auxw, aux::e + dir::count*i + dir::y);

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

            // Accumulate filter source while simultaneously scaling for FFT
            // alpha = dNx*dNz = 1 / chi scales for dealiasing
            o.diffwave_accumulate(0, 0, 1 / o.dgrid.chi(), fsrcw, i,
                                        1                , swave, i);

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
