//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

/** @file
 * @copydoc perfect.hpp
 */

#include "perfect.hpp"

#include <esio/esio.h>
#include <esio/error.h>
#include <gsl/gsl_errno.h>
#include <sys/file.h>

#include <suzerain/blas_et_al.hpp>
#include <suzerain/bl.h>
#include <suzerain/bspline.hpp>
#include <suzerain/channel.h>
#include <suzerain/coalescing_pool.hpp>
#include <suzerain/common.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/htstretch.h>
#include <suzerain/inorder.hpp>
#include <suzerain/largo_state.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/profile.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/rngstream.hpp>
#include <suzerain/samples.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/specification_noise.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/field.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/utility.hpp>

#include "definition_scenario.hpp"
#include "instantaneous.hpp"
#include "references.hpp"
#include "slowgrowth.hpp"

using boost::numeric_cast;
using std::size_t;

namespace suzerain {

namespace perfect {

std::vector<support::field>
default_fields()
{
    std::vector<support::field> retval(ndx::identifier.static_size);
    for (size_t i = 0; i < ndx::identifier.static_size; ++i) {
        retval[i].identifier   = ndx::identifier[i];
        retval[i].description += "Nondimensional ";
        retval[i].description += ndx::description[i];
        retval[i].location     = ndx::identifier[i];
    }
    return retval;
}

void
adjust_scenario(contiguous_state<4,complex_t> &swave,
                const definition_scenario& scenario,
                const specification_grid& grid,
                const pencil_grid& dgrid,
                const bsplineop& cop,
                const real_t old_Ma,
                const real_t old_gamma)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    bool quickreturn = true;
#pragma warning(push,disable:1572)
    if (old_Ma != scenario.Ma) {
#pragma warning(pop)
        INFO0("Changing state Mach number from "
              << old_Ma << " to " << scenario.Ma);
        quickreturn = false;
    }
#pragma warning(push,disable:1572)
    if (old_gamma != scenario.gamma) {
#pragma warning(pop)
        INFO0("Changing state ratio of specific heats from "
              << old_gamma << " to " << scenario.gamma);
        quickreturn = false;
    }
    if (quickreturn) {
        DEBUG0("Not rescaling energy as scenario parameters have not changed");
        return;
    }

    // Convert state to physical space collocation points
    operator_tools otool(grid, dgrid, cop);
    physical_view<state_count> sphys(dgrid, swave);
    for (size_t k = 0; k < state_count; ++k) {
        otool.bop_apply(0, 1.0, swave, k);
        dgrid.transform_wave_to_physical(&sphys.coeffRef(k,0));
    }

    // Adjust total energy by the necessary amount at every collocation point
    // This procedure is not the cheapest or least numerically noisy,
    // but it does re-use existing compute kernels in a readable way.
    INFO0("Holding density and temperature constant during changes");
    for (int offset = 0, j = dgrid.local_physical_start.y();
        j < dgrid.local_physical_end.y();
        ++j) {
        const int last_zxoffset = offset
                                + dgrid.local_physical_extent.z()
                                * dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            const real_t   e  (sphys(ndx::e,   offset));
            const Vector3r m  (sphys(ndx::mx,  offset),
                               sphys(ndx::my,  offset),
                               sphys(ndx::mz,  offset));
            const real_t   rho(sphys(ndx::rho, offset));

            // Compute temperature using old_gamma, old_Ma
            real_t p, T;
            rholut::p_T(scenario.alpha, scenario.beta, old_gamma, old_Ma,
                        rho, m, e, /*out*/ p, /*out*/ T);
            // Compute total energy from new gamma, Ma, rho, T
            rholut::p(scenario.gamma, rho, T, /*out*/ p);
            sphys(ndx::e, offset) = rholut::energy_internal(scenario.gamma, p)
                                  + rholut::energy_kinetic(scenario.Ma, rho, m);
        }
    }

    // Convert state back to wave space coefficients in X, Y, and Z
    // building FFT normalization constant into the mass matrix
    bsplineop_luz massluz(cop);
    const complex_t scale_factor = 1 / dgrid.chi();
    massluz.opform(1, &scale_factor, cop);
    massluz.factor();
    for (size_t i = 0; i < state_count; ++i) {
        dgrid.transform_physical_to_wave(&sphys.coeffRef(i, 0));
        otool.bop_solve(massluz, swave, i);
    }
}

void
add_noise(contiguous_state<4,complex_t> &state,
          const specification_noise& noise,
          const definition_scenario& scenario,
          const specification_grid& grid,
          const pencil_grid& dgrid,
          const bsplineop& cop,
          bspline &b)
{
    // FIXME Needs to be made aware of channel vs flat plate

    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  state.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // Ensure we were handed collocation-based operator matrices
    SUZERAIN_ENSURE(cop.get()->method == SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);

    const int Ny = grid.N.y();

#pragma warning(push,disable:1572)
    if (noise.percent == 0) {
#pragma warning(pop)
        DEBUG0("Zero noise added to velocity fields");
        return;
    }

    using inorder::wavenumber_abs;
    using inorder::wavenumber_max;
    using inorder::wavenumber_translatable;
    const real_t twopi = 2 * boost::math::constants::pi<real_t>();

    // Evaluate maximum fluctuation magnitude based on percentage of
    // centerline mean streamwise velocity and broadcast result.
    // Alright, alright... actually an approximate mean velocity.
    real_t maxfluct;
    if (dgrid.has_zero_zero_modes()) {
        complex_t momentum, density;
        const real_t centerline = grid.L.y() / 2;
        b.linear_combination(
                0, &state[ndx::mx][0][0][0], centerline, &momentum);
        INFO0("Centerline mean streamwise momentum at y = "
              << centerline << " is " << momentum);
        b.linear_combination(
                0, &state[ndx::rho][0][0][0], centerline, &density);
        INFO0("Centerline mean density at y = "
              << centerline << " is " << density);
        maxfluct = noise.percent / 100 * (abs(momentum) / abs(density));
    }
    SUZERAIN_MPICHKR(MPI_Bcast(&maxfluct, 1,
                mpi::datatype_of(maxfluct),
                dgrid.rank_zero_zero_modes, MPI_COMM_WORLD));
    INFO0("Adding velocity perturbations with maximum magnitude " << maxfluct);

    // Compute and display kxfrac_min, kxfrac_max constraints
#pragma warning(push,disable:2259)
    const int dkx_max = noise.kxfrac_max * wavenumber_max(grid.N.x());
    const int dkx_min = noise.kxfrac_min * wavenumber_max(grid.N.x());
#pragma warning(pop)
#pragma warning(push,disable:1572)
    if (noise.kxfrac_max != 1 || noise.kxfrac_min != 0) {
#pragma warning(pop)
        INFO0("Perturbations added only to absolute X wavenumbers in range ["
                << dkx_min << ":" << dkx_max << "]");
    }

    // Compute and display kzfrac_min, kzfrac_max constraints
#pragma warning(push,disable:2259)
    const int dkz_max = noise.kzfrac_max * wavenumber_max(grid.N.z());
    const int dkz_min = noise.kzfrac_min * wavenumber_max(grid.N.z());
#pragma warning(pop)
#pragma warning(push,disable:1572)
    if (noise.kzfrac_max != 1 || noise.kzfrac_min != 0) {
#pragma warning(pop)
    INFO0("Perturbations added only to absolute Z wavenumbers in range ["
            << dkz_min << ":" << dkz_max << "]");
    }

    // Form mass matrix to convert (wave, collocation point values, wave)
    // perturbations to (wave, coefficients, wave)
    bsplineop_luz massluz(cop);
    massluz.factor_mass(cop);

    // Set L'Ecuyer et al.'s rngstream seed.  Use a distinct Substream for each
    // wall-normal pencil to ensure process is a) repeatable despite changes in
    // processor count, b) easy to code, and c) embarrassingly parallel.
    rngstream rng;
    {
        array<unsigned long,6> seed;
        std::fill(seed.begin(), seed.end(), noise.seed);
        rng.SetSeed(seed.data());
    }
    rng.IncreasedPrecis(true);  // Use more bits of resolution

    // Approach is the following:
    //  0) Allocate storage for state and three additional scalar fields.
    //  1) Generate a random vector-valued field \tilde{A}.
    //     \tilde{A}'s x- and z-derivatives have zero mean by periodicity;
    //  2) Zero first two B-spline coefficients near walls so partial_y
    //     \tilde{A} vanishes at the wall
    //  3) This step intentionally left blank.
    //  4) Compute curl A in physical space and rescale so maximum
    //     pointwise norm of A is maxfluct.  curl A is solenoidal
    //     and now has the velocity perturbation properties we desire.
    //  5) Store curl A in physical space in the three scalar fields.
    //  6) Copy state into auxiliary state storage and bring to
    //     physical space.
    //  7) At each point compute velocity and internal energy.  Perturb
    //     velocity and compute new total energy using perturbed
    //     velocities.
    //  8) Bring perturbed state information back to wavespace.
    //  9) Overwrite state storage with the new perturbed state.

    //  0) Allocate storage for state and three additional scalar fields.
    scoped_ptr<contiguous_state<4,complex_t> > _s_ptr( // RAII
            support::allocate_padded_state<contiguous_state<4,complex_t> >(
                state_count + 3, dgrid));
    contiguous_state<4,complex_t> &s = *_s_ptr;               // Shorthand
    std::fill(s.range().begin(), s.range().end(), 0);        // Zero memory

    // 1) Generate a random vector-valued field \tilde{A}.
    // For each scalar component of \tilde{A}...
    for (size_t l = 0; l < 3; ++l) {

        for (int k = 0; k < grid.dN.z(); ++k) {
            if (!wavenumber_translatable(grid.N.z(), grid.dN.z(), k)) continue;


            for (int i = 0; i < grid.dN.x(); ++i) {
                if (!wavenumber_translatable(grid.N.x(), grid.dN.x(), i)) continue;

                // ...and advance rngstream to the (i, ., k) substream...
                // ...(necessary for processor-topology independence)...
                rng.ResetNextSubstream();

                // ...but only the rank holding the (i, ., k) pencil continues.
                if (   k <  dgrid.local_wave_start.z()
                    || k >= dgrid.local_wave_end.z()
                    || i <  dgrid.local_wave_start.x()
                    || i >= dgrid.local_wave_end.x()) continue;

                // Satisfy fluct_kzfrac_min, fluct_kzfrac_max constraints
                if (wavenumber_abs(grid.dN.z(), k) < dkz_min) continue;
                if (wavenumber_abs(grid.dN.z(), k) > dkz_max) continue;

                // Satisfy fluct_kxfrac_min, fluct_kxfrac_max constraints
                if (wavenumber_abs(grid.dN.x(), i) < dkx_min) continue;
                if (wavenumber_abs(grid.dN.x(), i) > dkx_max) continue;

                // Compute local indices for global (i, ., k).
                const int local_i = i - dgrid.local_wave_start.x();
                const int local_k = k - dgrid.local_wave_start.z();

                // ...generate coeffs with well-defined pseudorandom order.
                //  2) Zero first two B-spline coefficients near walls
                //     so partial_y \tilde{A} vanishes at the wall
                s[2*l][0][local_i][local_k] = 0;
                s[2*l][1][local_i][local_k] = 0;
                for (int j = 2; j < Ny - 2; ++j) {
                    const real_t magnitude = rng.RandU01();
                    const real_t phase     = rng.RandU01() * twopi;
                    s[2*l][j][local_i][local_k] = std::polar(magnitude, phase);
                }
                s[2*l][Ny - 2][local_i][local_k] = 0;
                s[2*l][Ny - 1][local_i][local_k] = 0;

            } // end X

        } // end Z

    } // end scalar components of A

    //  4) Compute curl A in physical space and rescale so maximum
    //     pointwise norm of A is maxfluct.  curl A is solenoidal
    //     and now has the velocity perturbation properties we desire.

    // Copy s[2l] into s[2l+1] as we need two copies to compute curl A
    for (size_t l = 0; l < 3; ++l) s[2*l+1] = s[2*l];

    // Prepare physical-space view of the wave-space storage
    physical_view<state_count+3> p(dgrid, s);

    // Initializing operator_base to access decomposition-ready utilities
    operator_tools otool(grid, dgrid, cop);

    // From Ax in s[0] compute \partial_y Ax
    otool.bop_apply(1, 1.0, s, 0);
    dgrid.transform_wave_to_physical(&p.coeffRef(0,0));

    // From Ax in s[1]  compute \partial_z Ax
    otool.bop_apply(0, 1.0, s, 1);
    otool.diffwave_apply(0, 1, 1.0, s, 1);
    dgrid.transform_wave_to_physical(&p.coeffRef(1,0));

    // From Ay in s[2] compute \partial_x Ay
    otool.bop_apply(0, 1.0, s, 2);
    otool.diffwave_apply(1, 0, 1.0, s, 2);
    dgrid.transform_wave_to_physical(&p.coeffRef(2,0));

    // From Ay in s[3] compute \partial_z Ay
    otool.bop_apply(0, 1.0, s, 3);
    otool.diffwave_apply(0, 1, 1.0, s, 3);
    dgrid.transform_wave_to_physical(&p.coeffRef(3,0));

    // From Az in s[4] compute \partial_x Az
    otool.bop_apply(0, 1.0, s, 4);
    otool.diffwave_apply(1, 0, 1.0, s, 4);
    dgrid.transform_wave_to_physical(&p.coeffRef(4,0));

    // From Az in s[5] compute \partial_y Az
    otool.bop_apply(1, 1.0, s, 5);
    dgrid.transform_wave_to_physical(&p.coeffRef(5,0));

    // Store curl A in s[{5,6,7}] and find global maximum magnitude of curl A
    real_t maxmagsquared = 0;
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Assert curl A is identically zero at the walls
                if (j == 0 || j == Ny - 1) {
#pragma warning(push,disable:1572)
                    assert(p(0, offset) == 0.0);
                    assert(p(1, offset) == 0.0);
                    assert(p(2, offset) == 0.0);
                    assert(p(3, offset) == 0.0);
                    assert(p(4, offset) == 0.0);
                    assert(p(5, offset) == 0.0);
#pragma warning(pop)
                }

                // The definition of curl gives the following components
                //   1) \partial_y A_z - \partial_z A_y
                //   2) \partial_z A_x - \partial_x A_z
                //   3) \partial_x A_y - \partial_y A_x
                // where the mean of \partial_x and \partial_z terms must be
                // zero by periodicity.  Components 1 and 3 may have nonzero
                // mean because they wall-normal derivatives contributions.
                const Vector3r curlA(p(5, offset) - p(3, offset),
                                     p(1, offset) - p(4, offset),
                                     p(2, offset) - p(0, offset));

                //  5) Store curl A in physical space in the 3 scalar fields.
                //
                // A nonzero mean in the x, y, and z directions is,
                // respectively, "corrected" by bulk forcing, the wall, and
                // viscous effects.  The spanwise viscous effects presumably
                // have the slowest timescale so rotate the components from 123
                // to 312 to reduce the simulation time before stationarity.
                // This rotation may introduce acoustic noise.
                p(state_count + 0, offset) = curlA.z();
                p(state_count + 1, offset) = curlA.x();
                p(state_count + 2, offset) = curlA.y();

                maxmagsquared = math::maxnan(
                        maxmagsquared, curlA.squaredNorm());

            } // end X

        } // end Z

    } // end Y
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, &maxmagsquared, 1,
                mpi::datatype<real_t>::value,
                MPI_MAX, MPI_COMM_WORLD));

    // Rescale curl A components so max ||curl A|| == maxfluct
    p.row(state_count + 0) *= (maxfluct / std::sqrt(maxmagsquared));
    p.row(state_count + 1) *= (maxfluct / std::sqrt(maxmagsquared));
    p.row(state_count + 2) *= (maxfluct / std::sqrt(maxmagsquared));

    //  6) Copy state into auxiliary state storage and bring to
    //     physical space.
    for (size_t i = 0; i < state_count; ++i) {
        s[i] = state[i];
        otool.bop_apply(0, 1.0, s, i);
        dgrid.transform_wave_to_physical(&p.coeffRef(i,0));
    }

    //  7) At each point compute velocity and internal energy.  Perturb
    //     velocity and compute new total energy using perturbed
    //     velocities.
    const real_t Ma = scenario.Ma;
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        for (int k = dgrid.local_physical_start.z();
            k < dgrid.local_physical_end.z();
            ++k) {

            for (int i = dgrid.local_physical_start.x();
                i < dgrid.local_physical_end.x();
                ++i, /* NB */ ++offset) {

                // Retrieve internal energy
                real_t         e(p(ndx::e,   offset));
                Vector3r       m(p(ndx::mx,  offset),
                                 p(ndx::my,  offset),
                                 p(ndx::mz,  offset));
                const real_t rho(p(ndx::rho, offset));
                const real_t e_int = rholut::energy_internal(Ma, rho, m, e);

                // Perturb momentum and compute updated total energy
                m.x() += rho * p(state_count + 0, offset);
                m.y() += rho * p(state_count + 1, offset);
                m.z() += rho * p(state_count + 2, offset);
                const real_t e_kin = rholut::energy_kinetic(Ma, rho, m);
                e = e_int + e_kin;

                // Store results back to state fields
                p(ndx::e,  offset) = e;
                p(ndx::mx, offset) = m.x();
                p(ndx::my, offset) = m.y();
                p(ndx::mz, offset) = m.z();

            } // end X

        } // end Z

    } // end Y

    //  8) Bring perturbed state information back to wavespace (no rho!)
    // Build FFT normalization constant into Y direction's mass matrix.
    const complex_t scale_factor = 1 / dgrid.chi();
    massluz.opform(1, &scale_factor, cop);
    massluz.factor();
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::e , 0));  // X, Z
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::mx, 0));  // X, Z
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::my, 0));  // X, Z
    dgrid.transform_physical_to_wave(&p.coeffRef(ndx::mz, 0));  // X, Z
    otool.bop_solve(massluz, s, ndx::e );                       // Y
    otool.bop_solve(massluz, s, ndx::mx);                       // Y
    otool.bop_solve(massluz, s, ndx::my);                       // Y
    otool.bop_solve(massluz, s, ndx::mz);                       // Y

    //  9) Overwrite state storage with the new perturbed state (not rho!)
    state[ndx::e ] = s[ndx::e ];
    state[ndx::mx] = s[ndx::mx];
    state[ndx::my] = s[ndx::my];
    state[ndx::mz] = s[ndx::mz];
}

// This looks like logic from navier_stokes.hpp but does not belong there.
// Reading through that file, especially the apply_operator implementation, is
// recommended before reviewing this logic.  This routine is definitely
// suboptimal but is expected to be invoked relatively infrequently.
std::auto_ptr<samples>
take_samples(const definition_scenario &scenario,
             const specification_largo &sg,
             const operator_tools &otool,
             const bspline &b,
             contiguous_state<4,complex_t> &swave,
             const real_t t)
{
    // State enters method as coefficients in X, Y, and Z directions
    SUZERAIN_TIMER_SCOPED("perfect::take_samples");

    // Shorthand for the operator_tools member(s) commonly used
    const pencil_grid &dgrid = otool.dgrid;

    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Shorthand
    const std::size_t Ny = swave.shape()[1];
    namespace acc = boost::accumulators;
    typedef contiguous_state<4,complex_t> state_type;

    // We need auxiliary scalar-field storage.  Prepare logical indices using a
    // struct for scoping (e.g. aux::rho_y).  Ordering will match usage below.
    struct aux { enum {
        e_y,   e_x,   e_z,
        mx_y,  mx_x,  mx_z,
        my_y,  my_x,  my_z,
        mz_y,  mz_x,  mz_z,
        rho_y, rho_x, rho_z,
        count // Sentry
    }; };

    // Obtain the auxiliary storage (likely from a pool to avoid fragmenting).
    // We assume no garbage values in the memory will impact us (for speed).
    scoped_ptr<state_type> _auxw_ptr(
            support::allocate_padded_state<state_type>(
                aux::count, dgrid)); // RAII
    state_type &auxw = *_auxw_ptr;                                 // Shorthand

    // Ensure state storage meets this routine's assumptions

    // Sanity check incoming swave's and auxw's shape and contiguity
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE((unsigned) swave.strides()[1] == 1u);
    SUZERAIN_ENSURE((unsigned) swave.strides()[2] == swave.shape()[1]);
    SUZERAIN_ENSURE((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);
    SUZERAIN_ENSURE(std::equal(swave.shape() + 1, swave.shape() + 4,
                               auxw.shape() + 1));
    SUZERAIN_ENSURE(std::equal(swave.strides() + 1, swave.strides() + 4,
                               auxw.strides() + 1));

    // Rank-specific details accumulated in ret to be MPI_Allreduce later
    std::auto_ptr<samples> ret(new samples(t, Ny));

    // Obtain samples available in wave-space from mean conserved state.
    // These coefficients are inherently averaged across the X-Z plane.
    if (dgrid.has_zero_zero_modes()) {
        ret->rho_E()        = Map<VectorXc>(swave[ndx::e  ].origin(), Ny).real();
        ret->rho_u().col(0) = Map<VectorXc>(swave[ndx::mx ].origin(), Ny).real();
        ret->rho_u().col(1) = Map<VectorXc>(swave[ndx::my ].origin(), Ny).real();
        ret->rho_u().col(2) = Map<VectorXc>(swave[ndx::mz ].origin(), Ny).real();
        ret->rho()          = Map<VectorXc>(swave[ndx::rho].origin(), Ny).real();
    }

    // Compute Y derivatives of total energy at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::e);
    otool.bop_accumulate(1,    1., swave, ndx::e, 0., auxw, aux::e_y);
    otool.bop_apply     (0,    1., swave, ndx::e);

    // Compute X- and Z- derivatives of total energy at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::e, 0., auxw, aux::e_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::e, 0., auxw, aux::e_z);

    // Compute Y derivatives of X momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::mx);
    otool.bop_accumulate(1,    1., swave, ndx::mx, 0., auxw, aux::mx_y);
    otool.bop_apply     (0,    1., swave, ndx::mx);

    // Compute X- and Z- derivatives of X momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::mx, 0., auxw, aux::mx_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::mx, 0., auxw, aux::mx_z);

    // Compute Y derivatives of Y momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::my);
    otool.bop_accumulate(1,    1., swave, ndx::my, 0., auxw, aux::my_y);
    otool.bop_apply     (0,    1., swave, ndx::my);

    // Compute X- and Z- derivatives of Y momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::my, 0., auxw, aux::my_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::my, 0., auxw, aux::my_z);

    // Compute Y derivatives of Z momentum at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::mz);
    otool.bop_accumulate(1,    1., swave, ndx::mz, 0., auxw, aux::mz_y);
    otool.bop_apply     (0,    1., swave, ndx::mz);

    // Compute X- and Z- derivatives of Z momentum at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::mz, 0., auxw, aux::mz_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::mz, 0., auxw, aux::mz_z);

    // Compute Y derivatives of density at collocation points
    // Zero wavenumbers present only for dealiasing along the way
    otool.zero_dealiasing_modes(   swave, ndx::rho);
    otool.bop_accumulate(1,    1., swave, ndx::rho, 0., auxw, aux::rho_y);
    otool.bop_apply     (0,    1., swave, ndx::rho);

    // Compute X- and Z- derivatives of density at collocation points
    // Zeros wavenumbers present only for dealiasing in the target storage
    otool.diffwave_accumulate(1, 0, 1., swave, ndx::rho,  0., auxw, aux::rho_x);
    otool.diffwave_accumulate(0, 1, 1., swave, ndx::rho,  0., auxw, aux::rho_z);

    // Initialize slow growth treatment if necessary
    slowgrowth sg_treater(sg, scenario.Ma);
    sg_treater.initialize();

    // With conserved state Fourier in X and Z but collocation in Y,
    // gather any RMS-like information necessary for slow growth forcing.
    sg_treater.gather_wavexz(otool, swave);

    // Collectively convert swave and auxw to physical space using parallel
    // FFTs. In physical space, we'll employ views to reshape the 4D row-major
    // (F, Y, Z, X) with contiguous (Y, Z, X) into a 2D (F, Y*Z*X) layout where
    // we know F a priori.  Reducing the dimensionality encourages linear
    // access and eases indexing overhead.
    physical_view<aux::count>  auxp(dgrid, auxw);
    physical_view<state_count> sphys(dgrid, swave);
    for (std::size_t i = 0; i < state_count; ++i) {
        dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }
    for (std::size_t i = 0; i < aux::count; ++i) {
        dgrid.transform_wave_to_physical(&auxp.coeffRef(i,0));
    }

    // Physical space is traversed linearly using a single offset 'offset'.
    // The two loop structure is used as positions x(i), and z(k) unnecessary.
    //
    // Two traversals occur:
    // (1) Computing mean profiles necessary for slow growth forcing.
    // (2) Gathering the full suite of desired samples, including slow growth.

    // Traversal:
    // (1) Computing mean profiles necessary for slow growth forcing.
    {
        instantaneous inst;
        collect_instantaneous(scenario, otool.grid, otool.dgrid, sphys, inst);
        sg_treater.gather_physical_cons(otool, inst);
        sg_treater.gather_physical_rqq (otool, inst);
    }

    // Traversal:
    // (2) Gathering the full suite of desired samples, including slow growth.
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Prepare any y-dependent slow growth computation
        sg_treater.inner_y(j, b.collocation_point(j));

        // Prepare logical indices using struct for scoping (e.g. ref::ux).
#define COMPONENT(quantity, component, offset, description) component,
        struct ref { enum {
            SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH(COMPONENT,
                                                 SUZERAIN_SAMPLES_PHYSICAL)
            count // Sentry
        }; };
#undef  COMPONENT

        // An array of summing_accumulator_type holds all running sums
        array<summing_accumulator_type, ref::count> acc;

        // Iterate across the j-th ZX plane
        const int last_zxoffset = offset
                                + dgrid.local_physical_extent.z()
                                * dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // Unpack total energy-related quantities
            const real_t e(sphys(ndx::e, offset));
            const Vector3r grad_e(auxp(aux::e_x, offset),
                                  auxp(aux::e_y, offset),
                                  auxp(aux::e_z, offset));

            // Unpack momentum-related quantities
            const Vector3r m(sphys(ndx::mx, offset),
                             sphys(ndx::my, offset),
                             sphys(ndx::mz, offset));
            const real_t div_m = auxp(aux::mx_x, offset)
                               + auxp(aux::my_y, offset)
                               + auxp(aux::mz_z, offset);
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

            // Unpack density-related quantities
            const real_t rho(sphys(ndx::rho, offset));
            const Vector3r grad_rho(auxp(aux::rho_x, offset),
                                    auxp(aux::rho_y, offset),
                                    auxp(aux::rho_z, offset));

            // Compute local quantities based upon state
            const real_t   rho2   = rho * rho;

            const Vector3r u      = rholut::u(
                                       rho, m);
            const real_t   div_u  = rholut::div_u(
                                        rho, grad_rho, m, div_m);
            const Matrix3r grad_u = rholut::grad_u(
                                        rho, grad_rho, m, grad_m);
            const Vector3r om       (grad_u(2,1) - grad_u(1,2),
                                     grad_u(0,2) - grad_u(2,0),
                                     grad_u(1,0) - grad_u(0,1));

            real_t p, T, mu, lambda;
            Vector3r grad_p, grad_T, grad_mu, grad_lambda;
            rholut::p_T_mu_lambda(
                scenario.alpha, scenario.beta, scenario.gamma, scenario.Ma,
                rho, grad_rho, m, grad_m, e, grad_e,
                p, grad_p, T, grad_T, mu, grad_mu, lambda, grad_lambda);
            const real_t a = std::sqrt(T);

            const Matrix3r tau = rholut::tau(
                                    mu, lambda, div_u, grad_u);
            const Vector3r tau_u = tau * u;

            // Accumulate quantities into sum_XXX using function syntax.
            acc[ref::E](e / rho);

            acc[ref::p](p);

            acc[ref::T](T);

            acc[ref::a](a);

            acc[ref::h0](e + p);

            acc[ref::H0]((e + p) / rho);

            acc[ref::ke](u.squaredNorm() / 2);

            acc[ref::mu](mu);

            acc[ref::nu](mu / rho);

            acc[ref::u](u.x());
            acc[ref::v](u.y());
            acc[ref::w](u.z());

            // Note code unit correction
            acc[ref::Ma](scenario.Ma * (u.x() / a));

            acc[ref::rho2](rho2);

            acc[ref::p2](p * p);

            acc[ref::T2](T * T);

            acc[ref::h02](std::pow(e + p, 2));

            acc[ref::H02](std::pow((e + p) / rho, 2));

            acc[ref::mu2](mu * mu);

            acc[ref::nu2](std::pow(mu / rho, 2));

            // Note code unit correction
            acc[ref::Ma2](std::pow(scenario.Ma * (u.x() / a), 2));

            acc[ref::symxx_grad_u]( grad_u(0,0)                   );
            acc[ref::symxy_grad_u]((grad_u(0,1) + grad_u(1,0)) / 2);
            acc[ref::symxz_grad_u]((grad_u(0,2) + grad_u(2,0)) / 2);
            acc[ref::symyy_grad_u]( grad_u(1,1)                   );
            acc[ref::symyz_grad_u]((grad_u(1,2) + grad_u(2,1)) / 2);
            acc[ref::symzz_grad_u]( grad_u(2,2)                   );

            acc[ref::symxx_rho_grad_u](rho *  grad_u(0,0)                   );
            acc[ref::symxy_rho_grad_u](rho * (grad_u(0,1) + grad_u(1,0)) / 2);
            acc[ref::symxz_rho_grad_u](rho * (grad_u(0,2) + grad_u(2,0)) / 2);
            acc[ref::symyy_rho_grad_u](rho *  grad_u(1,1)                   );
            acc[ref::symyz_rho_grad_u](rho * (grad_u(1,2) + grad_u(2,1)) / 2);
            acc[ref::symzz_rho_grad_u](rho *  grad_u(2,2)                   );

            acc[ref::symxx_rho2_grad_u](rho2 *  grad_u(0,0)                   );
            acc[ref::symxy_rho2_grad_u](rho2 * (grad_u(0,1) + grad_u(1,0)) / 2);
            acc[ref::symxz_rho2_grad_u](rho2 * (grad_u(0,2) + grad_u(2,0)) / 2);
            acc[ref::symyy_rho2_grad_u](rho2 *  grad_u(1,1)                   );
            acc[ref::symyz_rho2_grad_u](rho2 * (grad_u(1,2) + grad_u(2,1)) / 2);
            acc[ref::symzz_rho2_grad_u](rho2 *  grad_u(2,2)                   );

            acc[ref::gradx_T](grad_T.x());
            acc[ref::grady_T](grad_T.y());
            acc[ref::gradz_T](grad_T.z());

            acc[ref::rho_gradx_T](rho * grad_T.x());
            acc[ref::rho_grady_T](rho * grad_T.y());
            acc[ref::rho_gradz_T](rho * grad_T.z());

            acc[ref::rho2_gradx_T](rho2 * grad_T.x());
            acc[ref::rho2_grady_T](rho2 * grad_T.y());
            acc[ref::rho2_gradz_T](rho2 * grad_T.z());

            acc[ref::tau_colon_grad_u]((tau.transpose()*grad_u).trace());

            acc[ref::tauxx](tau(0,0));
            acc[ref::tauxy](tau(0,1));
            acc[ref::tauxz](tau(0,2));
            acc[ref::tauyy](tau(1,1));
            acc[ref::tauyz](tau(1,2));
            acc[ref::tauzz](tau(2,2));

            acc[ref::tauux](tau_u.x());
            acc[ref::tauuy](tau_u.y());
            acc[ref::tauuz](tau_u.z());

            acc[ref::p_div_u](p*div_u);

            acc[ref::rho_u_u](rho * u.x() * u.x());
            acc[ref::rho_u_v](rho * u.x() * u.y());
            acc[ref::rho_u_w](rho * u.x() * u.z());
            acc[ref::rho_v_v](rho * u.y() * u.y());
            acc[ref::rho_v_w](rho * u.y() * u.z());
            acc[ref::rho_w_w](rho * u.z() * u.z());

            acc[ref::rho2_u_u](rho2 * u.x() * u.x());
            acc[ref::rho2_u_v](rho2 * u.x() * u.y());
            acc[ref::rho2_u_w](rho2 * u.x() * u.z());
            acc[ref::rho2_v_v](rho2 * u.y() * u.y());
            acc[ref::rho2_v_w](rho2 * u.y() * u.z());
            acc[ref::rho2_w_w](rho2 * u.z() * u.z());

            acc[ref::rho_u_u_u](rho * u.x() * u.x() * u.x());
            acc[ref::rho_u_u_v](rho * u.x() * u.x() * u.y());
            acc[ref::rho_u_u_w](rho * u.x() * u.x() * u.z());
            acc[ref::rho_u_v_v](rho * u.x() * u.y() * u.y());
            acc[ref::rho_u_v_w](rho * u.x() * u.y() * u.z());
            acc[ref::rho_u_w_w](rho * u.x() * u.z() * u.z());
            acc[ref::rho_v_v_v](rho * u.y() * u.y() * u.y());
            acc[ref::rho_v_v_w](rho * u.y() * u.y() * u.z());
            acc[ref::rho_v_w_w](rho * u.y() * u.z() * u.z());
            acc[ref::rho_w_w_w](rho * u.z() * u.z() * u.z());

            acc[ref::rho2_u_u_u](rho2 * u.x() * u.x() * u.x());
            acc[ref::rho2_u_u_v](rho2 * u.x() * u.x() * u.y());
            acc[ref::rho2_u_u_w](rho2 * u.x() * u.x() * u.z());
            acc[ref::rho2_u_v_v](rho2 * u.x() * u.y() * u.y());
            acc[ref::rho2_u_v_w](rho2 * u.x() * u.y() * u.z());
            acc[ref::rho2_u_w_w](rho2 * u.x() * u.z() * u.z());
            acc[ref::rho2_v_v_v](rho2 * u.y() * u.y() * u.y());
            acc[ref::rho2_v_v_w](rho2 * u.y() * u.y() * u.z());
            acc[ref::rho2_v_w_w](rho2 * u.y() * u.z() * u.z());
            acc[ref::rho2_w_w_w](rho2 * u.z() * u.z() * u.z());

            acc[ref::rho_T_u](rho * T * u.x());
            acc[ref::rho_T_v](rho * T * u.y());
            acc[ref::rho_T_w](rho * T * u.z());

            acc[ref::rho2_T_u](rho2 * T * u.x());
            acc[ref::rho2_T_v](rho2 * T * u.y());
            acc[ref::rho2_T_w](rho2 * T * u.z());

            acc[ref::rho_E_u](e * u.x());
            acc[ref::rho_E_v](e * u.y());
            acc[ref::rho_E_w](e * u.z());

            acc[ref::rho2_E_u](e * m.x());
            acc[ref::rho2_E_v](e * m.y());
            acc[ref::rho2_E_w](e * m.z());

            acc[ref::rho_E_E](e * e / rho);

            acc[ref::rho2_E_E](e * e);

            acc[ref::rho_a](rho * a);

            acc[ref::rho2_a](rho2 * a);

            acc[ref::rho_mu](rho * mu);

            acc[ref::rho2_mu](rho2 * mu);

            acc[ref::mu_Sxx](mu * ( grad_u(0,0)                  - div_u / 3));
            acc[ref::mu_Sxy](mu * ((grad_u(0,1) + grad_u(1,0))/2            ));
            acc[ref::mu_Sxz](mu * ((grad_u(0,2) + grad_u(2,0))/2            ));
            acc[ref::mu_Syy](mu * ( grad_u(1,1)                  - div_u / 3));
            acc[ref::mu_Syz](mu * ((grad_u(1,2) + grad_u(2,1))/2            ));
            acc[ref::mu_Szz](mu * ( grad_u(2,2)                  - div_u / 3));

            acc[ref::mu_div_u](mu * div_u);

            acc[ref::mu_gradx_T](mu * grad_T.x());
            acc[ref::mu_grady_T](mu * grad_T.y());
            acc[ref::mu_gradz_T](mu * grad_T.z());

            acc[ref::u_u](u.x() * u.x());
            acc[ref::u_v](u.x() * u.y());
            acc[ref::u_w](u.x() * u.z());
            acc[ref::v_v](u.y() * u.y());
            acc[ref::v_w](u.y() * u.z());
            acc[ref::w_w](u.z() * u.z());

            acc[ref::u_u_u](u.x() * u.x() * u.x());
            acc[ref::u_u_v](u.x() * u.x() * u.y());
            acc[ref::u_u_w](u.x() * u.x() * u.z());
            acc[ref::u_v_v](u.x() * u.y() * u.y());
            acc[ref::u_v_w](u.x() * u.y() * u.z());
            acc[ref::u_w_w](u.x() * u.z() * u.z());
            acc[ref::v_v_v](u.y() * u.y() * u.y());
            acc[ref::v_v_w](u.y() * u.y() * u.z());
            acc[ref::v_w_w](u.y() * u.z() * u.z());
            acc[ref::w_w_w](u.z() * u.z() * u.z());

            acc[ref::T_u](T * u.x());
            acc[ref::T_v](T * u.y());
            acc[ref::T_w](T * u.z());

            acc[ref::a_u](a * u.x());
            acc[ref::a_v](a * u.y());
            acc[ref::a_w](a * u.z());

            acc[ref::omx](om.x());
            acc[ref::omy](om.y());
            acc[ref::omz](om.z());

            acc[ref::omx_omx](om.x() * om.x());
            acc[ref::omx_omy](om.x() * om.y());
            acc[ref::omx_omz](om.x() * om.z());
            acc[ref::omy_omy](om.y() * om.y());
            acc[ref::omy_omz](om.y() * om.z());
            acc[ref::omz_omz](om.z() * om.z());

            acc[ref::rho_omx](rho * om.x());
            acc[ref::rho_omy](rho * om.y());
            acc[ref::rho_omz](rho * om.z());

            acc[ref::rho2_omx](rho2 * om.x());
            acc[ref::rho2_omy](rho2 * om.y());
            acc[ref::rho2_omz](rho2 * om.z());

            acc[ref::rho_omx_omx](rho * om.x() * om.x());
            acc[ref::rho_omx_omy](rho * om.x() * om.y());
            acc[ref::rho_omx_omz](rho * om.x() * om.z());
            acc[ref::rho_omy_omy](rho * om.y() * om.y());
            acc[ref::rho_omy_omz](rho * om.y() * om.z());
            acc[ref::rho_omz_omz](rho * om.z() * om.z());

            acc[ref::rho2_omx_omx](rho2 * om.x() * om.x());
            acc[ref::rho2_omx_omy](rho2 * om.x() * om.y());
            acc[ref::rho2_omx_omz](rho2 * om.x() * om.z());
            acc[ref::rho2_omy_omy](rho2 * om.y() * om.y());
            acc[ref::rho2_omy_omz](rho2 * om.y() * om.z());
            acc[ref::rho2_omz_omz](rho2 * om.z() * om.z());

            // Compute any slow growth forcing applied to the flow
            largo_state state(e, m.x(), m.y(), m.z(), rho, p);
            largo_state slowgrowth_f;
            sg_treater.inner_xz(state, slowgrowth_f);

            // Gather basic slow growth forcing
            acc[ref::SrhoE      ](slowgrowth_f.e  );
            acc[ref::Srhou      ](slowgrowth_f.mx );
            acc[ref::Srhov      ](slowgrowth_f.my );
            acc[ref::Srhow      ](slowgrowth_f.mz );
            acc[ref::Srho       ](slowgrowth_f.rho);
            acc[ref::Srhou_dot_u](  slowgrowth_f.mx*u.x()
                                  + slowgrowth_f.my*u.y()
                                  + slowgrowth_f.mz*u.z());

            // Gather squared slow growth forcing to permit computing variances
            acc[ref::S2rhoE      ](std::pow(slowgrowth_f.e  , 2));
            acc[ref::S2rhou      ](std::pow(slowgrowth_f.mx , 2));
            acc[ref::S2rhov      ](std::pow(slowgrowth_f.my , 2));
            acc[ref::S2rhow      ](std::pow(slowgrowth_f.mz , 2));
            acc[ref::S2rho       ](std::pow(slowgrowth_f.rho, 2));
            acc[ref::S2rhou_dot_u](std::pow(  slowgrowth_f.mx*u.x()
                                            + slowgrowth_f.my*u.y()
                                            + slowgrowth_f.mz*u.z(), 2));
        } // end X // end Z

        // All accumulators should have seen a consistent number of samples
#ifndef NDEBUG
        const std::size_t miscreant = inconsistent_accumulation_count(acc);
        assert(!miscreant);
#endif

        // Extract y-specific sums into MPI-reduction-ready storage for y(j)
        using boost::accumulators::sum;
#define EXTRACT(quantity, component, offset, description) \
        ret->quantity()(j, offset) = sum(acc[ref::component]);
        SUZERAIN_SAMPLES_COMPONENTS_FOR_EACH(
                EXTRACT, SUZERAIN_SAMPLES_PHYSICAL)
#undef  EXTRACT

    } // end Y

    // Notice dgrid.rank_zero_zero_modes already contains "wave-sampled"
    // quantities while other ranks have zeros in those locations.

    // Allreduce to obtain global sums on every rank
    SUZERAIN_MPICHKR(MPI_Allreduce(
            MPI_IN_PLACE, ret->storage.data(), ret->storage.size(),
            mpi::datatype<samples::storage_type::Scalar>::value,
            MPI_SUM, MPI_COMM_WORLD));

    // Physical space sums, which are at collocation points, need to be
    // scaled by the dealiased extents and converted to coefficients.
    ret->physical() *= dgrid.chi();
    otool.masslu()->solve(ret->physical().outerSize(),
                          ret->physical().data(),
                          ret->physical().innerStride(),
                          ret->physical().outerStride());

    // Fill with NaNs those samples that were not computed by this method
    ret->implicit().setConstant(std::numeric_limits<real_t>::quiet_NaN());

    return ret;
}

// This logic is a trimmed down version of take_samples.
// See comments there.
std::auto_ptr<profile>
take_profile(const definition_scenario &scenario,
             const operator_tools &otool,
             contiguous_state<4,complex_t> &swave)
{
    // State enters method as coefficients in X, Y, and Z directions
    SUZERAIN_TIMER_SCOPED("perfect::take_profile");

    // Shorthand for the operator_tools member(s) commonly used
    const pencil_grid &dgrid = otool.dgrid;

    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };
    assert(static_cast<int>(ndx::e  ) < state_count);
    assert(static_cast<int>(ndx::mx ) < state_count);
    assert(static_cast<int>(ndx::my ) < state_count);
    assert(static_cast<int>(ndx::mz ) < state_count);
    assert(static_cast<int>(ndx::rho) < state_count);

    // Sanity check incoming swave's shape and contiguity
    const std::size_t Ny = swave.shape()[1];
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());
    SUZERAIN_ENSURE((unsigned) swave.strides()[1] == 1u);
    SUZERAIN_ENSURE((unsigned) swave.strides()[2] == swave.shape()[1]);
    SUZERAIN_ENSURE((unsigned) swave.strides()[3] == swave.shape()[1]*swave.shape()[2]);

    // Rank-specific details accumulated in ret to be MPI_Allreduce-d later
    std::auto_ptr<profile> ret(new profile(Ny));

    // Obtain samples available in wave-space from mean conserved state.
    // These coefficients are inherently averaged across the X-Z plane.
    if (dgrid.has_zero_zero_modes()) {
        ret->rho_u().col(0) = Map<VectorXc>(swave[ndx::mx ].origin(), Ny).real();
        ret->rho_u().col(1) = Map<VectorXc>(swave[ndx::my ].origin(), Ny).real();
        ret->rho()          = Map<VectorXc>(swave[ndx::rho].origin(), Ny).real();
    }

    // Transform to obtain physical space view of state on collocation points
    physical_view<state_count> sphys(dgrid, swave);
    for (std::size_t i = 0; i < state_count; ++i) {
        otool.zero_dealiasing_modes(swave, i);
        otool.bop_apply(0, 1., swave, i);
        dgrid.transform_wave_to_physical(&sphys.coeffRef(i,0));
    }

    // Physical space is traversed linearly using a single offset 'offset'.
    // The two loop structure is used as positions x(i), and z(k) unnecessary.
    for (int offset = 0, j = dgrid.local_physical_start.y();
         j < dgrid.local_physical_end.y();
         ++j) {

        // Type used to accumulate means versus wall-normal position
        typedef boost::accumulators::accumulator_set<
                real_t, boost::accumulators::stats<
                    boost::accumulators::tag::sum_kahan
                >
            > accumulator_type;

        // Accumulators for each mean quantity computed in physical space.
        // For example, quantity "foo" has accumulator "sum_foo".
        // Declared within Y loop so they are reset on each Y iteration.
#define DECLARE(r, data, tuple)                                              \
        accumulator_type BOOST_PP_CAT(sum_,BOOST_PP_TUPLE_ELEM(2, 0, tuple)) \
                [BOOST_PP_TUPLE_ELEM(2, 1, tuple)];
        BOOST_PP_SEQ_FOR_EACH(DECLARE,, SUZERAIN_PROFILE_PHYSICAL)
#undef DECLARE

        // Iterate across the j-th ZX plane
        const int last_zxoffset = offset
                                + dgrid.local_physical_extent.z()
                                * dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // Unpack conserved state
            const real_t   e(sphys(ndx::e,   offset));
            const Vector3r m(sphys(ndx::mx,  offset),
                             sphys(ndx::my,  offset),
                             sphys(ndx::mz,  offset));
            const real_t rho(sphys(ndx::rho, offset));

            // Compute quantities of interest
            const Vector3r u = rholut::u(rho, m);
            real_t p, T, mu, lambda;
            rholut::p_T_mu_lambda(scenario.alpha, scenario.beta,
                                  scenario.gamma, scenario.Ma,
                                  rho, m, e, p, T, mu, lambda);

            // Accumulate into sum_XXX using function syntax.
            sum_a [0](std::sqrt(T));
            sum_H0[0]((e + p) / rho);
            sum_ke[0](u.squaredNorm() / 2);
            sum_mu[0](mu);
            sum_T [0](T);
            sum_u [0](u.x());
            sum_u [1](u.y());

        } // end X // end Z

        // Move y-specific sums into MPI-reduction-ready storage for y(j) using
        // Eigen comma initialization syntax.  Yes, this is not stride 1.
#define EXTRACT_SUM(z, n, q) boost::accumulators::sum(sum_##q[n])
#define MOVE_SUM_INTO_TMP(r, data, tuple)                                  \
        ret->BOOST_PP_TUPLE_ELEM(2, 0, tuple)().row(j) <<                  \
             BOOST_PP_ENUM(BOOST_PP_TUPLE_ELEM(2, 1, tuple),               \
                           EXTRACT_SUM, BOOST_PP_TUPLE_ELEM(2, 0, tuple));
        BOOST_PP_SEQ_FOR_EACH(MOVE_SUM_INTO_TMP,, SUZERAIN_PROFILE_PHYSICAL)
#undef EXTRACT_SUM
#undef MOVE_SUM_INTO_TMP

    } // end Y

    // Notice dgrid.rank_zero_zero_modes already contains "wave-sampled"
    // profile while other ranks have zeros in those locations.

    // Allreduce to obtain global sums on every rank
    SUZERAIN_MPICHKR(MPI_Allreduce(
            MPI_IN_PLACE, ret->storage.data(), ret->storage.size(),
            mpi::datatype<profile::storage_type::Scalar>::value,
            MPI_SUM, MPI_COMM_WORLD));

    // Physical space sums, which are at collocation points, need to be
    // scaled by the dealiased extents and converted to coefficients.
    ret->physical() *= dgrid.chi();
    otool.masslu()->solve(ret->physical().outerSize(),
                          ret->physical().data(),
                          ret->physical().innerStride(),
                          ret->physical().outerStride());

    // BEWARE: The 4-argument take_profile(...) below relies on
    // the fact that swave has been converted to physical space
    // as of the end of this method.

    return ret;
}

std::auto_ptr<profile>
take_profile(const definition_scenario &scenario,
             const operator_tools &otool,
             contiguous_state<4,complex_t> &swave,
             instantaneous &inst)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };
    assert(static_cast<int>(ndx::e  ) < state_count);
    assert(static_cast<int>(ndx::mx ) < state_count);
    assert(static_cast<int>(ndx::my ) < state_count);
    assert(static_cast<int>(ndx::mz ) < state_count);
    assert(static_cast<int>(ndx::rho) < state_count);

    // Uses that above take_profile(...) destroys swave by conversion into sphys
    std::auto_ptr<profile> ret = take_profile(scenario, otool, swave);
    physical_view<state_count> sphys(otool.dgrid, swave);
    collect_instantaneous(scenario, otool.grid, otool.dgrid, sphys, inst);
    return ret;
}

void
collect_references(const definition_scenario &scenario,
                   const specification_grid& grid,
                   const pencil_grid &dgrid,
                   const physical_view<5> &sphys,
                   references &refs)
{
    SUZERAIN_TIMER_SCOPED("collect_references");

    // To avoid accumulating garbage, must zero y(j) not present on rank.
    // Clearing everything is expected to be a bit more performant.
    refs.set_zero(dgrid.global_physical_extent.y());

    // Sum reference quantities as a function of y(j) into refs
    for (int offset = 0, j = dgrid.local_physical_start.y();
        j < dgrid.local_physical_end.y();
        ++j) {

        // An array of summing_accumulator_type holds all running sums
        array<summing_accumulator_type, references::q::count> acc;

        // Redmine #2983 disables mu, lambda on the freestream boundary to
        // adhere to superset of Poinsot and Lele subsonic NRBC conditions
        const bool locallyviscous
            = !(grid.one_sided() && j+1 == dgrid.global_physical_extent.y());

        const int last_zxoffset = offset
                                + dgrid.local_physical_extent.z()
                                * dgrid.local_physical_extent.x();
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
                scenario.alpha, scenario.beta, scenario.gamma, scenario.Ma,
                rho, m, e, p, T, mu, lambda);

            // The linearization reference quantity implementation
            // requires a single global 1 / Re value, so unlike
            // the RHS implementation below setting 1 / Re = 0,
            // here we locally turn off viscosity as necessary.
            mu     *= static_cast<int>(locallyviscous);
            lambda *= static_cast<int>(locallyviscous);

            // Accumulate reference quantities into running sums...
            acc[references::q::rho](rho);
            acc[references::q::p  ](p);
            acc[references::q::p2 ](p*p);
            acc[references::q::T  ](T);
            acc[references::q::a  ](sqrt(T));

            // ...including simple velocity-related quantities...
            const Vector3r u = rholut::u(rho, m);
            acc[references::q::u ](u.x());
            acc[references::q::v ](u.y());
            acc[references::q::w ](u.z());
            acc[references::q::u2](u.squaredNorm());
            acc[references::q::uu](u.x()*u.x());
            acc[references::q::uv](u.x()*u.y());
            acc[references::q::uw](u.x()*u.z());
            acc[references::q::vv](u.y()*u.y());
            acc[references::q::vw](u.y()*u.z());
            acc[references::q::ww](u.z()*u.z());

            // ...including simple viscosity-related quantities...
            const real_t nu = mu / rho;
            acc[references::q::nu   ](nu);
            acc[references::q::nu_u ](nu*u.x());
            acc[references::q::nu_v ](nu*u.y());
            acc[references::q::nu_w ](nu*u.z());
            acc[references::q::nu_u2](nu*u.squaredNorm());
            acc[references::q::nu_uu](nu*u.x()*u.x());
            acc[references::q::nu_uv](nu*u.x()*u.y());
            acc[references::q::nu_uw](nu*u.x()*u.z());
            acc[references::q::nu_vv](nu*u.y()*u.y());
            acc[references::q::nu_vw](nu*u.y()*u.z());
            acc[references::q::nu_ww](nu*u.z()*u.z());

            // ...other, more complicated expressions...
            const Vector3r e_gradrho
                    = rholut::explicit_div_e_plus_p_u_refcoeff_grad_rho(
                            scenario.gamma, rho, m, e, p);
            acc[references::q::ex_gradrho](e_gradrho.x());
            acc[references::q::ey_gradrho](e_gradrho.y());
            acc[references::q::ez_gradrho](e_gradrho.z());

            acc[references::q::e_divm](
                    rholut::explicit_div_e_plus_p_u_refcoeff_div_m(
                        rho, e, p));

            acc[references::q::e_deltarho](
                    rholut::explicit_mu_div_grad_T_refcoeff_div_grad_rho(
                        scenario.gamma, mu, rho, e, p));

            // ...and, lastly, details needed for slow growth forcing.
            acc[references::q::rhou ](m.x()             );
            acc[references::q::rhov ](m.y()             );
            acc[references::q::rhow ](m.z()             );
            acc[references::q::rhoE ](e                 );
            acc[references::q::rhouu](m.x()* m.x() / rho);
            acc[references::q::rhouv](m.x()* m.y() / rho);
            acc[references::q::rhouw](m.x()* m.z() / rho);
            acc[references::q::rhovv](m.y()* m.y() / rho);
            acc[references::q::rhovw](m.y()* m.z() / rho);
            acc[references::q::rhoww](m.z()* m.z() / rho);
            acc[references::q::rhoEE](e    * e     / rho);

        } // end X // end Z

        // All accumulators should have seen a consistent number of samples
        // Failure usually indicates a coding indicator on add new quantities
#ifndef NDEBUG
        const std::size_t miscreant = inconsistent_accumulation_count(acc);
        assert(!miscreant);
#endif

        // Store sums into common block in preparation for MPI Allreduce
        for (int i = 0; i < static_cast<int>(references::q::count); ++i) {
            refs.row(i)[j] = boost::accumulators::sum(acc[i]);
        }

    } // end Y

    // Allreduce and scale refs sums to obtain means on all ranks
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, refs.data(), refs.size(),
            mpi::datatype<references::Scalar>::value, MPI_SUM, MPI_COMM_WORLD));
    refs *= dgrid.chi();
}

// This is a trimmed down version of collect_references.  See comments there.
void
collect_instantaneous(const definition_scenario &scenario,
                      const specification_grid &grid,
                      const pencil_grid &dgrid,
                      const physical_view<5> &sphys,
                      instantaneous &inst)
{
    SUZERAIN_TIMER_SCOPED("collect_instantaneous");

    // To avoid accumulating garbage, must zero y(j) not present on rank.
    // Clearing everything is expected to be a bit more performant.
    inst.set_zero(dgrid.global_physical_extent.y());

    // Sum reference quantities as a function of y(j) into inst*
    for (int offset = 0, j = dgrid.local_physical_start.y();
        j < dgrid.local_physical_end.y();
        ++j) {

        // An array of summing_accumulator_type holds all running sums
        array<summing_accumulator_type, instantaneous::q::count> acc;

        // Redmine #2983 disables mu, lambda on the freestream boundary to
        // adhere to superset of Poinsot and Lele subsonic NRBC conditions
        const bool locallyviscous
            = !(grid.one_sided() && j+1 == dgrid.global_physical_extent.y());

        const int last_zxoffset = offset
                                + dgrid.local_physical_extent.z()
                                * dgrid.local_physical_extent.x();
        for (; offset < last_zxoffset; ++offset) {

            // Unpack conserved state
            const real_t   e  (sphys(ndx::e,   offset));
            const Vector3r m  (sphys(ndx::mx,  offset),
                               sphys(ndx::my,  offset),
                               sphys(ndx::mz,  offset));
            const real_t   rho(sphys(ndx::rho, offset));

            // Compute quantities related to the equation of state
            real_t p;
            rholut::p(
                scenario.alpha, scenario.beta, scenario.gamma, scenario.Ma,
                rho, m, e, p);

            // Ask yourself if it should be used for viscous quantities...
            SUZERAIN_UNUSED(locallyviscous);

            // Accumulate reference quantities into running sums...
            acc[instantaneous::q::rho](rho);
            acc[instantaneous::q::p  ](p);
            acc[instantaneous::q::p2 ](p*p);

            // ...including simple velocity-related quantities...
            const Vector3r u = rholut::u(rho, m);
            acc[instantaneous::q::u ](u.x());
            acc[instantaneous::q::v ](u.y());
            acc[instantaneous::q::w ](u.z());
            acc[instantaneous::q::uu](u.x()*u.x());
            acc[instantaneous::q::uv](u.x()*u.y());
            acc[instantaneous::q::uw](u.x()*u.z());
            acc[instantaneous::q::vv](u.y()*u.y());
            acc[instantaneous::q::vw](u.y()*u.z());
            acc[instantaneous::q::ww](u.z()*u.z());

            // ...and, lastly, details needed for slow growth forcing.
            acc[instantaneous::q::rhou ](m.x()             );
            acc[instantaneous::q::rhov ](m.y()             );
            acc[instantaneous::q::rhow ](m.z()             );
            acc[instantaneous::q::rhoE ](e                 );
            acc[instantaneous::q::rhouu](m.x()* m.x() / rho);
            acc[instantaneous::q::rhouv](m.x()* m.y() / rho);
            acc[instantaneous::q::rhouw](m.x()* m.z() / rho);
            acc[instantaneous::q::rhovv](m.y()* m.y() / rho);
            acc[instantaneous::q::rhovw](m.y()* m.z() / rho);
            acc[instantaneous::q::rhoww](m.z()* m.z() / rho);
            acc[instantaneous::q::rhoEE](e    * e     / rho);

        } // end X // end Z

        // All accumulators should have seen a consistent number of samples
        // Failure usually indicates a coding indicator on add new quantities
#ifndef NDEBUG
        const std::size_t miscreant = inconsistent_accumulation_count(acc);
        assert(!miscreant);
#endif

        // Store sums into common block in preparation for MPI Allreduce
        for (int i = 0; i < static_cast<int>(instantaneous::q::count); ++i) {
            inst.col(i)[j] = boost::accumulators::sum(acc[i]);
        }

    } // end Y

    // Allreduce and scale inst sums to obtain means on all ranks
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE, inst.data(), inst.size(),
            mpi::datatype<instantaneous::Scalar>::value,
            MPI_SUM, MPI_COMM_WORLD));
    inst *= dgrid.chi();
}

void
compute_fluctuations(const bsplineop& cop,
                     const bsplineop_lu &masslu,
                     const instantaneous& inst,
                     ArrayX6r &reynolds,
                     ArrayX6r &favre,
                     ArrayX3r &prodterms)
{
    // Compute Reynolds-averaged velocity fluctuations
    reynolds.resize(inst.rows(), NoChange);
    reynolds.col(0) = inst.uu() - inst.u().square();
    reynolds.col(1) = inst.uv() - inst.u()*inst.v();
    reynolds.col(2) = inst.uw() - inst.u()*inst.w();
    reynolds.col(3) = inst.vv() - inst.v().square();
    reynolds.col(4) = inst.vw() - inst.v()*inst.w();
    reynolds.col(5) = inst.ww() - inst.w().square();

    // Prepare pointwise \bar{\rho} \widetilde{u_i'' u_j''}
    // which is NOT the final result for the method
    favre.resize(inst.rows(), NoChange);
    favre.col(0)=(inst.rhouu() - inst.rhou().square()   /inst.rho());
    favre.col(1)=(inst.rhouv() - inst.rhou()*inst.rhov()/inst.rho());
    favre.col(2)=(inst.rhouw() - inst.rhou()*inst.rhow()/inst.rho());
    favre.col(3)=(inst.rhovv() - inst.rhov().square()   /inst.rho());
    favre.col(4)=(inst.rhovw() - inst.rhov()*inst.rhow()/inst.rho());
    favre.col(5)=(inst.rhovw() - inst.rhow().square()   /inst.rho());

    // Compute the three separate contributions to production:
    //     - \bar{\rho} \widetilde{u''\otimes{}u''} : \nabla\tilde{u} =
    //         - bar_rho(y)*(    tilde_upp_vpp*tilde_u__y
    //                         + tilde_vpp_vpp*tilde_v__y
    //                         + tilde_vpp_wpp*tilde_w__y )
    //
    // Begin by forming pointwise \partial_y \tilde{u}_i
    prodterms.resize(inst.rows(), NoChange);
    prodterms.col(0) = inst.rhou();
    prodterms.col(1) = inst.rhov();
    prodterms.col(2) = inst.rhow();
    prodterms.colwise() /= inst.rho();
    masslu.solve(prodterms.cols(), prodterms.data(),
                 prodterms.innerStride(), prodterms.outerStride());
    cop.apply(1, prodterms.cols(), 1.0, prodterms.data(),
              prodterms.innerStride(), prodterms.outerStride());
    // End by scaling with the relevant fluctuation data.
    // Notice \bar{\rho} already present in contents of tmp.
    prodterms.col(0) *= - favre.col(1);
    prodterms.col(1) *= - favre.col(3);
    prodterms.col(2) *= - favre.col(4);

    // Divide \bar{\rho} \widetilde{u_i'' u_j''} by \bar{\rho}
    // to get the desired Favre-averaged fluctuations.
    favre.colwise() /= inst.rho();
}

void summarize_boundary_layer_nature(
        const profile &prof,
        const definition_scenario &scenario,
        const shared_ptr<specification_largo> &sg,
        const bsplineop_lu &masslu,
        bspline &b,
        suzerain_bl_local       &wall,
        suzerain_bl_viscous     &viscous,
        suzerain_bl_thicknesses &thick,
        suzerain_bl_local       &edge,
        suzerain_bl_local       &edge99,
        suzerain_bl_reynolds    &reynolds,
        suzerain_bl_qoi         &qoi,
        suzerain_bl_pg          &pg)
{
    // Repeated warnings from this method are annoying, but I hate to
    // summarily make them debug messages as they can be of genuine use
    // during postprocessing and to understand why NaNs arise.
    // TODO Permit library caller to specify warn once vs warn always
#   define WARNONCE0(l, m)                                              \
    do { static bool warn = true;                                       \
        if (warn) { WARN0(l, m); warn = false; } else { DEBUG0(l, m); } \
    } while (0);

    // Care taken to avoid any B-spline evaluation using NaN thicknesses
    // per Redmine #2940 as this can bring down a simulation from GSL_ERROR.
    using boost::math::isnan;

    // Prepare local state at the wall (y=0)
    // Uses B-spline 0th coefficient being wall value to reduce costs.
    std::fill_n(reinterpret_cast<double *>(&wall),
                sizeof(wall)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    wall.a     = prof.a()[0];
    wall.gamma = scenario.gamma;
    wall.mu    = prof.mu()[0];
    wall.Pr    = scenario.Pr;
    wall.rho   = prof.rho()[0];
    assert(&(wall.T) + 1 == &(wall.T__y)); // Next, compute both T and T__y
    b.linear_combination(1, prof.T().col(0).data(), 0.0, &(wall.T), 1);
    assert(&(wall.u) + 1 == &(wall.u__y)); // Next, compute both u and u__y
    b.linear_combination(1, prof.u().col(0).data(), 0.0, &(wall.u), 1);
    wall.v     = prof.u().col(1)[0];
    wall.y     = 0.0;

    // Compute viscous quantities based only on wall information
    if (const int err = suzerain_bl_compute_viscous(
                            scenario.Re,
                            &wall,
                            &viscous)) {
        WARNONCE0("profile",
                  "suzerain_bl_compute_viscous(...) returned "
                  << err << " (" << suzerain_strerror(err) << ")");
    }

    // Prepare inviscid baseflow coefficients if necessary
    ArrayX4r coeffs_inviscid; // Columns are H0, rhou, u, v
    if (sg && sg->formulation.enabled() && sg->baseflow) {
        coeffs_inviscid.resize(b.n(), NoChange);

        // Retrieve values on collocation points
        largo_state state;
        baseflow_interface &baseflow = *sg->baseflow;
        for (int j = 0; j < b.n(); ++j) {
            const real_t y_j = b.collocation_point(j);
            baseflow.conserved(y_j, state.as_is());
            baseflow.pressure (y_j, state.p);
            coeffs_inviscid(j, 0) = state.H0();
            coeffs_inviscid(j, 1) = state.mx;
            coeffs_inviscid(j, 2) = state.u();
            coeffs_inviscid(j, 3) = state.v();
        }

        // Convert to coefficients using mass matrix
        masslu.solve(coeffs_inviscid.outerSize(),
                     coeffs_inviscid.data(),
                     coeffs_inviscid.innerStride(),
                     coeffs_inviscid.outerStride());
    }

    // Compute boundary layer thicknesses, including delta
    if (0 == coeffs_inviscid.size()) {
        if (const int err = suzerain_bl_compute_thicknesses(
                                prof.H0().data(),
                                prof.ke().data(),
                                prof.rho_u().col(0).data(),
                                prof.u().col(0).data(),
                                &thick,
                                b.bw,
                                b.dbw)) {
            WARNONCE0("profile",
                      "suzerain_bl_compute_thicknesses(...) returned "
                      << err << " (" << suzerain_strerror(err) << ")");
        }
    } else {
        if (const int err = suzerain_bl_compute_thicknesses_baseflow(
                                prof.H0().data(),
                                prof.ke().data(),
                                prof.rho_u().col(0).data(),
                                prof.u().col(0).data(),
                                coeffs_inviscid.innerStride(),
                                coeffs_inviscid.col(0).data(),
                                coeffs_inviscid.col(1).data(),
                                coeffs_inviscid.col(2).data(),
                                coeffs_inviscid.col(3).data(),
                                &thick,
                                b.bw,
                                b.dbw)) {
            WARNONCE0("profile",
                      "suzerain_bl_compute_thicknesses_baseflow(...) returned "
                      << err << " (" << suzerain_strerror(err) << ")");
        }
    }

    // Evaluate state at the edge   (y=thick.delta  ) from B-spline coefficients
    // Evaluate state at the edge99 (y=thick.delta99) from B-spline coefficients
    // As this is futzy, use a loop to repeat the computations for each target
    {
        const double       y[2] = { thick.delta, thick.delta99 };
        suzerain_bl_local *e[2] = { &edge,       &edge99       };
        for (size_t i = 0; i < sizeof(y)/sizeof(y[0]); ++i) {
            std::fill_n(reinterpret_cast<double *>(e[i]),
                        sizeof(*e[i])/sizeof(double),
                        std::numeric_limits<double>::quiet_NaN());
            e[i]->gamma = scenario.gamma;
            e[i]->Pr    = scenario.Pr;
            e[i]->y     = y[i];
            if (SUZERAIN_UNLIKELY((isnan)(y[i]))) {
                // NOP as fill_n above ensured NaN propagated correctly
            } else {
                // Later, could more quickly evaluate basis once and then re-use that
                // result to repeatedly form the necessary linear combinations.
                b.linear_combination(0, prof.a().data(),        y[i], &e[i]->a);
                b.linear_combination(0, prof.mu().data(),       y[i], &e[i]->mu);
                b.linear_combination(0, prof.rho().data(),      y[i], &e[i]->rho);
                assert(&e[i]->T + 1 == &e[i]->T__y); // Next, both T and T__y
                b.linear_combination(1, prof.T().col(0).data(), y[i], &e[i]->T, 1);
                assert(&e[i]->u + 1 == &e[i]->u__y); // Next, both u and u__y
                b.linear_combination(1, prof.u().col(0).data(), y[i], &e[i]->u, 1);
                b.linear_combination(0, prof.u().col(1).data(), y[i], &e[i]->v);
            }
        }
    }

    // Compute Reynolds numbers
    if (0 == coeffs_inviscid.size()) {
        if (const int err = suzerain_bl_compute_reynolds(
                                scenario.Re,
                                &edge,
                                &edge99,
                                &thick,
                                &reynolds)) {
            WARNONCE0("profile",
                      "suzerain_bl_compute_reynolds(...) returned "
                      << err << " (" << suzerain_strerror(err) << ")");
        }
    } else {
        if (const int err = suzerain_bl_compute_reynolds_baseflow(
                                scenario.Re,
                                prof.H0().data(),
                                prof.ke().data(),
                                prof.rho_u().col(0).data(),
                                prof.u().col(0).data(),
                                1,
                                coeffs_inviscid.col(0).data(),
                                coeffs_inviscid.col(1).data(),
                                coeffs_inviscid.col(2).data(),
                                coeffs_inviscid.col(3).data(),
                                &edge,
                                &edge99,
                                &reynolds,
                                b.bw)) {
            WARNONCE0("profile",
                      "suzerain_bl_compute_reynolds_baseflow(...) returned "
                      << err << " (" << suzerain_strerror(err) << ")");
        }
    }

    // Compute general quantities of interest
    suzerain_bl_compute_qoi(scenario.Ma, scenario.Re,
                            &wall, &viscous, &edge99, &thick, &qoi);

    // Mean pressure and streamwise velocity gradients come from slow growth
    double edge99_p__x = 0, edge99_u__x = 0;
    if (sg && sg->formulation.enabled() && sg->baseflow) {
        const double delta99 = thick.delta99;
        if (SUZERAIN_UNLIKELY((isnan)(delta99))) {
            edge99_p__x = edge99_u__x = std::numeric_limits<real_t>::quiet_NaN();
        } else {
            largo_state base, dy, dx; // as_is()
            sg->baseflow->conserved(delta99,
                                    base.as_is(), dy.as_is(), dx.as_is());
            sg->baseflow->pressure (delta99,
                                    base.p, dy.p, dx.p);
            edge99_p__x = dx.p;                                       // Direct
            edge99_u__x = (dx.mx - dx.rho*base.mx/base.rho)/base.rho; // Chained
        }
    }
    if (const int err = suzerain_bl_compute_pg(
                            scenario.Ma,
                            scenario.Re,
                            &wall,
                            &viscous,
                            &edge99,
                            edge99_p__x,
                            edge99_u__x,
                            &thick,
                            &pg)) {
        WARNONCE0("profile",
                  "suzerain_bl_compute_pg(...) returned "
                  << err << " (" << suzerain_strerror(err) << ")");
    }

#   undef WARNONCE0
}

void summarize_channel_nature(
        const profile &prof,
        const definition_scenario &scenario,
        bspline &b,
        suzerain_channel_local   &wall,
        suzerain_channel_viscous &viscous,
        suzerain_channel_local   &center,
        suzerain_channel_qoi     &qoi)
{
    // Prepare local state at the wall (y=0)
    // Uses B-spline 0th coefficient being wall value to reduce costs.
    std::fill_n(reinterpret_cast<double *>(&wall),
                sizeof(wall)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    wall.y     = b.collocation_point(0);
    wall.a     = prof.a()[0];
    wall.gamma = scenario.gamma;
    wall.mu    = prof.mu()[0];
    wall.Pr    = scenario.Pr;
    wall.rho   = prof.rho()[0];
    assert(&(wall.T) + 1 == &(wall.T__y)); // Next, compute both T and T__y
    b.linear_combination(1, prof.T().col(0).data(), wall.y, &(wall.T), 1);
    assert(&(wall.u) + 1 == &(wall.u__y)); // Next, compute both u and u__y
    b.linear_combination(1, prof.u().col(0).data(), wall.y, &(wall.u), 1);
    wall.v     = prof.u().col(1)[0];

    // Compute viscous quantities based only on wall information
    suzerain_channel_compute_viscous(scenario.Re, &wall, &viscous);

    // Evaluate state at the centerline from B-spline basis
    std::fill_n(reinterpret_cast<double *>(&center),
                sizeof(center)/sizeof(double),
                std::numeric_limits<double>::quiet_NaN());
    center.y = b.collocation_point(b.n()-1) / 2;
    b.linear_combination(0, prof.a().data(),        center.y, &(center.a));
    center.gamma = scenario.gamma;
    b.linear_combination(0, prof.mu().data(),       center.y, &(center.mu));
    center.Pr    = scenario.Pr;
    b.linear_combination(0, prof.rho().data(),      center.y, &(center.rho));
    assert(&(center.T) + 1 == &(center.T__y)); // Next, compute both T and T__y
    b.linear_combination(1, prof.T().col(0).data(), center.y, &(center.T), 1);
    assert(&(center.u) + 1 == &(center.u__y)); // Next, compute both u and u__y
    b.linear_combination(1, prof.u().col(0).data(), center.y, &(center.u), 1);
    b.linear_combination(0, prof.u().col(1).data(), center.y, &(center.v));

    // Compute general quantities of interest
    suzerain_channel_compute_qoi(scenario.Ma, scenario.Re,
                                 &wall, &viscous, &center, &qoi);
}

} // namespace perfect

} // namespace suzerain
