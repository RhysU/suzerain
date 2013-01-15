//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
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

/** @file
 * @copydoc perfect.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include "perfect.hpp"

#include <esio/error.h>
#include <gsl/gsl_errno.h>
#include <sys/file.h>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/coalescing_pool.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/htstretch.h>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/rngstream.hpp>
#include <suzerain/shared_range.hpp>
#include <suzerain/support/field.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

using boost::numeric_cast;
using std::size_t;

namespace suzerain {

namespace perfect {

std::vector<support::field> default_fields()
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

void save_collocation_values(
        const esio_handle h,
        contiguous_state<4,complex_t>& swave,
        const scenario_definition& scenario,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  swave.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(swave.shape()[3]) == dgrid.local_wave_extent.z());

    // Convert coefficients into collocation point values
    // Transforms state from full-wave coefficients to full-physical points
    operator_tools otool(grid, dgrid, cop);
    for (size_t i = 0; i < swave.shape()[0]; ++i) {
        otool.bop_apply(0, 1, swave, i);
        otool.zero_dealiasing_modes(swave, i);
        dgrid.transform_wave_to_physical(
                reinterpret_cast<real_t *>(swave[i].origin()));
    }

    // Convert conserved rho{_E,_u,_v,_w,} into u, v, w, p, T
    physical_view<state_count> sphys(dgrid, swave);

    const real_t alpha = scenario.alpha;
    const real_t beta  = scenario.beta;
    const real_t gamma = scenario.gamma;
    const real_t Ma    = scenario.Ma;

    for (int o = 0; o < dgrid.local_physical_extent.prod(); ++o) {
        // Unpack conserved quantities from fields
        const real_t     e(sphys(ndx::e,   o));
        Vector3r         m(sphys(ndx::mx,  o),
                           sphys(ndx::my,  o),
                           sphys(ndx::mz,  o));
        const real_t   rho(sphys(ndx::rho, o));

        // Compute primitive quantities to be stored
        real_t p, T;
        rholut::p_T(alpha, beta, gamma, Ma, rho, m, e, p, T);
        m /= rho;

        // Pack primitive quantities back into fields (by position)
        sphys(0, o) = m.x(); // Now just X velocity
        sphys(1, o) = m.y(); // Now just Y velocity
        sphys(2, o) = m.z(); // Now just Z velocity
        sphys(3, o) = p;
        sphys(4, o) = T;
    }

    // HDF5 file storage locations and corresponding descriptions
    const array<const char *,state_count> prim_names = {{
        "u", "v", "w", "p", "T"
    }};
    const array<const char *,state_count> prim_descriptions = {{
        "streamwise velocity",
        "wall-normal velocity",
        "spanwise velocity",
        "pressure",
        "temperature",
    }};
    assert(prim_names.static_size == prim_descriptions.static_size);

    // Establish size of collective writes across all ranks and write fields
    esio_field_establish(h, grid.dN.y(), dgrid.local_physical_start.y(),
                                         dgrid.local_physical_extent.y(),
                            grid.dN.z(), dgrid.local_physical_start.z(),
                                         dgrid.local_physical_extent.z(),
                            grid.dN.x(), dgrid.local_physical_start.x(),
                                         dgrid.local_physical_extent.x());

    for (size_t i = 0; i < state_count; ++i) {

        std::string comment = "Nondimensional ";
        comment += prim_descriptions[i];
        comment += " stored row-major YZX on the 3D rectilinear grid defined"
                   " by taking the outer product of arrays"
                   " /collocation_points_y, /collocation_points_z, and"
                   " /collocation_points_z";

        esio_field_write(h, prim_names[i],
                reinterpret_cast<real_t *>(swave[i].origin()),
                0, 0, 0, comment.c_str());
    }
}

void load_collocation_values(
        const esio_handle h,
        contiguous_state<4,complex_t>& state,
        const scenario_definition& scenario,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop,
        const bspline& b)
{
    // We are only prepared to handle a fixed number of fields in this routine
    enum { state_count = 5 };

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(                  state.shape()[0]  == state_count);
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[1]) == dgrid.local_wave_extent.y());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[2]) == dgrid.local_wave_extent.x());
    SUZERAIN_ENSURE(numeric_cast<int>(state.shape()[3]) == dgrid.local_wave_extent.z());

    // This routine does no grid interpolation.  Yell loudly if necessary
    {
        // Check that restart file size matches runtime dealiased extents
        int cg, bg, ag;
        if (ESIO_SUCCESS != esio_field_size(h, "u", &cg, &bg, &ag)) {
            SUZERAIN_ERROR_VOID("Unable to find /u field size from restart",
                                SUZERAIN_EFAILED);
        }
        if (cg != grid.dN.y() || bg != grid.dN.z() || ag != grid.dN.x()) {
            ERROR0("Physical space restart fields have row-major YZX extents "
                   << "(" << cg << "," << bg << "," << ag << ")" << " but "
                   << "(" << grid.dN.y() << "," << grid.dN.z() << ","
                   << grid.dN.x() << ") are required");
            SUZERAIN_ERROR_VOID(
                    "Cannot interpolate during physical space restart",
                    SUZERAIN_EFAILED);
        }


        // Check that restart file specifies the same B-spline basis.
        // TODO Too restrictive!  Identical collocation points would be okay.
        // TODO Too restrictive?  Any floating point differences kill us.
        shared_ptr<bspline> Fb;
        shared_ptr<bsplineop> Fbop;
        support::load(h, Fb, Fbop);
        const double bsp_dist = b.distance_to(*Fb);
        const bool   bsp_same = bsp_dist < suzerain_bspline_distance_distinct;
        if (!bsp_same) {
            ERROR0("Physical restart has different wall-normal bases ("
                   << bsp_dist << ")");
            SUZERAIN_ERROR_VOID(
                    "Cannot interpolate during physical space restart",
                    SUZERAIN_EFAILED);
        }
    }

    // Establish size of collective reads across all ranks and read data
    physical_view<state_count> sphys(dgrid, state);
    esio_field_establish(h, grid.dN.y(), dgrid.local_physical_start.y(),
                                         dgrid.local_physical_extent.y(),
                            grid.dN.z(), dgrid.local_physical_start.z(),
                                         dgrid.local_physical_extent.z(),
                            grid.dN.x(), dgrid.local_physical_start.x(),
                                         dgrid.local_physical_extent.x());
    esio_field_read(h, "u", &sphys(0,0), 0, 0, 0);
    esio_field_read(h, "v", &sphys(1,0), 0, 0, 0);
    esio_field_read(h, "w", &sphys(2,0), 0, 0, 0);
    esio_field_read(h, "p", &sphys(3,0), 0, 0, 0);
    esio_field_read(h, "T", &sphys(4,0), 0, 0, 0);

    // Convert primitive u, v, w, p, and T into rho{_E,_u,_v,_w,}
    const real_t gamma = scenario.gamma;
    const real_t Ma    = scenario.Ma;

    for (int o = 0; o < dgrid.local_physical_extent.prod(); ++o) {
        // Unpack primitive quantities from fields (by position)
        Vector3r       m(sphys(0, o),  // Now just X velocity
                         sphys(1, o),  // Now just Y velocity
                         sphys(2, o)); // Now just Z velocity
        const real_t   p(sphys(3, o));
        const real_t   T(sphys(4, o));

        // Compute conserved quantities from primitive ones
        const real_t rho = gamma * p / T;   // Assumes EOS
        m               *= rho;             // Now m contains momentum
        const real_t e   = rholut::energy_kinetic(Ma, rho, m)
                         + rholut::energy_internal(gamma, p);

        // Pack conserved quantities into fields (by name)
        sphys(ndx::e,   o) = e;
        sphys(ndx::mx,  o) = m.x();
        sphys(ndx::my,  o) = m.y();
        sphys(ndx::mz,  o) = m.z();
        sphys(ndx::rho, o) = rho;
    }

    // Collectively convert physical state to wave space coefficients
    // Build FFT normalization constant into Y direction's mass matrix
    bsplineop_luz massluz(cop);
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.opform(1, &scale_factor, cop);
    massluz.factor();
    operator_tools otool(grid, dgrid, cop);
    for (size_t i = 0; i < state.shape()[0]; ++i) {
        dgrid.transform_physical_to_wave(&sphys.coeffRef(i, 0)); // X, Z
        otool.zero_dealiasing_modes(state, i);
        otool.bop_solve(massluz, state, i);                      // Y
    }
}

void load(const esio_handle h,
          contiguous_state<4,complex_t>& state,
          const scenario_definition& scenario,
          const grid_specification& grid,
          const pencil_grid& dgrid,
          const bsplineop& cop,
          const bspline& b)
{
    const std::vector<support::field> fields = default_fields();

    // Ensure state storage meets this routine's assumptions
    SUZERAIN_ENSURE(state.shape()[0]  == fields.size());

    // Check whether load_coefficients(...) should work
    bool trycoeffs = true;
    for (size_t i = 0; i < fields.size(); ++i) {
        int ncomponents = 0;
        switch (esio_field_sizev(h, fields[i].location.c_str(),
                                 0, 0, 0, &ncomponents))
        {
            case ESIO_SUCCESS:
                if (ncomponents != 2) {
                    DEBUG0("Field /" << fields[i].identifier
                                     << " looks fishy...");
                }
                break;
            case ESIO_NOTFOUND:
                DEBUG0("Field /" << fields[i].identifier
                                 << " not found in restart");
                trycoeffs = false;
                break;
            default:
                DEBUG0("Field /" << fields[i].identifier << " looks fishy...");
                break;
        }
    }

    // Dispatch to the appropriate restart loading logic
    DEBUG0("Started loading simulation fields");
    if (trycoeffs) {
        support::load_coefficients(h, fields, state, grid, dgrid, cop, b);
    } else {
        INFO0("Loading collocation-based, physical-space restart data");
        load_collocation_values(h, state, scenario, grid, dgrid, cop, b);
    }
    DEBUG0("Finished loading simulation fields");
}

void
adjust_scenario(contiguous_state<4,complex_t> &swave,
                const scenario_definition& scenario,
                const grid_specification& grid,
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
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
    massluz.opform(1, &scale_factor, cop);
    massluz.factor();
    for (size_t i = 0; i < state_count; ++i) {
        dgrid.transform_physical_to_wave(&sphys.coeffRef(i, 0));
        otool.bop_solve(massluz, swave, i);
    }
}

void
add_noise(contiguous_state<4,complex_t> &state,
          const noise_specification& noise,
          const scenario_definition& scenario,
          const grid_specification& grid,
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
                0, &state[ndx::mx][0][0][0], 1, &centerline, &momentum);
        INFO0("Centerline mean streamwise momentum at y = "
              << centerline << " is " << momentum);
        b.linear_combination(
                0, &state[ndx::rho][0][0][0], 1, &centerline, &density);
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
    const complex_t scale_factor = grid.dN.x() * grid.dN.z();
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

} // namespace perfect

} // namespace suzerain
