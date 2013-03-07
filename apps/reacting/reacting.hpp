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

#ifndef SUZERAIN_REACTING_HPP
#define SUZERAIN_REACTING_HPP

/** @file
 * Support logic for the reacting flow application.
 */

#include <esio/esio.h>
#ifdef HAVE_UNDERLING
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <underling/error.h>
#include <underling/underling.hpp>
#include <underling/underling_fftw.hpp>
#endif

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/noise_specification.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timestepper.hpp>

#include "scenario_definition.hpp"

namespace suzerain {

// Forward declaration
namespace support { class field; }

/**
 * Functionality used throughout the Suzerain reacting flow application.
 */
namespace reacting {

/** Return default nondimensional field information per \ref suzerain::ndx. */
std::vector<support::field> default_fields();

/**
 * Store the current simulation primitive state as collocation point values
 * into an open restart file.  Note that <tt>state</tt>'s contents are
 * destroyed.  Collocation point values required only for dealiasing purposes
 * <i>are</i> stored but only those points informed by non-dealiased state.
 * This method is less efficient and the restart data less flexible than that
 * produced by save_coefficients().
 */
void save_collocation_values(
        const esio_handle h,
        contiguous_state<4,complex_t>& swave,
        const scenario_definition& scenario,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop);

/**
 * Load the current simulation state from an open collocation point value
 * restart file.  Cannot handle interpolating onto a different grid.
 */
void load_collocation_values(
        const esio_handle h,
        contiguous_state<4,complex_t>& state,
        const scenario_definition& scenario,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        const bsplineop& cop,
        const bspline& b);

/**
 * Interrogate an open restart file and invoke either load_coefficients()
 * or load_collocation_values() as necessary.
 */
void load(const esio_handle h,
          contiguous_state<4,complex_t>& state,
          const scenario_definition& scenario,
          const grid_specification& grid,
          const pencil_grid& dgrid,
          const bsplineop& cop,
          const bspline& b);

/**
 * Hold temperature and density constant while changing the Mach number and
 * ratio of specific heats.  On input, \c state should contain total energy
 * fields using \c old_Ma and \c old_gamma.  On output \c state will contain
 * total energy fields using <tt>scenario.Ma</tt> and <tt>scenario.gamma</tt>.
 */
void
adjust_scenario(contiguous_state<4,complex_t> &swave,
                const scenario_definition& scenario,
                const grid_specification& grid,
                const pencil_grid& dgrid,
                const bsplineop& cop,
                const real_t old_Ma,
                const real_t old_gamma);

/**
 * Add random momentum field perturbations ("noise") according to
 * the provided noise_definition.
 */
void
add_noise(contiguous_state<4,complex_t> &state,
          const noise_specification& noise,
          //const scenario_definition& scenario,
          const grid_specification& grid,
          const pencil_grid& dgrid,
          const bsplineop& cop,
          bspline &b);

} // end namespace reacting

} // end namespace suzerain

#endif // SUZERAIN_REACTING_HPP
