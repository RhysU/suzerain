//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
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
// perfect.hpp: Support logic for the Suzerain perfect gas application
// $Id$

#ifndef SUZERAIN_PERFECT_HPP
#define SUZERAIN_PERFECT_HPP

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
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/inorder.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/support/time_definition.hpp>
#include <suzerain/timestepper.hpp>

#include "nsctpl_rholut.hpp"
#include "scenario_definition.hpp"

namespace suzerain {

/**
 * Functionality used throughout the Suzerain perfect gas application.
 */
namespace perfect {

/** Return default nondimensional field information per \ref suzerain::ndx */
std::vector<support::field> default_fields();

/** Manufactured solution employed throughout the channel code */
typedef nsctpl_rholut::manufactured_solution<real_t> manufactured_solution;

/**
 * Store manufactured solution parameters in a restart file.
 * Parameters are only stored when \c msoln evaluates as true.
 */
void store(const esio_handle h,
           const scenario_definition& scenario,
           const grid_specification& grid,
           const shared_ptr<manufactured_solution> & msoln);

/**
 * Load manufactured solution parameters from a restart file.
 * If the restart file contains active manufactured solution parameters, \c
 * msoln will be modified to contain an appropriate instance.  If it does not,
 * \c msoln will be reset.
 */
void load(const esio_handle h,
          const scenario_definition& scenario,
          const grid_specification& grid,
          shared_ptr<manufactured_solution>& msoln);

/**
 * Store the current simulation primitive state as collocation point values
 * into an open restart file.  Note that <tt>state</tt>'s contents are
 * destroyed.  Collocation point values required only for dealiasing purposes
 * <i>are</i> stored but only those points informed by non-dealiased state.
 * This method is less efficient and the restart data less flexible than that
 * produced by store_coefficients().
 */
void store_collocation_values(
        const esio_handle h,
        contiguous_state<4,complex_t>& swave,
        const scenario_definition& scenario,
        const grid_specification& grid,
        const pencil_grid& dgrid,
        bspline& b,
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
        bspline& b,
        const bsplineop& cop);

/**
 * Interrogate an open restart file and invoke either load_coefficients()
 * or load_collocation_values() as necessary.
 */
void load(const esio_handle h,
          contiguous_state<4,complex_t>& state,
          const scenario_definition& scenario,
          const grid_specification& grid,
          const pencil_grid& dgrid,
          bspline& b,
          const bsplineop& cop);

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
                bspline &b,
                const bsplineop& cop,
                const real_t old_Ma,
                const real_t old_gamma);

/** Options definitions for adding random noise to momentum fields */
class noise_definition : public support::definition_base {

public:

    /** Construct an instance with the given default values */
    explicit noise_definition(real_t fluct_percent = 0,
                              unsigned long fluct_seed = 12345);

    /**
     * Maximum fluctuation magnitude to add as a percentage
     * of centerline streamwise momentum.
     */
    real_t percent;

    /**
     * Fraction of the X direction wavenumbers in [0,1] below
     * which fluctuations will not be added.
     */
    real_t kxfrac_min;

    /**
     * Fraction of the X direction wavenumbers in [0,1] above
     * which fluctuations will not be added.
     */
    real_t kxfrac_max;

    /**
     * Fraction of the Z direction wavenumbers in [0,1] below
     * which fluctuations will not be added.
     */
    real_t kzfrac_min;

    /**
     * Fraction of the Z direction wavenumbers in [0,1] above
     * which fluctuations will not be added.
     */
    real_t kzfrac_max;

    /** rngstream generator seed (see L'Ecuyer et al. 2002) */
    unsigned long seed;

};

/**
 * Add random momentum field perturbations ("noise") according to
 * the provided noise_definition.
 */
void
add_noise(contiguous_state<4,complex_t> &state,
          const noise_definition& noisedef,
          const scenario_definition& scenario,
          const grid_specification& grid,
          const pencil_grid& dgrid,
          bspline &b,
          const bsplineop& cop);

/**
 * Accumulate the result of adding \c alpha times the manufactured solution \c
 * msoln times \c beta times the given wave-space state \c swave.  Setting
 * <tt>alpha=1</tt> and <tt>beta=0</tt> may be used to initialize a
 * manufactured solution field.  Setting <tt>alpha=-1</tt> and <tt>beta=1</tt>
 * may be used to compute error against the manufactured solution.  The
 * manufactured solution lives on \e only the non-dealiased, non-Nyquist modes.
 */
void accumulate_manufactured_solution(
        const real_t alpha,
        const manufactured_solution &msoln,
        const real_t beta,
        contiguous_state<4,complex_t> &swave,
        const grid_specification &grid,
        const pencil_grid &dgrid,
        bspline &b,
        const bsplineop &cop,
        const real_t simulation_time);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_HPP
