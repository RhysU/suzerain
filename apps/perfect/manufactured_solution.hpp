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

#ifndef SUZERAIN_PERFECT_MANUFACTURED_SOLUTION_HPP
#define SUZERAIN_PERFECT_MANUFACTURED_SOLUTION_HPP

/** @file
 * Method of manufactured solution details for the perfect gas application.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/support/definition_base.hpp>

#include "nsctpl_rholut.hpp"

namespace suzerain {

// Forward declarations
class grid_specification;
class pencil_grid;

namespace perfect {

// Forward declarations
class scenario_definition;

/**
 * Manufactured solution employed throughout the perfect gas code.
 *
 * By subclassing \ref support::definition_base we can pass an instance to
 * program_options::add_definition method to permit modifying parameters via
 * the command line.
 *
 * None of #alpha, #beta, #gamma, #Ma, #Pr, #Re or #Lx, #Ly, #Lz are modified
 * when used in this fashion.  The former group will need to be synced using
 * \ref match(const scenario_definition&) and the later group using \ref
 * match(const grid_specification&).
 */
class manufactured_solution
    : public nsctpl_rholut::manufactured_solution<real_t>,
      public support::definition_base
{
public:

    /**
     * Default caption used to collect together options per
     * support::program_options::add_definition.  Useful as a building block
     * when a non-default value is desired.
     */
    static const std::string default_caption;

    /** Default constructor using #default_caption. */
    manufactured_solution();

    /** Constructor permitting a non-default caption. */
    explicit manufactured_solution(const std::string& caption);

    /** Set #alpha, #beta, #gamma, #Ma, #Pr, and #Re to match \c scenario. */
    manufactured_solution& match(const scenario_definition& scenario);

    /** Set #Lx, #Ly, and #Lz to match \c grid. */
    manufactured_solution& match(const grid_specification& grid);

    /**
     * Set #rho, #u, #v, #w, and #T per nsctpl_rholut::isothermal_channel.
     *
     * Afterwards, \ref match(const scenario_definition&) and \ref match(const
     * grid_specification&) may need to be called.
     */
    manufactured_solution& isothermal_channel();

    /**
     * Set #rho, #u, #v, #w, and #T per nsctpl_rholut::isothermal_flat_plate.
     *
     * @copydetails isothermal_channel
     */
    manufactured_solution& isothermal_flat_plate();

    /** @copydoc support::definition_base::options_description() */
    boost::program_options::options_description options_description();

private:

    std::string caption;

};

/**
 * Save manufactured solution parameters in the given \c location in an
 * ESIO-based file.  Parameters are only stored when \c msoln evaluates as \c
 * true.
 */
void save(const esio_handle h,
          const shared_ptr<manufactured_solution> & msoln,
          const scenario_definition& scenario,
          const grid_specification& grid,
          const char *location = "manufactured_solution");

/**
 * Load manufactured solution parameters from the given \c location from an
 * ESIO-based file.
 *
 * If the file contains active manufactured solution parameters, \c msoln will
 * be modified to contain an appropriate instance.  Otherwise, \c msoln will be
 * reset and evaluate to \c false in a boolean context.
 */
void load(const esio_handle h,
          shared_ptr<manufactured_solution>& msoln,
          const scenario_definition& scenario,
          const grid_specification& grid,
          const char *location = "manufactured_solution");

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
        const bsplineop &cop,
        bspline &b,
        const real_t simulation_time);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_MANUFACTURED_SOLUTION_HPP
