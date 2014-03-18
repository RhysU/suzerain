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

#ifndef SUZERAIN_PERFECT_MANUFACTURED_SOLUTION_HPP
#define SUZERAIN_PERFECT_MANUFACTURED_SOLUTION_HPP

/** @file
 * Method of manufactured solution details for the perfect gas application.
 */

#include <suzerain/common.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/esio_fwd.hpp>
#include <suzerain/support/overridable.hpp>
#include <suzerain/support/populatable.hpp>

SUZERAIN_GCC_DIAG_OFF(unused-parameter);
SUZERAIN_GCC_DIAG_OFF(unused-variable);
#include "nsctpl_rholut.hpp"
SUZERAIN_GCC_DIAG_ON(unused-variable);
SUZERAIN_GCC_DIAG_ON(unused-parameter);

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;
class pencil_grid;
class specification_grid;

namespace perfect {

// Forward declarations
class definition_scenario;

/**
 * Manufactured solution employed throughout the perfect gas code.
 *
 * By subclassing \ref support::definition_base we can pass an instance to
 * program_options::add_definition method to permit modifying parameters via
 * the command line.
 *
 * None of #alpha, #beta, #gamma, #Ma, #Pr, #Re or #Lx, #Ly, #Lz are modified
 * when used in this fashion.  The former group will need to be synced using
 * \ref match(const definition_scenario&) and the later group using \ref
 * match(const specification_grid&).
 */
class manufactured_solution
    : public virtual support::definition_base
    , public virtual support::overridable<manufactured_solution>
    , public virtual support::populatable<manufactured_solution>
    , public nsctpl_rholut::manufactured_solution<real_t>
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
    manufactured_solution& match(const definition_scenario& scenario);

    /** Set #Lx, #Ly, and #Lz to match \c grid. */
    manufactured_solution& match(const specification_grid& grid);

    /**
     * Set #rho, #u, #v, #w, and #T per nsctpl_rholut::isothermal_channel.
     *
     * Afterwards, \ref match(const definition_scenario&) and \ref match(const
     * specification_grid&) may need to be called.
     */
    manufactured_solution& isothermal_channel();

    /**
     * Set #rho, #u, #v, #w, and #T per nsctpl_rholut::isothermal_flat_plate.
     *
     * @copydetails isothermal_channel
     */
    manufactured_solution& isothermal_flat_plate();

    /** @copydoc support::definition_base::options_description() */
    virtual boost::program_options::options_description options_description();

    /** @copydoc support::populatable::override */
    virtual void override(
            const manufactured_solution& that,
            const bool verbose = false);

    /** @copydoc support::populatable::populate */
    virtual void populate(
            const manufactured_solution& that,
            const bool verbose = false);

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
          const definition_scenario& scenario,
          const specification_grid& grid,
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
          const definition_scenario& scenario,
          const specification_grid& grid,
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
        const specification_grid &grid,
        const pencil_grid &dgrid,
        const bsplineop &cop,
        bspline &b,
        const real_t simulation_time);

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_MANUFACTURED_SOLUTION_HPP
