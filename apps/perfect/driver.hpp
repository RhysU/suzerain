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
// driver.hpp: A driver for nondimensional perfect gas simulations
// $Id$

#ifndef SUZERAIN_PERFECT_DRIVER_HPP
#define SUZERAIN_PERFECT_DRIVER_HPP

#include <suzerain/common.hpp>
#include <suzerain/support/driver_base.hpp>

#include "manufactured_solution.hpp"
#include "mean_quantities.hpp"
#include "nonlinear_operator_fwd.hpp"
#include "scenario_definition.hpp"

namespace suzerain {

namespace perfect {

/**
 * A driver class for managing a nondimensional perfect gas application using
 * details from \ref suzerain::rholut.  Instantiate from within \c main().
 */
class driver : public support::driver_base
{
    /** Provides simple access to the superclass type. */
    typedef support::driver_base super;

public:

    /** @copydoc driver_base::driver_base */
    driver(const std::string &application_synopsis,
           const std::string &description = "",
           const std::string &revstr = "");

    /**
     * @copydoc driver_base::initialize
     *
     * While #scenario is registered in #options during this method, #msoln is
     * not as not all perfect gas binaries will expose related options.
     */
    virtual std::vector<std::string> initialize(
            int argc,
            char **argv);

    /** Nondimensional scenario parameters used by physics routines. */
    shared_ptr<scenario_definition> scenario;

    /** Nondimensional manufactured solution optionally used by applications. */
    shared_ptr<manufactured_solution> msoln;

    /**
     * Data sharable between #L and #N to permit computing implicit forcing.
     * The user is responsible for wiring #common_block into #L and #N at
     * their construction time.
     */
    operator_common_block common_block;

    /**
     * Maintains instantaneously sampled wall-normal mean quantities.  Member
     * \c mean.t tracks the last time statistics were computed and is used as a
     * mechanism to avoid expensive recomputations.  Implicitly computed mean
     * quantities are tracked via \ref common_block.
     */
    mean_quantities mean;

    /**
     * When \c msoln is true, log messages containing the absolute
     * error in \c state_linear relative to the manufactured solution.
     * Destroys the contents of \c state_nonlinear during execution.
     */
    virtual void log_manufactured_solution_absolute_error(
            const std::string& timeprefix,
            const real_t t,
            const std::size_t nt);

    /**
     * Collectively compute statistics from #state_linear saving them into
     * #mean and destroying #state_nonlinear in the process.
     *
     * @param t Time to be saved in the statistics recorded in #mean.
     */
    virtual void compute_statistics(
            time_type t);

protected:

    /**
     * Beyond the inherited behavior, this method invokes
     * \ref log_manufactured_solution_absolute_error.
     *
     * @returns True or false per the superclass behavior.
     */
    virtual bool log_status_hook(
            const std::string& timeprefix,
            const time_type t,
            const step_type nt);

    /**
     * Beyond the inherited behavior, this method saves #scenario and #msoln.
     * @copydetails driver_base::save_metadata_hook
     */
    virtual void save_metadata_hook(
            const esio_handle esioh);

    /**
     * Beyond the inherited behavior, this method loads #scenario and #msoln.
     * @copydetails driver_base::save_metadata_hook
     */
    virtual void load_metadata_hook(
            const esio_handle esioh);

    /**
     * Beyond the inherited behavior, this method invokes \ref
     * compute_statistics and saves #mean.  Any currently cached results in
     * #mean will be reused <em>without recomputation</em> whenever <tt>mean.t
     * = controller->current_t()</tt>.
     * @copydetails driver_base::save_statistics_hook
     */
    virtual bool save_statistics_hook(
            const esio_handle esioh);

    /**
     * Beyond the inherited behavior, this method loads #mean.
     * @copydetails application_base::load_restart
     */
    virtual void load_statistics_hook(
            const esio_handle esioh);

    /**
     * The default interval is a fractional flow through time per #grid,
     * #bulk_rho_u, and #bulk_rho.
     */
    virtual void default_status_interval(time_type& t, step_type&);

    /** @copydoc default_status_interval */
    virtual void default_restart_interval(time_type& t, step_type&);

    /** @copydoc default_status_interval */
    virtual void default_statistics_interval(time_type& t, step_type&);

};

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_DRIVER_HPP
