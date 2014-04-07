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

#ifndef SUZERAIN_PERFECT_DRIVER_HPP
#define SUZERAIN_PERFECT_DRIVER_HPP

/** @file
 * A driver for nondimensional perfect gas simulations.
 */

#include <suzerain/support/definition_helm.hpp>
#include <suzerain/support/definition_isothermal.hpp>
#include <suzerain/support/definition_largo.hpp>
#include <suzerain/support/definition_radialflow.hpp>
#include <suzerain/support/driver_base.hpp>

#include "common_block.hpp"
#include "manufactured_solution.hpp"
#include "definition_scenario.hpp"

namespace suzerain {

// Forward declarations
class samples;

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
           const std::string &argument_synopsis = "",
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
    shared_ptr<definition_scenario> scenario;

    /** Nondimensional isothermal boundary condition parameters. */
    shared_ptr<support::definition_isothermal> isothermal;

    /** Slow growth parameters for use with the Largo library. */
    shared_ptr<support::definition_largo> sg;

    /** Radial flow baseflow specification for possible use with Largo. */
    shared_ptr<support::definition_radialflow> rad;

    /** Provides PID-based controls of Largo slow growth parameters.  */
    shared_ptr<support::definition_helm> helm;

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
    shared_ptr<samples> mean;

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
     * Log the linearization error present between \ref L and \ref N.
     *
     * More specifically, for \f$\partial_t u = \mathscr{L}u + Lu +
     * \left(N(u)-Lu\right)\f$ how well did \ref N compute a \f$-Lu\f$ term
     * matching with \ref L \f$+Lu\f$ contribution?  \ref L is assumed to also
     * compute a nontrivial \f$\mathscr{L}u\f$ which is independent of
     * reference values (e.g.  the divergence of the momentum within the mass
     * equation).
     *
     * If \ref L and \ref N "match" perfectly, the logged messages would be
     * identically zero.  Seeing more than acceptable (e.g. <tt>1e-10</tt>)
     * floating point error indicates something is amiss.
     *
     * The computation costs a significant fraction of a Runge--Kutta time step
     * but \e is low-storage-friendly.  It destroys the contents of all of \ref
     * common_block, \ref state_linear, and \ref state_nonlinear.  The
     * resultant message should be uniformly small regardless of \ref
     * state_linear but its particular value does depend on \ref state_linear.
     */
    virtual void log_linearization_error(
            const std::string& timeprefix,
            const real_t t,
            const std::size_t nt);

    /**
     * Log a variety of quantities of interest for a boundary layer or channel
     * flow simulation.
     *
     * The computation costs a significant fraction of a Runge--Kutta substep
     * but \e is low-storage-friendly.  It destroys the contents of \ref
     * state_nonlinear.
     */
    virtual void log_quantities_of_interest(
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

    /**
     * Collectively compute and save two-point spectra versus \f$k_x\f$ and
     * \f$k_z\f$ for \f$T\f$, \f$u\f$, \f$v\f$, \f$w\f$, \f$\rho\f$, and all
     * pairwise combinations thereof based on #state_linear destroying
     * #state_nonlinear in the process.
     */
    virtual void save_spectra_primitive(
            const esio_handle esioh);

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
     * Invokes \ref save_spectra_primitive after calling the superclass method.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool save_state_hook(
            const esio_handle esioh);

    /**
     * Beyond the inherited behavior, this method invokes \ref
     * compute_statistics and saves #mean.  Any currently cached results in
     * #mean will be reused <em>without recomputation</em> whenever
     * <tt>this->controller && mean.t == t</tt>.
     * @copydetails driver_base::save_statistics_hook
     */
    virtual bool save_statistics_hook(
            const esio_handle esioh,
            const time_type t);

    /**
     * Beyond the inherited behavior, this method loads #mean.
     * @copydetails application_base::load_restart
     */
    virtual void load_statistics_hook(
            const esio_handle esioh);

    /**
     * The default interval is a fractional flow through time per #grid and
     * possibly #scenario, #isothermal, and #sg.
     */
    virtual void default_restart_interval(time_type& t, step_type&);

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

} // end namespace perfect

} // end namespace suzerain

#endif // SUZERAIN_PERFECT_DRIVER_HPP
