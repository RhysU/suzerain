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
// driver_base.hpp: Application driver logic spanning multiple applications
// $Id$

#ifndef SUZERAIN_SUPPORT_DRIVER_BASE_HPP
#define SUZERAIN_SUPPORT_DRIVER_BASE_HPP

#include <suzerain/common.hpp>

#include <suzerain/restart_definition.hpp>
#include <suzerain/signal_definition.hpp>
#include <suzerain/statistics_definition.hpp>
#include <suzerain/timecontroller.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/timestepper.hpp>

#include "application_base.hpp"

namespace suzerain {

namespace support {

/** Provides well-scoped names related to driver signal handling */
namespace signal {

/** Potentially coincident actions to take due to POSIX signal receipt. */
enum action_type
{
    /** Output a status message. */
    log_status,

    /** Write a restart file */
    write_restart,

    /** Tear down the simulation (proactively due to --advance_wt limit). */
    teardown_reactive,

    /** Tear down the simulation (reactively due to an incoming signal). */
    teardown_proactive,

    /** Compute and write a statistics file. */
    write_statistics,

    /** A sentry giving the number of distinct actions possible. */
    count  // MUST BE LAST
};

/**
 * Type of global, atomic locations used to track incoming signal actions.
 * To be indexed using \ref action_type.
 */
typedef array<
        volatile sig_atomic_t, static_cast<std::size_t>(count)
    > volatile_received_type;

/**
 * Global, atomic locations used to track local POSIX signal receipt.  Beware
 * that usage by multiple clients, e.g. multiple \ref delta_t_allreducer
 * instances, may cause race conditions.
 */
extern volatile_received_type global_received;

/**
 * Type of possibly local, non-volatile locations used to track global receipt
 * of the same actions as \ref volatile_received_type.  To be indexed using
 * \ref action_type.
 */
typedef array<
        sig_atomic_t, static_cast<std::size_t>(count)
    > received_type;

} // end namespace signal

class delta_t_allreducer;

class field;

/**
 * An abstract driver base class for managing a Suzerain application.
 * Intended for time-varying, three-dimensional problems.
 * Instantiate from within \c main().
 *
 * Signal handling capabilities may misbehave if multiple instances
 * are executing within the same process.
 */
class driver_base : public application_base
{
    /** Provides simple access to the superclass type */
    typedef application_base super;

public:

    /** Type used to express simulation time quantities. */
    typedef timecontroller<real_t>::time_type time_type;

    /** Type used to express discrete simulation step quantities. */
    typedef timecontroller<real_t>::step_type step_type;

    /** @copydoc application_base::application_base */
    driver_base(const std::string &application_synopsis,
                const std::string &description = "",
                const std::string &revstr = "");

    /** Virtual destructor as appropriate for abstract base class */
    virtual ~driver_base();

    /** @copydoc application_base::initialize */
    virtual std::vector<std::string> initialize(int argc, char **argv);

    /**
     * The log4cxx configuration to use.  Files <tt>bulk.dat</tt>,
     * <tt>L2.mean.dat</tt>, and <tt>rms.fluct.dat</tt> collecting messages with
     * the names <tt>bulk</tt>, <tt>L2.mean</tt>, and <tt>rms.fluct</tt> have
     * been added.
     */
    virtual std::string log4cxx_config();

    /**
     * Information about the state fields.  Should be populated by subclasses
     * prior to \ref save_restart or \ref load_restart.
     */
    std::vector<support::field> fields;

    /** Operational details on saving restart files. */
    shared_ptr<restart_definition> restartdef;

    /** Operational details on saving statistics files. */
    shared_ptr<statistics_definition> statsdef;

    /** Operational details on advancing simulation time. */
    shared_ptr<time_definition> timedef;

    /** Controls the OS signals triggering various types of processing. */
    static signal_definition signaldef;

    /** Controls low storage method to be used for time advance. */
    shared_ptr<timestepper::lowstorage::method_interface<
            complex_t
        > > method; // FIXME Default in advance routine

    /**
     * Refers to a linear operator instance able to interoperate between
     * #state_linear and #state_nonlinear.  This is the interface to which
     * linear operators should be coded.
     */
    shared_ptr<timestepper::lowstorage::linear_operator<
                state_common_type, state_nonlinear_type
            > > L;

    /**
     * Refers to a nonlinear operator instance able to operate on
     * #state_nonlinear.  This is the interface to which nonlinear operators
     * should be coded.
     */
    shared_ptr<timestepper::nonlinear_operator<
                state_nonlinear_type
            > > N;

    /** Controls time advance, including callback processing. */
    shared_ptr<timecontroller<real_t> > tc; // FIXME

    /**
     * Did the previous time advance end in a predicted, controlled manner?
     */
    bool soft_teardown;

    /**
     * Build a fixed-width, human-friendly way to output the given simulation
     * time and time step number.  Care is taken to ensure sequential outputs
     * are minimally distinct.
     *
     * @param t  Simulation time to output
     * @param nt Simulation time step to output
     *
     * @return a string suitable for use in output status messages.
     */
    virtual std::string build_timeprefix(
            const time_type t,
            const step_type nt);

    /**
     * Routine to output status, generally called via the timecontroller. FIXME
     *
     * Invokes \ref log_status_bulk, \ref log_status_L2, \ref
     * log_status_boundary_state, and \ref log_status_hook.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool log_status(
            const time_type t,
            const step_type nt);

    /** Log messages containing mean L2 and RMS fluctuation information. */
    virtual void log_status_L2(
            const std::string& timeprefix,
            const char * const name_L2  = "L2.mean",
            const char * const name_rms = "rms.fluct");

    /** Log messages containing bulk quantities. */
    virtual void log_status_bulk(
            const std::string& timeprefix);

    /**
     * Log messages containing specific state quantities at the upper and lower
     * boundaries.  Density is reported as-is.  All other scalars are divided
     * by density.
     */
    virtual void log_status_specific_boundary_state(
            const std::string& timeprefix);

    // FIXME
    /**
     * Load the contents of a restart file into #state_nonlinear using the
     * current decomposition.  Subclasses should override this method
     * adding any desired functionality either before or after invoking
     * the superclass version.
     *
     * @param[in]  esioh An ESIO handle pointing to an open restart file.
     * @param[out] t     The simulation time stored in the restart file.
     */
    virtual void load_restart(
            esio_handle esioh,
            real_t &t);

    /**
     * Save time-independent metadata that should appear in all restart files.
     * Subclasses should not generally override this method but should instead
     * use \ref save_restart_metadata_hook.
     *
     * @see Member #restartdef to control the metadata file location.
     */
    virtual void save_restart_metadata();

    /**
     * Save state into a restart file.  The method \ref save_restart_metadata()
     * must have been called first as restart files contained cloned metadata.
     * Subclasses should not generally override this method but should instead
     * use \ref save_restart_hook.
     *
     * @param t  The simulation time to be stored in the restart file.
     * @param nt The time step number which is not stored in the restart file.
     *           No restart file is written when multiple invocations are
     *           performed in succession on the same \c nt.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @see Member #restartdef to control restart writing options.
     */
    virtual bool save_restart(
            const time_type t,
            const step_type nt);

    /**
     * Save statistics into a sample file.  The method \ref save_restart_metadata()
     * must have been called first as sample files contain cloned metadata.
     * Subclasses should not generally override this method but should instead
     * use \ref save_statistics_hook.
     *
     * @param t  The simulation time to be stored in the statistics file.
     * @param nt The time step number which is not stored in the file.
     *           No file is written when multiple invocations are
     *           performed in succession on the same \c nt.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @see Member #statsdef to control statistical sample writing options.
     */
    virtual bool save_statistics(
            const time_type t,
            const step_type nt);

    /**
     * Type maintaining the mean ratio of each delta_t_candidate to the minimum
     * delta_t_candidate selected during each time step.  Useful for
     * determining the relative restrictiveness of each stability criterion.
     */
    typedef std::vector<boost::accumulators::accumulator_set<
                real_t,
                boost::accumulators::stats<boost::accumulators::tag::mean>
        > > delta_t_ratios_type;

    /**
     * Maintains the mean ratio of each delta_t_candidate to the minimum
     * delta_t_candidate selected during each time step.
     */
    delta_t_ratios_type delta_t_ratios;

protected:

    /**
     * Extension point to permit adding arbitrary metadata to all restart and
     * statistics files during \ref save_restart_metadata.
     *
     * Subclasses should override this method adding or changing any desired
     * functionality either before or after invoking the superclass
     * implementation.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     */
    virtual void save_restart_metadata_hook(
            esio_handle esioh);

    /**
     * Extension point to permit adding arbitrary information into restart
     * files during \ref save_restart.  Subclasses should override this method
     * with the desired functionality.  Invoking the superclass method in the
     * override is optional.
     *
     * The default implementation saves the contents of #state_linear into a
     * restart file destroying #state_nonlinear in the process.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     */
    virtual bool save_restart_hook(
            esio_handle esioh);

    /**
     * Extension point to permit adding arbitrary information into statistical
     * sample files during \ref save_statistics.  Subclasses should override
     * this method with the desired functionality.  Invoking the superclass
     * method in the override is optional.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     */
    virtual bool save_statistics_hook(
            esio_handle esioh);

    /**
     * Hook permitting subclasses to output additional status information.
     * Invoked at the end of \ref log_status.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool log_status_hook(
            const std::string& timeprefix,
            const time_type t,
            const step_type nt);

    /**
     * Flag used to control whether \ref log_status_L2 shows headers.
     * The default implementation disables headers after the first invocation.
     */
    bool log_status_L2_show_header;

    /**
     * Flag used to control whether \ref log_status_bulk shows headers.
     * The default implementation disables headers after the first invocation.
     */
    bool log_status_bulk_show_header;

    /** Wall time elapsed during loading of state from the restart file */
    double wtime_load_state;

    /** Wall time at which we began time stepping */
    double wtime_advance_start;

    // FIXME
    /** Signal handler which mutates \c atomic_signal_received. */
    static void process_signal(const int sig);

    /** Tracks last time a status line was output */
    step_type last_status_nt;

    /** Tracks last time a restart file was written successfully */
    step_type last_restart_saved_nt;

    /** Tracks last time a statistics sample file was written successfully */
    step_type last_statistics_saved_nt;

    /**
     * Maintains if any signals were observed in \ref signal::global_received
     * during the last time step during a time advance.  To be indexed using
     * \ref signal::action_type.
     */
    signal::received_type signal_received;

private:

    /** Was a metadata file ever saved to disk? */
    bool metadata_created;

public:

    /**
     * Handles \c MPI_Allreduce calls necessary to obtain stable time steps on
     * all ranks as well as aggregating POSIX signals received on any rank.
     */
    const scoped_ptr<delta_t_allreducer> allreducer;

};

/**
 * A stateful functor that performs a MPI_Allreduce to determine the minimum
 * stable time step size across all ranks.  The same MPI Allreduce is used to
 * hide the cost of querying \ref signal::global_received across all ranks.
 */
class delta_t_allreducer : public timestepper::delta_t_reducer
{
    /** Provides simple access to the superclass type */
    typedef timestepper::delta_t_reducer super;

public:

    /** Construct an instance referring to external state. */
    delta_t_allreducer(
            const double& wtime_mpi_init,
            const double& wtime_fftw_planning,
            const shared_ptr<time_definition>& timedef,
            const double& wtime_load_state,
            const double& wtime_advance_start,
            const driver_base::step_type& last_status_nt,
            const driver_base::step_type& last_restart_saved_nt,
            driver_base::delta_t_ratios_type& delta_t_ratios,
            signal::received_type& signal_received)
        : wtime_mpi_init(wtime_mpi_init)
        , wtime_fftw_planning(wtime_fftw_planning)
        , timedef(timedef)
        , wtime_load_state(wtime_load_state)
        , wtime_advance_start(wtime_advance_start)
        , last_status_nt(last_status_nt)
        , last_restart_saved_nt(last_restart_saved_nt)
        , delta_t_ratios(delta_t_ratios)
        , signal_received(signal_received)
    {}

    /**
     * Returns the smallest entry in \c delta_t_candidates from any rank.  Must
     * be called collectively by all ranks in \c MPI_COMM_WORLD.  Mutates \ref
     * signal::global_received.
     */
    real_t operator()(const std::vector<real_t>& delta_t_candidates);

private:

    /** Reference to content maintained by a \ref application_base instance. */
    const double& wtime_mpi_init;

    /** Reference to content maintained by a \ref application_base instance. */
    const double& wtime_fftw_planning;

    /** Reference to content maintained by a \ref driver_base instance. */
    const shared_ptr<time_definition>& timedef;

    /** Reference to content maintained by a \ref driver_base instance. */
    const double& wtime_load_state;

    /** Reference to content maintained by a \ref driver_base instance. */
    const double& wtime_advance_start;

    /** Reference to content maintained by a \ref driver_base instance. */
    const driver_base::step_type& last_status_nt;

    /** Reference to content maintained by a \ref driver_base instance. */
    const driver_base::step_type& last_restart_saved_nt;

    /** Reference to content maintained by a \ref driver_base instance. */
    driver_base::delta_t_ratios_type& delta_t_ratios;

    /** Reference to content maintained by a \ref driver_base instance. */
    signal::received_type& signal_received;

    /**
     * Maintains coarse statistics on the wall time duration between calls.
     * This provides a rank-specific measure of the variability per time step
     * which includes all registered periodic callbacks (e.g. status and
     * restart processing).
     */
    boost::accumulators::accumulator_set<
            real_t,
            boost::accumulators::stats<boost::accumulators::tag::max,
                                       boost::accumulators::tag::mean,
                                       boost::accumulators::tag::variance>
        > period;
};

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_DRIVER_BASE_HPP
