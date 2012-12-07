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
#include <suzerain/support/application_base.hpp>
#include <suzerain/support/restart_definition.hpp>
#include <suzerain/support/signal_definition.hpp>
#include <suzerain/support/statistics_definition.hpp>
#include <suzerain/support/time_definition.hpp>
#include <suzerain/timecontroller.hpp>
#include <suzerain/timestepper.hpp>

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
 * instances by multiple \ref driver_base instances, may cause race conditions.
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
 * An driver base class for managing a Suzerain application.  Intended for
 * time-varying, three-dimensional problems.  Provides many hooks to permit
 * lifecycle-related behavior for long-running, parallel applications.
 * Instantiate from within \c main().
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

    /**
     * Controls the OS signals triggering various types of processing.
     * Must be static so that it can be queried within signal handler.
     */
    static signal_definition signaldef;

    /** Controls low storage method to be used for time advance. */
    shared_ptr<timestepper::lowstorage::method_interface<
            complex_t
        > > method;

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
    shared_ptr<timecontroller<real_t> > controller;

    /**
     * Ensure #method is valid for use by #controller.  That is, if
     * <tt>!method</tt> then \ref timestepper::lowstorage::smr91 is used per
     * #timedef.  Otherwise, #method is not modified.
     */
    virtual void prepare_method();

    /**
     * Use instance members to prepare a low-storage timecontroller instance in
     * #controller.  Any existing instance in #controller will be released.
     * All operational members, e.g. #L, #N, #timedef, #statsdef, should have
     * been initialized as desired prior to invocation.  The controller and all
     * appropriate periodic callbacks will be prepared.
     *
     * @param initial_t Initial simulation time to use.
     * @param chi       Time-independent scaling factor for #N on each substep.
     *                  Often used to hide FFT normalization costs.
     */
    virtual void prepare_controller(
            const time_type initial_t,
            const real_t chi);

    /**
     * Wrapper invoking prepare_controller(time_type,real_t) with
     * \c chi set for dealiased FFT normalization per #grid.
     *
     * @param initial_t Initial simulation time defaulting to zero.
     */
    void prepare_controller(
            const time_type initial_t = 0);

    /**
     * Use #controller to advance the simulation per #timedef.
     *
     * Notice <tt>controller</tt> advance routines could be invoked directly,
     * but they provide none of this additional postprocessing functionality.
     *
     * @param final_status    Invoke \ref log_status after advance completes?
     * @param final_restart   Invoke \ref save_restart after successful advance?
     * @param output_timers   Output performance details after advance completes?
     * @param output_stepping Output time step size metrics after advance completes?
     *
     * @return The wall time spent advancing simulation time on success.
     *         The negated wall time spend advancing the simulation on failure.
     */
    virtual double advance_controller(
            const bool final_status    = true,
            const bool final_restart   = true,
            const bool output_timers   = true,
            const bool output_stepping = true);

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
     * Routine to output status, generally called via the timecontroller.
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
     * Log messages containing state at the upper and lower boundaries.
     */
    virtual void log_status_boundary_state(
            const std::string& timeprefix);

    /**
     * Save time-independent metadata that must appear in all restart and
     * statistics files.  Though this logic will be invoked automatically the
     * first time \ref save_restart() or \ref save_statistics() is used, users
     * may wish to invoke the method explicitly during general initialization
     * prior to \ref advance_controller().
     *
     * Subclasses should not override this method but should instead use \ref
     * save_metadata_hook.
     *
     * @see Member #restartdef to control the metadata file location.
     */
    virtual void save_metadata();

    /**
     * Save state \e and statistics into a restart file.  Subclasses should not
     * override this method but should instead use \ref save_state_hook and
     * \ref save_statistics_hook.
     *
     * @param t  The simulation time to be stored in the restart file.
     * @param nt The time step number which is not stored in the restart file.
     *           No restart file is written when multiple invocations are
     *           performed in succession on the same \c nt.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @see Member #restartdef to control restart writing options,
     *      including the name of the file on disk.
     */
    virtual bool save_restart(
            const time_type t,
            const step_type nt);

    /**
     * Save statistics into a statistical sampling file.  Subclasses should not
     * override this method but should instead use \ref save_statistics_hook.
     *
     * @param t  The simulation time to be stored in the statistics file.
     * @param nt The time step number which is not stored in the file.
     *           No file is written when multiple invocations are
     *           performed in succession on the same \c nt.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @see Member #statsdef to control statistical sample writing options,
     *      including the name of the file on disk.
     */
    virtual bool save_statistics(
            const time_type t,
            const step_type nt);

    /**
     * Load time-independent metadata appearing in all restart and statistics
     * files.  Subclasses should not override this method but should instead
     * use \ref load_metadata_hook.
     *
     * @param[in]  esioh An ESIO handle pointing to an open, readable file.
     */
    virtual void load_metadata(
            const esio_handle esioh);

    /**
     * Load the contents of a restart file into #state_linear using the current
     * decomposition.  Subclasses should override this method adding any
     * desired functionality either before or after invoking the superclass
     * version.
     *
     * @param[in]  esioh An ESIO handle pointing to an open restart file.
     * @param[out] t     The simulation time stored in the restart file.
     */
    virtual void load_restart(
            const esio_handle esioh,
            real_t &t);

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
     * Hook permitting subclasses to output additional status information.
     * Invoked at the end of \ref log_status.
     *
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass implementation.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool log_status_hook(
            const std::string& timeprefix,
            const time_type t,
            const step_type nt);

    /**
     * Hook permitting saving arbitrary metadata to all restart and statistics
     * files during \ref save_metadata.
     *
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass implementation.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     */
    virtual void save_metadata_hook(
            const esio_handle esioh);

    /**
     * Hook permitting saving arbitrary information from \ref #state_nonlinear
     * into restart files during \ref save_restart.  Subclasses should override
     * this method with the desired functionality.  Invoking the superclass
     * method in the override is optional.
     *
     * The default implementation saves the contents of #state_nonlinear into
     * the provided ESIO handle file destroying #state_nonlinear in the
     * process.  When \c restartdef->physical is \c false, expansion
     * coefficients are written in all three directions using \ref
     * support::store_coefficients.  When it is \c true, values at collocation
     * points are written using support::store_collocation_Values.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool save_state_hook(
            const esio_handle esioh);

    /**
     * Hook permitting saving arbitrary information into statistical sample
     * files during \ref save_statistics as well as into restart files during
     * \ref save_restart.
     *
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass implementation.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     */
    virtual bool save_statistics_hook(
            const esio_handle esioh);

    /**
     * Hook permitting loading arbitrary metadata from all restart and
     * statistics files during \ref load_metadata.  In particular, changing the
     * number of scalar state fields in #fields should be done here when
     * necessary.
     *
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass implementation.
     *
     * @param esioh An ESIO handle pointing to an open, readable file.
     */
    virtual void load_metadata_hook(
            const esio_handle esioh);

    /**
     * Hook permitting loading information from restart files into
     * #state_nonlinear during \ref load_restart.  Subclasses should override
     * this method with the desired functionality.  Invoking the superclass
     * method in the override is optional.
     *
     * The default implementation into #state_nonlinear from the provided ESIO
     * handle file per #fields.  It can recover data stored using the \ref
     * save_state_hook implementation.
     *
     * @param esioh An ESIO handle pointing to an open, readable file.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual void load_state_hook(
            const esio_handle esioh);

    /**
     * Hook permitting loading arbitrary information from statistical sample
     * files during \ref load_statistics as well as from restart files during
     * \ref load_restart.
     *
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass implementation.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     */
    virtual void load_statistics_hook(
            const esio_handle esioh);

    /**
     * Compute default status output intervals.
     * Useful if the desired intervals depend on other scenario parameters.
     * Called by \ref prepare_controller(time_type,real_t).
     */
    virtual void default_status_interval(time_type&, step_type&)
    { /* NOP */ }

    /**
     * Compute default restart writing intervals.
     * \copydoc default_status_interval(time_type&,step_type&)
     */
    virtual void default_restart_interval(time_type&, step_type&)
    { /* NOP */ }

    /**
     * Compute default statistics writing intervals.
     * \copydoc default_status_interval(time_type&,step_type&)
     */
    virtual void default_statistics_interval(time_type&, step_type&)
    { /* NOP */ }

    /** Permits a subclass-specific default restart writing interval. */
    virtual time_type default_restart_interval()
    { return controller->forever_t(); }

    /** Permits a subclass-specific default statistics writing interval. */
    virtual time_type default_statistics_interval()
    { return controller->forever_t(); }

    /**
     * Did the previous time advance end in a predicted, controlled manner?
     */
    bool soft_teardown;

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
    double wtime_load_restart;

    /** Wall time at which we began time stepping */
    double wtime_advance_start;

    /** Tracks last time a status line was output */
    step_type last_status_nt;

    /** Tracks last time a restart file was written successfully */
    step_type last_restart_saved_nt;

    /** Tracks last time a statistics sample file was written successfully */
    step_type last_statistics_saved_nt;

private:

    /** Processes incoming signal actions via #controller callbacks. */
    bool process_any_signals_received(
            const time_type t,
            const step_type nt);

    /**
     * Maintains if any signals were observed in \ref signal::global_received
     * during the last time step during a time advance.  To be indexed using
     * \ref signal::action_type.
     */
    signal::received_type signal_received;

    /** Was a metadata file ever saved to disk? */
    bool metadata_saved;

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

    /**
     * Construct an instance querying and mutating external state.
     * Intended to be used in conjunction with a \ref driver_base instance.
     */
    delta_t_allreducer(
            const double& wtime_mpi_init,
            const double& wtime_fftw_planning,
            const shared_ptr<time_definition>& timedef,
            const double& wtime_load_restart,
            const double& wtime_advance_start,
            const driver_base::step_type& last_status_nt,
            const driver_base::step_type& last_restart_saved_nt,
            driver_base::delta_t_ratios_type& delta_t_ratios,
            signal::received_type& signal_received);

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
    const double& wtime_load_restart;

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
