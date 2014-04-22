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

#ifndef SUZERAIN_SUPPORT_DRIVER_BASE_HPP
#define SUZERAIN_SUPPORT_DRIVER_BASE_HPP

/** @file
 * Suzerain-based driver logic reusable across multiple applications.
 */

#include <suzerain/common.hpp>
#include <suzerain/extrema.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/support/application_base.hpp>
#include <suzerain/support/definition_restart.hpp>
#include <suzerain/support/definition_signal.hpp>
#include <suzerain/support/definition_statistics.hpp>
#include <suzerain/support/definition_time.hpp>
#include <suzerain/support/field.hpp>
#include <suzerain/timecontroller.hpp>

// Forward declarations
struct suzerain_bl_local;
struct suzerain_bl_pg;
struct suzerain_bl_qoi;
struct suzerain_bl_reynolds;
struct suzerain_bl_thicknesses;
struct suzerain_bl_viscous;
struct suzerain_channel_local;
struct suzerain_channel_qoi;
struct suzerain_channel_viscous;

namespace suzerain {

// Forward declarations
class operator_tools;
class summary;

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

    /** Tear down the simulation (reactively due to an incoming signal). */
    teardown_reactive,

    /** Tear down the simulation (proactively due to --advance_wt limit). */
    teardown_proactive,

    /** Halt the simulation (reactively due to an incoming signal). */
    halt_reactive,

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

/**
 * A driver base class for managing a Suzerain application.  Intended for
 * time-varying, three-dimensional problems.  Provides many hooks to permit
 * lifecycle-related behavior for long-running, parallel applications.
 * Instantiate from within \c main().
 *
 * Two separate types of file input/output are accommodated:
 * <ol>
 * <li>Restart files which contain metadata, state, and statistics.</li>
 * <li>Statistical sample files which contain metadata and statistics.</li>
 * </ol>
 * Though computing statistics and writing them via \ref save_statistics may
 * require an active parallel decomposition (i.e. \ref establish_decomposition
 * to have been invoked), loading statistics via \ref load_statistics should
 * not require a viable parallel pencil grid.
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
                const std::string &argument_synopsis = "",
                const std::string &description = "",
                const std::string &revstr = "");

    /** Virtual destructor as appropriate for abstract base class */
    virtual ~driver_base();

    /** @copydoc application_base::initialize */
    virtual std::vector<std::string> initialize(
            int argc,
            char **argv);

    /**
     * The log4cxx configuration to use.  Files <tt>bulk.dat</tt>,
     * <tt>L2.dat</tt>, and <tt>rms.dat</tt> collecting messages with the names
     * <tt>bc</tt>, <tt>bulk</tt>, <tt>L2</tt>, and <tt>rms</tt> have been
     * added.
     */
    virtual std::string log4cxx_config();

    /**
     * Information about the state fields.  Should be populated by subclasses
     * prior to \ref save_restart or \ref load_restart.
     */
    std::vector<support::field> fields;

    /** Operational details on saving restart files. */
    shared_ptr<definition_restart> restartdef;

    /** Operational details on saving statistics files. */
    shared_ptr<definition_statistics> statsdef;

    /** Operational details on advancing simulation time. */
    shared_ptr<definition_time> timedef;

    /**
     * Controls the OS signals triggering various types of processing.
     * Must be static so that it can be queried within signal handler.
     */
    static definition_signal signaldef;

    /** Controls low storage method to be used for time advance. */
    shared_ptr<lowstorage::method_interface<
            complex_t
        > > method;

    /**
     * Refers to a linear operator instance able to interoperate between
     * #state_linear and #state_nonlinear.  This is the interface to which
     * linear operators should be coded.
     */
    shared_ptr<lowstorage::linear_operator<
                state_common_type, state_nonlinear_type
            > > L;

    /**
     * Refers to a nonlinear operator instance able to operate on
     * #state_nonlinear.  This is the interface to which nonlinear operators
     * should be coded.
     */
    shared_ptr<lowstorage::operator_nonlinear<
                state_nonlinear_type
            > > N;

    /** Controls time advance, including callback processing. */
    shared_ptr<timecontroller<real_t> > controller;

    /** @copydoc application_base::reset */
    virtual void reset();

    /**
     * When possible, any operator_tools superclass of N is reused so that
     * callers may benefit from any cached factorizations.
     */
    virtual shared_ptr<operator_tools> obtain_operator_tools();

    /**
     * Ensure #method is valid for use by #controller.  That is, if
     * <tt>!method</tt> then \ref lowstorage::smr91 is used per
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
     *                  Often, zero is a perfectly find choice.
     * @param chi       Time-independent scaling factor for #N on each substep.
     *                  Often, \ref pencil_grid::chi() is used to hide FFT
     *                  normalization costs.
     */
    virtual void prepare_controller(
            const time_type initial_t,
            const real_t chi);

    /**
     * Use #controller to advance the simulation per #timedef.
     *
     * Notice <tt>controller</tt> advance routines could be invoked directly,
     * but they provide none of this additional postprocessing functionality.
     *
     * @param final_status     Invoke \ref log_status after advance completes?
     * @param final_statistics Invoke \ref save_statistics after successful
     *                         advance?  Default is \c false as any
     *                         \c final_restart will contain statistics too.
     * @param final_restart    Invoke \ref save_restart after successful
     *                         advance?
     * @param output_timers    Output performance details after advance
     *                         completes?
     * @param output_stepping  Output time step size metrics after advance
     *                         completes?
     *
     * @return The wall time spent advancing simulation time on success.
     *         The negated wall time spend advancing the simulation on failure.
     */
    virtual double advance_controller(
            const bool final_status     = true,
            const bool final_statistics = false,
            const bool final_restart    = true,
            const bool output_timers    = true,
            const bool output_stepping  = true);

    /**
     * Build a fixed-width, human-friendly way to output the given simulation
     * time and time step number.  Care is taken to ensure sequential outputs
     * are minimally distinct.
     *
     * Subclasses overriding this method should likely also override
     * \ref build_timeprefix_description.
     *
     * @param t  Simulation time to output
     * @param nt Simulation time step to output
     *
     * @return a string suitable for use in output status messages.
     */
    virtual std::string build_timeprefix(
            const time_type t,
            const step_type nt) const;

    /**
     * Build a fixed-width label matching the output of \ref build_timeprefix.
     * Additional left padding is required to match a given \c timeprefix
     * string.
     *
     * @return a string suitable for use in output status messages.
     */
    virtual std::string build_timeprefix_description(
            const char * describe_t  = "t",
            const char * describe_nt = "nt") const;

    /**
     * Compute digits required in mantissa to ensure time output is monotonic.
     * This process ensures multiple status lines are minimally distinct and
     * requires interrogating several details within #timedef.
     *
     * @return the maximum of <code>max(0, -floor(log10(dt)))</code>
     * across several pathological values of \c dt implied by #timedef.
     */
    virtual int build_timeprefix_mantissa_digits() const;

    /**
     * Routine to output status invoking \ref log_status_hook,
     * generally called via the timecontroller.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool log_status(
            const time_type t,
            const step_type nt);

    /**
     * Log messages containing mean L2 and RMS fluctuation information.
     *
     * RMS fluctuations are only logged when the grid can support them.
     */
    virtual void log_state_L2(
            const std::string& timeprefix,
            const char * const name_L2  = "state.L2",
            const char * const name_RMS = "state.RMS");

    /** Log messages containing bulk quantities. */
    virtual void log_state_bulk(
            const std::string& timeprefix,
            const char * const name_bulk = "state.bulk");

    /**
     * Log messages containing global minimum and maximum state values.
     *
     * Messages are logged only when the grid can support fluctuations.
     */
    virtual void log_state_extrema(
            const std::string& timeprefix,
            const char * const name_min  = "state.min",
            const char * const name_xmin = "state.xmin",
            const char * const name_ymin = "state.ymin",
            const char * const name_zmin = "state.zmin",
            const char * const name_max  = "state.max",
            const char * const name_xmax = "state.xmax",
            const char * const name_ymax = "state.ymax",
            const char * const name_zmax = "state.zmax",
            const char * const name_fneg = "state.fneg");

    /**
     * Log messages containing state at the upper and lower boundaries.
     */
    virtual void log_boundary_conditions(
            const std::string& timeprefix);

    /**
     * Log messages containing a wide variety of quantities of interest for a
     * boundary layer simulation.  Data to be logged must be provided in \c
     * viscous, \c thick, \c reynolds, \c qoi, and \c pg.
     *
     * This method produces well-formatted results given relevant computed
     * quantities.  Subclasses will likely provide an overload that gathers
     * these details in a simulation-dependent manner and then invoke this
     * method within a \ref log_status_hook() or \ref log_statistics_hook()
     * overload.
     *
     * @param wall        Local information from the wall.
     * @param viscous     Possibly computed by suzerain_bl_compute_viscous().
     * @param thick       Possibly computed by suzerain_bl_compute_thicknesses().
     * @param edge        Local state information from \f$y = \delta\f$.
     * @param edge99      Local state information from \f$y = \delta_{99}\f$.
     * @param reynolds    Possibly computed by suzerain_bl_compute_reynolds().
     * @param qoi         Possibly computed by suzerain_bl_compute_qoi().
     * @param pg          Possibly computed by suzerain_bl_compute_pg().
     *                    If \c NULL, no pressure gradient-related information
     *                    will be logged.
     * @param name_wall   Logging name for wall-related information.
     *                    If \c NULL, this message will not be logged.
     * @param name_visc   Logging name for viscous-related quantities.
     *                    If \c NULL, this message will not be logged.
     * @param name_thick  Logging name for boundary layer thickness details.
     *                    If \c NULL, this message will not be logged.
     * @param name_edge   Logging name for edge-related information.
     *                    If \c NULL, this message will not be logged.
     * @param name_edge99 Logging name for edge99-related information.
     *                    If \c NULL, this message will not be logged.
     * @param name_Re     Logging name for boundary layer Reynolds numbers.
     *                    If \c NULL, this message will not be logged.
     * @param name_qoi    Logging name for general quantities of interest.
     *                    If \c NULL, this message will not be logged.
     * @param name_pg     Logging name for pressure gradient-related information.
     *                    If \c NULL, this message will not be logged.
     */
    virtual void log_quantities_boundary_layer(
            const std::string& timeprefix,
            const suzerain_bl_local        * const wall,
            const suzerain_bl_viscous      * const viscous,
            const suzerain_bl_thicknesses  * const thick,
            const suzerain_bl_local        * const edge,
            const suzerain_bl_local        * const edge99,
            const suzerain_bl_reynolds     * const reynolds,
            const suzerain_bl_qoi          * const qoi,
            const suzerain_bl_pg           * const pg,
            const char * const name_wall   =  "bl.wall",
            const char * const name_visc   =  "bl.visc",
            const char * const name_thick  =  "bl.thick",
            const char * const name_edge   =  "bl.edge",
            const char * const name_edge99 =  "bl.edge99",
            const char * const name_Re     =  "bl.Re",
            const char * const name_qoi    =  "bl.qoi",
            const char * const name_pg     =  "bl.pg");

    /**
     * Log messages containing a wide variety of quantities of interest for a
     * channel flow simulation.  Data to be logged must be provided in \c wall,
     * \c viscous, \c center, and \c qoi.
     *
     * This method produces well-formatted results given relevant computed
     * quantities.  Subclasses will likely provide an overload that gathers
     * these details in a simulation-dependent manner and then invoke this
     * method within a \ref log_status_hook() or \ref log_statistics_hook()
     * overload.
     *
     * @param wall        Local information from the wall.
     * @param viscous     Possibly from suzerain_channel_compute_viscous().
     * @param center      Local information from the channel centerline.
     * @param qoi         Possibly from suzerain_channel_compute_qoi().
     * @param name_wall   Logging name for wall-related quantities.
     *                    If \c NULL, this message will not be logged.
     * @param name_visc   Logging name for viscous-related quantities.
     *                    If \c NULL, this message will not be logged.
     * @param name_center Logging name for centerline-related quantities.
     *                    If \c NULL, this message will not be logged.
     * @param name_qoi    Logging name for general quantities of interest.
     *                    If \c NULL, this message will not be logged.
     */
    virtual void log_quantities_channel(
            const std::string& timeprefix,
            const suzerain_channel_local   * const wall,
            const suzerain_channel_viscous * const viscous,
            const suzerain_channel_local   * const center,
            const suzerain_channel_qoi     * const qoi,
            const char * const name_wall   = "chan.wall",
            const char * const name_visc   = "chan.visc",
            const char * const name_center = "chan.center",
            const char * const name_qoi    = "chan.qoi");

    /**
     * Save time-independent metadata that must appear in all restart and
     * statistics files.  Though this logic will be invoked automatically the
     * first time \ref save_restart() or \ref save_statistics() is used, users
     * may wish to invoke the method explicitly during general initialization
     * prior to \ref advance_controller().
     *
     * @see Member #restartdef to control the metadata file location.
     */
    void save_metadata();

    /**
     * Save time-independent metadata that should appear in all restart and
     * statistics files.
     *
     * Subclasses should not override this method but should instead use \ref
     * save_metadata_hook.
     */
    virtual void save_metadata(
            const esio_handle esioh);

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
     * Save state \e and statistics into an automatically managed sequence of
     * restart files sharing common metadata.  Subclasses should not override
     * this method but should instead use \ref save_state_hook and \ref
     * save_statistics_hook.
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
     *      including the name of the file on disk and how many files
     *      in the sequence should be maintained.
     */
    virtual bool save_restart(
            const time_type t,
            const step_type nt);

    /**
     * Save state \e and statistics into a one-off restart file with one-off
     * metadata.  Subclasses should not override this method but should instead
     * use \ref save_state_hook and \ref save_statistics_hook.
     *
     * @param t         The simulation time to be stored in the restart file.
     * @param dstfile   The path to which the restart file will be written.
     * @param overwrite If \c true, clobber any existing file.
     *                  If \c false, handled per ESIO's \c esio_file_clone.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     *
     * @see Member #restartdef to set the location of temporary metadata files.
     */
    virtual bool save_restart(
            const time_type t,
            const std::string dstfile,
            const bool overwrite);

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
     * Load statistics from a statistical sampling file.  Subclasses should not
     * override this method but should instead use \ref load_statistics_hook.
     *
     * Notice that because \ref save_restart invokes \ref save_statistics_hook
     * and because both restart and statistics files share common metadata,
     * this method should also be invokable on an open restart file.
     *
     * @param[in]  esioh An ESIO handle pointing to an open statistics file.
     * @param[out] t     The simulation time stored in the statistics file.
     */
    virtual void load_statistics(
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
     * Hook permitting subclasses to output status information which is
     * invoked at the end of \ref log_status.
     *
     * Invokes \ref log_state_bulk, \ref log_state_L2, \ref log_state_extrema,
     * \ref log_boundary_conditions, and \ref log_status_hook.
     *
     * Subclasses should override this method adding any desired functionality
     * either before, after, or in lieu of invoking the superclass
     * implementation.
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
     * Hook permitting saving arbitrary information from \ref #state_nonlinear
     * into restart files during \ref save_restart.  Subclasses should override
     * this method with the desired functionality.  Invoking the superclass
     * method in the override is optional.
     *
     * The default implementation saves the contents of #state_nonlinear into
     * the provided ESIO handle file destroying #state_nonlinear in the
     * process.  When \c restartdef->physical is \c false, expansion
     * coefficients are written in all three directions using \ref
     * store_coefficients.  When it is \c true, values at collocation
     * points are written using store_collocation_Values.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool save_state_hook(
            const esio_handle esioh);

    /**
     * Hook permitting loading information from restart files into
     * #state_nonlinear during \ref load_restart.  Subclasses should override
     * this method with the desired functionality.  Invoking the superclass
     * method in the override is optional.
     *
     * The default implementation loads into #state_nonlinear from the provided
     * ESIO handle file per #fields.  It can recover data stored using the \ref
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
     * Hook permitting saving arbitrary information into statistical sample
     * files during \ref save_statistics as well as into restart files during
     * \ref save_restart.
     *
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass implementation.
     *
     * @param esioh An ESIO handle pointing to an open, writeable file.
     * @param t     The simulation time at which statistics are being requested.
     *              The time should not be saved, as this is the role of
     *              \ref save_statistics, but it may be used for logging
     *              purposes or to avoid re-computing statistics multiple
     *              times at one simulation time when #controller is used.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual bool save_statistics_hook(
            const esio_handle esioh,
            const time_type t);

    /**
     * Hook permitting loading arbitrary information from statistical sample
     * files during \ref load_statistics as well as from restart files during
     * \ref load_restart.
     *
     * Subclasses should override this method adding any desired functionality
     * either before or after invoking the superclass implementation.
     *
     * @param esioh An ESIO handle pointing to an open, writable file.
     *
     * @returns True if any active time advance should continue.
     *          False otherwise.
     */
    virtual void load_statistics_hook(
            const esio_handle esioh);

    /**
     * Compute default restart writing intervals.
     * Useful if the desired intervals depend on other scenario parameters.
     * Called by \ref prepare_controller(time_type,real_t).
     *
     * @param[out] dt Time between restarts measured in simulation time.
     * @param[out] nt Time between restarts measured in time steps.
     */
    virtual void default_restart_interval(time_type& dt, step_type& nt)
    {
        // No default behavior because neither dt nor nt is modified
        (void) dt;
        (void) nt;
    }

    /**
     * Compute default statistics writing intervals.
     *
     * By default, this is one-fourth the frequency of any non-zero, non-NaN
     * value returned by \ref default_restart_interval.
     *
     * \copydetails default_restart_interval(time_type&,step_type&)
     */
    virtual void default_statistics_interval(time_type& dt, step_type& nt);

    /**
     * Compute default status output intervals.
     *
     * By default, this is one-sixteenth the frequency of any non-zero, non-NaN
     * value returned by \ref default_statistics_interval.
     *
     * \copydetails default_restart_interval(time_type&,step_type&)
     */
    virtual void default_status_interval(time_type& dt, step_type& nt);

    /**
     * Atop commonly-defined infrastructure, a generic data summarization
     * application can be defined.  The \ref summary_run routine is designed by
     * to invoked on a subclass within \c main after the strings below are
     * passed to the constructor.
     * @{
     */

    /**
     * Provide to constructor as \c description when using \ref summary_run.
     */
    static const char * const summary_description;

    /**
     * Provide to constructor as \c argument_synopsis when using \ref
     * summary_run.
     */
    static const char * const summary_argument_synopsis;

    /**
     * A resource-owning mapping from simulation time to \ref summary results.
     * In-order traversal gives temporal evolution of summary statistics.
     */
    typedef boost::ptr_map<real_t, summary> summary_pool_type;

    /**
     * Run the generic data summarization application.
     * Do not invoke \ref initialize beforehand as this routine does so.
     *
     * On successful return, the instance will be populated with metadata
     * from either the target file specified on the command line or the
     * final restart file listed in \c argv.
     *
     * \param[in]  argc  Incoming arguments per <code>main(argc, ...)</code>
     * \param[in]  argv  Incoming arguments per <code>main(..., argv)</code>
     * \param[out] pool  Processed data summaries as a function of unique time.
     * \param[out] final When \c pool has non-zero size,
     *                   either mean or final profile information
     *                   on collocation points suitable for further analysis.
     *
     * \return Success or failure suitable for return to the OS.
     */
    int
    summary_run(int argc, char **argv, summary_pool_type& pool, summary& final);

    /**
     * @}
     */

    /**
     * Did the previous time advance end because of tear down receipt?
     */
    bool received_teardown;

    /**
     * Did the previous time advance end because of halt receipt?  This is like
     * #received_teardown except #advance_controller considers it a failure.
     */
    bool received_halt;

    /**
     * Tracking of whether headers have been shown for particular log messages.
     * Keys are logging-related header names, e.g. "state.RMS", with the value
     * indicating if the header has been output.
     */
    std::map<std::string,bool> header_shown;

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

    // FIXME: Carry extrema time? #3071
    /**
     * Maintains instantaneously sampled wall-normal extrema quantities.  
     * Member \c extrema.t tracks the last time statistics were computed 
     * and is used as a mechanism to avoid expensive recomputations.  
     */
    std::vector<field_extrema_xz>  extrema;

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

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

/**
 * A stateful functor that performs a MPI_Allreduce to determine the minimum
 * stable time step size across all ranks.  The same MPI Allreduce is used to
 * hide the cost of querying \ref signal::global_received across all ranks.
 */
class delta_t_allreducer : public lowstorage::delta_t_reducer
{
    /** Provides simple access to the superclass type */
    typedef lowstorage::delta_t_reducer super;

public:

    /**
     * Construct an instance querying and mutating external state.
     * Intended to be used in conjunction with a \ref driver_base instance.
     */
    delta_t_allreducer(
            const double& wtime_mpi_init,
            const double& wtime_fftw_planning,
            const shared_ptr<definition_time>& timedef,
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
    const shared_ptr<definition_time>& timedef;

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

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

} // end namespace support

} // end namespace suzerain

#endif // SUZERAIN_SUPPORT_DRIVER_BASE_HPP
