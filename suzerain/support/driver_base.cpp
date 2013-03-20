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
 * @copydoc driver_base.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/support/driver_base.hpp>

#include <suzerain/countof.h>
#include <suzerain/error.h>
#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/os.h>
#include <suzerain/state.hpp>
#include <suzerain/support/field.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/version.hpp>

// We want plain ol' signal(2) below but suzerain::support::signal interferes
// Prepare alias for the name via a function pointer so we may invoke it below
static sighandler_t (* const signal2)(int, sighandler_t) = &signal;

namespace suzerain {

namespace support {

namespace signal {

// Initialized to zero indicating no signals have been received
volatile_received_type global_received = {{/*0*/}};

} // end namespace signal

extern "C" { // So we may pass the address of process_signal(...) to signal(2)

static void driver_base_process_signal(const int sig)
{
    // Strictly speaking this handler performs too much work.  The design
    // choice was to have this extra work done on the (rare) signal receipt
    // rather than on the (frequent) polling of signal receipt status.

    using std::find;
    std::vector<int>::iterator end;

    // Determine if we should log status due to the signal
    end = driver_base::signaldef.status.end();
    if (find(driver_base::signaldef.status.begin(), end, sig) != end) {
        signal::global_received[signal::log_status] = sig;
    }

    // Determine if we should write a restart due to the signal
    end = driver_base::signaldef.restart.end();
    if (find(driver_base::signaldef.restart.begin(), end, sig) != end) {
        signal::global_received[signal::write_restart] = sig;
    }

    // Determine if we should tear down the simulation due to the signal
    end = driver_base::signaldef.teardown.end();
    if (find(driver_base::signaldef.teardown.begin(), end, sig) != end) {
        signal::global_received[signal::teardown_reactive] = sig;
    }

    // signal::global_received[signal::teardown_proactive] "detected"
    // within delta_t_allreducer::operator()(...).

    // Determine if we should compute and write statistics due to the signal
    end = driver_base::signaldef.statistics.end();
    if (find(driver_base::signaldef.statistics.begin(), end, sig) != end) {
        signal::global_received[signal::write_statistics] = sig;
    }
}

} // end extern "C"

// Static member initialization
signal_definition driver_base::signaldef;

driver_base::driver_base(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &revstr)
    : application_base(application_synopsis,
                       argument_synopsis,
                       description,
                       revstr)
    , restartdef(make_shared<restart_definition>(
                /* metadata    */ "metadata.h5.XXXXXX",
                /* uncommitted */ "uncommitted.h5.XXXXXX",
                /* destination */ "restart#.h5",
                /* retain      */ 1,
                /* dt          */ 0,
                /* nt          */ 0))
    , statsdef(make_shared<statistics_definition>(
                /* destination */ "sample#.h5"))
    , timedef(make_shared<time_definition>(
                /* advance_dt  */ 0,
                /* advance_nt  */ 0,
                /* advance_wt  */ 0,
                /* status_dt   */ 0,
                /* status_nt   */ 0,
                /* min_dt      */ 1e-8,
                /* max_dt      */ 1))
    , method()
    , L()
    , N()
    , controller()
    , soft_teardown(true)
    , log_status_L2_header_shown(false)
    , log_status_bulk_header_shown(false)
    , wtime_load_restart(std::numeric_limits<double>::quiet_NaN())
    , wtime_advance_start(std::numeric_limits<double>::quiet_NaN())
    , last_status_nt(std::numeric_limits<step_type>::max())
    , last_restart_saved_nt(std::numeric_limits<step_type>::max())
    , last_statistics_saved_nt(std::numeric_limits<step_type>::max())
    , metadata_saved(false)
    , allreducer(new delta_t_allreducer(this->wtime_mpi_init,
                                        this->wtime_fftw_planning,
                                        this->timedef,
                                        this->wtime_load_restart,
                                        this->wtime_advance_start,
                                        this->last_status_nt,
                                        this->last_restart_saved_nt,
                                        this->delta_t_ratios,
                                        this->signal_received))
    , who("driver")
{
    std::fill(signal_received.begin(), signal_received.end(), 0);
}

std::string
driver_base::log4cxx_config()
{
    return super::log4cxx_config() + // Append to the default configuration
        "\n"
        "## Collect \"bulk\" messages into bulk.dat mimicking LOG file behavior\n"
        "log4j.logger.bulk=INHERITED, BULK\n"
        "log4j.appender.BULK=${log4j.appender.LOG}\n"
        "log4j.appender.BULK.filename=bulk.dat\n"
        "log4j.appender.BULK.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.BULK.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.BULK.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
        "\n"
        "## Collect \"L2.mean\" messages into L2.mean.dat mimicking LOG file behavior\n"
        "log4j.logger.L2.mean=INHERITED, L2MEAN\n"
        "log4j.appender.L2MEAN=${log4j.appender.LOG}\n"
        "log4j.appender.L2MEAN.filename=L2.mean.dat\n"
        "log4j.appender.L2MEAN.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.L2MEAN.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.L2MEAN.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
        "\n"
        "## Collect \"rms.fluct\" messages into rms.fluct.dat mimicking LOG file behavior\n"
        "log4j.logger.rms.fluct=INHERITED, RMSFLUCT\n"
        "log4j.appender.RMSFLUCT=${log4j.appender.LOG}\n"
        "log4j.appender.RMSFLUCT.filename=rms.fluct.dat\n"
        "log4j.appender.RMSFLUCT.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.RMSFLUCT.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.RMSFLUCT.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
    ;
}

std::vector<std::string>
driver_base::initialize(int argc, char **argv)
{
    // Add problem definitions to options
    options.add_definition(*restartdef);
    options.add_definition(*statsdef  );
    options.add_definition(*timedef   );
    options.add_definition( signaldef );

    // Process incoming arguments by invoking superclass method
    std::vector<std::string> positional = super::initialize(argc, argv);

    // Generating unique file names as needed using mkstemp(3)
    // See Redmine ticket #2385 towards a nicer long-term solution
    {
        // Pack a temporary buffer with the three file name templates
        array<size_t,5> pos = {{ 0,
                                 restartdef->metadata.length()    + 1,
                                 restartdef->uncommitted.length() + 1,
                                 restartdef->destination.length() + 1,
                                 statsdef  ->destination.length() + 1 }};
        std::partial_sum(pos.begin(), pos.end(), pos.begin());
        scoped_array<char> buf(new char[pos[4]]);
        strcpy(&buf[pos[0]], restartdef->metadata.c_str());
        strcpy(&buf[pos[1]], restartdef->uncommitted.c_str());
        strcpy(&buf[pos[2]], restartdef->destination.c_str());
        strcpy(&buf[pos[3]], statsdef  ->destination.c_str());

        // Generate unique files to be overwritten and/or just file names.
        // File generation relies on template semantics of mkstemp(3).
        // Error checking kept minimal as failures here should not be fatal.
        // This is not particularly robust but it should serve our needs.
        if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0) {
            if (boost::ends_with(restartdef->metadata, "XXXXXX")) {
                close(mkstemp(&buf[pos[0]]));  // Clobbered later...
            }
            if (boost::ends_with(restartdef->uncommitted, "XXXXXX")) {
                close(mkstemp(&buf[pos[1]]));  // Possibly clobbered later...
                unlink(&buf[pos[1]]);          // ...so remove any evidence
            }
            if (boost::ends_with(restartdef->destination, "XXXXXX")) {
                close(mkstemp(&buf[pos[2]]));  // Not clobbered later...
                unlink(&buf[pos[2]]);          // ...so remove any evidence
            }
            if (boost::ends_with(statsdef->destination, "XXXXXX")) {
                close(mkstemp(&buf[pos[3]]));  // Not clobbered later...
                unlink(&buf[pos[3]]);          // ...so remove any evidence
            }
        }

        // Broadcast any generated names to all ranks and unpack values
        SUZERAIN_MPICHKR(MPI_Bcast(buf.get(), pos[4],
                         mpi::datatype<char>(), 0, MPI_COMM_WORLD));
        restartdef->metadata    = &buf[pos[0]];
        restartdef->uncommitted = &buf[pos[1]];
        restartdef->destination = &buf[pos[2]];
        statsdef  ->destination = &buf[pos[3]];
    }

    return positional;
}

driver_base::~driver_base()
{
    // Attempt to remove any lingering metadata file
    if (metadata_saved) {
        if (mpi::comm_rank(MPI_COMM_WORLD) == 0) {
            if (0 == unlink(restartdef->metadata.c_str())) {
                DEBUG("Cleaned up temporary file "
                      << restartdef->metadata);
            } else {
                WARN("Error cleaning up temporary file "
                     << restartdef->metadata);
            }
        }
    }

    // Preserve restartdef->uncommitted as it may help post mortem debugging as
    // during operation committing a file removes the uncommitted temporary
}

void
driver_base::prepare_method()
{
    if (!method) {
        INFO0(who, "Preparing SMR91 timestepping scheme using evmagfactor "
              << timedef->evmagfactor);
        this->method.reset(new timestepper::lowstorage::method<
                    timestepper::lowstorage::smr91,
                    state_common_type::element
                >(timedef->evmagfactor));
    }
}

void
driver_base::prepare_controller(
        const driver_base::time_type initial_t,
        const real_t chi)
{
    using std::fmod;

    INFO0(who, "Preparing controller at initial_t " << initial_t
          << " for step size in [" << timedef->min_dt
          << ", " << timedef->max_dt << "]");

    // Instantiate the low-storage timecontroller
    if (!method) prepare_method();    // Only call if necessary
    SUZERAIN_ENSURE(method);          // Defend against segfaults...
    SUZERAIN_ENSURE(L);               // ...if everything not initialized
    SUZERAIN_ENSURE(N);
    SUZERAIN_ENSURE(state_linear);
    SUZERAIN_ENSURE(state_nonlinear);
    SUZERAIN_ENSURE(restartdef);
    SUZERAIN_ENSURE(statsdef);
    SUZERAIN_ENSURE(timedef);
    controller.reset(make_lowstorage_timecontroller(
                *method, *allreducer, *L, chi, *N,
                *state_linear, *state_nonlinear,
                initial_t, timedef->min_dt, timedef->max_dt));

    // Register status callbacks status_{dt,nt} as requested.
    // However, the first callback occurs modulo the requested dt value.
    // to permit sampling rate continuity across shutdown/restart.
    // If no non-default, non-zero values were provided, permit override.
    if (    options.variables()["status_dt"].defaulted()
         && options.variables()["status_nt"].defaulted()
         && !timedef->status_dt
         && !timedef->status_nt) {
        default_status_interval(timedef->status_dt, timedef->status_nt);
    }
    controller->add_periodic_callback(
            (timedef->status_dt ?   timedef->status_dt
                                  - fmod(initial_t, timedef->status_dt)
                                : controller->forever_t()),
            (timedef->status_dt ? timedef->status_dt
                                : controller->forever_t()),
            (timedef->status_nt ? timedef->status_nt
                                : controller->forever_nt()),
            boost::bind(&driver_base::log_status, this, _1, _2));

    // Register restart-writing callbacks restart_{dt,nt} as requested.
    // Again, the first callback occurs modulo the requested dt value.
    // If no non-default, non-zero values were provided, permit override.
    if (    options.variables()["restart_dt"].defaulted()
         && options.variables()["restart_nt"].defaulted()
         && !restartdef->dt
         && !restartdef->nt) {
        default_restart_interval(restartdef->dt, restartdef->nt);
    }
    controller->add_periodic_callback(
            (restartdef->dt ? restartdef->dt - fmod(initial_t, restartdef->dt)
                            : controller->forever_t()),
            (restartdef->dt ? restartdef->dt : controller->forever_t()),
            (restartdef->nt ? restartdef->nt : controller->forever_nt()),
            boost::bind(&driver_base::save_restart, this, _1, _2));

    // Register statistics-related callbacks per statistics_{dt,nt}.
    // Again, the first callback occurs modulo the requested dt value.
    // If no non-default, non-zero values were provided, be sensible.
    if (   options.variables()["statistics_dt"].defaulted()
        && options.variables()["statistics_nt"].defaulted()
        && !statsdef->dt
        && !statsdef->nt) {
        default_statistics_interval(statsdef->dt, statsdef->nt);
    }
    controller->add_periodic_callback(
            (statsdef->dt ? statsdef->dt - fmod(initial_t, statsdef->dt)
                          : controller->forever_t()),
            (statsdef->dt ? statsdef->dt : controller->forever_t()),
            (statsdef->nt ? statsdef->nt : controller->forever_nt()),
            boost::bind(&driver_base::save_statistics, this, _1, _2));

    // Register any necessary signal handling logic once per unique signal
    {
        logging::logger_type log_signal = logging::get_logger("signal");

        // Obtain a set of signal numbers which we need to register
        std::vector<int> s;
        s.insert(s.end(), signaldef.status.begin(),
                          signaldef.status.end());
        s.insert(s.end(), signaldef.restart.begin(),
                          signaldef.restart.end());
        s.insert(s.end(), signaldef.statistics.begin(),
                          signaldef.statistics.end());
        s.insert(s.end(), signaldef.teardown.begin(),
                          signaldef.teardown.end());
        sort(s.begin(), s.end());
        s.erase(unique(s.begin(), s.end()), s.end());

        // Register the signal handler for each of these signals
        typedef std::vector<int>::const_iterator const_iterator;
        for (const_iterator i = s.begin(); i != s.end(); ++i) {
            const char * name = suzerain_signal_name(*i);
            if (SIG_ERR != signal2(*i, &driver_base_process_signal)) {
                if (name) {
                    DEBUG0(log_signal, "Registered handler for " << name);
                } else {
                    DEBUG0(log_signal, "Registered handler for " << *i);
                }
            } else {
                if (name) {
                    WARN0(log_signal,
                          "Unable to register handler for " << name);
                } else {
                    WARN0(log_signal,
                          "Unable to register handler for " << *i);
                }
            }
        }

        // Iff we registered any handlers, process signal receipt in stepper.
        // Notice signal receipt include --advance_wt calling us a pumpkin.
        // We can afford this every time step because of delta_t_allreducer.
        if (s.size() > 0 || timedef->advance_wt > 0) {
            controller->add_periodic_callback(
                    controller->forever_t(), 1,
                    boost::bind(&driver_base::process_any_signals_received,
                                this, _1, _2));
        }
    }
}

void
driver_base::prepare_controller(
        const driver_base::time_type initial_t)
{
    // Nonlinear scaling factor (N_x N_z)^(-1) from write up section 2.1
    SUZERAIN_ENSURE(grid);
    SUZERAIN_ENSURE(grid->dN.x());
    SUZERAIN_ENSURE(grid->dN.z());
    return prepare_controller(initial_t,
                              real_t(1) / (grid->dN.x() * grid->dN.z()));
}

double
driver_base::advance_controller(
            const bool final_status,
            const bool final_restart,
            const bool output_timers,
            const bool output_stepping)
{
    SUZERAIN_ENSURE(controller);

    // Advance time according to timedef->advance_dt, advance_nt criteria
    if (output_timers) {
#ifdef SUZERAIN_HAVE_GRVY
        grvy_log_setlevel(GRVY_ERROR);  // Suppress GRVY resolution warnings
        grvy_timer_reset();
#endif
    }
    bool success = true;
    const time_type t_initial  = controller->current_t();
    const step_type nt_initial = controller->current_nt();
    soft_teardown = false;
    wtime_advance_start = MPI_Wtime();
#pragma warning(push,disable:1572)
    switch ((!!timedef->advance_dt << 1) + !!timedef->advance_nt) {
#pragma warning(pop)
        case 3:
            INFO0(who, "Advancing simulation by at most "
                  << timedef->advance_dt << " units of physical time");
            INFO0(who, "Advancing simulation by at most "
                  << timedef->advance_nt << " discrete time steps");
            success = controller->advance(t_initial + timedef->advance_dt,
                                          timedef->advance_nt);
            break;
        case 2:
            INFO0(who, "Advancing simulation by at most "
                  << timedef->advance_dt << " units of physical time");
            success = controller->advance(t_initial + timedef->advance_dt);
            break;
        case 1:
            INFO0(who, "Advancing simulation by at most "
                  << timedef->advance_nt << " discrete time steps");
            success = controller->step(timedef->advance_nt);
            break;
        case 0:
            if (options.variables()["advance_nt"].defaulted()) {
                INFO0(who, "Advancing simulation until terminated by a signal");
                success = controller->advance();
            } else {
                INFO0(who, "Simulation will not be advanced in time");
            }
            break;
        default:
            FATAL0(who, "Sanity error in advance_controller");
            return EXIT_FAILURE;
    }
    const double wtime_advance_end = MPI_Wtime();
    const step_type nsteps = controller->current_nt() - nt_initial;
    if (output_timers) {
#ifdef SUZERAIN_HAVE_GRVY
        grvy_timer_finalize();
        grvy_log_setlevel(GRVY_INFO);   // Re-enable GRVY warnings
#endif
    }
    if (soft_teardown) {
        INFO0(who, "Time controller stopped in reaction to tear down signal");
        success = true; // ...treat like successful advance
    } else if (!success && controller->current_dt() < controller->min_dt()) {
        WARN0(who, "Time controller halted because step "
              << controller->current_dt() << " was smaller than min_dt "
              << controller->min_dt() );
    } else if (!success) {
        WARN0(who, "Time controller halted unexpectedly");
    }

    // Output status if it was not just output during time advancement
    if (   final_status
        && last_status_nt != controller->current_nt()) {
        log_status(controller->current_t(), controller->current_nt());
    }

    // Save a final restart if one was not just saved during time advancement
    if (   final_restart
        && success
        && last_restart_saved_nt != controller->current_nt()) {
        INFO0(who, "Saving restart file after time controller finished");
        save_restart(controller->current_t(), controller->current_nt());
    }

    // Postprocess GRVY timer information on time advance now that restart is
    // safely on disk.  Reduces likelihood that GRVY hiccups cause data loss.
    // Only summarize when time advance took long enough to be interesting
    const real_t wtime_advance = wtime_advance_end - wtime_advance_start;
    if (   output_timers
        && nsteps > 0
        && dgrid->has_zero_zero_modes()
        && wtime_advance > 5 /*seconds*/) {

#ifdef SUZERAIN_HAVE_GRVY
        static const char header[]
            = "GRVY timings from MPI rank with zero-zero modes:";

        // GRVY uses printf so futzing required to employ logging subsystem
        // Error handling is unsophisticated and ugly.  Quel dommage, but
        // some attempt is made to fail back to printf if tmpfile fails.

        int saved_stdout = -1;
        FILE * const tmp = tmpfile();
        if (!tmp) {
            WARN("Could not open temporary file to read "
                 << header << " " << strerror(errno));
            INFO(header);
        } else {
            fflush(stdout);
            saved_stdout = dup(STDOUT_FILENO);
            dup2(fileno(tmp), STDOUT_FILENO);
            puts(header);
        }

        grvy_timer_summarize();

        if (tmp) {
            fflush(stdout);
            fflush(tmp);
            const long len = ftell(tmp);
            char * const buf = (char *) calloc(len + 1, 1);
            rewind(tmp);
            if (!buf) {
                WARN("Could not allocate buffer to read "
                     << header << " " << strerror(errno));
            } else {
                fread(buf, sizeof(buf[0]), len, tmp);
            }
            dup2(saved_stdout, STDOUT_FILENO);
            close(saved_stdout);
            if (buf) INFO(buf);
            free(buf);
            fclose(tmp);
        }
#endif

        // Admit what overhead we're neglecting in the following calculations
        INFO0(who, "Advancement rate calculations ignore "
              << MPI_Wtime() - wtime_mpi_init - wtime_advance
              << " seconds of fixed overhead");

        // Advance rate measured in a (mostly) problem-size-agnostic metric
        INFO0(who, "Advancing at "
              << wtime_advance / nsteps / grid->N.prod()
              << " wall seconds per time step per grid point");

        // Advance rate measured in a problem-size-dependent metric
        INFO0(who, "Advancing at " << wtime_advance / nsteps
              << " wall seconds per time step");

        // Advance rate measured in nondimensional simulation time units
        INFO0(who, "Advancing at "
              << wtime_advance / (controller->current_t() - t_initial)
              << " wall seconds per simulation time unit");
    }

    // Output details on time advancement (whenever advancement occurred)
    if (output_stepping && nsteps > 0) {
        using std::numeric_limits;
        std::ostringstream msg;
        msg.precision(static_cast<int>(numeric_limits<real_t>::digits10*0.75));
        msg << "Advanced simulation from t_initial = " << t_initial
            << " to t_final = " << controller->current_t()
            << " in " << nsteps << " steps";
        INFO0(who, msg.str());
        msg.str("");
        msg.precision(static_cast<int>(numeric_limits<real_t>::digits10*0.50));
        msg << "Min/avg/max/std of delta_t: "
            << controller->taken_min()  << ", "
            << controller->taken_mean() << ", "
            << controller->taken_max()  << ", "
            << controller->taken_stddev();
        INFO0(who, msg.str());
        msg.str("");
        msg << "Mean delta_t criteria versus minimum criterion: ";
        const size_t n = delta_t_ratios.size();
        for (size_t i = 0; i < n; ++i) {
            namespace acc = boost::accumulators;
            msg << acc::mean(delta_t_ratios[i]);
            if (i < n-1) msg << ", ";
        }
        INFO0(who, msg.str());
    }

    return success ? wtime_advance : -wtime_advance;
}

/**
 * Width constant used in build_timeprefix and friends for \c nt.
 * If more than ten million steps have elapsed, bully for you.
 */
static const std::streamsize timeprefix_setw_nt = 7;

std::string
driver_base::build_timeprefix(
        const driver_base::time_type t,
        const driver_base::step_type nt)
{
    // Precision computations ensure multiple status lines minimally distinct
    const std::streamsize np = build_timeprefix_mantissa_digits();

    // Build string using the computed precision information
    std::ostringstream oss;
    if (np > 0) {
        oss.setf(std::ios::fixed, std::ios::floatfield);
        const std::streamsize oldprec = oss.precision(np);
        oss << t;
        oss.precision(oldprec);
        oss.unsetf(std::ios::fixed);
    } else {
        oss << t;
    }
    oss << ' ' << std::setw(timeprefix_setw_nt) << nt;

    return oss.str();
}

std::string
driver_base::build_timeprefix_description(
            const char * describe_t,
            const char * describe_nt)
{
    // Build string with inter-label spacing matching build_timeprefix output
    std::ostringstream oss;
    oss << describe_t
        << ' '
        << std::setw(timeprefix_setw_nt) << describe_nt;
    return oss.str();
}

int
driver_base::build_timeprefix_mantissa_digits()
{
    using std::max;
    using std::floor;
    using std::log10;

    real_t n = 0;

    // Restart, statistics, and status frequencies measured in simulation time
    if (restartdef->dt > 0) {
        n = max(n, -floor(log10(restartdef->dt)));
    }
    if (statsdef->dt > 0) {
        n = max(n, -floor(log10(statsdef->dt)));
    }
    if (timedef->status_dt > 0) {
        n = max(n, -floor(log10(timedef->status_dt)));
    }

    // Restart, statistics, and status frequencies measured in
    // time steps checking against minimum step sizes
    if (restartdef->nt > 0) {
        n = max(n, -floor(log10(timedef->min_dt * restartdef->nt)));
    }
    if (statsdef->nt > 0) {
        n = max(n, -floor(log10(timedef->min_dt * statsdef->nt)));
    }
    if (timedef->status_nt > 0) {
        n = max(n, -floor(log10(timedef->min_dt * timedef->status_nt)));
    }

    // Maximum simulation duration measured as time
    if (timedef->advance_dt > 0) {
        n = max(n, -floor(log10(timedef->advance_dt)));
    }

    // Maximum time step sizes as a human-friendly service to the user
    if (timedef->max_dt > 0) {
        n = max(n, -floor(log10(timedef->max_dt)));
    }

    // Pad the result by one decimal place so the user realizes
    // that the time steps are not exactly what is output as
    // some drift will likely be noticeable in this final digit
    n += 1;

    // While n was real-valued to please the type system,
    // the result mathematically had to be an integer.
    return static_cast<int>(n);
}

bool
driver_base::log_status(
        const driver_base::time_type t,
        const driver_base::step_type nt)
{
    // Notice collective operations are never inside logging macros!

    // Defensively avoid multiple invocations with no intervening changes
    if (last_status_nt == nt) {
        DEBUG0(who, "Cowardly refusing to repeatedly show status at nt = "
               << nt);
        return true;
    }

    // Common message prefix used across all status-related routines
    const std::string timeprefix(build_timeprefix(t, nt));

    SUZERAIN_TIMER_SCOPED("driver_base::log_status");

    // Log information about the various quantities of interest
    log_status_bulk(timeprefix);
    log_status_L2(timeprefix);
    log_status_boundary_state(timeprefix);

    // Permit subclasses to dump arbitrary status information.  E.g. MMS error
    const bool retval = log_status_hook(timeprefix, t, nt);

    last_status_nt = nt; // Maintain last status time step

    return retval;
}

void
driver_base::log_status_L2(
        const std::string& timeprefix,
        const char * const name_L2,
        const char * const name_rms)
{
    // Avoid computational cost when logging is disabled
    logging::logger_type log_L2  = logging::get_logger(name_L2);
    logging::logger_type log_rms = logging::get_logger(name_rms);
    if (!INFO0_ENABLED(log_L2) && !INFO0_ENABLED(log_rms)) return;

    // Show headers only on first invocation
    std::ostringstream msg;
    if (!log_status_L2_header_shown) {
        msg << std::setw(timeprefix.size()) << build_timeprefix_description();
        for (size_t k = 0; k < fields.size(); ++k)
            msg << ' ' << std::setw(fullprec<>::width) << fields[k].identifier;
        INFO0(log_L2, msg.str());
        INFO0(log_rms, msg.str());
        msg.str("");
        log_status_L2_header_shown = true;
    }

    // Collective computation of the L_2 norms
    state_nonlinear->assign(*state_linear);
    const std::vector<field_L2> result
        = compute_field_L2(*state_nonlinear, *grid, *dgrid, *gop);

    // Build and log L2 of mean conserved state
    msg << timeprefix;
    for (size_t k = 0; k < result.size(); ++k) {
        msg << ' ' << fullprec<>(result[k].mean());
    }
    INFO0(log_L2, msg.str());

    // Build and log root-mean-squared-fluctuations of conserved state
    // RMS fluctuations are a scaling factor away from L2 fluctuations
    const real_t rms_coeff = 1/std::sqrt(grid->L.x()*grid->L.y()*grid->L.z());
    msg.str("");
    msg << timeprefix;
    for (size_t k = 0; k < result.size(); ++k) {
        msg << ' ' << fullprec<>(rms_coeff*result[k].fluctuating());
    }
    INFO0(log_rms, msg.str());
}

void
driver_base::log_status_bulk(
        const std::string& timeprefix)
{
    // Only continue on the rank housing the zero-zero modes...
    if (!dgrid->has_zero_zero_modes()) return;

    // ...and when logging is enabled.  Notice INFO not INFO0 is used.
    logging::logger_type bulk_state = logging::get_logger("bulk.state");
    if (!INFO_ENABLED(bulk_state)) return;

    // Show headers only on first invocation
    std::ostringstream msg;
    if (!log_status_bulk_header_shown) {
        msg << std::setw(timeprefix.size()) << build_timeprefix_description();
        for (size_t k = 0; k < fields.size(); ++k)
            msg << ' ' << std::setw(fullprec<>::width) << fields[k].identifier;
        INFO0(bulk_state, msg.str());
        msg.str("");
        log_status_bulk_header_shown = true;
    }

    // Compute operator for finding bulk quantities from coefficients
    VectorXr bulkcoeff(b->n());
    b->integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid->L.y();

    // Prepare the status message and log it
    msg << timeprefix;
    for (size_t k = 0; k < state_linear->shape()[0]; ++k) {
        Map<VectorXc> mean(
                (*state_linear)[k].origin(), state_linear->shape()[1]);
        msg << ' ' << fullprec<>(bulkcoeff.dot(mean.real()));
    }
    INFO(bulk_state, msg.str());
}

void
driver_base::log_status_boundary_state(
        const std::string& timeprefix)
{
    // Only continue on the rank housing the zero-zero modes.
    if (!dgrid->has_zero_zero_modes()) return;

    logging::logger_type nick[2] = { logging::get_logger("bc.lower"),
                                     logging::get_logger("bc.upper")  };

    // Indices at the lower and upper walls.  Use that bc collocation point
    // values are nothing but the first and last B-spline coefficient values.
    size_t bc[2] = { 0, state_linear->shape()[1] - 1 };

    for (size_t l = 0; l < SUZERAIN_COUNTOF(bc); ++l) {

        // Avoid computational cost when logging is disabled
        if (!DEBUG_ENABLED(nick[l])) continue;

        std::ostringstream msg;
        msg << timeprefix;

        for (size_t k = 0; k < fields.size(); ++k) {
            real_t val = ((*state_linear)[k][bc[l]][0][0]).real();
            msg << ' ' << fullprec<>(val);
        }
        DEBUG(nick[l], msg.str());
    }
}

void
driver_base::save_metadata()
{
    SUZERAIN_TIMER_SCOPED("driver_base::save_metadata");

    DEBUG0(who, "Saving metadata temporary file: " << restartdef->metadata);

    esio_handle esioh = esio_handle_initialize(MPI_COMM_WORLD);
    esio_file_create(esioh, restartdef->metadata.c_str(), 1 /* overwrite */);

    // Per Ticket #2595 we cannot ignore string-based attributes set on "/" in
    // version 1.8.x of h5diff.  This makes ignoring version information
    // problematic for regression tests.  As a workaround, we write a dataset
    // whose path can be excluded by h5diff version 1.8.6 or later.
    static const char ignorable_location[] = "metadata_generated";

    // Store time(2) result as floating point to avoid sizeof(time_t) issues
    const real_t now = static_cast<real_t>(std::time(NULL));
    // Though all ranks invoked time(2), only the value from rank zero is saved
    int procid;
    esio_handle_comm_rank(esioh, &procid);
    esio_line_establish(esioh, 1, 0, (procid == 0 ? 1 : 0));
    esio_line_write(esioh, ignorable_location, &now, 0,
            "Time since the Epoch when metadata saved according to time(2)");

    // Off that entry we may hang any number of ignorable attributes
    esio_string_set(esioh, ignorable_location, "revision", revstr.c_str());

    // TODO Broadcast who from rank zero and save
    // TODO Broadcast uname details from rank zero and save

    // Save known metadata and then invoke subclass extension point
    save_grid_and_operators(esioh);
    if (timedef) save(esioh, *timedef);
    save_metadata_hook(esioh);

    esio_file_close(esioh);
    esio_handle_finalize(esioh);

    metadata_saved = true;
}

void
driver_base::load_metadata(
        const esio_handle esioh)
{
    SUZERAIN_TIMER_SCOPED("driver_base::load_metadata");

    load_grid_and_operators(esioh);

    SUZERAIN_ENSURE(timedef);
    load(esioh, *timedef);

    // Invoke subclass extension point
    load_metadata_hook(esioh);
}

// The file created below must match that made by other save_restart overloads!
bool
driver_base::save_restart(
        const driver_base::time_type t,
        const driver_base::step_type nt)
{
    // Defensively avoid multiple invocations with no intervening changes
    if (last_restart_saved_nt == nt) {
        DEBUG0(who, "Cowardly refusing to save multiple restarts at nt = "
               << nt);
        return true;
    }

    if (!metadata_saved) save_metadata();
    const std::string timeprefix(build_timeprefix(t, nt));

    // Time only after prerequisites are satisfied
    SUZERAIN_TIMER_SCOPED("driver_base::save_restart");

    const double starttime = MPI_Wtime();
    INFO0(who, timeprefix << " Starting to save restart");
    esio_handle esioh = esio_handle_initialize(MPI_COMM_WORLD);

    DEBUG0(who, "Cloning " << restartdef->metadata
           << " to " << restartdef->uncommitted);
    esio_file_clone(esioh, restartdef->metadata.c_str(),
                    restartdef->uncommitted.c_str(), 1 /*overwrite*/);
    save_time(esioh, t);

    // Invoke subclass extension points for both restart AND statistics
    state_nonlinear->assign(*state_linear);
    const bool continue_advancing =    save_state_hook(esioh)
                                    && save_statistics_hook(esioh, t);

    DEBUG0(who, "Committing " << restartdef->uncommitted
           << " as a restart file using template " << restartdef->destination);
    esio_file_close_restart(esioh, restartdef->destination.c_str(),
                            restartdef->retain);
    esio_handle_finalize(esioh);

    const double elapsed = MPI_Wtime() - starttime;
    INFO0(who, timeprefix << " Committed restart in " << elapsed << " seconds");

    last_restart_saved_nt = nt; // Maintain last successful restart time step

    return continue_advancing;
}

// The file created below must match that made by other save_restart overloads!
bool
driver_base::save_restart(
        const driver_base::time_type t,
        const std::string dstfile,
        const bool overwrite)
{
    const double starttime = MPI_Wtime();
    INFO0(who, "Starting to save restart file " << dstfile);
    SUZERAIN_TIMER_SCOPED("driver_base::save_restart (one-off)");

    // Save metadata ensuring previously set values are not reused
    metadata_saved = false;
    save_metadata();

    esio_handle esioh = esio_handle_initialize(MPI_COMM_WORLD);
    DEBUG0(who, "Cloning " << restartdef->metadata << " to " << dstfile);
    esio_file_clone(esioh, restartdef->metadata.c_str(),
                    dstfile.c_str(), overwrite);
    save_time(esioh, t);

    // Invoke subclass extension points for both restart AND statistics
    state_nonlinear->assign(*state_linear);
    const bool continue_advancing =    save_state_hook(esioh)
                                    && save_statistics_hook(esioh, t);

    esio_file_close(esioh);
    esio_handle_finalize(esioh);

    const double elapsed = MPI_Wtime() - starttime;
    INFO0(who, "Committed restart file " << dstfile
          << " in "<< elapsed << " seconds");

    return continue_advancing;
}

void
driver_base::load_restart(
        const esio_handle esioh,
        real_t& t)
{
    SUZERAIN_TIMER_SCOPED("driver_base::load_restart");

    const double begin = MPI_Wtime();

    load_metadata(esioh);
    load_time(esioh, t);
    establish_decomposition();
    establish_state_storage(fields.size(), fields.size());
    load_state_hook(esioh);                  // Invoke subclass extension point
    state_linear->assign(*state_nonlinear);  // Copy into state_linear
    load_statistics_hook(esioh);             // Invoke subclass extension point

    wtime_load_restart = MPI_Wtime() - begin;
}

bool
driver_base::save_statistics(
        const driver_base::time_type t,
        const driver_base::step_type nt)
{
    // Defensively avoid multiple invocations with no intervening changes
    if (last_statistics_saved_nt == nt) {
        DEBUG0(who, "Cowardly refusing to save multiple samples at nt = "
               << nt);
        return true;
    }

    if (!metadata_saved) save_metadata();
    const std::string timeprefix(build_timeprefix(t, nt));

    // Time only after prerequisites are satisfied
    SUZERAIN_TIMER_SCOPED("driver_base::save_statistics");

    const double starttime = MPI_Wtime();
    DEBUG0(who, timeprefix << " Started to store statistics");
    esio_handle esioh = esio_handle_initialize(MPI_COMM_WORLD);

    // We use restartdef->{metadata,uncommitted} for statistics too.
    DEBUG0(who, "Cloning " << restartdef->metadata
           << " to " << restartdef->uncommitted);
    esio_file_clone(esioh, restartdef->metadata.c_str(),
                    restartdef->uncommitted.c_str(), 1 /*overwrite*/);
    save_time(esioh, t);

    // Invoke subclass extension point
    const bool continue_advancing = save_statistics_hook(esioh, t);

    DEBUG0(who, "Committing " << restartdef->uncommitted
           << " as a statistics file using template "
           << statsdef->destination);
    esio_file_close_restart(esioh, statsdef->destination.c_str(),
                            statsdef->retain);
    esio_handle_finalize(esioh);

    const double elapsed = MPI_Wtime() - starttime;
    INFO0(who, timeprefix << " Committed statistics in "
          << elapsed << " seconds");

    last_statistics_saved_nt = nt; // Maintain last successful statistics nt

    return continue_advancing;
}

void
driver_base::load_statistics(
        const esio_handle esioh,
        real_t& t)
{
    SUZERAIN_TIMER_SCOPED("driver_base::load_statistics");

    load_metadata(esioh);
    // Per documentation, establish_decomposition() not invoked!
    load_time(esioh, t);
    load_statistics_hook(esioh);  // Invoke subclass extension point
}

bool
driver_base::log_status_hook(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    SUZERAIN_UNUSED(timeprefix);
    SUZERAIN_UNUSED(t);
    SUZERAIN_UNUSED(nt);

    return true;
}

void
driver_base::save_metadata_hook(
        const esio_handle esioh)
{
    SUZERAIN_UNUSED(esioh);

    // For example:
    //     perfect::store(h, scenario);
    //     perfect::store(h, scenario, grid, msoln);
}

void
driver_base::load_metadata_hook(
        const esio_handle esioh)
{
    SUZERAIN_UNUSED(esioh);

    // For example:
    //     perfect::load(h, scenario);
    //     perfect::load(h, scenario, grid, msoln);
}

bool
driver_base::save_state_hook(
        const esio_handle esioh)
{
    SUZERAIN_ENSURE(restartdef);
    SUZERAIN_ENSURE(state_linear);
    SUZERAIN_ENSURE(state_nonlinear);

    // Save either coefficients or collocation values, as requested
    // The former is a destructive operation clobbering *state_nonlinear
    if (restartdef->physical) {
        save_collocation_values(
                esioh, fields, *state_nonlinear, *grid, *dgrid, *cop);
    } else {
        save_coefficients(
                esioh, fields, *state_nonlinear, *grid, *dgrid);
    }

    return true;
}

void
driver_base::load_state_hook(
        const esio_handle esioh)
{
    SUZERAIN_ENSURE(grid);
    SUZERAIN_ENSURE(dgrid);

    // Check if load_coefficients(), load_collocation_values() should work
    bool allreal = true, allcplx = true;
    for (size_t i = 0; i < fields.size(); ++i) {
        int ncomponents = 0;
        const int status = esio_field_sizev(esioh, fields[i].location.c_str(),
                                            0, 0, 0, &ncomponents);
        if (status == ESIO_SUCCESS) {
            allreal = allreal && (ncomponents == 1);
            allcplx = allcplx && (ncomponents == 2);
            if (ncomponents == 0 || ncomponents > 2) {
                WARN0(who, "Field /" << fields[i].location
                      << " looks fishy; ncomponents = " << ncomponents);
            }
        } else {
            WARN0(who, "Field /" << fields[i].location
                  << " not found in restart");
        }
    }

    // Ensure our state storage is up to the task at hand
    establish_state_storage(fields.size(), fields.size());

    // Dispatch to the appropriate state-loading logic
    if (allcplx) {
        load_coefficients(
                esioh, fields, *state_nonlinear, *grid, *dgrid, *cop, *b);
    } else if (allreal) {
        load_collocation_values(
                esioh, fields, *state_nonlinear, *grid, *dgrid, *cop, *b);
    } else {
        SUZERAIN_ERROR_VOID(
                "Unable to load state from file", SUZERAIN_EFAILED);
    }
}

bool
driver_base::save_statistics_hook(
        const esio_handle esioh,
        const driver_base::time_type t)
{
    SUZERAIN_UNUSED(esioh);
    SUZERAIN_UNUSED(t);

    // For example:
    //     perfect::store(esioh, samples);

    return true;
}

void
driver_base::load_statistics_hook(
        const esio_handle esioh)
{
    SUZERAIN_UNUSED(esioh);

    // For example:
    //     perfect::load(esioh, samples);
}

bool
driver_base::process_any_signals_received(
        const time_type t,
        const step_type nt)
{
    // The signal logger is looked up whenever logging occurs, rather
    // than once at the beginning of the routine via something like
    //   logging::logger_type log_signal = logging::get_logger("signal");
    // to avoid overhead in the common case when nothing is logged.
    static const char log_signal[] = "signal";

    // this->allreducer performs the Allreduce necessary to get local status
    // from signal::global_received into actions within this->signal_received.

    // Keep advancing time unless keep_advancing is set false
    bool keep_advancing = true;

    if (signal_received[signal::log_status]) {
        INFO0(log_signal,
              "Outputting simulation status due to receipt of "
              << suzerain_signal_name(signal_received[signal::log_status]));
        keep_advancing = keep_advancing && log_status(t, nt);
    }

    if (signal_received[signal::write_restart]) {
        INFO0(log_signal,
              "Writing restart file due to receipt of "
              << suzerain_signal_name(signal_received[signal::write_restart]));
        keep_advancing = keep_advancing && save_restart(t, nt);
    }

    if (signal_received[signal::teardown_reactive]) {
        const char * const name = suzerain_signal_name(
                signal_received[signal::teardown_reactive]);
        INFO0(log_signal,
              "Initiating teardown due to receipt of " << name);
        soft_teardown  = true;
        keep_advancing = false;
        switch (signal_received[signal::teardown_reactive]) {
            case SIGINT:
            case SIGTERM:
                INFO0(log_signal,
                      "Receipt of another " << name <<
                      " will forcibly terminate program");
                signal2(signal_received[signal::teardown_reactive], SIG_DFL);
                break;
        }
    }

    if (signal_received[signal::teardown_proactive]) {
        INFO0(log_signal,
              "Initiating proactive teardown because of wall time constraint");
        soft_teardown  = true;
        keep_advancing = false;
    }

    if (signal_received[signal::write_statistics]) {
        const char * const name = suzerain_signal_name(
                signal_received[signal::write_statistics]);
        INFO0(log_signal,
              "Computing and writing statistics due to receipt of " << name);
        keep_advancing = keep_advancing && save_statistics(t, nt);
    }

    // Clear signal_received to defensively avoid stale data bugs.
    // These would only be problematic if process_any_signals_received
    // was run multiple times in between delta_t_allreducer invocations.
    signal_received.assign(0);

    return keep_advancing;
}

delta_t_allreducer::delta_t_allreducer(
        const double& wtime_mpi_init,
        const double& wtime_fftw_planning,
        const shared_ptr<time_definition>& timedef,
        const double& wtime_load_restart,
        const double& wtime_advance_start,
        const driver_base::step_type& last_status_nt,
        const driver_base::step_type& last_restart_saved_nt,
        driver_base::delta_t_ratios_type& delta_t_ratios,
        signal::received_type& signal_received)
    : wtime_mpi_init(wtime_mpi_init)
    , wtime_fftw_planning(wtime_fftw_planning)
    , timedef(timedef)
    , wtime_load_restart(wtime_load_restart)
    , wtime_advance_start(wtime_advance_start)
    , last_status_nt(last_status_nt)
    , last_restart_saved_nt(last_restart_saved_nt)
    , delta_t_ratios(delta_t_ratios)
    , signal_received(signal_received)
{
}

real_t delta_t_allreducer::operator()(
        const std::vector<real_t>& delta_t_candidates)
{
    // Copy incoming candidates so we may mutate them
    std::vector<real_t> candidates(delta_t_candidates);

    // Take atomic snapshot of and then clear signal::global_received
    signal_received = signal::global_received;
    static const signal::volatile_received_type::value_type zero = 0;
    signal::global_received.assign(zero);

    if (DEBUG_ENABLED()) {
        const int rank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
        for (std::size_t i = 0; i < signal::received_type::static_size; ++i) {
            if (signal_received[i]) {
                DEBUG("Received signal number " << signal_received[i]
                      << " on rank " << rank);
            }
        }
    }

    // When possible, obtain time step period statistics and a projected
    // wall time for when we could complete the next time step and dump a
    // restart file.  If the projection is after --advance_wt, register
    // teardown.  Logic here is key to proactive soft_teardown success.
    static double wtime_last = std::numeric_limits<double>::quiet_NaN();
    if (timedef->advance_wt > 0 && (boost::math::isfinite)(wtime_last)) {

        // Accumulate time step period statistics
        const double wtime = MPI_Wtime();
        period(wtime - wtime_last);

        // Find a pessimistic time for the next time step completion...
        // (that is, finish current step *and* finish another one)
        namespace acc = boost::accumulators;
        double wtime_projected = wtime
            + 2*(acc::mean(period) + 3*std::sqrt(acc::variance(period)));
        // ...to which we add a pessimistic estimate for dumping a restart
        if (last_restart_saved_nt == std::numeric_limits<size_t>::max()) {
            wtime_projected += 2*wtime_load_restart;  // Load as surrogate
        } else {
            wtime_projected += (acc::max)(period);  // Includes dumps
        }
        // ...to which we add an estimate of other finalization costs
        wtime_projected += 2*(wtime_advance_start - wtime_mpi_init
                                                  - wtime_fftw_planning);

        // Raise a "signal" if we suspect we cannot teardown quickly enough
        signal_received[signal::teardown_proactive]
            = (wtime_projected >= wtime_mpi_init + timedef->advance_wt);
        if (DEBUG_ENABLED() && signal_received[signal::teardown_proactive]) {
            const int rank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
            DEBUG("Rank " << rank << " projects delaying teardown until "
                  << wtime_projected - wtime_mpi_init
                  << " elapsed seconds is dangerous.");
        }
    }
    wtime_last = MPI_Wtime();

    // Push a negated version of signal_received onto end of candidates
    candidates.reserve(candidates.size() + signal::received_type::static_size);
    std::transform(signal_received.begin(), signal_received.end(),
                   std::back_inserter(candidates),
                   std::negate<signal::received_type::value_type>());

    // Allreduce so each rank knows the minimum of all candidates
    assert(candidates.size() > 0);
    SUZERAIN_MPICHKR(MPI_Allreduce(MPI_IN_PLACE,
                                   &candidates.front(), candidates.size(),
                                   suzerain::mpi::datatype<real_t>::value,
                                   MPI_MIN, MPI_COMM_WORLD));

    // Negate the signal_received details once again to get a logical MAX
    // stored within signal_received.  Erase temporaries from candidates.
    std::transform(candidates.begin() + delta_t_candidates.size(),
                   candidates.begin() + delta_t_candidates.size()
                                      + signal::received_type::static_size,
                   signal_received.begin(),
                   std::negate<signal::received_type::value_type>());
    candidates.erase(candidates.begin() + delta_t_candidates.size(),
                     candidates.begin() + delta_t_candidates.size()
                                        + signal::received_type::static_size);

    // Delegate finding-the-minimum work on each rank to superclass
    // whose logic should enforce requirement that min(NaN,x) == NaN
    const real_t delta_t = super::operator()(candidates);

    // Update delta_t_ratios using the just chosen delta_t
    // isnan used to avoid a NaN from destroying all accumulator data
    if (!(boost::math::isnan)(delta_t)) {
        const size_t n = candidates.size();
        if (SUZERAIN_UNLIKELY(delta_t_ratios.size() < n)) {
            delta_t_ratios.resize(n);
        }
        for (size_t i = 0; i < n; ++i) {
            delta_t_ratios[i](candidates[i] / delta_t);
        }
    }

    return delta_t;
}

void
driver_base::default_statistics_interval(
        time_type& dt,
        step_type& nt)
{
    // Possibly compute the default restart writing interval
    using std::numeric_limits;
    time_type tmp_dt = numeric_limits<time_type>::quiet_NaN();
    step_type tmp_nt = 0;
    default_restart_interval(tmp_dt, tmp_nt);

    // Output statistics one-fourth as often as we restart
    if (boost::math::isnormal(tmp_dt)) {
        dt = tmp_dt / 4;
    }
    if (tmp_nt) {
        using std::min;
        nt = min(tmp_nt / 4, static_cast<step_type>(1));
    }
}

void
driver_base::default_status_interval(
        time_type& dt,
        step_type& nt)
{
    // Possibly compute the default statistics writing interval
    using std::numeric_limits;
    time_type tmp_dt = numeric_limits<time_type>::quiet_NaN();
    step_type tmp_nt = 0;
    default_statistics_interval(tmp_dt, tmp_nt);

    // Output status one-fourth as often as we write statistics
    if (boost::math::isnormal(tmp_dt)) {
        dt = tmp_dt / 4;
    }
    if (tmp_nt) {
        using std::min;
        nt = min(tmp_nt / 4, static_cast<step_type>(1));
    }
}

} // end namespace support

} // end namespace suzerain
