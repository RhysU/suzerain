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

/** @file
 * @copydoc driver_base.hpp
 */

#include <suzerain/support/driver_base.hpp>

#include <esio/esio.h>
#include <esio/error.h>

#include <suzerain/bl.h>
#include <suzerain/channel.h>
#include <suzerain/countof.h>
#include <suzerain/error.h>
#include <suzerain/extrema.hpp>
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
    //
    // Determine if we should halt the simulation due to the signal
    end = driver_base::signaldef.halt.end();
    if (find(driver_base::signaldef.halt.begin(), end, sig) != end) {
        signal::global_received[signal::halt_reactive] = sig;
    }

    // Determine if we should compute and write statistics due to the signal
    end = driver_base::signaldef.statistics.end();
    if (find(driver_base::signaldef.statistics.begin(), end, sig) != end) {
        signal::global_received[signal::write_statistics] = sig;
    }
}

} // end extern "C"

// Static member initialization
definition_signal driver_base::signaldef;

driver_base::driver_base(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &revstr)
    : application_base(application_synopsis,
                       argument_synopsis,
                       description,
                       revstr)
    , restartdef(make_shared<definition_restart>(
                /* metadata    */ "metadata.h5.XXXXXX",
                /* uncommitted */ "uncommitted.h5.XXXXXX",
                /* destination */ "restart#####.h5"))
    , statsdef(make_shared<definition_statistics>(
                /* destination */ "sample#####.h5"))
    , timedef(make_shared<definition_time>(
                /* advance_dt   */ 0,
                /* advance_nt   */ 0,
                /* advance_wt   */ 0,
                /* status_dt    */ 0,
                /* status_nt    */ 0,
                /* status_final */ true,
                /* min_dt       */ 1e-8,
                /* max_dt       */ 1))
    , method()
    , L()
    , N()
    , controller()
    , received_teardown(true)
    , received_halt(false)
    , header_shown()
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
    using std::fill;
    fill(signal_received.begin(), signal_received.end(), 0);
}

std::string
driver_base::log4cxx_config()
{
    return super::log4cxx_config() + // Append to the default configuration
        "\n"
        "## Collect \"bc\" messages into only bc.dat mimicking LOG file behavior\n"
        "log4j.logger.bc=INHERITED, BC\n"
        "log4j.additivity.bc=false\n"
        "log4j.appender.BC=${log4j.appender.LOG}\n"
        "log4j.appender.BC.filename=bc.dat\n"
        "log4j.appender.BC.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.BC.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.BC.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
        "\n"
        "## Collect \"bl\"   messages into only qoi.dat mimicking LOG file behavior\n"
        "## Collect \"chan\" messages into only qoi.dat mimicking LOG file behavior\n"
        "log4j.logger.bl=INHERITED, QOI\n"
        "log4j.additivity.bl=false\n"
        "log4j.logger.chan=INHERITED, QOI\n"
        "log4j.additivity.chan=false\n"
        "log4j.appender.QOI=${log4j.appender.LOG}\n"
        "log4j.appender.QOI.filename=qoi.dat\n"
        "log4j.appender.QOI.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.QOI.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.QOI.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
        "\n"
        "## Collect \"state\" messages into only state.dat mimicking LOG file behavior\n"
        "log4j.logger.state=INHERITED, STATE\n"
        "log4j.additivity.state=false\n"
        "log4j.appender.STATE=${log4j.appender.LOG}\n"
        "log4j.appender.STATE.filename=state.dat\n"
        "log4j.appender.STATE.append=${log4j.appender.LOG.append}\n"
        "log4j.appender.STATE.layout=${log4j.appender.LOG.layout}\n"
        "log4j.appender.STATE.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
        "\n"
        "## Additionally ensure \"state.L2\" messages reach the CONSOLE and LOG\n"
        "log4j.logger.state.L2=INHERITED, CONSOLE, LOG\n"
        "\n"
        "## Additionally ensure \"state.RMS\" messages reach the CONSOLE and LOG\n"
        "log4j.logger.state.RMS=INHERITED, CONSOLE, LOG\n"
    ;
}

std::vector<std::string>
driver_base::initialize(int argc, char **argv)
{
    // Only add groups of options when non-trivial at initialization
    if (restartdef) options.add_definition(*restartdef);
    if (statsdef  ) options.add_definition(*statsdef  );
    if (timedef   ) {
        options.add_definition(*timedef   ); // Only care about signals,
        options.add_definition( signaldef ); // when we will advance time
    }

    // Process incoming arguments by invoking superclass method
    std::vector<std::string> positional = super::initialize(argc, argv);

    // Generating unique file names as needed using mkstemp(3)
    // See Redmine ticket #2385 towards a nicer long-term solution
    //
    // Awful, awful checks permit empty {restart,stats}def
    // while maintaining a single broadcast operation.
    {
        // Pack a temporary buffer with the four file name templates
        array<size_t,5> pos = {{
            0,
            restartdef ? restartdef->metadata.length()    + 1 : 1,
            restartdef ? restartdef->uncommitted.length() + 1 : 1,
            restartdef ? restartdef->destination.length() + 1 : 1,
            statsdef   ? statsdef  ->destination.length() + 1 : 1
        }};
        std::partial_sum(pos.begin(), pos.end(), pos.begin());
        std::vector<char> buf(pos[4], '\0');
        if (restartdef) {
            strcpy(&buf[pos[0]], restartdef->metadata.c_str());
            strcpy(&buf[pos[1]], restartdef->uncommitted.c_str());
            strcpy(&buf[pos[2]], restartdef->destination.c_str());
        }
        if (statsdef) {
            strcpy(&buf[pos[3]], statsdef  ->destination.c_str());
        }

        // Generate unique files to be overwritten and/or just file names.
        // File generation relies on template semantics of mkstemp(3).
        // Error checking kept minimal as failures here should not be fatal.
        // This is not particularly robust but it should serve our needs.
        if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0) {
            if (restartdef) {
                if (boost::ends_with(restartdef->metadata, "XXXXXX")) {
                    close(mkstemp(&buf[pos[0]]));  // Clobbered later...
                }
                if (boost::ends_with(restartdef->uncommitted, "XXXXXX")) {
                    close(mkstemp(&buf[pos[1]]));  // Possibly clobbered later.
                    unlink(&buf[pos[1]]);          // ...so remove any evidence
                }
                if (boost::ends_with(restartdef->destination, "XXXXXX")) {
                    close(mkstemp(&buf[pos[2]]));  // Not clobbered later...
                    unlink(&buf[pos[2]]);          // ...so remove any evidence
                }
            }
            if (statsdef) {
                if (boost::ends_with(statsdef->destination, "XXXXXX")) {
                    close(mkstemp(&buf[pos[3]]));  // Not clobbered later...
                    unlink(&buf[pos[3]]);          // ...so remove any evidence
                }
            }
        }

        // Broadcast any generated names to all ranks and unpack values
        SUZERAIN_MPICHKR(MPI_Bcast(&buf[0], pos[4],
                         mpi::datatype<char>(), 0, MPI_COMM_WORLD));
        if (restartdef) {
            restartdef->metadata    = &buf[pos[0]];
            restartdef->uncommitted = &buf[pos[1]];
            restartdef->destination = &buf[pos[2]];
        }
        if (statsdef) {
            statsdef  ->destination = &buf[pos[3]];
        }
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
        this->method.reset(new lowstorage::method<
                    lowstorage::smr91,
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
    controller.reset(make_controller(
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
        s.insert(s.end(), signaldef.halt.begin(),
                          signaldef.halt.end());
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

double
driver_base::advance_controller(
            const bool final_status,
            const bool final_statistics,
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
    received_teardown = false;
    received_halt = false;
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
    if (received_teardown) {
        INFO0(who, "Time controller stopped in reaction to tear down request");
        success = true; // ...treat like successful advance
    } else if (received_halt) {
        INFO0(who, "Time controller halted because of incoming signal");
        success = false; // ...treat like failed advance
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

    // Save a final statistics if not just saved during time advancement
    if (   final_statistics
        && success
        && last_statistics_saved_nt != controller->current_nt()) {
        INFO0(who, "Saving statistics file after time controller finished");
        save_statistics(controller->current_t(), controller->current_nt());
    }

    // Save a final restart if not just saved during time advancement
    // Smart caching by applications will prevent pointless recomputation here
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
        msg << "Min/avg/max/std of finite delta_t: "
            << controller->taken_min()  << ", "
            << controller->taken_mean() << ", "
            << controller->taken_max()  << ", "
            << controller->taken_stddev();
        INFO0(who, msg.str());
        msg.str("");
        msg << "Mean finite delta_t criteria versus minimum criterion: ";
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
        const driver_base::step_type nt) const
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
            const char * describe_nt) const
{
    // Build string with inter-label spacing matching build_timeprefix output
    std::ostringstream oss;
    oss << describe_t
        << ' '
        << std::setw(timeprefix_setw_nt) << describe_nt;
    return oss.str();
}

int
driver_base::build_timeprefix_mantissa_digits() const
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

    // Permit subclasses to dump arbitrary status information
    // (for example, logging method of manufactured solution error)
    // either before, after, or in lieu of driver_base::log_status_hook.
    const bool retval = log_status_hook(timeprefix, t, nt);

    last_status_nt = nt; // Maintain last status time step

    return retval;
}

// Often-reused logic to show a "t nt field[0]::identifier..." header once
template <class Logger>
static void
maybe_timeprefix_fields_identifiers(driver_base& db,
                                    const std::string& timeprefix,
                                    Logger& log,
                                    bool& shown)
{
    if (!shown) {
        std::ostringstream msg;
        msg << std::setw(timeprefix.size()) << db.build_timeprefix_description();
        for (size_t k = 0; k < db.fields.size(); ++k) {
            msg << ' ' << std::setw(fullprec<>::width) << db.fields[k].identifier;
        }
        INFO0(log, msg.str());
        shown = true;
    }
}

void
driver_base::log_state_L2(
        const std::string& timeprefix,
        const char * const name_L2,
        const char * const name_RMS)
{
    // RMS fluctuations only make sense to log when either X or Z is nontrivial
    const bool nontrivial_rms_possible = grid->N.x() * grid->N.z() > 1;

    // Show headers only on first invocation
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_L2,
                                        header_shown[name_L2]);
    if (nontrivial_rms_possible) {
        maybe_timeprefix_fields_identifiers(*this, timeprefix, name_RMS,
                                            header_shown[name_RMS]);
    }

    // Collective computation of the $L^2_{xyz}$ norms
    state_nonlinear->assign_from(*state_linear);
    const std::vector<field_L2xyz> result
        = compute_field_L2xyz(*state_nonlinear, *grid, *dgrid, *gop);

    // Build and log L2 of mean conserved state
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t k = 0; k < result.size(); ++k) {
        msg << ' ' << fullprec<>(result[k].mean);
    }
    INFO0(name_L2, msg.str());

    // Build and log root-mean-squared-fluctuations of conserved state
    // RMS fluctuations are a scaling factor away from L2 fluctuations
    if (nontrivial_rms_possible) {
        real_t RMS_coeff = 1/std::sqrt(grid->L.x()*grid->L.y()*grid->L.z());
        msg.str("");
        msg << timeprefix;
        for (size_t k = 0; k < result.size(); ++k) {
            msg << ' ' << fullprec<>(RMS_coeff*result[k].fluctuating);
        }
        INFO0(name_RMS, msg.str());
    }
}

void
driver_base::log_state_bulk(
        const std::string& timeprefix,
        const char * const name_bulk)
{
    // Only continue on the rank housing the zero-zero modes...
    //...which currently must also be rank zero for logging purposes.
    if (dgrid->has_zero_zero_modes()) {
        if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) != 0)
            SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    } else {
        return;
    }

    // Show headers only on first invocation
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_bulk,
                                        header_shown[name_bulk]);

    // Compute operator for finding bulk quantities from coefficients
    VectorXr bulkcoeff(b->n());
    b->integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid->L.y();

    // Prepare the status message and log it
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t k = 0; k < state_linear->shape()[0]; ++k) {
        Map<VectorXc> mean(
                (*state_linear)[k].origin(), state_linear->shape()[1]);
        msg << ' ' << fullprec<>(bulkcoeff.dot(mean.real()));
    }
    INFO0(name_bulk, msg.str());
}

// FIXME Redmine #3056 load/save extrema from restart files
// Otherwise, we now require undesirable parallel FFTs to get status!
// FIXME Lazily recompute extrema only when necessary
// Otherwise, saving restarts and statistics files forces recomputation
void
driver_base::log_state_extrema(
        const std::string& timeprefix,
        const char * const name_min  ,
        const char * const name_xmin ,
        const char * const name_ymin ,
        const char * const name_zmin ,
        const char * const name_max  ,
        const char * const name_xmax ,
        const char * const name_ymax ,
        const char * const name_zmax ,
        const char * const name_fneg )
{
    // Fluctuations only make sense to log when either X or Z is nontrivial
    const bool nontrivial_possible = grid->N.x() * grid->N.z() > 1;
    if (!nontrivial_possible) return;

    // Show headers only on first invocation
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_min,
                                        header_shown[name_min]);
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_xmin,
                                        header_shown[name_xmin]);
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_ymin,
                                        header_shown[name_ymin]);
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_zmin,
                                        header_shown[name_zmin]);

    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_max,
                                        header_shown[name_max]);
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_xmax,
                                        header_shown[name_xmax]);
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_ymax,
                                        header_shown[name_ymax]);
    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_zmax,
                                        header_shown[name_zmax]);

    maybe_timeprefix_fields_identifiers(*this, timeprefix, name_fneg,
                                        header_shown[name_fneg]);

    // Collective computation of the global extrema
    state_nonlinear->assign_from(*state_linear);
    const std::vector<field_extrema_xyz> result
        = compute_field_extrema_xyz(*state_nonlinear, *grid, *dgrid, *cop);

    // Build and log global minimum and maximum values
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(result[f].min);
    }
    INFO0(name_min, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(grid->x(result[f].imin));
    }
    INFO0(name_xmin, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(b->collocation_point(result[f].jmin));
    }
    INFO0(name_ymin, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(grid->z(result[f].kmin));
    }
    INFO0(name_zmin, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(result[f].max);
    }
    INFO0(name_max, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(grid->x(result[f].imax));
    }
    INFO0(name_xmax, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(b->collocation_point(result[f].jmax));
    }
    INFO0(name_ymax, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(grid->z(result[f].kmax));
    }
    INFO0(name_zmax, msg.str());
    msg.str("");
    msg << timeprefix;
    for (size_t f = 0; f < result.size(); ++f) {
        msg << ' ' << fullprec<>(result[f].fneg);
    }
    INFO0(name_fneg, msg.str());
    msg.str("");
}

void
driver_base::log_boundary_conditions(
        const std::string& timeprefix)
{
    // Only continue on the rank housing the zero-zero modes...
    //...which currently must also be rank zero for logging purposes.
    if (dgrid->has_zero_zero_modes()) {
        if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) != 0)
            SUZERAIN_ERROR_VOID_UNIMPLEMENTED();
    } else {
        return;
    }

    // Named loggers for lower/upper state at the 0th, 1st, and 2nd derivatives
    static const char * nick[2][3] = {
        { "bc.lower.d0", "bc.lower.d1", "bc.lower.d2" },
        { "bc.upper.d0", "bc.upper.d1", "bc.upper.d2" }
    };

    // Indices at the lower and upper walls
    size_t bc[2] = { 0, state_linear->shape()[1] - 1 };
    assert(SUZERAIN_COUNTOF(bc) == SUZERAIN_COUNTOF(nick));

    // Buffer for state evaluation and derivatives
    Array3Xc values(values.rows(), fields.size());
    assert(values.rows() == SUZERAIN_COUNTOF(nick[0]));  // Consistency

    for (size_t l = 0; l < SUZERAIN_COUNTOF(bc); ++l) {

        // Evaluate 0th, 1st, and 2nd derivatives of all state.  B-spline
        // recursion used instead of matrix operators for memory locality.
        SUZERAIN_ENSURE((unsigned) state_linear->strides()[1] == 1u);
        for (size_t k = 0; k < fields.size(); ++k) {
            b->linear_combination(values.rows() - 1,
                                  (*state_linear)[k].origin(),
                                  b->collocation_point(bc[l]),
                                  values.col(k).data(), 1u);
        }

        // Show mean boundary state on every invocation, possibly with headers
        for (size_t m = 0; m < (unsigned) values.rows(); ++m) {
            maybe_timeprefix_fields_identifiers(
                    *this, timeprefix, nick[l][m], header_shown[nick[l][m]]);
            std::ostringstream msg;
            msg << timeprefix;
            for (size_t k = 0; k < fields.size(); ++k) {
                using std::abs;
                msg << ' ' << fullprec<>(abs(values(m,k)));
            }
            INFO0(nick[l][m], msg.str());
        }
    }
}

// A building block called from both log_{channel,boundary_layer}_quantities
template< class Local >
void log_quantities_local_helper(
        const driver_base&  driver,
        const std::string&  timeprefix,
        const Local * const local,
        const char  * const name,
        bool&               header_shown)
{
    using std::setw;

    if (name) {
        std::ostringstream msg;
        if (!header_shown) {
            header_shown = true;
            msg << setw(timeprefix.size())
                << driver.build_timeprefix_description()
                << ' ' << setw(fullprec<>::width) << "a"
                << ' ' << setw(fullprec<>::width) << "gamma"
                << ' ' << setw(fullprec<>::width) << "mu"
                << ' ' << setw(fullprec<>::width) << "Pr"
                << ' ' << setw(fullprec<>::width) << "rho"
                << ' ' << setw(fullprec<>::width) << "T"
                << ' ' << setw(fullprec<>::width) << "u"
                << ' ' << setw(fullprec<>::width) << "v";
            INFO0(name, msg.str());
            msg.str("");
        }
        msg << timeprefix
            << ' ' << fullprec<>(local->a)
            << ' ' << fullprec<>(local->gamma)
            << ' ' << fullprec<>(local->mu)
            << ' ' << fullprec<>(local->Pr)
            << ' ' << fullprec<>(local->rho)
            << ' ' << fullprec<>(local->T)
            << ' ' << fullprec<>(local->u)
            << ' ' << fullprec<>(local->v);
        INFO0(name, msg.str());
    }
}

void driver_base::log_quantities_boundary_layer(
        const std::string& timeprefix,
        const suzerain_bl_local       * const wall,
        const suzerain_bl_viscous     * const viscous,
        const suzerain_bl_thicknesses * const thick,
        const suzerain_bl_local       * const edge,
        const suzerain_bl_local       * const edge99,
        const suzerain_bl_reynolds    * const reynolds,
        const suzerain_bl_qoi         * const qoi,
        const suzerain_bl_pg          * const pg,
        const char * const name_wall,
        const char * const name_visc,
        const char * const name_thick,
        const char * const name_edge,
        const char * const name_edge99,
        const char * const name_Re,
        const char * const name_qoi,
        const char * const name_pg)
{
    using std::setw;

    log_quantities_local_helper(*this, timeprefix, wall, name_wall,
                                header_shown[name_wall]);

    if (const char * const name = name_visc) {   // Evil, but helps avoid typos
        std::ostringstream msg;
        bool& log_header_shown = header_shown[name];
        if (!log_header_shown) {
            log_header_shown = true;
            msg << setw(timeprefix.size()) << build_timeprefix_description()
                << ' ' << setw(fullprec<>::width) << "cf"
                << ' ' << setw(fullprec<>::width) << "delta_nu"
                << ' ' << setw(fullprec<>::width) << "tau_w"
                << ' ' << setw(fullprec<>::width) << "u_tau"
                << ' ' << setw(fullprec<>::width) << "Bq"
                << ' ' << setw(fullprec<>::width) << "v_wallplus";
            INFO0(name, msg.str());
            msg.str("");
        }
        msg << timeprefix
            << ' ' << fullprec<>(qoi->cf)
            << ' ' << fullprec<>(viscous->delta_nu)
            << ' ' << fullprec<>(viscous->tau_w)
            << ' ' << fullprec<>(viscous->u_tau)
            << ' ' << fullprec<>(qoi->Bq)
            << ' ' << fullprec<>(qoi->v_wallplus);
        INFO0(name, msg.str());
    }

    if (const char * const name = name_thick) {
        std::ostringstream msg;
        bool& log_header_shown = header_shown[name];
        if (!log_header_shown) {
            log_header_shown = true;
            msg << setw(timeprefix.size()) << build_timeprefix_description()
                << ' ' << setw(fullprec<>::width) << "delta"
                << ' ' << setw(fullprec<>::width) << "delta1"
                << ' ' << setw(fullprec<>::width) << "delta2"
                << ' ' << setw(fullprec<>::width) << "delta3"
                << ' ' << setw(fullprec<>::width) << "deltaH0"
                << ' ' << setw(fullprec<>::width) << "delta99";
            INFO0(name, msg.str());
            msg.str("");
        }
        msg << timeprefix
            << ' ' << fullprec<>(thick->delta)
            << ' ' << fullprec<>(thick->delta1)
            << ' ' << fullprec<>(thick->delta2)
            << ' ' << fullprec<>(thick->delta3)
            << ' ' << fullprec<>(thick->deltaH0)
            << ' ' << fullprec<>(thick->delta99);
        INFO0(name, msg.str());
    }

    log_quantities_local_helper(*this, timeprefix, edge, name_edge,
                                header_shown[name_edge]);

    log_quantities_local_helper(*this, timeprefix, edge99, name_edge99,
                                header_shown[name_edge99]);

    if (const char * const name = name_Re) {
        std::ostringstream msg;
        bool& log_header_shown = header_shown[name];
        if (!log_header_shown) {
            log_header_shown = true;
            msg << setw(timeprefix.size()) << build_timeprefix_description()
                << ' ' << setw(fullprec<>::width) << "Re_delta"
                << ' ' << setw(fullprec<>::width) << "Re_delta1"
                << ' ' << setw(fullprec<>::width) << "Re_delta2"
                << ' ' << setw(fullprec<>::width) << "Re_delta3"
                << ' ' << setw(fullprec<>::width) << "Re_deltaH0"
                << ' ' << setw(fullprec<>::width) << "Re_delta99";
            INFO0(name, msg.str());
            msg.str("");
        }
        msg << timeprefix
            << ' ' << fullprec<>(reynolds->delta)
            << ' ' << fullprec<>(reynolds->delta1)
            << ' ' << fullprec<>(reynolds->delta2)
            << ' ' << fullprec<>(reynolds->delta3)
            << ' ' << fullprec<>(reynolds->deltaH0)
            << ' ' << fullprec<>(reynolds->delta99);
        INFO0(name, msg.str());
    }

    if (const char * const name = name_qoi) {
        std::ostringstream msg;
        bool& log_header_shown = header_shown[name];
        if (!log_header_shown) {
            log_header_shown = true;
            msg << setw(timeprefix.size()) << build_timeprefix_description()
                << ' ' << setw(fullprec<>::width) << "Ma_e"
                << ' ' << setw(fullprec<>::width) << "Ma_tau"
                << ' ' << setw(fullprec<>::width) << "ratio_rho"
                << ' ' << setw(fullprec<>::width) << "ratio_nu"
                << ' ' << setw(fullprec<>::width) << "ratio_T"
                << ' ' << setw(fullprec<>::width) << "Re_tau";
            INFO0(name, msg.str());
            msg.str("");
        }
        msg << timeprefix
            << ' ' << fullprec<>(qoi->Ma_e)
            << ' ' << fullprec<>(qoi->Ma_tau)
            << ' ' << fullprec<>(qoi->ratio_rho)
            << ' ' << fullprec<>(qoi->ratio_nu)
            << ' ' << fullprec<>(qoi->ratio_T)
            << ' ' << fullprec<>(thick->delta99 / viscous->delta_nu);
        INFO0(name, msg.str());
    }

    if (pg) {
        if (const char * const name = name_pg) {
            std::ostringstream msg;
            bool& log_header_shown = header_shown[name];
            if (!log_header_shown) {
                log_header_shown = true;
                msg << setw(timeprefix.size()) << build_timeprefix_description()
                    << ' ' << setw(fullprec<>::width) << "Clauser"
                    << ' ' << setw(fullprec<>::width) << "Lambda_n"
                    << ' ' << setw(fullprec<>::width) << "Launder_e"
                    << ' ' << setw(fullprec<>::width) << "Launder_w"
                    << ' ' << setw(fullprec<>::width) << "Pohlhausen"
                    << ' ' << setw(fullprec<>::width) << "p_ex";
                INFO0(name, msg.str());
                msg.str("");
            }
            msg << timeprefix
                << ' ' << fullprec<>(pg->Clauser)
                << ' ' << fullprec<>(pg->Lambda_n)
                << ' ' << fullprec<>(pg->Launder_e)
                << ' ' << fullprec<>(pg->Launder_w)
                << ' ' << fullprec<>(pg->Pohlhausen)
                << ' ' << fullprec<>(pg->p_ex);
            INFO0(name, msg.str());
        }
    }
}

void driver_base::log_quantities_channel(
        const std::string& timeprefix,
        const suzerain_channel_local   * const wall,
        const suzerain_channel_viscous * const viscous,
        const suzerain_channel_local   * const center,
        const suzerain_channel_qoi     * const qoi,
        const char * const name_wall,
        const char * const name_visc,
        const char * const name_center,
        const char * const name_qoi)
{
    using std::setw;          // Brevity

    log_quantities_local_helper(*this, timeprefix, wall, name_wall,
                                header_shown[name_wall]);

    if (const char * const name = name_visc) {
        std::ostringstream msg;
        bool& log_header_shown = header_shown[name];
        if (!log_header_shown) {
            log_header_shown = true;
            msg << setw(timeprefix.size()) << build_timeprefix_description()
                << ' ' << setw(fullprec<>::width) << "cf"
                << ' ' << setw(fullprec<>::width) << "delta_nu"
                << ' ' << setw(fullprec<>::width) << "tau_w"
                << ' ' << setw(fullprec<>::width) << "u_tau"
                << ' ' << setw(fullprec<>::width) << "Bq"
                << ' ' << setw(fullprec<>::width) << "v_wallplus";
            INFO0(name, msg.str());
            msg.str("");
        }
        msg << timeprefix
            << ' ' << fullprec<>(qoi->cf)
            << ' ' << fullprec<>(viscous->delta_nu)
            << ' ' << fullprec<>(viscous->tau_w)
            << ' ' << fullprec<>(viscous->u_tau)
            << ' ' << fullprec<>(qoi->Bq)
            << ' ' << fullprec<>(qoi->v_wallplus);
        INFO0(name, msg.str());
    }

    log_quantities_local_helper(*this, timeprefix, center, name_center,
                                header_shown[name_center]);

    if (const char * const name = name_qoi) {
        std::ostringstream msg;
        bool& log_header_shown = header_shown[name];
        if (!log_header_shown) {
            log_header_shown = true;
            msg << setw(timeprefix.size()) << build_timeprefix_description()
                << ' ' << setw(fullprec<>::width) << "Ma_c"
                << ' ' << setw(fullprec<>::width) << "Ma_tau"
                << ' ' << setw(fullprec<>::width) << "Pr_w"
                << ' ' << setw(fullprec<>::width) << "Re_c"
                << ' ' << setw(fullprec<>::width) << "Re_tau";
            INFO0(name, msg.str());
            msg.str("");
        }
        msg << timeprefix
            << ' ' << fullprec<>(qoi->Ma_c)
            << ' ' << fullprec<>(qoi->Ma_tau)
            << ' ' << fullprec<>(qoi->Pr_w)
            << ' ' << fullprec<>(qoi->Re_c)
            << ' ' << fullprec<>(qoi->Re_tau);
        INFO0(name, msg.str());
    }
}

void
driver_base::save_metadata()
{
    SUZERAIN_TIMER_SCOPED("driver_base::save_metadata");

    DEBUG0(who, "Saving metadata temporary file: " << restartdef->metadata);

    esio_handle esioh = esio_handle_initialize(MPI_COMM_WORLD);
    esio_file_create(esioh, restartdef->metadata.c_str(), 1 /* overwrite */);

    save_metadata(esioh);

    esio_file_close(esioh);
    esio_handle_finalize(esioh);

    metadata_saved = true;
}

void
driver_base::save_metadata(
        const esio_handle esioh)
{
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

    // Save metadata known by us
    save_grid_and_operators(esioh);
    if (timedef) timedef->save(esioh);

    // Invoke subclass extension point
    return save_metadata_hook(esioh);
}

void
driver_base::load_metadata(
        const esio_handle esioh)
{
    SUZERAIN_TIMER_SCOPED("driver_base::load_metadata");

    load_grid_and_operators(esioh);

    if (!timedef) {
        timedef = make_shared<definition_time>(
                std::numeric_limits<real_t>::quiet_NaN());
    }
    timedef->load(esioh);

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
    state_nonlinear->assign_from(*state_linear);
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
    state_nonlinear->assign_from(*state_linear);
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
    load_state_hook(esioh);                       // Invoke extension point
    state_linear->assign_from(*state_nonlinear);  // Copy into state_linear
    load_statistics_hook(esioh);                  // Invoke extension point

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
    SUZERAIN_UNUSED(t);
    SUZERAIN_UNUSED(nt);

    // Log information about physics-independent quantities of interest
    log_state_bulk(timeprefix);
    log_state_L2(timeprefix);
    log_state_extrema(timeprefix);
    log_boundary_conditions(timeprefix);

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

    // FIXME: choose between computed and cached extrema stats #3071
    state_nonlinear->assign_from(*state_linear);
    extrema = 
        compute_field_extrema_xz(*state_nonlinear, *grid, *dgrid, *cop);

    save_extrema(esioh, fields, extrema, *grid);

    return true;
}

void
driver_base::load_statistics_hook(
        const esio_handle esioh)
{
    SUZERAIN_UNUSED(esioh);

    // For example:
    //     perfect::load(esioh, samples);
    //
    // FIXME: Add load_extrema method #3071
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

    {
        signal::received_type::value_type signum; // Assigned in "if"
        if (    (signum = signal_received[signal::teardown_reactive])
            ||  (signum = signal_received[signal::halt_reactive    ])) {
            const char * const name = suzerain_signal_name(signum);
            INFO0(log_signal, "Tearing down due to receipt of " << name);
            received_teardown = signal_received[signal::teardown_reactive];
            received_halt     = signal_received[signal::halt_reactive    ];
            keep_advancing    = false;
            switch (signum) {
            case SIGINT:
            case SIGTERM:
                INFO0(log_signal,
                      "Receipt of another " << name <<
                      " will forcibly terminate program");
                signal2(signum, SIG_DFL);
                break;
            }
        }
    }

    if (signal_received[signal::teardown_proactive]) {
        INFO0(log_signal,
              "Initiating proactive teardown because of wall time constraint");
        received_teardown = true;
        keep_advancing    = false;
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

    // Output status one-sixteenth as often as we write statistics
    if (boost::math::isnormal(tmp_dt)) {
        dt = tmp_dt / 16;
    }
    if (tmp_nt) {
        using std::min;
        nt = min(tmp_nt / 16, static_cast<step_type>(1));
    }
}

delta_t_allreducer::delta_t_allreducer(
        const double& wtime_mpi_init,
        const double& wtime_fftw_planning,
        const shared_ptr<definition_time>& timedef,
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
    , who("reduce")
{
}

real_t delta_t_allreducer::operator()(
        const std::vector<real_t>& delta_t_candidates)
{
    // Perform one-time lookup of named logger, shadowing this->who.
    logging::logger_type who = logging::get_logger(this->who);

    // Squawk informatively on any non-positive delta_t_candidates values
    {
        const std::size_t N = delta_t_candidates.size();
        for (std::size_t i = 0; i < N; ++i) {
            if (SUZERAIN_LIKELY(delta_t_candidates[i] > 0)) {
                // NOP:  Written weirdly to make NaNs trigger the else branch
            } else {
                const int rank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
                WARN(who, "Non-positive delta_t_candidates[" << i << "] = "
                           << delta_t_candidates[i] << " on rank " << rank);
            }
        }
    }

    // Copy incoming candidates so we may mutate them
    std::vector<real_t> candidates(delta_t_candidates);

    // Take atomic snapshot of and then clear signal::global_received
    signal_received = signal::global_received;
    static const signal::volatile_received_type::value_type zero = 0;
    signal::global_received.assign(zero);
    if (DEBUG_ENABLED(who)) {
        for (std::size_t i = 0; i < signal::received_type::static_size; ++i) {
            if (SUZERAIN_UNLIKELY(signal_received[i])) {
                const int rank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
                DEBUG(who, "Received signal number " << signal_received[i]
                           << " on rank " << rank);
            }
        }
    }

    // When possible, obtain time step period statistics and a projected
    // wall time for when we could complete the next time step and dump a
    // restart file.  If the projection is after --advance_wt, register
    // teardown.  Logic here is key to proactive soft teardown success.
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
        if (DEBUG_ENABLED(who) && signal_received[signal::teardown_proactive]) {
            const int rank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
            DEBUG(who, "Rank " << rank << " projects delaying teardown until "
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

} // end namespace support

} // end namespace suzerain
