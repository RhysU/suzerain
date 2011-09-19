//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
//
// Please see http://pecos.ices.utexas.edu for more information.
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
//
// channel_explicit.cpp: A fully explicit channel calculation using Suzerain
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#define EIGEN_DEFAULT_IO_FORMAT \
        Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ";", "", "", "[", "]")

#include <suzerain/common.hpp>
#include <gsl/gsl_errno.h>
#include <esio/esio.h>
#include <esio/error.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/error.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/os.h>
#include <suzerain/pencil.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/signal_definition.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/utility.hpp>

#include "logging.hpp"
#include "precision.hpp"
#include "channel.hpp"
#include "explicit_op.hpp"

#pragma warning(disable:383 1572)

using boost::array;
using boost::make_shared;
using boost::numeric_cast;
using boost::scoped_ptr;
using boost::shared_ptr;
using std::numeric_limits;

// Explicit timestepping scheme uses only complex_t 4D ContiguousState
// State indices range over (scalar field, Y, X, Z) in wave space
typedef suzerain::ContiguousState<4,complex_t> state_type;

// Global scenario parameters initialized in main().  These are declared const
// to avoid accidental modification but have their const-ness const_cast away
// where necessary to load settings.
using suzerain::problem::ScenarioDefinition;
using suzerain::problem::GridDefinition;
using suzerain::problem::RestartDefinition;
using suzerain::problem::TimeDefinition;
using channel::NoiseDefinition;
using suzerain::problem::SignalDefinition;
static const ScenarioDefinition<real_t> scenario;
static const GridDefinition grid;
static const RestartDefinition restart(
        /* metadata     */ "metadata.h5.XXXXXX",
        /* uncommitted  */ "uncommitted.h5.XXXXXX",
        /* desttemplate */ "restart#.h5",
        /* retain       */ 1,
        /* dt           */ 0,
        /* nt           */ 0);
static const TimeDefinition<real_t> timedef(
        /* advance_dt                */ 0,
        /* advance_nt                */ 0,
        /* advance_wt                */ 0,
        /* status_dt                 */ 0,
        /* status_nt                 */ 0,
        /* min_dt                    */ 1e-8,
        /* max_dt                    */ 0,
        /* evmagfactor per Venugopal */ 0.72);
static const NoiseDefinition  noisedef;
static const SignalDefinition sigdef;

// Global details initialized in main()
static shared_ptr<      suzerain::bspline>              b;
static shared_ptr<      suzerain::bsplineop>            bop;    // Collocation
static shared_ptr<      suzerain::bsplineop>            gop;    // Galerkin L2
static shared_ptr<      suzerain::bsplineop_luz>        bopluz;
static shared_ptr<const suzerain::pencil_grid>          dgrid;
static shared_ptr<      channel::manufactured_solution> msoln;

// State details specific to this rank initialized in main()
static shared_ptr<state_type> state_linear;
static shared_ptr<state_type> state_nonlinear;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void)
{
    if (esioh) esio_handle_finalize(esioh);
}

/**
 * <tt>atexit</tt> callback to remove the metadata file.  Do \i NOT remove
 * restart.uncommitted as it may help post mortem debugging.
 */
static void atexit_metadata(void) {
    if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0) {
        if (0 == unlink(restart.metadata.c_str())) {
            DEBUG("Cleaned up temporary file " << restart.metadata);
        } else {
            WARN("Error cleaning up temporary file " << restart.metadata);
        }
    }
}

// Formatting parameters used within append_real()
static const int append_real_prec  = numeric_limits<real_t>::digits10 * 0.75;
static const int append_real_width = append_real_prec + 5;

/** Provides nice formatting of real-valued quantities for status lines */
template<class CharT, class Traits, class Number>
static std::basic_ostream<CharT,Traits>& append_real(
        std::basic_ostream<CharT,Traits>& os,
        Number value)
{
    // Magic "2" is the width of a sign and a decimal point
    static const real_t fixedmax
        = std::pow(real_t(10), append_real_width - append_real_prec - 2);

    // Magic 3 is the width of a sign, leading zero, and decimal point
    static const real_t fixedmin
        = std::pow(real_t(10), -(append_real_width - append_real_prec - 3));

    // Format in fixed or scientific form as appropriate in given width
    // Care taken to not perturb observable ostream state after function call
    std::ios::fmtflags savedflags;
    std::streamsize savedprec;
    if (value >= fixedmin && value <= fixedmax) {
        savedflags = os.setf(std::ios::fixed      | std::ios::right,
                             std::ios::floatfield | std::ios::adjustfield);
        savedprec = os.precision(append_real_prec);
    } else {
        savedflags = os.setf(std::ios::scientific | std::ios::right,
                             std::ios::floatfield | std::ios::adjustfield);
        savedprec = os.precision(append_real_width - 9);
    }
    os << std::setw(append_real_width) << static_cast<real_t>(value);
    os.precision(savedprec);
    os.setf(savedflags);

    return os;
}

/** Log messages containing mean and fluctuating L2 information */
static void information_L2(const std::string& timeprefix)
{
    // Avoid computational cost when logging is disabled
    logging::logger_type L2_mean  = logging::get_logger("L2.mean");
    logging::logger_type L2_fluct = logging::get_logger("L2.fluct");
    if (!INFO0_ENABLED(L2_mean) && !INFO0_ENABLED(L2_fluct)) return;

    // Collective computation of the L_2 norms
    const array<channel::L2,channel::field::count> L2
        = channel::field_L2(*state_linear, scenario, grid, *dgrid, *gop);

    // Build and log L2 of mean conserved state
    std::ostringstream msg;
    msg << timeprefix;
    for (std::size_t k = 0; k < L2.size(); ++k) {
        append_real(msg << ' ', L2[k].mean());
    }
    INFO0(L2_mean, msg.str());

    // Build and log L2 of fluctuating conserved state
    msg.str("");
    msg << timeprefix;
    for (std::size_t k = 0; k < L2.size(); ++k) {
        append_real(msg << ' ', L2[k].fluctuating());
    }
    INFO0(L2_fluct, msg.str());
}

/** Build a message containing bulk quantities (intended for root rank only) */
static void information_bulk(const std::string& timeprefix)
{
    // Avoid computational cost when logging is disabled
    logging::logger_type bulk_state = logging::get_logger("bulk.state");
    if (!INFO0_ENABLED(bulk_state)) return;

    // Compute operator for finding bulk quantities from coefficients
    Eigen::VectorXr bulkcoeff(b->n());
    b->integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= scenario.Ly;

    // Prepare the status message and log it
    std::ostringstream msg;
    msg << timeprefix;
    for (std::size_t k = 0; k < state_linear->shape()[0]; ++k) {
        Eigen::Map<Eigen::VectorXc> mean(
                (*state_linear)[k].origin(), state_linear->shape()[1]);
        append_real(msg << ' ', bulkcoeff.dot(mean.real()));
    }
    INFO0(bulk_state, msg.str());
}

/** Build a message containing specific state quantities at the wall */
static void information_specific_wall_state(const std::string& timeprefix)
{
    namespace ndx = channel::field::ndx;

    logging::logger_type nick[2] = { logging::get_logger("wall.lower"),
                                     logging::get_logger("wall.upper")  };

    // Indices at the lower and upper walls.  Use that wall collocation point
    // values are nothing but the first and last B-spline coefficient values.
    std::size_t wall[2] = { 0, state_linear->shape()[1] - 1 };

    // Message lists rho, u, v, w, and total energy at walls
    for (std::size_t l = 0; l < sizeof(wall)/sizeof(wall[0]); ++l) {

        // Avoid computational cost when logging is disabled
        if (!DEBUG0_ENABLED(nick[l])) continue;

        std::ostringstream msg;
        msg << timeprefix;

        const real_t rho = ((*state_linear)[ndx::rho][wall[l]][0][0]).real();
        append_real(msg << ' ', rho);
        assert(ndx::rho == 0);
        for (std::size_t k = 1; k < channel::field::count; ++k) {
            append_real(msg << ' ' ,
                        ((*state_linear)[k][wall[l]][0][0]).real() / rho);
        }
        DEBUG0(nick[l], msg.str());
    }
}

/**
 * Build a message for the absolute error versus a manufactured solution.
 * Uses state_nonlinear as a scratch space and is not cheap.
 */
static void information_manufactured_solution_absolute_error(
        const std::string& timeprefix,
        const real_t simulation_time)
{
    // Avoid computational cost when logging is disabled
    logging::logger_type mms_abserr = logging::get_logger("mms.abserr");
    if (!INFO0_ENABLED(mms_abserr)) return;

    // Compute L2 of error of state against manufactured solution
    assert(msoln);
    state_nonlinear->assign(*state_linear);
    channel::accumulate_manufactured_solution(
            1, *msoln, -1, *state_nonlinear,
            scenario, grid, *dgrid, *b, *bop, simulation_time);
    const array<channel::L2,channel::field::count> L2
        = channel::field_L2(*state_nonlinear, scenario, grid, *dgrid, *gop);

    // Output absolute global errors for each field
    std::ostringstream msg;
    msg << timeprefix;
    for (std::size_t k = 0; k < channel::field::count; ++k) {
        append_real(msg << ' ', L2[k].total());
    }
    INFO0(mms_abserr, msg.str());
}

/** Tracks last time we output a status line */
static std::size_t last_status_nt = numeric_limits<std::size_t>::max();

/** Routine to output status.  Signature for TimeController use. */
static bool log_status(real_t t, std::size_t nt)
{
    // Notice collective operations are never inside logging macros!

    using std::max;
    using std::floor;
    using std::log10;

    // Defensively avoid multiple invocations with no intervening changes
    if (last_status_nt == nt) {
        DEBUG0("Cowardly refusing to repeatedly show status at nt = " << nt);
        return true;
    }

    // Build time- and timestep-specific status prefix.
    // Precision computations ensure multiple status lines minimally distinct
    std::ostringstream oss;
    real_t np = 0;
    if (timedef.status_dt > 0) {
        np = max(np, -floor(log10(timedef.status_dt)));
    }
    if (timedef.status_nt > 0) {
        np = max(np, -floor(log10(timedef.min_dt * timedef.status_nt)) + 1);
    }
    if (np > 0) {
        oss.setf(std::ios::fixed, std::ios::floatfield);
        const std::streamsize oldprec = oss.precision(np);
        oss << t;
        oss.precision(oldprec);
        oss.unsetf(std::ios::fixed);
    } else {
        oss << t;
    }
    oss << ' ' << std::setw(7) << nt;
    const std::string timeprefix = oss.str();

    // Log information about the various quantities of interest
    information_bulk(timeprefix);
    information_L2(timeprefix);
    information_specific_wall_state(timeprefix);

    // Log errors versus any manufactured solution in use
    if (msoln) {
        information_manufactured_solution_absolute_error(timeprefix, t);
    }

    last_status_nt = nt; // Maintain last status time step

    return true;
}

/** Tracks last time a restart file was written successfully */
static std::size_t last_restart_saved_nt = numeric_limits<std::size_t>::max();

/** Routine to store a restart file.  Signature for TimeController use. */
static bool save_restart(real_t t, std::size_t nt)
{
    // Defensively avoid multiple invocations with no intervening changes
    if (last_restart_saved_nt == nt) {
        DEBUG0("Cowardly refusing to save multiple restarts at nt = " << nt);
        return true;
    }

    const double starttime = MPI_Wtime();
    DEBUG0("Started to store restart at t = " << t << " and nt = " << nt);

    DEBUG0("Cloning " << restart.metadata << " to " << restart.uncommitted);
    esio_file_clone(esioh, restart.metadata.c_str(),
                    restart.uncommitted.c_str(), 1 /*overwrite*/);
    channel::store_time(esioh, t);

    if (restart.physical) {
        DEBUG0("Storing primitive collocation point values into "
               << restart.uncommitted);
        channel::store_collocation_values(
                esioh, *state_linear, *state_nonlinear,
                scenario, grid, *dgrid, *b, *bop);
    } else {
        DEBUG0("Storing conserved coefficients into " << restart.uncommitted);
        channel::store_coefficients(
                esioh, *state_linear, *state_nonlinear, scenario, grid, *dgrid);
    }

    DEBUG0("Committing " << restart.uncommitted
           << " as a restart file using template " << restart.desttemplate);
    esio_file_close_restart(
            esioh, restart.desttemplate.c_str(), restart.retain);

    const double elapsed = MPI_Wtime() - starttime;
    INFO0("Successfully wrote restart at t = " << t << " for nt = " << nt
          << " in " << elapsed << " seconds");

    last_restart_saved_nt = nt; // Maintain last successful restart time step

    return true; // Continue time advancement
}

/**
 * Type of atomic locations used to track local receipt of the following
 * signal-based actions:
 *
 * \li \c 0 Output a status message
 * \li \c 1 Write a restart file
 * \li \c 2 Tear down the simulation (reactively  due to an incoming signal)
 * \li \c 3 Tear down the simulation (proactively due to --advance_wt limit)
 */
typedef array<sig_atomic_t,4> atomic_signal_received_t;

/** Atomic locations used to track local signal receipt. */
static atomic_signal_received_t atomic_signal_received = {{ /*0*/ }};

/** Signal handler which mutates \c atomic_signal_received. */
static void process_signal(int sig)
{

    // Strictly speaking this handler performs too much work.  The design
    // choice was to have this extra work done on the (rare) signal receipt
    // rather than on the (frequent) polling of signal receipt status.

    std::vector<int>::const_iterator end;

    // Determine if we should output status due to the signal
    end = sigdef.status.end();
    if (std::find(sigdef.status.begin(), end, sig) != end) {
        atomic_signal_received[0] = sig;
    }

    // Determine if we should write a restart due to the signal
    end = sigdef.restart.end();
    if (std::find(sigdef.restart.begin(), end, sig) != end) {
        atomic_signal_received[1] = sig;
    }

    // Determine if we should tear down the simulation due to the signal
    end = sigdef.teardown.end();
    if (std::find(sigdef.teardown.begin(), end, sig) != end) {
        atomic_signal_received[2] = sig;
    }
}

/**
 * Type of non-atomic locations used to track global receipt of the
 * same actions as atomic_signal_received;
 */
typedef array<int,atomic_signal_received_t::static_size> signal_received_t;

/** Non-atomic locations used to track global signal receipt. */
static signal_received_t signal_received = {{ /*0*/ }};

/** Flag used to indicate early time advancement stopping is legit. */
static bool soft_teardown = false;

/**
 * Routine to check for incoming signals on any rank.
 * Signature for TimeController use.
 */
static bool process_any_signals_received(real_t t, std::size_t nt)
{
    // DeltaTAllreducer performs the Allreduce necessary to get local status
    // from atomic_signal_received into global status in signal_received.

    // Keep advancing time unless keep_advancing is set false
    bool keep_advancing = true;

    if (signal_received[0]) {
        INFO0("Outputting simulation status due to receipt of "
              << suzerain_signal_name(signal_received[0]));
        keep_advancing = keep_advancing && log_status(t, nt);
    }

    if (signal_received[1]) {
        INFO0("Writing restart file due to receipt of "
              << suzerain_signal_name(signal_received[1]));
        keep_advancing = keep_advancing && save_restart(t, nt);
    }

    if (signal_received[2]) {
        const char * const name = suzerain_signal_name(signal_received[2]);
        INFO0("Initiating teardown due to receipt of " << name);
        soft_teardown  = true;
        keep_advancing = false;
        switch (signal_received[2]) {
            case SIGINT:
            case SIGTERM:
                INFO0("Receipt of another " << name <<
                      " will forcibly terminate program");
                signal(signal_received[2], SIG_DFL);
                break;
        }
    }

    if (signal_received[3]) {
        INFO0("Initiating proactive teardown because of wall time constraint");
        soft_teardown  = true;
        keep_advancing = false;
    }

    // Clear signal_received to defensively avoid stale data bugs.
    // These would only be problematic if process_any_signals_received
    // was run multiple times in between DeltaTAllreducer invocations.
    signal_received.assign(0);

    return keep_advancing;
}

/** Wall time at which MPI_Init completed */
static double wtime_mpi_init;

/** Wall time elapsed during loading of state from the restart file */
static double wtime_load_state;

/** Wall time at which we began time stepping */
static double wtime_advance_start;

/**
 * A stateful functor that performs an MPI Allreduce to determine the minimum
 * stable time step size across all ranks.  The same MPI Allreduce is used to
 * hide the cost of querying atomic_signal_received across all ranks.
 */
static class DeltaTAllreducer {

private:

    // Maintains coarse statistics on the wall time duration between calls.
    // This provides a rank-specific measure of the variability per time step
    // which includes all registered periodic callbacks (e.g. status and
    // restart processing).
    boost::accumulators::accumulator_set<
            real_t,
            boost::accumulators::stats<boost::accumulators::tag::max,
                                       boost::accumulators::tag::mean,
                                       boost::accumulators::tag::variance>
        > period;

public:

    real_t operator()(const std::vector<real_t>& delta_t_candidates) {

        // Copy incoming candidates so we may mutate them
        std::vector<real_t> candidates(delta_t_candidates);

        // Take atomic snapshop of and then clear atomic_signal_received
        signal_received = atomic_signal_received;
        atomic_signal_received.assign(0);

        if (DEBUG_ENABLED()) {
            const int rank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
            for (int i = 0; i < signal_received_t::static_size; ++i) {
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
        if (timedef.advance_wt > 0 && (boost::math::isfinite)(wtime_last)) {

            // Accumulate time step period statistics
            const double wtime = MPI_Wtime();
            period(wtime - wtime_last);

            // Find a pessimistic time for the next time step completion...
            // (that is, finish current step *and* finish another one)
            namespace acc = boost::accumulators;
            double wtime_projected = wtime
                + 2*(acc::mean(period) + 3*std::sqrt(acc::variance(period)));
            // ...to which we add a pessimistic estimate for dumping a restart
            if (last_restart_saved_nt == numeric_limits<std::size_t>::max()) {
                wtime_projected += 2*wtime_load_state;  // Load as surrogate
            } else {
                wtime_projected += (acc::max)(period);  // Includes dumps
            }
            // ...to which we add an estimate of other finalization costs
            wtime_projected += 2*(wtime_advance_start - wtime_mpi_init);

            // Raise a "signal" if we suspect we cannot teardown quickly enough
            signal_received[3] = (wtime_projected >=   wtime_mpi_init
                                                     + timedef.advance_wt);
            if (DEBUG_ENABLED() && signal_received[3]) {
                const int rank = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
                DEBUG("Rank " << rank << " projects delaying teardown until "
                      << wtime_projected - wtime_mpi_init
                      << " elapsed seconds is dangerous.");
            }
        }
        wtime_last = MPI_Wtime();

        // Push a negated version of signal_received onto end of candidates
        candidates.reserve(candidates.size() + signal_received_t::static_size);
        std::transform(signal_received.begin(), signal_received.end(),
                       std::back_inserter(candidates),
                       std::negate<signal_received_t::value_type>());

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
                                          + signal_received_t::static_size,
                       signal_received.begin(),
                       std::negate<signal_received_t::value_type>());
        candidates.erase(candidates.begin() + delta_t_candidates.size(),
                         candidates.begin() + delta_t_candidates.size()
                                            + signal_received_t::static_size);

        // Delegate remaining work on each rank to the usual DeltaTReducer
        return suzerain::timestepper::DeltaTReducer()(candidates);
    }

} delta_t_allreducer;

/**
 * Default log4cxx configuration (differs from channel::log4cxx_config).
 * <tt>${FOO}</tt> syntax may be used to pick up environment variables in addition to properties
 * See http://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PatternLayout.html
 * and https://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PropertyConfigurator.html
 */
static const char log4cxx_config[] =
    "## Set root logger level (e.g. TRACE, DEBUG, INFO) and one or more root appenders\n"
    "# Output INFO or higher messages on CONSOLE and in LOG\n"
    "log4j.rootLogger=INFO, CONSOLE, LOG\n"
    "\n"
    "## Configure output to CONSOLE\n"
    "log4j.appender.CONSOLE=org.apache.log4j.ConsoleAppender\n"
    "log4j.appender.CONSOLE.layout=org.apache.log4j.PatternLayout\n"
    "log4j.appender.CONSOLE.layout.ConversionPattern=%-5p %8r %-10c %m%n\n"
    "\n"
    "## Configure output to LOG file with formatting mimicking CONSOLE\n"
    "log4j.appender.LOG=org.apache.log4j.FileAppender\n"
    "log4j.appender.LOG.append=true\n"
    "log4j.appender.LOG.filename=log.dat\n"
    "log4j.appender.LOG.layout=${log4j.appender.CONSOLE.layout}\n"
    "log4j.appender.LOG.layout.ConversionPattern=${log4j.appender.CONSOLE.layout.ConversionPattern}\n"
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
    "## Collect \"L2.fluct\" messages into L2.fluct.dat mimicking LOG file behavior\n"
    "log4j.logger.L2.fluct=INHERITED, L2FLUCT\n"
    "log4j.appender.L2FLUCT=${log4j.appender.LOG}\n"
    "log4j.appender.L2FLUCT.filename=L2.fluct.dat\n"
    "log4j.appender.L2FLUCT.append=${log4j.appender.LOG.append}\n"
    "log4j.appender.L2FLUCT.layout=${log4j.appender.LOG.layout}\n"
    "log4j.appender.L2FLUCT.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
;

/** Main driver logic */
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                          // Initialize MPI
    wtime_mpi_init = MPI_Wtime();                    // Record MPI_Init time
    atexit((void (*) ()) MPI_Finalize);              // Finalize MPI at exit
    logging::initialize(MPI_COMM_WORLD,              // Initialize logging
                        log4cxx_config);
    esioh = esio_handle_initialize(MPI_COMM_WORLD);  // Initialize ESIO
    atexit(&atexit_esio);                            // Finalize ESIO

    // Record invocation for posterity and to aid in debugging
    {
        std::ostringstream os;
        std::copy(argv, argv+argc, std::ostream_iterator<const char *>(os," "));
        INFO0("Invocation: " << os.str());
    }

    // Hook error handling into logging infrastructure
    gsl_set_error_handler(
            &channel::mpi_abort_on_error_handler_gsl);
    suzerain_set_error_handler(
            &channel::mpi_abort_on_error_handler_suzerain);
    esio_set_error_handler(
            &channel::mpi_abort_on_error_handler_esio);

    DEBUG0("Processing command line arguments and response files");
    std::string restart_file;
    bool default_advance_nt;
    {
        suzerain::ProgramOptions options(
                "Suzerain-based explicit compressible channel simulation",
                "RESTART-FILE");
        options.add_definition(
                const_cast<ScenarioDefinition<real_t>&>(scenario));
        options.add_definition(
                const_cast<GridDefinition&>(grid));
        options.add_definition(
                const_cast<RestartDefinition&>(restart));
        options.add_definition(
                const_cast<TimeDefinition<real_t>&>(timedef));
        options.add_definition(
                const_cast<NoiseDefinition&>(noisedef));
        options.add_definition(
                const_cast<SignalDefinition&>(sigdef));

        std::vector<std::string> positional = options.process(argc, argv);

        if (positional.size() != 1) {
            FATAL0("Exactly one restart file name must be specified");
            return EXIT_FAILURE;
        }
        restart_file = positional[0];

        default_advance_nt = options.variables()["advance_nt"].defaulted();
    }

    INFO0("Loading details from restart file: " << restart_file);
    esio_file_open(esioh, restart_file.c_str(), 0 /* read-only */);
    channel::load(esioh, const_cast<ScenarioDefinition<real_t>&>(scenario));
    channel::load(esioh, const_cast<GridDefinition&>(grid));
    channel::load(esioh, scenario, msoln);
    esio_file_close(esioh);

    if (msoln) {
        INFO0("Restart file prescribes a manufactured solution");
        if (boost::math::isnormal(scenario.bulk_rho)) {
            WARN0("Manufactured solution incompatible with bulk_rho = "
                  << scenario.bulk_rho);
        }
        if (boost::math::isnormal(scenario.bulk_rhou)) {
            WARN0("Manufactured solution incompatible with bulk_rhou = "
                  << scenario.bulk_rhou);
        }
    }

    INFO0("Using B-splines of order " << (grid.k - 1)
          << " on [0, " << scenario.Ly << "] with "
          << grid.N.y() << " DOF stretched per htdelta " << grid.htdelta);
    channel::create(grid.N.y(), grid.k, 0.0, scenario.Ly,
                    grid.htdelta, b, bop);
    assert(b->k() == grid.k);
    assert(b->n() == grid.N.y());
    bopluz = make_shared<suzerain::bsplineop_luz>(*bop);
    bopluz->form_mass(*bop);
    gop.reset(new suzerain::bsplineop(*b, 0, SUZERAIN_BSPLINEOP_GALERKIN_L2));

    // Compute and display a couple of discretization quality metrics.
    {
        // Temporarily work in real-valued quantities as it is a bit simpler.
        suzerain::bsplineop_lu boplu(*bop);
        boplu.form_mass(*bop);

        // Compute and display discrete conservation error magnitude
        Eigen::MatrixXXr mat = Eigen::MatrixXXr::Identity(b->n(),b->n());
        boplu.solve(b->n(), mat.data(), 1, b->n());         // M^-1
        bop->apply(1, b->n(), 1.0, mat.data(), 1, b->n());  // D*M^-1
        boplu.solve(b->n(), mat.data(), 1, b->n());         // M^-1*D*M^-1
        Eigen::VectorXr vec(b->n());
        b->integration_coefficients(0, vec.data());
        vec = vec.transpose() * mat;                        // w^{T}*M^-1*D*M^-1
        vec.head<1>()[0] -= -1;                             // Exact head
        vec.tail<1>()[0] -=  1;                             // Exact tail
        double relerr = vec.norm() / std::sqrt(2);          // Exact 2-norm
        INFO0("B-spline discrete conservation relative error near "
              << relerr * 100 << "%");

        // Compute and display condition number
        double rcond;
        boplu.rcond(&rcond);
        INFO0("B-spline mass matrix has condition number near "
              << (1 / rcond));
    }

    // Generating unique file names as needed using mkstemp(3)
    {
        // Pack a temporary buffer with the three file name templates
        array<std::size_t,4> pos = {{
                0,
                restart.metadata.length()     + 1,
                restart.uncommitted.length()  + 1,
                restart.desttemplate.length() + 1
        }};
        std::partial_sum(pos.begin(), pos.end(), pos.begin());
        boost::scoped_array<char> buf(new char[pos[3]]);
        strcpy(&buf[pos[0]], restart.metadata.c_str());
        strcpy(&buf[pos[1]], restart.uncommitted.c_str());
        strcpy(&buf[pos[2]], restart.desttemplate.c_str());

        // Generate unique files to be overwritten and/or just file names.
        // File generation relies on template semantics of mkstemp(3).
        // Error checking kept minimal as failures here should not be fatal.
        // This is not particularly robust but it should serve our needs.
        if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0) {
            if (boost::ends_with(restart.metadata, "XXXXXX")) {
                close(mkstemp(&buf[pos[0]]));  // Clobbered later...
            }
            if (boost::ends_with(restart.uncommitted, "XXXXXX")) {
                close(mkstemp(&buf[pos[1]]));  // Possibly clobbered later...
                unlink(&buf[pos[1]]);          // ...so remove any evidence
            }
            if (boost::ends_with(restart.desttemplate, "XXXXXX")) {
                close(mkstemp(&buf[pos[2]]));  // Not clobbered later...
                unlink(&buf[pos[2]]);          // ...so remove any evidence
            }
        }

        // Broadcast any generated names to all ranks and unpack values
        SUZERAIN_MPICHKQ(MPI_Bcast(buf.get(), pos[3],
                         suzerain::mpi::datatype<char>(), 0,
                         MPI_COMM_WORLD));
        const_cast<RestartDefinition&>(restart).metadata     = &buf[pos[0]];
        const_cast<RestartDefinition&>(restart).uncommitted  = &buf[pos[1]];
        const_cast<RestartDefinition&>(restart).desttemplate = &buf[pos[2]];
    }

    DEBUG0("Saving metadata temporary file: " << restart.metadata);
    {
        atexit(&atexit_metadata); // Delete any lingering metadata file at exit
        esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
        esio_file_create(h, restart.metadata.c_str(), 1 /* overwrite */);
        channel::store(h, scenario);
        channel::store(h, grid, scenario.Lx, scenario.Lz);
        channel::store(h, b, bop, gop);
        channel::store(h, scenario, msoln);
        esio_file_close(h);
        esio_handle_finalize(h);
    }

    // Display global degree of freedom information
    INFO0("Global number of unknowns:         " << (  grid.N.prod()
                                                    * channel::field::count));
    INFO0("Grid degrees of freedom    (GDOF): " << grid.N.prod());
    INFO0("GDOF by direction           (XYZ): " << grid.N);
    INFO0("Dealiased GDOF by direction (XYZ): " << grid.dN);
    INFO0("Number of MPI ranks:               "
          << suzerain::mpi::comm_size(MPI_COMM_WORLD));

    // Initialize pencil_grid which handles P3DFFT setup/teardown RAII
    dgrid = make_shared<suzerain::pencil_grid>(grid.dN, grid.P);
    assert((grid.dN == dgrid->global_physical_extent).all());
    INFO0("Rank grid used for decomposition: " << dgrid->processor_grid);
    DEBUG("Local wave start      (XYZ): " << dgrid->local_wave_start);
    DEBUG("Local wave end        (XYZ): " << dgrid->local_wave_end);
    DEBUG("Local wave extent     (XYZ): " << dgrid->local_wave_extent);
    DEBUG("Local physical start  (XYZ): " << dgrid->local_physical_start);
    DEBUG("Local physical end    (XYZ): " << dgrid->local_physical_end);
    DEBUG("Local physical extent (XYZ): " << dgrid->local_physical_extent);

    // Create state storage for linear operator
    // TODO Have state_linear only store non-dealiased state
    state_linear = make_shared<state_type>(
            suzerain::to_yxz(channel::field::count, dgrid->local_wave_extent));

    // Load restart information into state_linear, including simulation time
    esio_file_open(esioh, restart_file.c_str(), 0 /* read-only */);
    real_t initial_t;
    channel::load_time(esioh, initial_t);
    {
        const double begin = MPI_Wtime();
        channel::load(esioh, *state_linear, scenario, grid, *dgrid, *b, *bop);
        wtime_load_state = MPI_Wtime() - begin;
    }
    esio_file_close(esioh);

    // If requested, add noise to the momentum fields at startup (expensive).
    channel::add_noise(*state_linear, noisedef,
                       scenario, grid, *dgrid, *b, *bop);

    // Create the state storage for nonlinear operator with appropriate padding
    state_nonlinear.reset(channel::allocate_padded_state<state_type>(
                channel::field::count, *dgrid));

    // Dump some state shape and stride information for debugging purposes
    DEBUG("Linear state shape      (FYXZ): "
          << suzerain::multi_array::shape_array(*state_linear));
    DEBUG("Nonlinear state shape   (FYXZ): "
          << suzerain::multi_array::shape_array(*state_nonlinear));
    DEBUG("Linear state strides    (FYXZ): "
          << suzerain::multi_array::strides_array(*state_linear));
    DEBUG("Nonlinear state strides (FYXZ): "
          << suzerain::multi_array::strides_array(*state_nonlinear));

    // Prepare generic timestepping handles
    using suzerain::timestepper::lowstorage::ILowStorageMethod;
    scoped_ptr<ILowStorageMethod<complex_t> > m;

    using suzerain::timestepper::lowstorage::ILinearOperator;
    scoped_ptr<ILinearOperator<state_type,state_type> > L;

    using suzerain::timestepper::INonlinearOperator;
    scoped_ptr<INonlinearOperator<state_type> > N;

    using suzerain::timestepper::TimeController;
    scoped_ptr<TimeController<real_t> > tc;

    // Prepare timestepping details specific to the chosen operators.
    // Nonlinear scaling factor (N_x N_z)^(-1) from write up section 2.1
    // (Spatial discretization) is modified for dealiasing and included here.
    m.reset(new suzerain::timestepper::lowstorage::SMR91Method<complex_t>(
                timedef.evmagfactor));
    L.reset(new channel::BsplineMassOperatorIsothermal(
                scenario, grid, *dgrid, *b, *bop));
    N.reset(new channel::NonlinearOperator(
                scenario, grid, *dgrid, *b, *bop, msoln));
    tc.reset(make_LowStorageTimeController(
                *m, delta_t_allreducer,
                *L, real_t(1)/(grid.dN.x()*grid.dN.z()), *N,
                *state_linear, *state_nonlinear,
                initial_t, timedef.min_dt, timedef.max_dt));

    // Register status callbacks status_{dt,nt} as requested.
    // When neither is provided, default and update timedef with new value.
    {
        TimeController<real_t>::time_type dt;
        if (timedef.status_dt) {
            dt = timedef.status_dt;
        } else if (timedef.advance_dt) {
            dt = timedef.advance_dt / restart.retain / 5;
            const_cast<real_t &>(timedef.status_dt) = dt;
        } else {
            dt = tc->forever_t();
        }

        TimeController<real_t>::step_type nt;
        if (timedef.status_nt) {
            nt = timedef.status_nt;
        } else if (timedef.advance_nt) {
            nt = timedef.advance_nt / restart.retain / 5;
            nt = std::max<TimeController<real_t>::step_type>(1, nt);
            const_cast<int &>(timedef.status_nt) = nt;
        } else {
            nt = tc->forever_nt();
        }

        tc->add_periodic_callback(dt, nt, &log_status);
    }

    // Register restart-writing callbacks restart_{dt,nt} as requested.
    // When neither is provided, default and update restart with new value.
    {
        TimeController<real_t>::time_type dt;
        if (restart.dt) {
            dt = restart.dt;
        } else if (timedef.advance_dt) {
            dt = timedef.advance_dt / restart.retain;
            const_cast<real_t &>(restart.dt) = dt;
        } else {
            dt = tc->forever_t();
        }

        TimeController<real_t>::step_type nt;
        if (restart.nt) {
            nt = restart.nt;
        } else if (timedef.advance_nt) {
            nt = timedef.advance_nt / restart.retain;
            nt = std::max<TimeController<real_t>::step_type>(1, nt);
            const_cast<int &>(restart.nt) = nt;
        } else {
            nt = tc->forever_nt();
        }

        tc->add_periodic_callback(dt, nt, &save_restart);
    }

    // Register any necessary signal handling logic once per unique signal
    {
        // Obtain a set of signal numbers which we need to register
        std::vector<int> s;
        s.insert(s.end(), sigdef.status.begin(),   sigdef.status.end());
        s.insert(s.end(), sigdef.restart.begin(),  sigdef.restart.end());
        s.insert(s.end(), sigdef.teardown.begin(), sigdef.teardown.end());
        std::sort(s.begin(), s.end());
        s.erase(std::unique(s.begin(), s.end()), s.end());

        // Register the signal handler for each of these signals
        typedef std::vector<int>::const_iterator const_iterator;
        for (const_iterator i = s.begin(); i != s.end(); ++i) {
            const char * name = suzerain_signal_name(*i);
            if (SIG_ERR != signal(*i, process_signal)) {
                if (name) {
                    DEBUG0("Registered signal handler for " << name);
                } else {
                    DEBUG0("Registered signal handler for " << *i);
                }
            } else {
                if (name) {
                    WARN0("Unable to register signal handler for " << name);
                } else {
                    WARN0("Unable to register signal handler for " << *i);
                }
            }
        }

        // Iff we registered any handlers, process signal receipt in stepper.
        // Notice signal receipt include --advance_wt calling us a pumpkin.
        // We can afford this every time step because of DeltaTAllreducer.
        if (s.size() > 0 || timedef.advance_wt > 0) {
            tc->add_periodic_callback(tc->forever_t(), 1,
                                      process_any_signals_received);
        }
    }

    // Advance time according to advance_dt, advance_nt criteria
    wtime_advance_start = MPI_Wtime();
    bool advance_success = true;
    switch ((!!timedef.advance_dt << 1) + !!timedef.advance_nt) {
        case 3:
            INFO0("Advancing simulation by at most " << timedef.advance_dt
                   << " units of physical time");
            INFO0("Advancing simulation by at most " << timedef.advance_nt
                   << " discrete time steps");
            advance_success = tc->advance(initial_t + timedef.advance_dt,
                                          timedef.advance_nt);
            break;
        case 2:
            INFO0("Advancing simulation by " << timedef.advance_dt
                  << " units of physical time");
            advance_success = tc->advance(initial_t + timedef.advance_dt);
            break;
        case 1:
            INFO0("Advancing simulation " << timedef.advance_nt
                   << " discrete time steps");
            advance_success = tc->step(timedef.advance_nt);
            break;
        case 0:
            if (!default_advance_nt) {
                INFO0("Simulation will not be advanced");
            } else {
                INFO0("Advancing simulation until terminated by a signal");
                advance_success = tc->advance();
            }
            break;
        default:
            FATAL0("Sanity error in time control");
            return EXIT_FAILURE;
    }
    if (soft_teardown) {
        INFO0("TimeController stopped advancing due to teardown signal");
        advance_success = true; // ...treat like successful advance
    } else if (!advance_success && tc->current_dt() < tc->min_dt()) {
        WARN0("TimeController halted because step " << tc->current_dt()
              << " was smaller than min_dt " << tc->min_dt() );
    } else if (!advance_success) {
        WARN0("TimeController halted unexpectedly");
    }
    const double wtime_advance_end = MPI_Wtime();

    // Output status if it was not just output during time advancement
    if (last_status_nt != tc->current_nt()) {
        log_status(tc->current_t(), tc->current_nt());
    }

    // Save a final restart if one was not just saved during time advancement
    if (advance_success && last_restart_saved_nt != tc->current_nt()) {
        INFO0("Saving final restart file");
        save_restart(tc->current_t(), tc->current_nt());
    }

    // Output details on time advancement
    {
        std::ostringstream msg;
        msg.precision(static_cast<int>(
                    numeric_limits<real_t>::digits10 * 0.75));
        msg << "Advanced simulation from t_initial = " << initial_t
            << " to t_final = " << tc->current_t()
            << " in " << tc->current_nt() << " steps";
        INFO0(msg.str());
    }
    INFO0("Min/mean/max/stddev of delta_t: "
          << tc->taken_min()  << ", "
          << tc->taken_mean() << ", "
          << tc->taken_max()  << ", "
          << tc->taken_stddev());

    // Output simulation advancement rate using flow through time language
    const real_t flowthroughs = (tc->current_t() - initial_t)
                              * (scenario.Ly / scenario.Lx);
    if (flowthroughs > 0) {
        INFO0("Simulation advance corresponds to "
              << flowthroughs << " flow throughs");
        INFO0("Simulation advancing at wall time per flow through of "
              << (wtime_advance_end - wtime_advance_start) / flowthroughs
              << " seconds");
        INFO0("Advancement rate calculation ignores "
              << (MPI_Wtime() - wtime_mpi_init
                  - (wtime_advance_end - wtime_advance_start))
              << " seconds of fixed overhead");
    }

    return advance_success ? EXIT_SUCCESS : EXIT_FAILURE;
}
