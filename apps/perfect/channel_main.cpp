//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
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
 * A perfect gas channel simulation driver.
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#define EIGEN_DEFAULT_IO_FORMAT \
        Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", ";", "", "", "[", "]")

#include <esio/error.h>
#include <esio/esio.h>
#include <gsl/gsl_errno.h>
#ifdef HAVE_UNDERLING
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <underling/underling.h>
#include <underling/error.h>
#endif

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/countof.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi_datatype.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/os.h>
#include <suzerain/pencil.hpp>
#include <suzerain/pre_gsl.h>
#include <suzerain/support/fftw_definition.hpp>
#include <suzerain/support/field.hpp>
#include <suzerain/support/grid_definition.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/noise_definition.hpp>
#include <suzerain/support/program_options.hpp>
#include <suzerain/support/restart_definition.hpp>
#include <suzerain/support/signal_definition.hpp>
#include <suzerain/support/statistics_definition.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/support/time_definition.hpp>
#include <suzerain/support/time_definition.hpp>
#include <suzerain/timers.h>
#include <suzerain/utility.hpp>
#include <suzerain/version.hpp>
#include <suzerain/zgbsv_specification.hpp>

#include "manufactured_solution.hpp"
#include "mean_quantities.hpp"
#include "perfect.hpp"

#include "channel_treatment.hpp"
#include "explicit_operator.hpp"
#include "hybrid_operator.hpp"

// Provided by channel_main_svnrev.{c,h} to speed recompilation
#pragma warning(push,disable:1419)
extern "C" const char revstr[];
#pragma warning(pop)

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
using boost::numeric_cast;
using std::numeric_limits;
using std::size_t;
using suzerain::array;
using suzerain::complex_t;
using suzerain::fullprec;
using suzerain::real_t;
using suzerain::scoped_ptr;
using suzerain::shared_ptr;
namespace logging = suzerain::support::logging;
namespace ndx     = suzerain::ndx;
namespace perfect = suzerain::perfect;
namespace support = suzerain::support;

// FIXME Generalize as part of Redmine ticket #2480
// We are only prepared to deal with 5 equations
static const std::vector<support::field> fields = perfect::default_fields();

// Timestepping scheme uses only complex_t 4D contiguous_state
// State indices range over (scalar field, Y, X, Z) in wave space
typedef suzerain::interleaved_state<4,complex_t> linear_state_type;
typedef suzerain::contiguous_state<4,complex_t>  nonlinear_state_type;

// Global scenario parameters initialized in main().  These are declared const
// to avoid accidental modification but have their const-ness const_cast away
// where necessary to load settings.
using suzerain::perfect::scenario_definition;
using suzerain::support::fftw_definition;
using suzerain::support::grid_definition;
using suzerain::support::noise_definition;
using suzerain::support::restart_definition;
using suzerain::support::signal_definition;
using suzerain::support::statistics_definition;
using suzerain::support::time_definition;

static const scenario_definition scenario;
static const grid_definition grid;
static const fftw_definition fftwdef(
        suzerain::fftw::measure, suzerain::fftw::estimate);
static const restart_definition restart(
        /* metadata    */ "metadata.h5.XXXXXX",
        /* uncommitted */ "uncommitted.h5.XXXXXX",
        /* destination */ "restart#.h5",
        /* retain      */ 1,
        /* dt          */ 0,
        /* nt          */ 0);
static const statistics_definition statsdef(
        /* destination */ "sample#.h5");
static const time_definition timedef(
        /* advance_dt                */ 0,
        /* advance_nt                */ 0,
        /* advance_wt                */ 0,
        /* status_dt                 */ 0,
        /* status_nt                 */ 0,
        /* min_dt                    */ 1e-8,
        /* max_dt                    */ 1);
static const noise_definition  noisedef;
static const signal_definition sigdef;

// Global details initialized in main()
static shared_ptr<      suzerain::bspline>              b;
static shared_ptr<      suzerain::bsplineop>            cop;    // Collocation
static shared_ptr<      suzerain::bsplineop>            gop;    // Galerkin L2
static shared_ptr<const suzerain::pencil_grid>          dgrid;
static shared_ptr<      perfect::manufactured_solution> msoln;

/** <tt>atexit</tt> callback to ensure we finalize underling. */
static void atexit_underling(void) {
    dgrid.reset();        // Runs pencil_grid destructors
#ifdef HAVE_UNDERLING
    underling_cleanup();  // Cleans up the library
#endif
}

// State details specific to this rank initialized in main()
static shared_ptr<linear_state_type>    state_linear;
static shared_ptr<nonlinear_state_type> state_nonlinear;

// Common storage shared between the linear and nonlinear operators
// which also includes instantaneous mean quantity statistics
static perfect::operator_common_block common_block;

// The last collection of mean quantities sample obtained
static perfect::mean_quantities sample;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void)
{
    if (esioh) esio_handle_finalize(esioh);
}

/**
 * <tt>atexit</tt> callback to remove the metadata file.  Do \em NOT remove
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

/** Log messages containing mean L2 and RMS fluctuation information */
static void information_L2(const std::string& prefix,
                           const char * const name_L2  = "L2.mean",
                           const char * const name_rms = "rms.fluct")
{
    // Avoid computational cost when logging is disabled
    logging::logger_type log_L2  = logging::get_logger(name_L2);
    logging::logger_type log_rms = logging::get_logger(name_rms);
    if (!INFO0_ENABLED(log_L2) && !INFO0_ENABLED(log_rms)) return;

    // Show headers only on first invocation
    std::ostringstream msg;
    static bool show_header = true;
    if (show_header) {
        msg << prefix;
        for (size_t k = 0; k < fields.size(); ++k)
            msg << ' ' << std::setw(fullprec<>::width) << fields[k].identifier;
        INFO0(log_L2, msg.str());
        INFO0(log_rms, msg.str());
        msg.str("");
        show_header = false;
    }

    // Collective computation of the L_2 norms
    state_nonlinear->assign(*state_linear);
    const std::vector<suzerain::field_L2> L2
        = suzerain::compute_field_L2(*state_nonlinear, grid, *dgrid, *gop);

    // Build and log L2 of mean conserved state
    msg << prefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(L2[k].mean());
    }
    INFO0(log_L2, msg.str());

    // Build and log root-mean-squared-fluctuations of conserved state
    // RMS fluctuations are a scaling factor away from L2 fluctuations
    const real_t rms_coeff = 1/std::sqrt(grid.L.x()*grid.L.y()*grid.L.z());
    msg.str("");
    msg << prefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(rms_coeff*L2[k].fluctuating());
    }
    INFO0(log_rms, msg.str());
}

/** Build a message containing bulk quantities */
static void information_bulk(const std::string& prefix)
{
    // Only continue on the rank housing the zero-zero modes...
    if (!dgrid->has_zero_zero_modes()) return;

    // ...and when logging is enabled.  Notice INFO not INFO0 is used.
    logging::logger_type bulk_state = logging::get_logger("bulk.state");
    if (!INFO_ENABLED(bulk_state)) return;

    // Show headers only on first invocation
    std::ostringstream msg;
    static bool show_header = true;
    if (show_header) {
        msg << prefix;
        for (size_t k = 0; k < fields.size(); ++k)
            msg << ' ' << std::setw(fullprec<>::width) << fields[k].identifier;
        INFO0(bulk_state, msg.str());
        msg.str("");
        show_header = false;
    }

    // Compute operator for finding bulk quantities from coefficients
    suzerain::VectorXr bulkcoeff(b->n());
    b->integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= grid.L.y();

    // Prepare the status message and log it
    msg << prefix;
    for (size_t k = 0; k < state_linear->shape()[0]; ++k) {
        suzerain::Map<suzerain::VectorXc> mean(
                (*state_linear)[k].origin(), state_linear->shape()[1]);
        msg << ' ' << fullprec<>(bulkcoeff.dot(mean.real()));
    }
    INFO(bulk_state, msg.str());
}

/** Build a message containing specific state quantities at the wall */
static void information_specific_wall_state(const std::string& prefix)
{
    // Only continue on the rank housing the zero-zero modes.
    if (!dgrid->has_zero_zero_modes()) return;

    logging::logger_type nick[2] = { logging::get_logger("wall.lower"),
                                     logging::get_logger("wall.upper")  };

    // Indices at the lower and upper walls.  Use that wall collocation point
    // values are nothing but the first and last B-spline coefficient values.
    size_t wall[2] = { 0, state_linear->shape()[1] - 1 };

    // Message lists rho, u, v, w, and total energy at walls
    for (size_t l = 0; l < SUZERAIN_COUNTOF(wall); ++l) {

        // Avoid computational cost when logging is disabled
        if (!DEBUG_ENABLED(nick[l])) continue;

        std::ostringstream msg;
        msg << prefix;

        const real_t rho = ((*state_linear)[ndx::rho][wall[l]][0][0]).real();
        for (size_t k = 0; k < fields.size(); ++k) {
            real_t val = (k == ndx::rho)
                       ? rho
                       : ((*state_linear)[k][wall[l]][0][0]).real() / rho;
            msg << ' ' << fullprec<>(val);
        }
        DEBUG(nick[l], msg.str());
    }
}

/**
 * Build a message for the absolute error versus a manufactured solution.
 * Uses state_nonlinear as a scratch space and is not cheap.
 */
static void information_manufactured_solution_absolute_error(
        const std::string& prefix,
        const real_t simulation_time)
{
    // Avoid computational cost when logging is disabled
    logging::logger_type mms_abserr = logging::get_logger("mms.abserr");
    if (!INFO0_ENABLED(mms_abserr)) return;

    // Compute L2 of error of state against manufactured solution
    assert(msoln);
    state_nonlinear->assign(*state_linear);
    perfect::accumulate_manufactured_solution(
            1, *msoln, -1, *state_nonlinear,
            grid, *dgrid, *cop, *b, simulation_time);
    const std::vector<suzerain::field_L2> L2
        = suzerain::compute_field_L2(*state_nonlinear, grid, *dgrid, *gop);

    // Output absolute global errors for each field
    std::ostringstream msg;
    msg << prefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(L2[k].total());
    }
    INFO0(mms_abserr, msg.str());
}

/** Tracks last time we output a status line */
static size_t last_status_nt = numeric_limits<size_t>::max();

/** Routine to output status.  Signature for timecontroller use. */
static bool log_status(real_t t, size_t nt)
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

    SUZERAIN_TIMER_SCOPED("log_status");

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

static void sample_statistics(real_t t)
{
    using perfect::mean_quantities;

    // Defensively avoid multiple invocations with no intervening changes
#pragma warning(push,disable:1572)
    if (sample.t == t) {
#pragma warning(pop)
        DEBUG0("Cowardly refusing to re-sample statistics at t = " << t);
        return;
    }

    SUZERAIN_TIMER_SCOPED("sample_statistics");

    const double starttime = MPI_Wtime();

    // Obtain mean samples from instantaneous fields
    state_nonlinear->assign(*state_linear);
    sample = perfect::sample_mean_quantities(
            scenario, grid, *dgrid, *cop, *state_nonlinear, t);

    // Obtain mean quantities computed via implicit forcing (when possible)
    if (common_block.means.rows() == sample.storage.rows()) {
       sample.f().col(0) = common_block.f();  // Only streamwise momentum...
       sample.f().rightCols<2>().setZero();   // ...not wall-normal, spanwise
       sample.f_dot_u()      = common_block.f_dot_u();
       sample.qb()           = common_block.qb();
       sample.CrhoE()        = common_block.CrhoE();
       sample.Crhou().col(0) = common_block.Crhou();
       sample.Crhou().col(1) = common_block.Crhov();
       sample.Crhou().col(2) = common_block.Crhow();
       sample.Crho()         = common_block.Crho();
       sample.Crhou_dot_u()  = common_block.Crhou_dot_u();
    } else {
        WARN0("Could not obtain mean quantities set by implicit forcing");
        sample.f          ().setConstant(numeric_limits<real_t>::quiet_NaN());
        sample.f_dot_u    ().setConstant(numeric_limits<real_t>::quiet_NaN());
        sample.qb         ().setConstant(numeric_limits<real_t>::quiet_NaN());
        sample.CrhoE      ().setConstant(numeric_limits<real_t>::quiet_NaN());
        sample.Crhou      ().setConstant(numeric_limits<real_t>::quiet_NaN());
        sample.Crho       ().setConstant(numeric_limits<real_t>::quiet_NaN());
        sample.Crhou_dot_u().setConstant(numeric_limits<real_t>::quiet_NaN());
    }

    // Convert implicit values on collocation points to coefficients
    suzerain::bsplineop_lu mass(*cop);
    mass.opform_mass(*cop);
    mass.factor();
    mass.solve(mean_quantities::nscalars::implicit,
               sample.storage.middleCols<mean_quantities::nscalars::implicit>(
                   mean_quantities::nscalars::physical).data(),
               sample.storage.innerStride(), sample.storage.outerStride());

    const double elapsed = MPI_Wtime() - starttime;
    INFO0("Computed statistics at t = " << t
          << " in " << elapsed << " seconds");
}

/** Tracks last time a restart file was written successfully */
static size_t last_restart_saved_nt = numeric_limits<size_t>::max();

/** Routine to store a restart file.  Signature for timecontroller use. */
static bool save_restart(real_t t, size_t nt)
{
    // Defensively avoid multiple invocations with no intervening changes
    if (last_restart_saved_nt == nt) {
        DEBUG0("Cowardly refusing to save multiple restarts at nt = " << nt);
        return true;
    }

    SUZERAIN_TIMER_SCOPED("save_restart");

    const double starttime = MPI_Wtime();
    DEBUG0("Started to store restart at t = " << t << " and nt = " << nt);

    DEBUG0("Cloning " << restart.metadata << " to " << restart.uncommitted);
    esio_file_clone(esioh, restart.metadata.c_str(),
                    restart.uncommitted.c_str(), 1 /*overwrite*/);
    support::save_time(esioh, t);

    // Copy state into state_nonlinear for possibly destructive processing
    state_nonlinear->assign(*state_linear);
    if (restart.physical) {
        DEBUG0("Storing primitive collocation point values into "
               << restart.uncommitted);
        perfect::save_collocation_values(
                esioh, *state_nonlinear, scenario, grid, *dgrid, *cop);
    } else {
        DEBUG0("Storing conserved coefficients into " << restart.uncommitted);
        support::save_coefficients(
                esioh, fields, *state_nonlinear, grid, *dgrid);
    }

    // Include statistics in the restart file
    sample_statistics(t);
    perfect::save(esioh, sample);

    DEBUG0("Committing " << restart.uncommitted
           << " as a restart file using template " << restart.destination);
    esio_file_close_restart(
            esioh, restart.destination.c_str(), restart.retain);

    const double elapsed = MPI_Wtime() - starttime;
    INFO0("Successfully wrote restart at t = " << t << " for nt = " << nt
          << " in " << elapsed << " seconds");

    last_restart_saved_nt = nt; // Maintain last successful restart time step

    return true; // Continue time advancement
}

/** Routine to write a statistics file.  Signature for timecontroller use. */
static bool save_statistics(real_t t, size_t nt)
{
    SUZERAIN_TIMER_SCOPED("save_statistics");

    const double starttime = MPI_Wtime();
    DEBUG0("Started to save statistics at t = " << t << " and nt = " << nt);

    // We use restart.{metadata,uncommitted} for statistics too.
    DEBUG0("Cloning " << restart.metadata << " to " << restart.uncommitted);
    esio_file_clone(esioh, restart.metadata.c_str(),
                    restart.uncommitted.c_str(), 1 /*overwrite*/);
    support::save_time(esioh, t);

    // Compute statistics and write to file
    sample_statistics(t);
    perfect::save(esioh, sample);

    DEBUG0("Committing " << restart.uncommitted
           << " as a statistics file using template " << statsdef.destination);
    esio_file_close_restart(
            esioh, statsdef.destination.c_str(), statsdef.retain);

    const double elapsed = MPI_Wtime() - starttime;
    INFO0("Successfully wrote statistics at t = " << t << " for nt = " << nt
          << " in " << elapsed << " seconds");

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
 * \li \c 4 Compute and write a statistics file
 */
typedef array<sig_atomic_t,5> atomic_signal_received_t;

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

    // atomic_signal_received[3] handled outside this routine

    // Determine if we should compute and write statistics due to the signal
    end = sigdef.statistics.end();
    if (std::find(sigdef.statistics.begin(), end, sig) != end) {
        atomic_signal_received[4] = sig;
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
 * Signature for timecontroller use.
 */
static bool process_any_signals_received(real_t t, size_t nt)
{
    // delta_t_allreducer performs the Allreduce necessary to get local status
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

    if (signal_received[4]) {
        INFO0("Computing and writing statistics due to receipt of "
              << suzerain_signal_name(signal_received[4]));
        keep_advancing = keep_advancing && save_statistics(t, nt);
    }

    // Clear signal_received to defensively avoid stale data bugs.
    // These would only be problematic if process_any_signals_received
    // was run multiple times in between delta_t_allreducer invocations.
    signal_received.assign(0);

    return keep_advancing;
}

/** Wall time at which MPI_Init completed */
static double wtime_mpi_init;

/** Wall time elapsed during FFTW planning */
static double wtime_fftw_planning;

/** Wall time elapsed during loading of state from the restart file */
static double wtime_load_state;

/** Wall time at which we began time stepping */
static double wtime_advance_start;

/**
 * A stateful functor that performs an MPI Allreduce to determine the minimum
 * stable time step size across all ranks.  The same MPI Allreduce is used to
 * hide the cost of querying atomic_signal_received across all ranks.
 */
static class delta_t_allreducer
{

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

    // Maintains the mean ratio of each delta_t_candidate to the minimum
    // delta_t_candidate selected on each operator() invocation.  Useful
    // for determining the relative restrictiveness of each criterion.
    std::vector<boost::accumulators::accumulator_set<
                real_t,
                boost::accumulators::stats<boost::accumulators::tag::mean>
        > > normalized_ratios;

    // Provides small default capacity for normalized_ratios
    delta_t_allreducer() : normalized_ratios(2) {}

    real_t operator()(const std::vector<real_t>& delta_t_candidates)
    {

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
            if (last_restart_saved_nt == numeric_limits<size_t>::max()) {
                wtime_projected += 2*wtime_load_state;  // Load as surrogate
            } else {
                wtime_projected += (acc::max)(period);  // Includes dumps
            }
            // ...to which we add an estimate of other finalization costs
            wtime_projected += 2*(wtime_advance_start - wtime_mpi_init
                                                      - wtime_fftw_planning);

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

        // Delegate finding-the-minimum work on each rank to delta_t_reducer
        // delta_t_reducer logic enforced requirement that min(NaN,x) == NaN
        const real_t delta_t
                = suzerain::timestepper::delta_t_reducer()(candidates);

        // Update normalized_ratios using the just chosen delta_t
        // isnan used to avoid a NaN from destroying all accumulator data
        if (!(boost::math::isnan)(delta_t)) {
            const size_t n = candidates.size();
            if (SUZERAIN_UNLIKELY(normalized_ratios.size() < n)) {
                normalized_ratios.resize(n);
            }
            for (size_t i = 0; i < n; ++i) {
                normalized_ratios[i](candidates[i] / delta_t);
            }
        }

        return delta_t;
    }

} delta_t_allreducer;

/**
 * Default log4cxx configuration (differs from support::log4cxx_config).
 * <tt>${FOO}</tt> syntax may be used to pick up environment variables in addition to properties
 * See http://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PatternLayout.html
 * and https://logging.apache.org/log4j/1.2/apidocs/org/apache/log4j/PropertyConfigurator.html
 */
static const char log4cxx_config[] =
    "## Set root logger level (e.g. TRACE, DEBUG, INFO) and one or more root appenders\n"
    "## Output INFO or higher messages on CONSOLE and in LOG\n"
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
    "## Collect \"rms.fluct\" messages into rms.fluct.dat mimicking LOG file behavior\n"
    "log4j.logger.rms.fluct=INHERITED, RMSFLUCT\n"
    "log4j.appender.RMSFLUCT=${log4j.appender.LOG}\n"
    "log4j.appender.RMSFLUCT.filename=rms.fluct.dat\n"
    "log4j.appender.RMSFLUCT.append=${log4j.appender.LOG.append}\n"
    "log4j.appender.RMSFLUCT.layout=${log4j.appender.LOG.layout}\n"
    "log4j.appender.RMSFLUCT.layout.ConversionPattern=${log4j.appender.LOG.layout.ConversionPattern}\n"
;

/** Main driver logic */
int main(int argc, char **argv)
{
#ifdef SUZERAIN_HAVE_GRVY
    grvy_timer_init("channel");                      // Initialize GRVY Timers
#endif
    MPI_Init(&argc, &argv);                          // Initialize MPI...
    wtime_mpi_init = MPI_Wtime();                    // Record MPI_Init time
    atexit((void (*) ()) MPI_Finalize);              // ...finalize at exit
    logging::initialize(MPI_COMM_WORLD,              // Initialize logging
                        log4cxx_config);
#ifdef HAVE_UNDERLING
    underling_init(&argc, &argv, 0);                 // Initialize underling...
#endif
    atexit(atexit_underling);                        // ...finalize at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD);  // Initialize ESIO
    atexit(&atexit_esio);                            // ...finalize at exit

    // Hook error handling into logging infrastructure
    gsl_set_error_handler(
            &support::mpi_abort_on_error_handler_gsl);
    suzerain_set_error_handler(
            &support::mpi_abort_on_error_handler_suzerain);
    esio_set_error_handler(
            &support::mpi_abort_on_error_handler_esio);
#ifdef HAVE_UNDERLING
    underling_set_error_handler(
            &support::mpi_abort_on_error_handler_underling);
#endif

    DEBUG0("Processing command line arguments and response files");
    std::string restart_file;
    std::string solver_spec(static_cast<std::string>(suzerain::zgbsv_specification()));
    bool use_explicit  = false;
    bool use_implicit  = false;
    bool use_yang11    = false;
    bool use_smr91     = false;
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
    bool use_p3dfft    = false;
    bool use_underling = false;
#endif
    bool default_advance_nt;
    bool default_statistics;
    {
        suzerain::support::program_options options(
                "Suzerain-based explicit compressible channel simulation",
                "RESTART-FILE", /* TODO description */ "", revstr);
        options.add_definition(const_cast<scenario_definition  &>(scenario));
        options.add_definition(const_cast<grid_definition      &>(grid    ));
        options.add_definition(const_cast<fftw_definition      &>(fftwdef ));
        options.add_definition(const_cast<restart_definition   &>(restart ));
        options.add_definition(const_cast<statistics_definition&>(statsdef));
        options.add_definition(const_cast<time_definition      &>(timedef ));
        options.add_definition(const_cast<noise_definition     &>(noisedef));
        options.add_definition(const_cast<signal_definition    &>(sigdef  ));

        options.add_options()
            ("explicit", "Use purely explicit operators")
            ("implicit", "Use hybrid implicit/explicit operators")
            ("solver",   boost::program_options::value<std::string>(&solver_spec)
                             ->default_value(solver_spec),
                         "Use the specified algorithm for any implicit solves")
            ("smr91",    "Advance time per Spalart, Moser, and Rogers 1991")
            ("yang11",   "Advance time per Shan Yang's 2011 thesis")
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
            ("p3dfft",    "Use P3DFFT for MPI-parallel FFTs")
            ("underling", "Use underling for MPI-parallel FFTs")
#endif
            ;
        std::vector<std::string> positional = options.process(argc, argv);

        // Select type of timestepping operators to use (default implicit)
        options.conflicting_options("implicit", "explicit");
        if (options.variables().count("explicit")) {
            use_explicit = true;
        } else {
            use_implicit = true;
        }

        // Select type of timestepping schemes to use (default smr91)
        options.conflicting_options("smr91", "yang11");
        if (options.variables().count("yang11")) {
            use_yang11 = true;
        } else {
            use_smr91 = true;
        }

#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
        // Select pencil decomposition and FFT library to use (default p3dfft)
        options.conflicting_options("p3dfft", "underling");
        if (options.variables().count("underling")) {
            use_underling = true;
        } else {
            use_p3dfft = true;
        }
#endif

        // Record build and invocation for posterity and to aid in debugging
        std::ostringstream os;
        std::copy(argv, argv+argc, std::ostream_iterator<const char *>(os," "));
        INFO0("Invocation: " << os.str());
        INFO0("Build:      " << suzerain::version("", revstr));

        switch (options.verbose()) {
            case 0:                   break;
            case 1:  DEBUG0_ENABLE(); break;
            default: TRACE0_ENABLE(); break;
        }
        switch (options.verbose_all()) {
            case 0:                   break;
            case 1:  DEBUG_ENABLE();  break;
            default: TRACE_ENABLE();  break;
        }

        if (positional.size() != 1) {
            FATAL0("Exactly one restart file name must be specified");
            return EXIT_FAILURE;
        }
        restart_file = positional[0];

        default_advance_nt = options.variables()["advance_nt"].defaulted();
        default_statistics =  options.variables()["statistics_dt"].defaulted()
                           && options.variables()["statistics_nt"].defaulted();
    }

    INFO0("Loading details from restart file: " << restart_file);
    esio_file_open(esioh, restart_file.c_str(), 0 /* read-only */);
    // Mach and gamma are "pushed" and "popped" around loading the restart file
    // to permit fixing temperature and density when changing the scenario.
    real_t restart_Ma, restart_gamma;
    {
        real_t cli_Ma = scenario.Ma, cli_gamma = scenario.gamma;
        const_cast<scenario_definition&>(scenario).Ma
            = const_cast<scenario_definition&>(scenario).gamma
            = numeric_limits<real_t>::quiet_NaN();
        const_cast<scenario_definition&>(scenario).load(esioh);
        restart_Ma    = scenario.Ma;
        restart_gamma = scenario.gamma;
        const_cast<scenario_definition&>(scenario).Ma
                = ((boost::math::isnan)(cli_Ma)) ? restart_Ma : cli_Ma;
        const_cast<scenario_definition&>(scenario).gamma
                = ((boost::math::isnan)(cli_gamma)) ? restart_gamma : cli_gamma;
    }
    const_cast<grid_definition&>(grid).load(esioh);
    support::load(esioh, const_cast<time_definition&>(timedef));
    perfect::load(esioh, msoln, scenario, grid);
    esio_file_close(esioh);

    if (msoln) {
        INFO0("Restart file prescribes a manufactured solution");
        if (boost::math::isnormal(scenario.bulk_rho)) {
            WARN0("Manufactured solution incompatible with bulk_rho = "
                  << scenario.bulk_rho);
        }
        if (boost::math::isnormal(scenario.bulk_rho_u)) {
            WARN0("Manufactured solution incompatible with bulk_rho_u = "
                  << scenario.bulk_rho_u);
        }
    }

    // Modify IEEE settings after startup complete as startup relies on NaNs
    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(suzerain::mpi::comm_rank(MPI_COMM_WORLD));

    support::create(grid.N.y(), grid.k, 0.0, grid.L.y(), grid.htdelta, b, cop);
    assert(b->k() == grid.k);
    assert(b->n() == grid.N.y());
    gop.reset(new suzerain::bsplineop(*b, 0, SUZERAIN_BSPLINEOP_GALERKIN_L2));

    // Compute and display a couple of discretization quality metrics.
    {
        // Temporarily work in real-valued quantities as it is a bit simpler.
        suzerain::bsplineop_lu boplu(*cop);
        boplu.opform_mass(*cop);
        double norm;
        boplu.opnorm(norm);
        boplu.factor();

        // Compute and display discrete conservation error magnitude
        suzerain::MatrixXXr mat = suzerain::MatrixXXr::Identity(b->n(),b->n());
        boplu.solve(b->n(), mat.data(), 1, b->n());         // M^-1
        cop->apply(1, b->n(), 1.0, mat.data(), 1, b->n());  // D*M^-1
        boplu.solve(b->n(), mat.data(), 1, b->n());         // M^-1*D*M^-1
        suzerain::VectorXr vec(b->n());
        b->integration_coefficients(0, vec.data());
        vec = vec.transpose() * mat;                        // w^{T}*M^-1*D*M^-1
        vec.head<1>()[0] -= -1;                             // Exact head
        vec.tail<1>()[0] -=  1;                             // Exact tail
        double relerr = vec.norm() / std::sqrt(2);          // Exact 2-norm
        INFO0("B-spline discrete conservation relative error near "
              << relerr * 100 << "%");

        // Compute and display condition number
        double rcond;
        boplu.rcond(norm, rcond);
        INFO0("B-spline mass matrix has condition number near "
              << (1 / rcond));
    }

    // Generating unique file names as needed using mkstemp(3)
    {
        // Pack a temporary buffer with the three file name templates
        array<size_t,5> pos = {{ 0,
                                 restart.metadata.length()     + 1,
                                 restart.uncommitted.length()  + 1,
                                 restart.destination.length()  + 1,
                                 statsdef.destination.length() + 1 }};
        std::partial_sum(pos.begin(), pos.end(), pos.begin());
        suzerain::scoped_array<char> buf(new char[pos[4]]);
        strcpy(&buf[pos[0]], restart.metadata.c_str());
        strcpy(&buf[pos[1]], restart.uncommitted.c_str());
        strcpy(&buf[pos[2]], restart.destination.c_str());
        strcpy(&buf[pos[3]], statsdef.destination.c_str());

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
            if (boost::ends_with(restart.destination, "XXXXXX")) {
                close(mkstemp(&buf[pos[2]]));  // Not clobbered later...
                unlink(&buf[pos[2]]);          // ...so remove any evidence
            }
            if (boost::ends_with(statsdef.destination, "XXXXXX")) {
                close(mkstemp(&buf[pos[3]]));  // Not clobbered later...
                unlink(&buf[pos[3]]);          // ...so remove any evidence
            }
        }

        // Broadcast any generated names to all ranks and unpack values
        SUZERAIN_MPICHKQ(MPI_Bcast(buf.get(), pos[4],
                         suzerain::mpi::datatype<char>(), 0,
                         MPI_COMM_WORLD));
        const_cast<restart_definition&>(restart).metadata        = &buf[pos[0]];
        const_cast<restart_definition&>(restart).uncommitted     = &buf[pos[1]];
        const_cast<restart_definition&>(restart).destination     = &buf[pos[2]];
        const_cast<statistics_definition&>(statsdef).destination = &buf[pos[3]];
    }

    DEBUG0("Saving metadata temporary file: " << restart.metadata);
    {
        atexit(&atexit_metadata); // Delete any lingering metadata file at exit
        esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
        esio_file_create(h, restart.metadata.c_str(), 1 /* overwrite */);
        esio_string_set(h, "/", "generated_by",
                        (std::string("channel ") + revstr).c_str()); // Ticket #2595
        scenario.save(h);
        grid.save(h);
        support::save(h, b, cop, gop);
        support::save(h, timedef);
        perfect::save(h, msoln, scenario, grid);
        esio_file_close(h);
        esio_handle_finalize(h);
    }

    // Display global degree of freedom information
    INFO0("Global number of unknowns:         " << (  grid.N.prod()
                                                    * fields.size()));
    INFO0("Grid degrees of freedom    (GDOF): " << grid.N.prod());
    INFO0("GDOF by direction           (XYZ): " << grid.N);
    INFO0("Dealiased GDOF by direction (XYZ): " << grid.dN);
    INFO0("Number of MPI ranks:               "
          << suzerain::mpi::comm_size(MPI_COMM_WORLD));

    // Initialize pencil_grid which handles P3DFFT setup/teardown RAII
    INFO0("Preparing MPI transpose and Fourier transform execution plans...");
    {
        const double begin = MPI_Wtime();
        fftw_set_timelimit(fftwdef.plan_timelimit);
        support::wisdom_broadcast(fftwdef.plan_wisdom);
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
        if (use_p3dfft) {
            dgrid = suzerain::make_shared<suzerain::pencil_grid_p3dfft>(
                    grid.dN, grid.P, fftwdef.rigor_fft, fftwdef.rigor_mpi);
        } else if (use_underling) {
            dgrid = suzerain::make_shared<suzerain::pencil_grid_underling>(
                    grid.dN, grid.P, fftwdef.rigor_fft, fftwdef.rigor_mpi);
        } else {
#endif
            dgrid = suzerain::make_shared<suzerain::pencil_grid_default>(
                    grid.dN, grid.P, fftwdef.rigor_fft, fftwdef.rigor_mpi);
#if defined(SUZERAIN_HAVE_P3DFFT) && defined(SUZERAIN_HAVE_UNDERLING)
        }
#endif
        support::wisdom_gather(fftwdef.plan_wisdom);
        wtime_fftw_planning = MPI_Wtime() - begin;
    }
    INFO0("MPI transpose and Fourier transform planning by "
          << dgrid->implementation() << " took "
          << wtime_fftw_planning << " seconds");
    assert((grid.dN == dgrid->global_physical_extent).all());
    INFO0("Rank grid used for decomposition: " << dgrid->processor_grid);
    INFO0("Zero-zero modes located on MPI_COMM_WORLD rank "
           << dgrid->rank_zero_zero_modes);
    { // Display normalized workloads metrics relative to zero-zero workload
        real_t sendbuf[4];
        if (dgrid->has_zero_zero_modes()) {
            sendbuf[0] = (real_t) dgrid->local_wave_extent.prod();
            sendbuf[1] = (real_t) dgrid->local_physical_extent.prod();
        }
        SUZERAIN_MPICHKQ(MPI_Bcast(sendbuf, 2,
                         suzerain::mpi::datatype<real_t>(), 0,
                         MPI_COMM_WORLD));
        sendbuf[0] = dgrid->local_wave_extent.prod()     / sendbuf[0];
        sendbuf[1] = dgrid->local_physical_extent.prod() / sendbuf[1];
        sendbuf[2] = -sendbuf[0];
        sendbuf[3] = -sendbuf[1];

        real_t recvbuf[4];
        SUZERAIN_MPICHKQ(MPI_Reduce(sendbuf, recvbuf, 2,
                         suzerain::mpi::datatype<real_t>(),
                         MPI_SUM, 0, MPI_COMM_WORLD));

        const std::size_t nranks = suzerain::mpi::comm_size(MPI_COMM_WORLD);
        const real_t mean_w = recvbuf[0] / nranks;
        const real_t mean_p = recvbuf[1] / nranks;
        SUZERAIN_MPICHKQ(MPI_Reduce(sendbuf, recvbuf, 4,
                         suzerain::mpi::datatype<real_t>(),
                         MPI_MIN, 0, MPI_COMM_WORLD));
        INFO0("Wave space zero-zero normalized workloads     (min/mean/max): "
              << recvbuf[0] << ", " << mean_w << ", " << -recvbuf[2]);
        INFO0("Physical space zero-zero normalized workloads (min/mean/max): "
              << recvbuf[1] << ", " << mean_p << ", " << -recvbuf[3]);
    }
    DEBUG("Local wave start      (XYZ): " << dgrid->local_wave_start);
    DEBUG("Local wave end        (XYZ): " << dgrid->local_wave_end);
    DEBUG("Local wave extent     (XYZ): " << dgrid->local_wave_extent);
    DEBUG("Local physical start  (XYZ): " << dgrid->local_physical_start);
    DEBUG("Local physical end    (XYZ): " << dgrid->local_physical_end);
    DEBUG("Local physical extent (XYZ): " << dgrid->local_physical_extent);

    // Create the state storage for nonlinear operator with appropriate padding
    state_nonlinear.reset(support::allocate_padded_state<nonlinear_state_type>(
                fields.size(), *dgrid));

    // Dump some state shape and stride information for debugging purposes
    DEBUG("Nonlinear state shape   (FYXZ): "
          << suzerain::multi_array::shape_array(*state_nonlinear));
    DEBUG("Nonlinear state strides (FYXZ): "
          << suzerain::multi_array::strides_array(*state_nonlinear));

    // Load restart information into state_nonlinear, including simulation time
    esio_file_open(esioh, restart_file.c_str(), 0 /* read-only */);
    real_t initial_t;
    support::load_time(esioh, initial_t);
    {
        const double begin = MPI_Wtime();
        sample.t = initial_t;         // For idempotent --advance_nt=0...
        perfect::load(esioh, sample); // ...when no grid rescaling employed
        perfect::load(esioh, *state_nonlinear,
                      scenario, grid, *dgrid, *cop, *b);
        wtime_load_state = MPI_Wtime() - begin;
    }
    esio_file_close(esioh);

    // If necessary, adjust total energy to account for scenario changes
    perfect::adjust_scenario(*state_nonlinear, scenario, grid, *dgrid, *cop,
                             restart_Ma, restart_gamma);

    // If requested, add noise to the momentum fields at startup (expensive).
    perfect::add_noise(*state_nonlinear, noisedef,
                       scenario, grid, *dgrid, *cop, *b);

    // Create state storage for linear operator usage
    state_linear = suzerain::make_shared<linear_state_type>(
            suzerain::to_yxz(fields.size(), dgrid->local_wave_extent));

    // Dump some state shape and stride information for debugging purposes
    DEBUG("Linear state shape      (FYXZ): "
          << suzerain::multi_array::shape_array(*state_linear));
    DEBUG("Linear state strides    (FYXZ): "
          << suzerain::multi_array::strides_array(*state_linear));

    // Copy (possibly perturbed) state from state_nonlinear into state_linear
    state_linear->assign(*state_nonlinear);

    // Prepare chosen time stepping scheme and required operators
    //
    // Operators are managed via shared_ptrs for both type erasure
    // (to avoid depending on specific operator details) and to
    // facilitate runtime selection of operators.
    //
    // The linear state type chosen for linear_operator is a superclass of both
    // contiguous_state and interleaved_state to allow swapping one for another
    // if so desired.  However, this is unlikely to be useful in conjunction
    // with hybrid implicit/explicit operators.
    shared_ptr<suzerain::timestepper::lowstorage::method_interface<
            complex_t
        > > m;
    if (use_smr91) {
        m.reset(new suzerain::timestepper::lowstorage::method<
                    suzerain::timestepper::lowstorage::smr91,
                    complex_t
                >(timedef.evmagfactor));
    } else if (use_yang11) {
        m.reset(new suzerain::timestepper::lowstorage::method<
                    suzerain::timestepper::lowstorage::yang11,
                    complex_t
                >(timedef.evmagfactor));
    } else {
        FATAL0("Sanity error in timestepping scheme selection");
        return EXIT_FAILURE;
    }
    shared_ptr<suzerain::timestepper::lowstorage::linear_operator<
            suzerain::multi_array::ref<complex_t,4>,
            nonlinear_state_type
        > > L;
    shared_ptr<suzerain::timestepper::nonlinear_operator<
            nonlinear_state_type
        > > N;

    using perfect::channel_treatment;
    if (use_explicit) {
        INFO0("Initializing explicit timestepping operators");
        L.reset(new channel_treatment<perfect::isothermal_mass_operator>(
                    scenario, grid, *dgrid, *cop, *b, common_block));
        N.reset(new perfect::explicit_nonlinear_operator(
                scenario, grid, *dgrid, *cop, *b, common_block, msoln));
    } else if (use_implicit) {
        INFO0("Initializing hybrid implicit/explicit timestepping operators");
        L.reset(new channel_treatment<perfect::isothermal_hybrid_linear_operator>(
                    solver_spec, scenario,
                    grid, *dgrid, *cop, *b, common_block));
        N.reset(new perfect::hybrid_nonlinear_operator(
                scenario, grid, *dgrid, *cop, *b, common_block, msoln));
    } else {
        FATAL0("Sanity error in operator selection");
        return EXIT_FAILURE;
    }

    // Prepare timecontroller for managing the time advance.
    // Nonlinear scaling factor chi = (N_x N_z)^(-1) from write up section 2.1
    // (Spatial discretization) accounts for dealiasing and included here.
    using suzerain::timecontroller;
    const real_t chi = real_t(1)/(grid.dN.x()*grid.dN.z());
    scoped_ptr<timecontroller<real_t> > tc(make_lowstorage_timecontroller(
                *m, delta_t_allreducer, *L, chi, *N,
                *state_linear, *state_nonlinear,
                initial_t, timedef.min_dt, timedef.max_dt));

    // Register status callbacks status_{dt,nt} as requested.
    tc->add_periodic_callback(
            (timedef.status_dt ? timedef.status_dt : tc->forever_t()),
            (timedef.status_nt ? timedef.status_nt : tc->forever_nt()),
            &log_status);

    // Register restart-writing callbacks restart_{dt,nt} as requested.
    tc->add_periodic_callback(
            (restart.dt ? restart.dt : tc->forever_t()),
            (restart.nt ? restart.nt : tc->forever_nt()),
            &save_restart);

    // Register statistics-related callbacks per statistics_{dt,nt}.
    // If no non-default, non-zero values were provided, be sensible.
    if (default_statistics && !statsdef.dt && !statsdef.nt) {
        const real_t flowthrough_time
                = grid.L.x()/(scenario.bulk_rho_u/scenario.bulk_rho);
        if (boost::math::isnormal(flowthrough_time)) {
            const_cast<real_t &>(statsdef.dt) = flowthrough_time / 4;
        }
    }
    tc->add_periodic_callback(
            (statsdef.dt ? statsdef.dt : tc->forever_t()),
            (statsdef.nt ? statsdef.nt : tc->forever_nt()),
            &save_statistics);

    // Register any necessary signal handling logic once per unique signal
    {
        // Obtain a set of signal numbers which we need to register
        std::vector<int> s;
        s.insert(s.end(), sigdef.status.begin(),     sigdef.status.end());
        s.insert(s.end(), sigdef.restart.begin(),    sigdef.restart.end());
        s.insert(s.end(), sigdef.statistics.begin(), sigdef.statistics.end());
        s.insert(s.end(), sigdef.teardown.begin(),   sigdef.teardown.end());
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
        // We can afford this every time step because of delta_t_allreducer.
        if (s.size() > 0 || timedef.advance_wt > 0) {
            tc->add_periodic_callback(tc->forever_t(), 1,
                                      process_any_signals_received);
        }
    }

    // Advance time according to advance_dt, advance_nt criteria
#ifdef SUZERAIN_HAVE_GRVY
    grvy_log_setlevel(GRVY_ERROR);  // Suppress GRVY timer resolution warnings
    grvy_timer_reset();
#endif
    wtime_advance_start = MPI_Wtime();
    bool advance_success = true;
    switch ((!!timedef.advance_dt << 1) + !!timedef.advance_nt) {
        case 3:
            INFO0("Advancing simulation by at most " << timedef.advance_dt
                   << " units of physical time using evmagfactor "
                   << timedef.evmagfactor);
            INFO0("Advancing simulation by at most " << timedef.advance_nt
                   << " discrete time steps using evmagfactor "
                   << timedef.evmagfactor);
            advance_success = tc->advance(initial_t + timedef.advance_dt,
                                          timedef.advance_nt);
            break;
        case 2:
            INFO0("Advancing simulation by " << timedef.advance_dt
                  << " units of physical time using evmagfactor "
                  << timedef.evmagfactor);
            advance_success = tc->advance(initial_t + timedef.advance_dt);
            break;
        case 1:
            INFO0("Advancing simulation " << timedef.advance_nt
                   << " discrete time steps using evmagfactor "
                   << timedef.evmagfactor);
            advance_success = tc->step(timedef.advance_nt);
            break;
        case 0:
            if (!default_advance_nt) {
                INFO0("Simulation will not be advanced");
            } else {
                INFO0("Advancing simulation until terminated by a signal"
                      " using evmagfactor " << timedef.evmagfactor);
                advance_success = tc->advance();
            }
            break;
        default:
            FATAL0("Sanity error in time control");
            return EXIT_FAILURE;
    }
    const double wtime_advance_end = MPI_Wtime();
#ifdef SUZERAIN_HAVE_GRVY
    grvy_timer_finalize();
    grvy_log_setlevel(GRVY_INFO);   // Re-enable GRVY warnings
#endif
    if (soft_teardown) {
        INFO0("timecontroller stopped advancing due to teardown signal");
        advance_success = true; // ...treat like successful advance
    } else if (!advance_success && tc->current_dt() < tc->min_dt()) {
        WARN0("timecontroller halted because step " << tc->current_dt()
              << " was smaller than min_dt " << tc->min_dt() );
    } else if (!advance_success) {
        WARN0("timecontroller halted unexpectedly");
    }

    // Output status if it was not just output during time advancement
    if (last_status_nt != tc->current_nt()) {
        log_status(tc->current_t(), tc->current_nt());
    }

    // Save a final restart if one was not just saved during time advancement
    if (advance_success && last_restart_saved_nt != tc->current_nt()) {
        INFO0("Saving final restart file");
        save_restart(tc->current_t(), tc->current_nt());
    }

    // Postprocess GRVY timer information on time advance now that restart is
    // safely on disk.  Reduces likelihood that GRVY hiccups cause data loss.
#ifdef SUZERAIN_HAVE_GRVY
    if (tc->current_nt() && dgrid->has_zero_zero_modes()) {
        // Only summarize when time advance took long enough to be interesting
        if (wtime_advance_end - wtime_advance_start > 5 /*seconds*/) {

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
        }
    }
#endif

    // Whenever we advanced the simulation, find the linearization
    // error present between the two chosen hybrid operator implementations.
    //
    // More specifically, for \partial_t u = \mathscr{L}u + Lu + (N(u)-Lu) did
    // the nonlinear_operator compute a "-Lu" matching with the
    // linear_operator's "+Lu" contribution?  linear_operator is assumed to
    // also compute \mathscr{L}u which is independent of reference values (e.g.
    // the divergence of the momentum within the mass equation).
    //
    // Procedure performed whenever we advance time to (a) ensure the test is
    // run often but (b) not costly when we wish to merely down-sample restart
    // files or convert to/from physical space (e.g. --advance_nt=0).  In those
    // circumstances the user has already paid the memory overhead for
    // nonlinear operator applications.  The incremental runtime overhead is
    // therefore no memory and computation smaller than a complete Runge-Kutta
    // step.  Done at shutdown as the check destroys the state information so
    // it may be performed within low storage memory requirements.
    //
    // Success requires that the following preconditions all hold:
    //    i) At substep zero, nonlinear_operator computes N(u) - Lu
    //       including the necessary reference quantities.
    //   ii) The reference quantities are stored in common_block
    //       and may be expressly zeroed.  Zeroing them nukes linearized
    //       contributions from I{Nonlinear,Linear}Operator but not
    //       linear contributions from linear_operator.
    //  iii) On non-zero substeps for zero reference values,
    //       nonlinear_operator computes only N(u) and linear_operator
    //       computes only \left(M+\varphi\mathscr{L}\right)u.
    //   iv) At the beginning of this process, state_linear contains
    //       valid information and adheres to boundary conditions.
    //    v) The explicit_nonlinear_operator application requires an
    //       auxiliary scaling factor "chi" to account for Fourier
    //       transform normalization needs.
    //   vi) Using suzerain::diffwave::apply zeros wavenumbers
    //       used only for dealiasing purposes-- this data is unimportant
    //       from the perspective of measuring actual linearization error.
    if (advance_success && tc->current_nt()) {
        const double starttime = MPI_Wtime();
        state_nonlinear->assign(*state_linear);
        common_block.set_zero(grid.dN.y());  // Defensive
        N->apply_operator(tc->current_t(), *state_nonlinear,
                m->evmaxmag_real(), m->evmaxmag_imag(), /*substep*/0);
        L->accumulate_mass_plus_scaled_operator(
                1., *state_linear, chi,  *state_nonlinear, *m, 0, /*substep*/0);
        common_block.set_zero(grid.dN.y());  // Zero reference quantities
        L->accumulate_mass_plus_scaled_operator(
                1., *state_linear, -1., *state_nonlinear, *m, 0, /*substep*/0);
        for (size_t k = 0; k < fields.size(); ++k) {
            suzerain::diffwave::apply(0, 0, 1., (*state_nonlinear)[k].origin(),
                grid.L.x(), grid.L.z(), dgrid->global_wave_extent.y(),
                grid.N.x(), grid.dN.x(),
                dgrid->local_wave_start.x(), dgrid->local_wave_end.x(),
                grid.N.z(), grid.dN.z(),
                dgrid->local_wave_start.z(), dgrid->local_wave_end.z());
        }
        state_nonlinear->exchange(*state_linear);
        N->apply_operator(tc->current_t(), *state_nonlinear,
                m->evmaxmag_real(), m->evmaxmag_imag(), /*substep*/1);
        for (size_t k = 0; k < fields.size(); ++k) {
            suzerain::diffwave::apply(0, 0, 1., (*state_nonlinear)[k].origin(),
                grid.L.x(), grid.L.z(), dgrid->global_wave_extent.y(),
                grid.N.x(), grid.dN.x(),
                dgrid->local_wave_start.x(), dgrid->local_wave_end.x(),
                grid.N.z(), grid.dN.z(),
                dgrid->local_wave_start.z(), dgrid->local_wave_end.z());
        }
        state_nonlinear->add_scaled(1/chi, *state_linear);
        const std::vector<suzerain::field_L2> L2
            = suzerain::compute_field_L2(*state_nonlinear, grid, *dgrid, *gop);
        const double elapsed = MPI_Wtime() - starttime;
        DEBUG0("Computed linearization error in " << elapsed << " seconds");

        // When operators "match", state_nonlinear should be identically zero.
        // Seeing more than floating point error indicates something is amiss.
        std::ostringstream msg;
        msg << "Linearization error";
        for (size_t k = 0; k < L2.size(); ++k) {
            msg << ' ' << fullprec<>(L2[k].total());
        }
        INFO0("lin.abserr", msg.str());
    }
    // Beware, state_linear and state_nonlinear now contain garbage!

    // Output details on time advancement (whenever advancement occurred)
    if (tc->current_nt()) {
        std::ostringstream msg;
        msg.precision(static_cast<int>(numeric_limits<real_t>::digits10*0.75));
        msg << "Advanced simulation from t_initial = " << initial_t
            << " to t_final = " << tc->current_t()
            << " in " << tc->current_nt() << " steps";
        INFO0(msg.str());
        msg.str("");
        msg.precision(static_cast<int>(numeric_limits<real_t>::digits10*0.50));
        msg << "Min/mean/max/stddev of delta_t: "
            << tc->taken_min()  << ", "
            << tc->taken_mean() << ", "
            << tc->taken_max()  << ", "
            << tc->taken_stddev();
        INFO0(msg.str());
        msg.str("");
        msg << "Mean delta_t criteria versus minimum criterion: ";
        const size_t n = delta_t_allreducer.normalized_ratios.size();
        for (size_t i = 0; i < n; ++i) {
            namespace acc = boost::accumulators;
            msg << acc::mean(delta_t_allreducer.normalized_ratios[i]);
            if (i < n-1) msg << ", ";
        }
        INFO0(msg.str());
    }

    // Output simulation advance rates when we advanced the simulation
    if (tc->current_nt()) {
        const real_t wtime_advance = (wtime_advance_end - wtime_advance_start);

        // Advance rate measured in a (mostly) problem-size-agnostic metric
        INFO0("Advancing at " << wtime_advance/tc->current_nt()/grid.N.prod()
                              << " wall seconds per time step per grid point");

        // Advance rate measured in a problem-size-dependent metric
        INFO0("Advancing at " << wtime_advance/tc->current_nt()
                              << " wall seconds per time step");

        // Advance rate measured in nondimensional simulation time units
        INFO0("Advancing at " << wtime_advance/(tc->current_t() - initial_t)
                              << " wall seconds per simulation time unit");

        // Advance rate measured in flow through based on bulk velocity
        // (where bulk velocity is estimated from bulk momentum and density)
        const real_t flowthrough_time
                = grid.L.x()/(scenario.bulk_rho_u/scenario.bulk_rho);
        const real_t flowthroughs
                = (tc->current_t() - initial_t)/flowthrough_time;
        if (boost::math::isnormal(flowthrough_time)) {
            INFO0("Advancing at " << wtime_advance/flowthroughs
                                  << " wall seconds per flow through");
        }

        // Admit what overhead we've neglected in those calculations
        INFO0("Advancement rate calculations ignore "
                              << MPI_Wtime() - wtime_mpi_init - wtime_advance
                              << " seconds of fixed overhead");
    }

    // Avoid destruction order fiasco from support::allocate_padded_state
    state_nonlinear.reset();
    state_linear.reset();

    return advance_success ? EXIT_SUCCESS : EXIT_FAILURE;
}
