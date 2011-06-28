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
#include <suzerain/pencil.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/restart_definition.hpp>
#include <suzerain/time_definition.hpp>
#include <suzerain/utility.hpp>

#include "logger.hpp"
#include "precision.hpp"
#include "channel.hpp"
#include "explicit_op.hpp"

#pragma warning(disable:383 1572)

using boost::make_shared;
using boost::numeric_cast;
using boost::scoped_ptr;
using boost::shared_ptr;
using std::numeric_limits;

// Explicit timestepping scheme uses only complex_t 4D NoninterleavedState
// State indices range over (scalar field, Y, X, Z) in wave space
typedef suzerain::NoninterleavedState<4,complex_t> state_type;

// Global scenario parameters initialized in main().  These are declared const
// to avoid accidental modification but have their const-ness const_cast away
// where necessary to load settings.
using suzerain::problem::ScenarioDefinition;
using suzerain::problem::GridDefinition;
using suzerain::problem::RestartDefinition;
using suzerain::problem::TimeDefinition;
static const ScenarioDefinition<real_t> scenario;
static const GridDefinition grid;
static const RestartDefinition<> restart(
        /* metadata     */ "metadata.h5",
        /* uncommitted  */ "uncommitted.h5",
        /* desttemplate */ "restart#.h5",
        /* retain       */ 1,
        /* restart_dt   */ 0,
        /* restart_nt   */ 0);
static const TimeDefinition<real_t> timedef(
        /* advance_dt                */ 0,
        /* advance_nt                */ 0,
        /* status_dt                 */ 0,
        /* status_nt                 */ 0,
        /* min_dt                    */ 0,
        /* max_dt                    */ 0,
        /* evmagfactor per Venugopal */ 0.72);

// Global details initialized in main()
static shared_ptr<      suzerain::bspline>       b;
static shared_ptr<      suzerain::bsplineop>     bop;    // Collocation
static shared_ptr<      suzerain::bsplineop>     gop;    // Galerkin L2
static shared_ptr<      suzerain::bsplineop_luz> bopluz;
static shared_ptr<const suzerain::pencil_grid>   dgrid;
static shared_ptr<nsctpl_rholut::manufactured_solution<real_t> > msoln;

// State details specific to this rank initialized in main()
static shared_ptr<state_type> state_linear;
static shared_ptr<state_type> state_nonlinear;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

/** <tt>atexit</tt> callback to remove the metadata file. */
static void atexit_metadata(void) {
    if (suzerain::mpi::comm_rank(MPI_COMM_WORLD) == 0) {
        if (0 == unlink(restart.metadata().c_str())) {
            DEBUG("Cleaned up temporary file " << restart.metadata());
        } else {
            WARN("Error cleaning up temporary file " << restart.metadata());
        }
    }
}

/** Build a message containing mean and fluctuating L2 information */
static std::string information_L2() {

    // Collective computation of the L_2 norms
    const boost::array<channel::L2,channel::field::count> L2
        = channel::field_L2(*state_linear, scenario, grid, *dgrid, *gop);

    // Prepare the status message
    std::ostringstream msg;
    msg << "mean|fluct L2 =";
    msg.precision(static_cast<int>(numeric_limits<real_t>::digits10 * 0.75));
    for (std::size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << L2[k].mean() << ' ' << L2[k].fluctuating();
    }

    return msg.str();
}

/** Build a message containing bulk quantities (intended for root rank only) */
static std::string information_bulk() {

    // Compute operator for finding bulk quantities from coefficients
    Eigen::VectorXr bulkcoeff(b->n());
    b->integration_coefficients(0, bulkcoeff.data());
    bulkcoeff /= scenario.Ly;

    // Prepare the status message
    std::ostringstream msg;
    msg << "bulk state =";
    msg.precision(static_cast<int>(numeric_limits<real_t>::digits10 * 0.75));
    for (std::size_t k = 0; k < state_linear->shape()[0]; ++k) {
        Eigen::Map<Eigen::VectorXc> mean(
                (*state_linear)[k].origin(), state_linear->shape()[1]);
        msg << ' ' << bulkcoeff.dot(mean.real());
    }

    return msg.str();
}

/** Build a message containing specific state quantities at the wall */
static std::string information_specific_wall_state() {
    namespace ndx = channel::field::ndx;

    // Indices at the lower and upper walls.  Use that wall collocation point
    // values are nothing but the first and last B-spline coefficient values.
    std::size_t wall[2] = { 0, state_linear->shape()[1] - 1 };
    real_t      rho[2]  = {
        ((*state_linear)[ndx::rho][wall[0]][0][0]).real(),
        ((*state_linear)[ndx::rho][wall[1]][0][0]).real()
    };

    // Prepare the status message:
    // Message lists lower then upper u, v, w, and internal energy at walls
    std::ostringstream msg;
    msg << "specific wall state = ";
    msg.precision(std::numeric_limits<real_t>::digits10 / 2);
    assert(ndx::rho == 0);
    for (std::size_t k = ndx::rho; k < channel::field::count; ++k) {
        for (std::size_t l = 0; l < sizeof(wall)/sizeof(wall[0]); ++l) {
            msg << ' ' << ((*state_linear)[k][wall[l]][0][0]).real() / rho[l];
        }
    }

    return msg.str();
}

/**
 * Build a message for the relative error versus a manufactured solution.
 * Uses state_nonlinear as a scratch space and is seriously collective
 * and not cheap.
 */
static std::string information_manufactured_solution_relative_error(
        const real_t simulation_time)
{
    assert(msoln);

    // First, compute L2 of error of state against manufactured solution
    state_nonlinear->assign(*state_linear);
    channel::accumulate_manufactured_solution(
            1, *msoln, -1, *state_nonlinear,
            scenario, grid, *dgrid, *b, *bop, simulation_time);
    const boost::array<channel::L2,channel::field::count> L2_error
        = channel::field_L2(*state_nonlinear, scenario, grid, *dgrid, *gop);

    // Next, compute L2 of manufactured solution
    channel::accumulate_manufactured_solution(
            1, *msoln, 0, *state_nonlinear,
            scenario, grid, *dgrid, *b, *bop, simulation_time);
    const boost::array<channel::L2,channel::field::count> L2_msoln
        = channel::field_L2(*state_nonlinear, scenario, grid, *dgrid, *gop);

    // Last, compute and output relative global errors for each field
    std::ostringstream msg;
    msg << "relative MMS error = ";
    msg.precision(static_cast<int>(numeric_limits<real_t>::digits10));
    for (std::size_t k = 0; k < channel::field::count; ++k) {
        msg << ' ' << L2_error[k].total() / L2_msoln[k].total();
    }

    return msg.str();
}

/** Tracks last time we output a status line */
static std::size_t last_status_nt = numeric_limits<std::size_t>::max();

/** Routine to output status.  Signature for TimeController use. */
static bool log_status(real_t t, std::size_t nt) {

    // Save resources by returning early when no status necessary
    if (!INFO_ENABLED) return true;

    // Build time- and timestep-specific status prefix.
    // Precision computations ensure multiple status lines minimally distinct
    std::ostringstream timeprefix;
    timeprefix << "t = ";
    const real_t nplaces = (timedef.status_dt > 0)
                         ? -std::floor(std::log10(timedef.status_dt))
                         : 0;
    if (nplaces > 0) {
        timeprefix.setf(std::ios::fixed,std::ios::floatfield);
        const std::streamsize oldprec = timeprefix.precision(nplaces);
        timeprefix << t;
        timeprefix.precision(oldprec);
        timeprefix.unsetf(std::ios::fixed);
    } else {
        timeprefix << t;
    }
    timeprefix << ", nt = " << nt << ", ";

    // On root only, compute and show bulk state quantities
    INFO0(timeprefix.str() << information_bulk());

    // Collectively compute and log L2 mean and fluctuating information
    const std::string msg_l2 = information_L2() /* collective! expensive! */;
    INFO0(timeprefix.str() << msg_l2);

    // On root only, compute and show specific state at the walls
    DEBUG0(timeprefix.str() << information_specific_wall_state());

    // If using manufactured solution, compute and log relative error
    if (msoln) {
        const std::string msg_relerr /* collective! expensive! */
              = information_manufactured_solution_relative_error(t);
        INFO0(timeprefix.str() << msg_relerr);
    }

    last_status_nt = nt; // Maintain last status time step

    return true;
}

/** Tracks last time a restart file was written successfully */
static std::size_t last_restart_saved_nt = numeric_limits<std::size_t>::max();

/** Routine to store a restart file.  Signature for TimeController use. */
static bool save_restart(real_t t, std::size_t nt)
{
    esio_file_clone(esioh, restart.metadata().c_str(),
                    restart.uncommitted().c_str(), 1 /*overwrite*/);

    DEBUG0("Started to store restart at t = " << t << " and nt = " << nt);
    channel::store_time(esioh, t);

    DEBUG0("Started to store simulation fields");
    channel::store(esioh, *state_linear, grid, *dgrid);

    DEBUG0("Started to commit restart file");
    esio_file_close_restart(esioh,
                            restart.desttemplate().c_str(),
                            restart.retain());

    INFO0("Successfully wrote restart at t = " << t << " for nt = " << nt);

    last_restart_saved_nt = nt; // Maintain last successful restart time step

    return true; // Continue time advancement
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    const double wtime_mpi_init = MPI_Wtime();      // Record MPI_Init time
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // Finalize ESIO at exit

    // Obtain some basic MPI environment details
    const int nranks = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    // Establish MPI-savvy, rank-dependent logging names
    name_logger_within_comm_world();

    // Hook error handling into logging infrastructure
    gsl_set_error_handler(
            &channel::mpi_abort_on_error_handler_gsl);
    suzerain_set_error_handler(
            &channel::mpi_abort_on_error_handler_suzerain);
    esio_set_error_handler(
            &channel::mpi_abort_on_error_handler_esio);

    DEBUG0("Processing command line arguments and response files");
    std::string restart_file;
    bool default_advance_dt;
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
                const_cast<RestartDefinition<>&>(restart));
        options.add_definition(
                const_cast<TimeDefinition<real_t>&>(timedef));
        std::vector<std::string> positional = options.process(argc, argv);

        if (positional.size() != 1) {
            FATAL0("Exactly one restart file name must be specified");
            return EXIT_FAILURE;
        }
        restart_file = positional[0];

        default_advance_dt = options.variables()["advance_dt"].defaulted();
        default_advance_nt = options.variables()["advance_nt"].defaulted();
    }

    if (default_advance_dt && default_advance_nt) {
        FATAL0("Either --advance_dt or --advance_nt is required");
        return EXIT_FAILURE;
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

    INFO0("Saving metadata temporary file: " << restart.metadata());
    {
        esio_handle h = esio_handle_initialize(MPI_COMM_WORLD);
        esio_file_create(h, restart.metadata().c_str(), 1 /* overwrite */);
        channel::store(h, scenario);
        channel::store(h, grid, scenario.Lx, scenario.Lz);
        channel::store(h, b, bop, gop);
        channel::store(h, scenario, msoln);
        esio_file_close(h);
        esio_handle_finalize(h);
        atexit(&atexit_metadata); // Delete lingering metadata file at exit
    }

    // Display global degree of freedom information
    INFO0("Global degrees of freedom  (DOF): " << grid.N.prod());
    INFO0("DOF by direction           (XYZ): " << grid.N);
    INFO0("Dealiased DOF by direction (XYZ): " << grid.dN);
    INFO0("Number of MPI ranks:              " << nranks);

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
    channel::load(esioh, *state_linear, grid, *dgrid, *b, *bop);
    esio_file_close(esioh);

    // Create the state storage for nonlinear operator with appropriate padding
    // to allow P3DFFTification.  Must clear to avoid lingering NaN issues.
    {
        using suzerain::to_yxz;
        using suzerain::prepend;
        using suzerain::strides_cm;
        state_nonlinear = make_shared<state_type>(
                to_yxz(channel::field::count, dgrid->local_wave_extent),
                prepend(dgrid->local_wave_storage(), strides_cm(
                        to_yxz(dgrid->local_wave_extent)))
                );
    }
    suzerain::multi_array::fill(*state_nonlinear, 0);

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
                *m, *L, real_t(1)/(grid.dN.x()*grid.dN.z()), *N,
                *state_linear, *state_nonlinear,
                initial_t, timedef.min_dt, timedef.max_dt));

    // Register status callbacks status_{dt,nt} as requested
    // When either is not provided, default to a reasonable behavior
    {
        TimeController<real_t>::time_type dt;
        if (timedef.status_dt) {
            dt = timedef.status_dt;
        } else if (timedef.advance_dt) {
            dt = timedef.advance_dt / restart.retain() / 5;
        } else {
            dt = tc->forever_t();
        }

        TimeController<real_t>::step_type nt;
        if (timedef.status_nt) {
            nt = timedef.status_nt;
        } else if (timedef.advance_nt) {
            nt = timedef.advance_nt / restart.retain() / 5;
            nt = std::max<TimeController<real_t>::step_type>(1, nt);
        } else {
            nt = tc->forever_nt();
        }

        tc->add_periodic_callback(dt, nt, &log_status);
    }

    // Register restart-writing callbacks restart_{dt,nt} as requested
    // When either is not provided, default to a reasonable behavior
    {
        TimeController<real_t>::time_type dt;
        if (restart.restart_dt()) {
            dt = restart.restart_dt();
        } else if (timedef.advance_dt) {
            dt = timedef.advance_dt / restart.retain();
        } else {
            dt = tc->forever_t();
        }

        TimeController<real_t>::step_type nt;
        if (restart.restart_nt()) {
            nt = restart.restart_nt();
        } else if (timedef.advance_nt) {
            nt = timedef.advance_nt / restart.retain();
            nt = std::max<TimeController<real_t>::step_type>(1, nt);
        } else {
            nt = tc->forever_nt();
        }

        tc->add_periodic_callback(dt, nt, &save_restart);
    }

    // Advance time according to advance_dt, advance_nt criteria
    const double wtime_advance_start = MPI_Wtime();
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
            if (!default_advance_dt) {
                INFO0("Advancing simulation until forcibly terminated");
                advance_success = tc->advance();
            } else if (!default_advance_nt) {
                WARN0("Simulation will not be advanced");
            } else {
                FATAL0("Sanity error in time control");
                return EXIT_FAILURE;
            }
            break;
        default:
            FATAL0("Sanity error in time control");
            return EXIT_FAILURE;
    }
    if (!advance_success) {
        WARN0("TimeController stopped advancing time unexpectedly");
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

    // Output statistics on time advancement
    INFO0("Advanced simulation from t_initial = " << initial_t
          << " to t_final = " << tc->current_t()
          << " in " << tc->current_nt() << " steps");
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
