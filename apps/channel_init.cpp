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
// channel_ex.cpp: Initialize restart files for use with Suzerain
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <log4cxx/logger.h>
#include <esio/esio.h>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/utility.hpp>

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
namespace sz = ::suzerain;
namespace pb = ::suzerain::problem;
namespace po = ::boost::program_options;
using boost::make_shared;
using boost::math::constants::pi;
using boost::shared_ptr;
using std::numeric_limits;

// Introduce scalar- and complex-valued typedefs
// Currently only real_t == double is supported by many, many components
typedef double               real_t;
typedef std::complex<real_t> complex_t;

// Global scenario parameters initialized in main()
static pb::ScenarioDefinition<real_t> def_scenario(
        /* default_Re    */ 100,
        /* default_Pr    */ real_t(7)/real_t(10),
        /* default_gamma */ real_t(14)/real_t(10),
        /* default_beta  */ real_t(2)/real_t(3),
        /* default_Lx    */ 4*pi<real_t>(),
        /* default_Ly    */ 2,
        /* default_Lz    */ 4*pi<real_t>()/3);
static pb::GridDefinition<real_t> def_grid(
        /* default_Nx    */ 1,
        /* default_DAFx  */ real_t(3)/real_t(2),
        /* default_Ny    */ 16,
        /* default_k     */ 6,
        /* default_Nz    */ 1,
        /* default_DAFz  */ real_t(3)/real_t(2));

// Global grid-details initialized in main()
static shared_ptr<sz::bspline>     bspw;
static shared_ptr<sz::bspline_luz> bspluzw;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // Finalize ESIO at exit

    // Obtain some basic MPI environment details.
    const int nproc  = sz::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = sz::mpi::comm_rank(MPI_COMM_WORLD);
    log4cxx::LoggerPtr log = log4cxx::Logger::getLogger(
            sz::mpi::comm_rank_identifier(MPI_COMM_WORLD));

    // Log only warnings and above from ranks 1 and higher when not debugging
    if (procid > 0 && !log->isDebugEnabled()) {
        log->setLevel(log4cxx::Level::getWarn());
    }

    // Ensure that we're running in a single processor environment
    if (nproc > 1) {
        LOG4CXX_FATAL(log, argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // Process incoming program arguments from command line, input files
    std::string restart_file;
    real_t htdelta;
    {
        sz::ProgramOptions options(
                "Suzerain-based compressible channel initialization");
        options.add_definition(def_scenario);
        options.add_definition(def_grid);
        using ::suzerain::validation::ensure_positive;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_positive(ensure_positive<real_t>);
        options.add_options()
            ("create", po::value<std::string>(&restart_file)
                ->default_value("restart0.h5"),
             "Name of new restart file to create"),
            ("htdelta", po::value<real_t>(&htdelta)
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"htdelta"))
                ->default_value(7),
             "Hyperbolic tangent stretching parameter for breakpoints")
        ;
        options.process(argc, argv);
    }

    LOG4CXX_INFO(log, "Creating new restart file " << restart_file);
    esio_file_create(esioh, restart_file.c_str(), false /* no clobber */);
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing basic scenario parameters");
    esio_attribute_write_double(esioh, "Re",    def_scenario.Re());
    esio_attribute_write_double(esioh, "Pr",    def_scenario.Pr());
    esio_attribute_write_double(esioh, "gamma", def_scenario.gamma());
    esio_attribute_write_double(esioh, "beta",  def_scenario.beta());
    esio_attribute_write_double(esioh, "Lx",    def_scenario.Lx());
    esio_attribute_write_double(esioh, "Ly",    def_scenario.Ly());
    esio_attribute_write_double(esioh, "Lz",    def_scenario.Lz());
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Finding B-spline basis of uniform order "
                      << (def_grid.k() - 1) << " on [0, Ly] with "
                      << def_grid.Ny() << " DOF");
    {
        boost::scoped_array<real_t> buf(new real_t[def_scenario.Ny()]);

        // Compute and store breakpoint locations
        const int nbreak = def_grid.Ny() + 2 - k;
        sz::math::linspace(0.0, 1.0, nbreak, buf.begin()); // Uniform [0, 1]
        for (int i = 0; i < nbreak; ++i) {                 // Stretch 'em out
            buf[i] = def_scenario.Ly()
                   * suzerain_htstretch2(htdelta, 1.0, buf[i]);
        }
        esio_attribute_writev_double(
                esioh, "breakpoints", buf.begin(), nbreak);

        // Generate the B-spline workspace based on order and breakpoints
        bspw = make_shared<sz::bspline>(def_grid.k(), 2, nbreak, breakpoints);
        assert(static_cast<unsigned>(bspw->ndof()) == def_grid.Ny());

        // Store collocation points to restart file
        bspw->collocation_points(buf.begin(), 1);
        esio_attribute_writev_double(
                esioh, "colpoints", buf.begin(), bspw->ndof());
    }
    esio_file_flush(esioh);


    const int nbreak = def_grid.Ny() + 2 - def_grid.k();
    real_t *breakpoints = (real_t *) sz::blas::malloc(nbreak*sizeof(real_t));
    assert(breakpoints);
    sz::math::linspace(0.0, 1.0, nbreak, breakpoints); // Uniform [0, 1]
    for (int i = 0; i < nbreak; ++i) {                 // Stretch 'em out
        breakpoints[i] = def_scenario.Ly()
                       * suzerain_htstretch2(htdelta, 1.0, breakpoints[i]);
    }
    esio_attribute_writev_double(esioh, "breakpoints", breakpoints, nbreak);
    esio_file_flush(esioh);

    bspw = make_shared<sz::bspline>(def_grid.k(), 2, nbreak, breakpoints);
    assert(static_cast<unsigned>(bspw->ndof()) == def_grid.Ny());
    sz::blas::free(breakpoints);

    // Initialize B-spline workspace to find coeffs from collocation points
    bspluzw = make_shared<sz::bspline_luz>(*bspw);
    bspluzw->form_mass(*bspw);

    // Initialize pencil_grid which handles P3DFFT setup/teardown RAII
    pg = make_shared<sz::pencil_grid>(def_grid.dealiased_extents(),
                                      def_grid.processor_grid());
    LOG4CXX_INFO(log, "Processor count: " << nproc);
    LOG4CXX_INFO(log, "Processor grid used: " << pg->processor_grid());
    LOG4CXX_DEBUG(log, "Local dealiased wave start  (XYZ): "
                       << pg->local_wave_start());
    LOG4CXX_DEBUG(log, "Local dealiased wave end    (XYZ): "
                       << pg->local_wave_end());
    LOG4CXX_DEBUG(log, "Local dealiased wave extent (XYZ): "
                       << pg->local_wave_extent());

    // Compute how much non-dealiased XYZ state is local to this rank
    // Additional munging necessary X direction has (Nx/2+1) complex values
    const boost::array<sz::pencil_grid::index,3> state_start
        = pg->local_wave_start();
    const boost::array<sz::pencil_grid::index,3> state_end = {{
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[0]/2+1,
                                             pg->local_wave_end()[0]),
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[1],
                                             pg->local_wave_end()[1]),
        std::min<sz::pencil_grid::size_type>(def_grid.global_extents()[2],
                                             pg->local_wave_end()[2])
    }};
    const boost::array<sz::pencil_grid::index,3> state_extent = {{
        std::max<sz::pencil_grid::index>(state_end[0] - state_start[0], 0),
        std::max<sz::pencil_grid::index>(state_end[1] - state_start[1], 0),
        std::max<sz::pencil_grid::index>(state_end[2] - state_start[2], 0)
    }};
    LOG4CXX_DEBUG(log, "Local state wave start  (XYZ): " << state_start);
    LOG4CXX_DEBUG(log, "Local state wave end    (XYZ): " << state_end);
    LOG4CXX_DEBUG(log, "Local state wave extent (XYZ): " << state_extent);

    // Create the state storage for the linear and nonlinear operators
    // with appropriate padding to allow nonlinear state to be P3DFFTified
    state_type state_linear(sz::to_yxz(5, state_extent));
    state_type state_nonlinear(
            sz::to_yxz(5, state_extent),
            sz::prepend(pg->local_wave_storage(),
                        sz::strides_cm(sz::to_yxz(pg->local_wave_extent()))));
    if (log->isDebugEnabled()) {
        boost::array<sz::pencil_grid::index,4> strides;
        std::copy(state_linear.strides(),
                  state_linear.strides() + 4, strides.begin());
        LOG4CXX_DEBUG(log, "Linear state strides    (FYXZ): " << strides);
        std::copy(state_nonlinear.strides(),
                  state_nonlinear.strides() + 4, strides.begin());
        LOG4CXX_DEBUG(log, "Nonlinear state strides (FYXZ): " << strides);
    }

    // Instantiate the operators and timestepping details
    // See write up section 2.1 (Spatial Discretization) for coefficient origin
    const sz::timestepper::lowstorage::SMR91Method<complex_t> smr91;
    MassOperator L(   def_scenario.Lx() * def_scenario.Lz()
                    * def_grid.Nx()     * def_grid.Nz());
    NonlinearOperator N;

    // Take a timestep
    sz::timestepper::lowstorage::step(smr91, L, N, state_linear, state_nonlinear);
}
