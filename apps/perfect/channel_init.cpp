//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 Rhys Ulerich
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

/** @file
 * Initialize perfect gas, channel restart files.
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <boost/math/special_functions/gamma.hpp>
#include <esio/error.h>
#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/pre_gsl.h>
#include <suzerain/support/definition_base.hpp>
#include <suzerain/support/field.hpp>
#include <suzerain/support/grid_definition.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/program_options.hpp>
#include <suzerain/support/support.hpp>
#include <suzerain/support/time_definition.hpp>
#include <suzerain/utility.hpp>
#include <suzerain/validation.hpp>
#include <suzerain/version.hpp>

#include "manufactured_solution.hpp"
#include "mean_quantities.hpp"
#include "perfect.hpp"

// Provided by channel_init_svnrev.{c,h} to speed recompilation
#pragma warning(push,disable:1419)
extern "C" const char revstr[];
#pragma warning(pop)

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
using boost::math::constants::pi;
using boost::numeric_cast;
using std::numeric_limits;
using suzerain::complex_t;
using suzerain::make_shared;
using suzerain::real_t;
namespace ndx     = suzerain::ndx;
namespace perfect = suzerain::perfect;
namespace support = suzerain::support;

// FIXME Generalize as part of Redmine ticket #2480
// We are only prepared to deal with 5 equations
static const std::vector<support::field> fields = perfect::default_fields();

// Global parameters initialized in main()
static suzerain::perfect::scenario_definition scenario(
        /* Re         */ 100,
        /* Ma         */ 1.5,
        /* Pr         */ 0.7,
        /* bulk_rho   */ 1,
        /* bulk_rho_u */ 1,
        /* alpha      */ 0,
        /* beta       */ 2.0 /3.0,
        /* gamma      */ 1.4);
static suzerain::support::grid_definition grid(
        /* Lx      */ 4 * pi<real_t>(),
        /* Nx      */ 1,
        /* DAFx    */ 1.5,
        /* Ly      */ 2,
        /* Ny      */ 32,
        /* k       */ 8,
        /* htdelta */ 3,
        /* Lz      */ 4.0 * pi<real_t>() / 3.0,
        /* Nz      */ 1,
        /* DAFz    */ 1.5);
static suzerain::support::time_definition timedef(
        /* evmagfactor per Venugopal */ 0.72);
static suzerain::shared_ptr<const suzerain::pencil_grid> dgrid;
static suzerain::shared_ptr<perfect::manufactured_solution> msoln(
            new perfect::manufactured_solution(
                  perfect::manufactured_solution::default_caption
                + " (active only when --mms supplied)"));

/** <tt>atexit</tt> callback to ensure we finalize underling. */
static void atexit_underling(void) {
    dgrid.reset();        // Runs pencil_grid destructors
#ifdef HAVE_UNDERLING
    underling_cleanup();  // Cleans up the library
#endif
}

// Global B-spline related-details initialized in main()
static suzerain::shared_ptr<suzerain::bspline>       b;
static suzerain::shared_ptr<suzerain::bsplineop>     cop;    // Collocation
static suzerain::shared_ptr<suzerain::bsplineop>     gop;    // Galerkin L2
static suzerain::shared_ptr<suzerain::bsplineop_luz> bopluz;

// Explicit timestepping scheme uses only complex_t 4D contiguous_state
// State indices range over (scalar field, Y, X, Z) in wave space
typedef suzerain::contiguous_state<4,complex_t> state_type;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI...
    atexit((void (*) ()) MPI_Finalize);             // ...finalize at exit
    suzerain::support::logging::initialize(         // Initialize logging
            MPI_COMM_WORLD, support::log4cxx_config);
#ifdef HAVE_UNDERLING
    underling_init(&argc, &argv, 0);                // Initialize underling...
#endif
    atexit(atexit_underling);                       // ...finalize at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // ...finalize at exit

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

    // Process incoming program arguments from command line, input files
    std::string restart_file;
    bool   clobber = false;
    real_t mms     = -1;
    real_t npower  = 1;
    {
        suzerain::support::program_options options(
                "Suzerain-based compressible channel initialization",
                "RESTART-FILE", /* TODO description */ "", revstr);

        namespace po = boost::program_options;
        std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_nonnegative(
                    suzerain::validation::ensure_nonnegative<real_t>);

        options.add_definition(scenario);
        options.add_definition(grid);
        options.add_definition(timedef);
        options.add_definition(msoln->isothermal_channel());
        options.add_options()
            ("clobber", "Overwrite an existing restart file?")
            ("npower",
             boost::program_options::value(&npower)->default_value(npower),
             "Power n in (0, 1] used to control the flatness of the"
             " \"parabolic\" streamwise velocity profile (y*(L-y))^n.")
            ("mms",
             boost::program_options::value(&mms)
                ->notifier(std::bind2nd(ptr_fun_ensure_nonnegative, "mms")),
             "If given, prepare a manufactured solution at the specified time.")
        ;

        std::vector<std::string> positional = options.process(argc, argv);
        msoln->match(scenario);
        msoln->match(grid);

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

        clobber = options.variables().count("clobber");
    }

    if (npower <= 0 || npower > 1) {
        FATAL("npower in (0,1] required");
        return EXIT_FAILURE;
    }

    if (grid.k < 4 /* cubics */) {
        FATAL("k >= 4 required for two non-trivial spatial derivatives");
        return EXIT_FAILURE;
    }

    // Initialization done under assumptions bulk_rho == 1 && bulk_rho_u == 1
    if (scenario.bulk_rho_u != 1) {
        WARN0("Forcing bulk streamwise momentum to be one");
        scenario.bulk_rho_u = 1;
    }
    if (scenario.bulk_rho != 1) {
        WARN0("Forcing bulk density to be one");
        scenario.bulk_rho = 1;
    }

    if (mms >= 0) {
        INFO0("Manufactured solution will be initialized at t = " << mms);
        INFO0("Disabling bulk_rho and bulk_rho_u constraints"
              " due to manufactured solution use");
        scenario.bulk_rho   = numeric_limits<real_t>::quiet_NaN();
        scenario.bulk_rho_u = numeric_limits<real_t>::quiet_NaN();
    } else {
        msoln.reset();
    }

    // Modify IEEE settings after startup complete as startup relies on NaNs
    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(suzerain::mpi::comm_rank(MPI_COMM_WORLD));

    support::create(grid.N.y(), grid.k, 0.0, grid.L.y(), grid.htdelta, b, cop);
    gop.reset(new suzerain::bsplineop(*b, 0, SUZERAIN_BSPLINEOP_GALERKIN_L2));

    INFO0("Creating new restart file " << restart_file);
    esio_file_create(esioh, restart_file.c_str(), clobber);
    esio_string_set(esioh, "/", "generated_by",
                    (std::string("channel ") + revstr).c_str()); // Ticket #2595
    perfect::save(esioh, scenario);
    support::save(esioh, grid);
    support::save(esioh, b, cop, gop);
    support::save(esioh, timedef);
    perfect::save(esioh, msoln, scenario, grid);
    esio_file_flush(esioh);

    INFO0("Initializing B-spline workspaces");
    bopluz = make_shared<suzerain::bsplineop_luz>(*cop);
    bopluz->factor_mass(*cop);

    INFO0("Initializing pencil_grid to obtain parallel decomposition details");
    dgrid = make_shared<suzerain::pencil_grid_default>(grid.dN, grid.P);
    assert((grid.dN == dgrid->global_physical_extent).all());


    INFO0("Allocating storage for the distributed state fields");
    state_type swave(suzerain::to_yxz(fields.size(), dgrid->local_wave_extent));
    state_type stemp(suzerain::to_yxz(1, dgrid->local_wave_extent));

    INFO0("Initializing data on collocation points values in physical space");
    if (mms >= 0) {

        // Use a canned manufactured solution routine for initialization
        perfect::accumulate_manufactured_solution(
                1, *msoln, 0, swave, grid, *dgrid, *cop,*b, mms);

    } else {

        // Use a simple parabolic velocity profile

        // Initializing operator_base to access decomposition-ready utilities
        suzerain::operator_base o(grid, *dgrid, *cop, *b);

        // State viewed as a 2D Eigen::Map ordered (F, Y*Z*X).
        suzerain::physical_view<> sphys(*dgrid, swave);

        // Find normalization required to have (y*(L-y))^npower integrate to one
        real_t factor;
        if (npower == 1) {
            // Mathematica: (Integrate[(x (L-x)),{x,0,L}]/L)^(-1)
            factor = 6 / std::pow(grid.L.y(), 2);
        } else {
            // Mathematica: (Integrate[(x (L - x))^n, {x, 0, L}]/L)^(-1)
            //      -  (Gamma[-n] Gamma[3/2+n])
            //       / (2^(-1-2 n) L^(1+2 n) \[Pi]^(3/2) Csc[n \[Pi]])
            using boost::math::constants::pi;
            const real_t num1   = std::sin(npower * pi<real_t>());
            const real_t num2   = boost::math::tgamma(-npower);
            const real_t num3   = boost::math::tgamma(real_t(3)/2 + npower);
            const real_t denom1 = std::pow(         2, -1-2*npower);
            const real_t denom2 = std::pow(grid.L.y(),  1+2*npower);
            const real_t denom3 = std::pow(pi<real_t>(),  real_t(3)/2);
            factor = - (num1 * num2 * num3 * grid.L.y())
                   /   (denom1 * denom2 * denom3);
        }

        // Physical space is traversed linearly using a single offset 'offset'.
        // The three loop structure is present to provide the global absolute
        // positions x(i), y(j), and z(k) where necessary.
        for (int offset = 0, j = dgrid->local_physical_start.y();
             j < dgrid->local_physical_end.y();
             ++j) {

            const real_t y = o.y(j);

            for (int k = dgrid->local_physical_start.z();
                k < dgrid->local_physical_end.z();
                ++k) {

                for (int i = dgrid->local_physical_start.x();
                    i < dgrid->local_physical_end.x();
                    ++i, /* NB */ ++offset) {

                    // Initialize primitive state
                    const real_t rho = 1;
                    const real_t u   = factor
                                     * std::pow(y * (grid.L.y() - y), npower);
                    const real_t v   = 0;
                    const real_t w   = 0;
                    const real_t T   = 1;

                    // Compute and store the conserved state from primitives
                    const real_t e = T / (scenario.gamma*(scenario.gamma - 1))
                                   +  (scenario.Ma*scenario.Ma/2)
                                     *(u*u + v*v + w*w);
                    sphys(ndx::e,   offset) = rho * e;
                    sphys(ndx::mx,  offset) = rho * u;
                    sphys(ndx::my,  offset) = rho * v;
                    sphys(ndx::mz,  offset) = rho * w;
                    sphys(ndx::rho, offset) = rho;

                } // end X

            } // end Z

        } // end Y

        // Build FFT normalization constant into Y direction's mass matrix
        suzerain::bsplineop_luz massluz(*cop);
        const complex_t scale_factor = grid.dN.x() * grid.dN.z();
        massluz.opform(1, &scale_factor, *cop);
        massluz.factor();

        for (std::size_t i = 0; i < swave.shape()[0]; ++i) {
            dgrid->transform_physical_to_wave(&sphys.coeffRef(i, 0));  // X, Z
            o.bop_solve(massluz, swave, i);                            // Y
        }

    }

    INFO0("Writing state fields to restart file");
    support::save_coefficients(esioh, fields, swave, grid, *dgrid);
    esio_file_flush(esioh);

    real_t t;
    if (mms < 0) {
        INFO0("Storing new simulation time of zero");
        t = 0;
    } else {
        INFO0("Storing simulation time to match manufactured solution");
        t = mms;
    }
    support::save_time(esioh, t);
    esio_file_flush(esioh);

    INFO0("Computing mean quantities from state fields");
    perfect::mean_quantities sample = perfect::sample_mean_quantities(
            scenario, grid, *dgrid, *cop, swave, t);

    INFO0("Writing mean quantities to restart file");
    perfect::save(esioh, sample);
    esio_file_flush(esioh);

    INFO0("Closing newly initialized restart file");
    esio_file_close(esioh);
}
