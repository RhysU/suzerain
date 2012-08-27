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
// channel_ex.cpp: Initialize restart files for use with Suzerain
// $Id$

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <boost/math/special_functions/gamma.hpp>
#include <esio/error.h>
#include <esio/esio.h>
#include <suzerain/error.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pre_gsl.h>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/utility.hpp>
#include <suzerain/version.hpp>

#include "../logging.hpp"
#include "../precision.hpp"
#include "../support.hpp"

// Provided by channel_init_svnrev.{c,h} to speed recompilation
extern "C" const char revstr[];

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
using boost::make_shared;
using boost::math::constants::pi;
using boost::numeric_cast;
using boost::shared_ptr;
using std::numeric_limits;

// Global parameters initialized in main()
using suzerain::problem::ScenarioDefinition;
using suzerain::problem::GridDefinition;
using suzerain::problem::TimeDefinition;
static ScenarioDefinition<real_t> scenario(
        /* Re        */ "100",
        /* Ma        */ "1.5",
        /* Pr        */ "0.7",
        /* bulk_rho  */ "1",
        /* bulk_rhou */ "1",
        /* alpha     */ "0",
        /* beta      */ "2/3",
        /* gamma     */ "1.4",
        /* Lx        */ "4*pi",
        /* Ly        */ "2",
        /* Lz        */ "4*pi/3");
static GridDefinition grid(
        /* Nx      */ 1,
        /* DAFx    */ 1.5,
        /* Ny      */ 32,
        /* k       */ 8,
        /* htdelta */ 3,
        /* Nz      */ 1,
        /* DAFz    */ 1.5);
static TimeDefinition<real_t> timedef(
        /* evmagfactor per Venugopal */ "0.72");
static shared_ptr<const suzerain::pencil_grid> dgrid;
static shared_ptr<channel::manufactured_solution> msoln(
            new channel::manufactured_solution);

/** <tt>atexit</tt> callback to ensure we finalize underling. */
static void atexit_underling(void) {
    dgrid.reset();        // Runs pencil_grid destructors
#ifdef HAVE_UNDERLING
    underling_cleanup();  // Cleans up the library
#endif
}

// Global B-spline related-details initialized in main()
static shared_ptr<suzerain::bspline>       b;
static shared_ptr<suzerain::bsplineop>     bop;    // Collocation
static shared_ptr<suzerain::bsplineop>     gop;    // Galerkin L2
static shared_ptr<suzerain::bsplineop_luz> bopluz;

// Explicit timestepping scheme uses only complex_t 4D ContiguousState
// State indices range over (scalar field, Y, X, Z) in wave space
typedef suzerain::ContiguousState<4,complex_t> state_type;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

/** Options definitions for tweaking the manufactured solution */
class MSDefinition : public suzerain::problem::IDefinition {

public:

    MSDefinition(channel::manufactured_solution &ms)
        : IDefinition("Manufactured solution parameters"
                      " (active only when --mms supplied)")
    {
        ms.rho.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects density field",     _1, _2));
        ms.u.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects X velocity field",  _1, _2));
        ms.v.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects Y velocity field",  _1, _2));
        ms.w.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects Z velocity field",  _1, _2));
        ms.T.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects temperature field", _1, _2));
    }

private:

    static void option_adder(
            boost::program_options::options_description_easy_init ei,
            const char *desc,
            const std::string &name,
            real_t &v)
    {
        ei(name.c_str(),
           boost::program_options::value(&v)->default_value(v),
           desc);
    }
};

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI...
    atexit((void (*) ()) MPI_Finalize);             // ...finalize at exit
    logging::initialize(MPI_COMM_WORLD,             // Initialize logging
                        channel::log4cxx_config);
#ifdef HAVE_UNDERLING
    underling_init(&argc, &argv, 0);                // Initialize underling...
#endif
    atexit(atexit_underling);                       // ...finalize at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // ...finalize at exit

    // Hook error handling into logging infrastructure
    gsl_set_error_handler(
            &channel::mpi_abort_on_error_handler_gsl);
    suzerain_set_error_handler(
            &channel::mpi_abort_on_error_handler_suzerain);
    esio_set_error_handler(
            &channel::mpi_abort_on_error_handler_esio);
#ifdef HAVE_UNDERLING
    underling_set_error_handler(
            &channel::mpi_abort_on_error_handler_underling);
#endif

    // Process incoming program arguments from command line, input files
    std::string restart_file;
    bool   clobber = false;
    real_t mms     = -1;
    real_t npower  = 1;
    {
        suzerain::ProgramOptions options(
                "Suzerain-based compressible channel initialization",
                "RESTART-FILE", /* TODO description */ "", revstr);

        namespace po = ::boost::program_options;
        using ::suzerain::validation::ensure_nonnegative;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_nonnegative(ensure_nonnegative<real_t>);

        isothermal_channel(*msoln);
        MSDefinition msdef(*msoln);

        options.add_definition(scenario);
        options.add_definition(grid);
        options.add_definition(timedef);
        options.add_definition(msdef);
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

    // Initialization done under assumptions bulk_rho == 1 && bulk_rhou == 1
    if (scenario.bulk_rhou != 1) {
        WARN0("Forcing bulk streamwise momentum to be one");
        scenario.bulk_rhou = 1;
    }
    if (scenario.bulk_rho != 1) {
        WARN0("Forcing bulk density to be one");
        scenario.bulk_rho = 1;
    }

    if (mms >= 0) {
        INFO0("Manufactured solution will be initialized at t = " << mms);
        msoln->alpha = scenario.alpha;
        msoln->beta  = scenario.beta;
        msoln->gamma = scenario.gamma;
        msoln->Ma    = scenario.Ma;
        msoln->Re    = scenario.Re;
        msoln->Pr    = scenario.Pr;
        msoln->Lx    = scenario.Lx;
        msoln->Ly    = scenario.Ly;
        msoln->Lz    = scenario.Lz;

        INFO0("Disabling bulk_rho and bulk_rhou constraints"
              " due to manufactured solution use");
        scenario.bulk_rho  = numeric_limits<real_t>::quiet_NaN();
        scenario.bulk_rhou = numeric_limits<real_t>::quiet_NaN();
    } else {
        msoln.reset();
    }

    // Modify IEEE settings after startup complete as startup relies on NaNs
    DEBUG0("Establishing floating point environment from GSL_IEEE_MODE");
    mpi_gsl_ieee_env_setup(suzerain::mpi::comm_rank(MPI_COMM_WORLD));

    channel::create(grid.N.y(), grid.k, 0.0, scenario.Ly, grid.htdelta, b, bop);
    gop.reset(new suzerain::bsplineop(*b, 0, SUZERAIN_BSPLINEOP_GALERKIN_L2));

    INFO0("Creating new restart file " << restart_file);
    esio_file_create(esioh, restart_file.c_str(), clobber);
    esio_string_set(esioh, "/", "generated_by",
                    (std::string("channel_init ") + revstr).c_str());
    channel::store(esioh, scenario);
    channel::store(esioh, grid, scenario.Lx, scenario.Lz);
    channel::store(esioh, b, bop, gop);
    channel::store(esioh, timedef);
    channel::store(esioh, scenario, msoln);
    esio_file_flush(esioh);

    INFO0("Initializing B-spline workspaces");
    bopluz = make_shared<suzerain::bsplineop_luz>(*bop);
    bopluz->factor_mass(*bop);

    INFO0("Initializing pencil_grid to obtain parallel decomposition details");
    dgrid = make_shared<suzerain::pencil_grid_default>(grid.dN, grid.P);
    assert((grid.dN == dgrid->global_physical_extent).all());


    INFO0("Allocating storage for the distributed state fields");
    state_type swave(suzerain::to_yxz(
                channel::field::count, dgrid->local_wave_extent));
    state_type stemp(suzerain::to_yxz(1, dgrid->local_wave_extent));

    INFO0("Initializing data on collocation points values in physical space");
    if (mms >= 0) {

        // Use a canned manufactured solution routine for initialization
        channel::accumulate_manufactured_solution(
                1, *msoln, 0, swave, scenario, grid, *dgrid, *b, *bop, mms);

    } else {

        // Use a simple parabolic velocity profile

        // Initializing OperatorBase to access decomposition-ready utilities
        suzerain::OperatorBase<real_t> obase(scenario, grid, *dgrid, *b, *bop);

        // State viewed as a 2D Eigen::Map ordered (F, Y*Z*X).
        channel::physical_view<channel::field::count>::type sphys
            = channel::physical_view<channel::field::count>::create(*dgrid, swave);

        // Find normalization required to have (y*(L-y))^npower integrate to one
        real_t factor;
        if (npower == 1) {
            // Mathematica: (Integrate[(x (L-x)),{x,0,L}]/L)^(-1)
            factor = 6 / std::pow(scenario.Ly, 2);
        } else {
            // Mathematica: (Integrate[(x (L - x))^n, {x, 0, L}]/L)^(-1)
            //      -  (Gamma[-n] Gamma[3/2+n])
            //       / (2^(-1-2 n) L^(1+2 n) \[Pi]^(3/2) Csc[n \[Pi]])
            using boost::math::constants::pi;
            const real_t num1   = std::sin(npower * pi<real_t>());
            const real_t num2   = boost::math::tgamma(-npower);
            const real_t num3   = boost::math::tgamma(real_t(3)/2 + npower);
            const real_t denom1 = std::pow(           2, -1-2*npower);
            const real_t denom2 = std::pow( scenario.Ly,  1+2*npower);
            const real_t denom3 = std::pow(pi<real_t>(),  real_t(3)/2);
            factor = - (num1 * num2 * num3 * scenario.Ly)
                   /   (denom1 * denom2 * denom3);
        }

        // Physical space is traversed linearly using a single offset 'offset'.
        // The three loop structure is present to provide the global absolute
        // positions x(i), y(j), and z(k) where necessary.
        size_t offset = 0;
        for (int j = dgrid->local_physical_start.y();
             j < dgrid->local_physical_end.y();
             ++j) {

            const real_t y = obase.y(j);

            for (int k = dgrid->local_physical_start.z();
                k < dgrid->local_physical_end.z();
                ++k) {

                for (int i = dgrid->local_physical_start.x();
                    i < dgrid->local_physical_end.x();
                    ++i, /* NB */ ++offset) {

                    // Initialize primitive state
                    const real_t rho = 1;
                    const real_t u   = factor
                                     * std::pow(y * (scenario.Ly - y), npower);
                    const real_t v   = 0;
                    const real_t w   = 0;
                    const real_t T   = 1;

                    // Compute and store the conserved state from primitives
                    const real_t e = T / (scenario.gamma*(scenario.gamma - 1))
                                   + (scenario.Ma*scenario.Ma/2)*(u*u + v*v + w*w);
                    namespace ndx = channel::field::ndx;
                    sphys(channel::field::ndx::rho,  offset) = rho;
                    sphys(channel::field::ndx::rhou, offset) = rho * u;
                    sphys(channel::field::ndx::rhov, offset) = rho * v;
                    sphys(channel::field::ndx::rhow, offset) = rho * w;
                    sphys(channel::field::ndx::rhoe, offset) = rho * e;

                } // end X

            } // end Z

        } // end Y

        // Build FFT normalization constant into Y direction's mass matrix
        suzerain::bsplineop_luz massluz(*bop);
        const complex_t scale_factor = grid.dN.x() * grid.dN.z();
        massluz.opform(1, &scale_factor, *bop);
        massluz.factor();

        for (std::size_t i = 0; i < channel::field::count; ++i) {
            dgrid->transform_physical_to_wave(&sphys.coeffRef(i, 0));  // X, Z
            obase.bop_solve(massluz, swave, i);                        // Y
        }

    }

    INFO0("Writing state fields to restart file");
    channel::store_coefficients(esioh, swave, scenario, grid, *dgrid);
    esio_file_flush(esioh);

    real_t t;
    if (mms < 0) {
        INFO0("Storing new simulation time of zero");
        t = 0;
    } else {
        INFO0("Storing simulation time to match manufactured solution");
        t = mms;
    }
    channel::store_time(esioh, t);
    esio_file_flush(esioh);

    INFO0("Computing mean quantities from state fields");
    channel::mean samples = channel::sample_mean_quantities(
            scenario, grid, *dgrid, *b, *bop, swave, t);

    INFO0("Writing mean quantities to restart file");
    channel::store(esioh, samples);
    esio_file_flush(esioh);

    INFO0("Closing newly initialized restart file");
    esio_file_close(esioh);
}
