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
#include <esio/esio.h>
#include <esio/error.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <suzerain/error.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/rholt.hpp>
#include <suzerain/svehla.h>
#include <suzerain/utility.hpp>

#include "logger.hpp"
#include "precision.hpp"
#include "channel.hpp"

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
static ScenarioDefinition<real_t> scenario(
        /* Re        */ 100,
        /* Pr        */ real_t(7)/real_t(10),
        /* Ma        */ real_t(115)/real_t(100),
        /* bulk_rho  */ 1,
        /* bulk_rhou */ real_t(115)/real_t(100),
        /* alpha     */ real_t(0),
        /* beta      */ real_t(2)/real_t(3),
        /* gamma     */ real_t(14)/real_t(10),
        /* Lx        */ 4*pi<real_t>(),
        /* Ly        */ 2,
        /* Lz        */ 4*pi<real_t>()/3);
static GridDefinition grid(
        /* Nx      */ 1,
        /* DAFx    */ 1.5,
        /* Ny      */ 16,
        /* k       */ 6,
        /* htdelta */ 3,
        /* Nz      */ 1,
        /* DAFz    */ 1.5);
static shared_ptr<const suzerain::pencil_grid> dgrid;

// Global B-spline related-details initialized in main()
static shared_ptr<suzerain::bspline>       b;
static shared_ptr<suzerain::bsplineop>     bop;    // Collocation
static shared_ptr<suzerain::bsplineop>     gop;    // Galerkin L2
static shared_ptr<suzerain::bsplineop_luz> bopluz;

// Explicit timestepping scheme uses only complex_t 4D NoninterleavedState
// State indices range over (scalar field, Y, X, Z) in wave space
typedef suzerain::NoninterleavedState<4,complex_t> state_type;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

/** Description of equation to find wall temperature */
static const char tsolver_desc[]
    = "mu(T) Re / (p / R / T) / L - Ma ( gamma R T )**(1/2) == 0";

/** Struct used for iterative solver to find wall temperature */
struct tsolver {
    double Re, L, gamma, R, Ma, p_wall;
    gsl_spline *s;
};

/** Function for iterative solver to find wall temperature */
static double f_tsolver(double T_wall, void *params) {
    tsolver *p            = (tsolver *) params;
    const double rho_wall = p->p_wall / (p->R * T_wall);
    const double mu       = gsl_spline_eval(p->s, T_wall, NULL);
    const double lhs      = mu * p->Re / (rho_wall * p->L);
    const double rhs      = p->Ma * sqrt(p->gamma * p->R * T_wall);
    return lhs - rhs;
}

/** Struct used for evaluating momentum and total energy profiles */
struct mesolver {
    double Ma, L, gamma, R, T_wall;
};

/**
 * Function for evaluating the momentum profile given rho* = 1.  Complex-valued
 * result to play nicely with suzerain::bspline::zfind_interpolation_rhs().
 */
static void zf_msolver(double y, void *params, double z[2]) {
    mesolver *p = (mesolver *) params;
    z[0] = p->Ma * 6 * y * (p->L - y) / (p->L * p->L);
    z[1] = 0;
}

/**
 * Function for evaluating the total energy profile given rho*, T* = 1.
 * Complex-valued result to play nicely with
 * suzerain::bspline::zfind_interpolation_rhs().
 */
static void zf_esolver(double y, void *params, double z[2]) {
    mesolver *p = (mesolver *) params;
    double m[2];
    zf_msolver(y, p, m);
    z[0] = 1 / (p->gamma * (p->gamma - 1)) + 0.5 * m[0] * m[0];
    z[1] = 0;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // Finalize ESIO at exit

    // Establish MPI-savvy, rank-dependent logging names
    name_logger_within_comm_world();

    // Hook error handling into logging infrastructure
    gsl_set_error_handler(
            &channel::mpi_abort_on_error_handler_gsl);
    suzerain_set_error_handler(
            &channel::mpi_abort_on_error_handler_suzerain);
    esio_set_error_handler(
            &channel::mpi_abort_on_error_handler_esio);

    // Process incoming program arguments from command line, input files
    std::string restart_file;
    bool clobber  = false;
    real_t M      = SUZERAIN_SVEHLA_AIR_M;
    real_t p_wall = GSL_CONST_MKSA_STD_ATMOSPHERE / 100000;
    {
        suzerain::ProgramOptions options(
                "Suzerain-based compressible channel initialization",
                "[RESTART-FILE]");
        namespace po = ::boost::program_options;

        // Cast away const so options processing can modify settings
        options.add_definition(scenario);
        options.add_definition(grid);

        using ::suzerain::validation::ensure_positive;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_positive(ensure_positive<real_t>);
        using ::suzerain::validation::ensure_nonnegative;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_nonnegative(ensure_nonnegative<real_t>);

        options.add_options()
            ("clobber", "Overwrite an existing restart file?")
            ("M", po::value<real_t>(&M)
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"M"))
                ->default_value(M),
             "Molar mass of species in grams per mole")
            ("p_wall", po::value<real_t>(&p_wall)
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"p_wall"))
                ->default_value(p_wall),
             "Pressure in N/m^2 used to obtain reference wall density")
        ;
        std::vector<std::string> positional = options.process(argc, argv);

        if (positional.size() != 1) {
            FATAL0("Exactly one restart file name must be specified");
            return EXIT_FAILURE;
        }
        restart_file = positional[0];

        clobber = options.variables().count("clobber");
    }
    const real_t R = GSL_CONST_MKSA_MOLAR_GAS * GSL_CONST_NUM_KILO / M;

    if (grid.k < 4 /* cubics */) {
        FATAL("k >= 4 required for two non-trivial spatial derivatives");
        return EXIT_FAILURE;
    }

    // Initialization done under assumptions bulk_rho == 1 && Ma == bulk_rhou
    if (scenario.bulk_rhou != scenario.Ma) {
        WARN0("Forcing bulk streamwise momentum to match the Mach number");
        scenario.bulk_rhou = scenario.Ma;
    }
    if (scenario.bulk_rho != 1) {
        WARN0("Forcing bulk density to be one");
        scenario.bulk_rho = 1;
    }

    INFO0("Creating B-spline basis of order " << (grid.k - 1)
          << " on [0, " << scenario.Ly << "] with "
          << grid.N.y() << " DOF stretched per htdelta " << grid.htdelta);
    channel::create(grid.N.y(), grid.k, 0.0, scenario.Ly, grid.htdelta, b, bop);
    gop.reset(new suzerain::bsplineop(*b, 0, SUZERAIN_BSPLINEOP_GALERKIN_L2));

    INFO0("Creating new restart file " << restart_file);
    esio_file_create(esioh, restart_file.c_str(), clobber);
    channel::store(esioh, scenario);
    channel::store(esioh, grid, scenario.Lx, scenario.Lz);
    channel::store(esioh, b, bop, gop);
    esio_file_flush(esioh);

    INFO0("Computing derived, dimensional reference parameters");
    real_t T_wall, rho_wall;
    {
        if (scenario.gamma != 1.4) {
            WARN0("Using air viscosity values for non-air!");
        }

        tsolver p;
        p.Re     = scenario.Re;
        p.L      = scenario.Ly;
        p.gamma  = scenario.gamma;
        p.R      = R;
        p.Ma     = scenario.Ma;
        p.p_wall = p_wall;
        p.s = suzerain_svehla_air_mu_vs_T();
        assert(p.s);

        INFO0("Solving " << tsolver_desc << " for T_wall");
        INFO0("Using Re = " << p.Re
              << ", L = " << p.L
              << ", gamma = " << p.gamma
              << ", R = " << p.R
              << ", Ma = " << p.Ma
              << ", p_wall = " << p.p_wall
              << ", and mu(T) from Svehla 1962");

        gsl_function F;
        F.function = &f_tsolver;
        F.params   = &p;

        gsl_root_fsolver * solver
            = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        assert(solver);
        real_t low  = 100;  // (degrees K) Low side of Svehla's data
        real_t high = 5000; // (degrees K) High side of Svehla's data
        {
            const real_t low_residual = GSL_FN_EVAL(&F,low);
            const real_t high_residual = GSL_FN_EVAL(&F,high);
            if ((low_residual >= 0) == (high_residual >= 0)) {
                FATAL0("Likely no solution in interval ["
                       << low <<"," << high << "]!");
            }
        }
        gsl_root_fsolver_set(solver, &F, low, high);

        const real_t tol = numeric_limits<real_t>::epsilon() * 1e4;
        int status = GSL_CONTINUE;
        for (int iter = 0; iter < 100 && status == GSL_CONTINUE; ++iter) {
            gsl_root_fsolver_iterate(solver);
            T_wall = gsl_root_fsolver_root(solver);
            low  = gsl_root_fsolver_x_lower(solver);
            high = gsl_root_fsolver_x_upper(solver);
            status = gsl_root_test_interval(low, high, tol, 0);
        }
        INFO0(std::setprecision(numeric_limits<real_t>::digits10)
              << "Wall temperature "
              << T_wall
              << " K gives residual "
              << std::scientific
              << GSL_FN_EVAL(&F,T_wall));

        rho_wall = p_wall / R / T_wall;
        INFO0("Wall density is " << rho_wall << " kg/m^3");

        gsl_root_fsolver_free(solver);
        gsl_spline_free(p.s);
    }

    // Initialize B-splines to find coefficients from collocation points
    bopluz = make_shared<suzerain::bsplineop_luz>(*bop);
    bopluz->form_mass(*bop);

    // Initialize pencil_grid to obtain parallel decomposition details
    dgrid = make_shared<suzerain::pencil_grid>(grid.dN, grid.P);
    assert((grid.dN == dgrid->global_physical_extent).all());

    // Allocate and clear storage for the distributed state
    state_type state(suzerain::to_yxz(
                channel::field::count, dgrid->local_wave_extent));
    suzerain::multi_array::fill(state, 0);

    INFO0("Computing mean, nondimensional profiles");
    if (dgrid->local_wave_start.x() == 0 && dgrid->local_wave_start.z() == 0) {
        namespace ndx = channel::field::ndx;

        // Create 1D mean views from 4D state storage
        const boost::multi_array_types::index_range all;
        boost::array_view_gen<state_type,1>::type mean_rho
                = state[boost::indices[ndx::rho][all][0][0]];
        boost::array_view_gen<state_type,1>::type mean_rhou
                = state[boost::indices[ndx::rhou][all][0][0]];
        boost::array_view_gen<state_type,1>::type mean_rhov
                = state[boost::indices[ndx::rhov][all][0][0]];
        boost::array_view_gen<state_type,1>::type mean_rhow
                = state[boost::indices[ndx::rhow][all][0][0]];
        boost::array_view_gen<state_type,1>::type mean_rhoe
                = state[boost::indices[ndx::rhoe][all][0][0]];

        // Prepare parameter struct for evaluation routines
        mesolver params;
        params.Ma     = scenario.Ma;
        params.L      = scenario.Ly;
        params.gamma  = scenario.gamma;
        params.R      = R;
        params.T_wall = T_wall;

        // Nondimensional density is the constant one
        suzerain::multi_array::fill(mean_rho, 1);

        // Find streamwise momentum as a function of wall-normal position
        const suzerain_zfunction zF_rhou = { &zf_msolver, &params };
        bop->interpolation_rhs(&zF_rhou, &mean_rhou[0], *b);
        bopluz->solve(1, &mean_rhou[0], 1, grid.N.y());

        // Nondimensional wall-normal momentum is zero
        suzerain::multi_array::fill(mean_rhov, 0);

        // Nondimensional spanwise momentum is zero
        suzerain::multi_array::fill(mean_rhow, 0);

        // Find total energy as a function of wall-normal position
        const suzerain_zfunction zF_rhoe = { &zf_esolver, &params };
        bop->interpolation_rhs(&zF_rhoe, &mean_rhoe[0], *b);
        bopluz->solve(1, &mean_rhoe[0], 1, grid.N.y());
    }

    INFO0("Writing state fields to restart file");
    channel::store(esioh, state, grid, *dgrid);
    esio_file_flush(esioh);

    // Store new simulation time of zero
    channel::store_time(esioh, 0);
    esio_file_flush(esioh);

    INFO0("Closing newly initialized restart file");
    esio_file_close(esioh);
}
