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
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/htstretch.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/scenario_definition.hpp>
#include <suzerain/svehla.h>
#include <suzerain/utility.hpp>

#include "logger.hpp"
#include "channel_common.hpp"

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
using boost::make_shared;
using boost::math::constants::pi;
using boost::numeric_cast;
using boost::scoped_array;
using boost::shared_ptr;
using std::numeric_limits;

// Global scenario parameters initialized in main()
using suzerain::problem::ScenarioDefinition;
using suzerain::problem::GridDefinition;
static const ScenarioDefinition<real_t> scenario(
        /* Re    */ 100,
        /* Pr    */ real_t(7)/real_t(10),
        /* gamma */ real_t(14)/real_t(10),
        /* beta  */ real_t(2)/real_t(3),
        /* Lx    */ 4*pi<real_t>(),
        /* Ly    */ 2,
        /* Lz    */ 4*pi<real_t>()/3);
static const GridDefinition<real_t> grid(
        /* Nx    */ 1,
        /* DAFx  */ real_t(3)/real_t(2),
        /* Ny    */ 16,
        /* k     */ 6,
        /* Nz    */ 1,
        /* DAFz  */ real_t(3)/real_t(2));

// Global B-spline -details initialized in main()
static shared_ptr<suzerain::bspline>     bspw;
static shared_ptr<suzerain::bspline_luz> bspluzw;

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

/** Function for evaluating the momentum profile given rho* = 1 */
static double f_msolver(double y, void *params) {
    mesolver *p = (mesolver *) params;
    return p->Ma * 6 * y * (p->L - y) / (p->L * p->L);
}

/** Function for evaluating the total energy profile given rho*, T* = 1 */
static double f_esolver(double y, void *params) {
    mesolver *p = (mesolver *) params;
    const double m = f_msolver(y, p);
    return 1 / (p->gamma * (p->gamma - 1)) + 0.5 * m * m;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit
    esioh = esio_handle_initialize(MPI_COMM_WORLD); // Initialize ESIO
    atexit(&atexit_esio);                           // Finalize ESIO at exit

    // Obtain some basic MPI environment details.
    const int nranks = suzerain::mpi::comm_size(MPI_COMM_WORLD);

    // Establish MPI-savvy, rank-dependent logging names
    name_logger_within_comm_world();

    // Ensure that we're running in a single processor environment
    if (nranks > 1) {
        FATAL(argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // Process incoming program arguments from command line, input files
    real_t Ma     = 1.15;
    real_t M      = SUZERAIN_SVEHLA_AIR_M;
    real_t p_wall = GSL_CONST_MKSA_STD_ATMOSPHERE / 100000;
    std::string restart_file;
    real_t htdelta;
    suzerain::ProgramOptions options(
            "Suzerain-based compressible channel initialization");
    {
        namespace po = ::boost::program_options;

        // Cast away const so options processing can modify settings
        options.add_definition(
                const_cast<ScenarioDefinition<real_t>& >(scenario));
        options.add_definition(
                const_cast<GridDefinition<real_t>& >(grid));

        using ::suzerain::validation::ensure_positive;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_positive(ensure_positive<real_t>);
        using ::suzerain::validation::ensure_nonnegative;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_nonnegative(ensure_nonnegative<real_t>);

        options.add_options()
            ("create", po::value<std::string>(&restart_file)
                ->default_value("restart0.h5"),
             "Name of new restart file to create")
            ("clobber", "Overwrite existing restart file?")
            ("Ma", po::value<real_t>(&Ma)
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"Ma"))
                ->default_value(Ma),
             "Mach number based on bulk velocity and wall sound speed")
            ("M", po::value<real_t>(&M)
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"M"))
                ->default_value(M),
             "Molar mass of species in grams per mole")
            ("p_wall", po::value<real_t>(&p_wall)
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"p_wall"))
                ->default_value(p_wall),
             "Pressure in N/m^2 used to obtain reference wall density")
            ("htdelta", po::value<real_t>(&htdelta)
                ->notifier(std::bind2nd(ptr_fun_ensure_nonnegative,"htdelta"))
                ->default_value(7),
             "Hyperbolic tangent stretching parameter for breakpoints")
        ;
        options.process(argc, argv);
    }
    const real_t R = GSL_CONST_MKSA_MOLAR_GAS * GSL_CONST_NUM_KILO / M;

    if (grid.k < 4 /* cubics */) {
        FATAL("k >= 4 required to compute two non-trivial spatial derivatives");
        return EXIT_FAILURE;
    }

    if (grid.Nx != 1) {
        FATAL(argv[0] << " can only handle Nx == 1");
        return EXIT_FAILURE;
    }

    if (grid.Nz != 1) {
        FATAL(argv[0] << " can only handle Nz == 1");
        return EXIT_FAILURE;
    }

    INFO("Creating B-spline basis of uniform order "
         << (grid.k - 1) << " on [0, Ly] with "
         << grid.Ny << " DOF");
    {
        scoped_array<real_t> buf(new real_t[grid.Ny]);

        // Compute breakpoint locations
        const int nbreak = grid.Ny + 2 - grid.k;
        suzerain::math::linspace(
                0.0, scenario.Ly, nbreak, buf.get()); // Uniform
        if (htdelta == 0.0) {
            INFO("Breakpoints distributed uniformly");
        } else {
            INFO("Breakpoints stretched with hyperbolic tangent delta = "
                 << htdelta);
            for (int i = 0; i < nbreak; ++i) {        // Stretch
                buf[i] = scenario.Ly
                       * suzerain_htstretch2(htdelta, scenario.Ly, buf[i]);
            }
        }

        // Generate the B-spline workspace based on order and breakpoints
        // Maximum non-trivial derivative operators included
        bspw = make_shared<suzerain::bspline>(
                grid.k, grid.k - 2, nbreak, buf.get());
        assert(static_cast<unsigned>(bspw->ndof()) == grid.Ny);
    }

    INFO("Creating new restart file " << restart_file);
    esio_file_create(esioh, restart_file.c_str(),
                     options.variables().count("clobber"));
    store(esioh, scenario);
    store(esioh, grid, scenario.Lx, scenario.Lz);
    store(esioh, bspw);
    esio_file_flush(esioh);

    INFO("Computing derived, dimensional reference parameters");
    real_t T_wall, rho_wall;
    {
        if (scenario.gamma != 1.4) {
            WARN("Using air viscosity values for non-air!");
        }

        tsolver p;
        p.Re     = scenario.Re;
        p.L      = scenario.Ly;
        p.gamma  = scenario.gamma;
        p.R      = R;
        p.Ma     = Ma;
        p.p_wall = p_wall;
        p.s = suzerain_svehla_air_mu_vs_T();
        assert(p.s);

        INFO("Solving " << tsolver_desc << " for T_wall");
        INFO("Using Re = " << p.Re
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
                FATAL("Likely no solution in interval ["
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
        INFO(std::setprecision(numeric_limits<real_t>::digits10)
             << "Wall temperature "
             << T_wall
             << " K gives residual "
             << std::scientific
             << GSL_FN_EVAL(&F,T_wall));

        rho_wall = p_wall / R / T_wall;
        INFO("Wall density is " << rho_wall << " kg/m^3");

        gsl_root_fsolver_free(solver);
        gsl_spline_free(p.s);
    }

    // Initialize B-spline workspace to find coeffs from collocation points
    bspluzw = make_shared<suzerain::bspline_luz>(*bspw);
    bspluzw->form_mass(*bspw);

    INFO("Computing nondimensional mean profiles for restart");
    {
        const int Ny = numeric_cast<int>(grid.Ny);
        esio_field_establish(esioh, 1, 0, 1, 1, 0, 1, Ny, 0, Ny);
        scoped_array<complex_t> buf(new complex_t[Ny]);

        // Nondimensional spanwise and wall-normal velocities are zero
        std::fill_n(buf.get(), Ny, complex_t(0,0));
        complex_field_write(
                esioh, "rhov", buf.get(), 0, 0, 0, field_descriptions[2]);
        complex_field_write(
                esioh, "rhow", buf.get(), 0, 0, 0, field_descriptions[3]);

        // Nondimensional density is the constant one
        std::fill_n(buf.get(), Ny, complex_t(1,0));
        complex_field_write(
                esioh, "rho", buf.get(), 0, 0, 0, field_descriptions[0]);

        // Set up to evaluate Y momentum and total energy profile coefficients
        scoped_array<double> rhs(new double[Ny]);
        mesolver params;
        params.Ma     = Ma;
        params.L      = scenario.Ly;
        params.gamma  = scenario.gamma;
        params.R      = R;
        params.T_wall = T_wall;

        suzerain_function F;
        F.params = &params;

        // Find Y momentum coefficients
        F.function = &f_msolver;
        bspw->find_interpolation_problem_rhs(&F, rhs.get());
        for (int i = 0; i < Ny; ++i) buf[i] = rhs[i];
        bspluzw->solve(1, buf.get(), 1, Ny);
        complex_field_write(
                esioh, "rhou", buf.get(), 0, 0, 0, field_descriptions[1]);

        // Find total energy coefficients
        F.function = &f_esolver;
        bspw->find_interpolation_problem_rhs(&F, rhs.get());
        for (int i = 0; i < Ny; ++i) buf[i] = rhs[i];
        bspluzw->solve(1, buf.get(), 1, Ny);
        complex_field_write(
                esioh, "rhoe", buf.get(), 0, 0, 0, field_descriptions[4]);

    }
    esio_file_flush(esioh);

    // Store new simulation zero time
    store_time(esioh, 0);
    esio_file_flush(esioh);

    INFO("Closing newly initialized restart file");
    esio_file_close(esioh);
}
