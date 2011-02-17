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
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <suzerain/bspline.hpp>
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

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
namespace sz = ::suzerain;
namespace pb = ::suzerain::problem;
namespace po = ::boost::program_options;
using boost::make_shared;
using boost::math::constants::pi;
using boost::numeric_cast;
using boost::scoped_array;
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

// Global B-spline -details initialized in main()
static shared_ptr<sz::bspline>     bspw;
static shared_ptr<sz::bspline_luz> bspluzw;

/** Global handle for ESIO operations across MPI_COMM_WORLD. */
static esio_handle esioh = NULL;

/** <tt>atexit</tt> callback to ensure we finalize esioh. */
static void atexit_esio(void) {
    if (esioh) esio_handle_finalize(esioh);
}

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
    return rhs - lhs;
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
    real_t Ma     = 1.15;
    real_t M      = SUZERAIN_SVEHLA_AIR_M;
    real_t p_wall = GSL_CONST_MKSA_STD_ATMOSPHERE;
    std::string restart_file;
    real_t htdelta;
    sz::ProgramOptions options(
            "Suzerain-based compressible channel initialization");
    {
        options.add_definition(def_scenario);
        options.add_definition(def_grid);
        using ::suzerain::validation::ensure_positive;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_positive(ensure_positive<real_t>);
        options.add_options()
            ("create", po::value<std::string>(&restart_file)
                ->default_value("restart0.h5"),
             "Name of new restart file to create")
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
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"htdelta"))
                ->default_value(7),
             "Hyperbolic tangent stretching parameter for breakpoints")
        ;
        options.process(argc, argv);
    }
    const real_t R = GSL_CONST_MKSA_MOLAR_GAS * GSL_CONST_NUM_KILO / M;



    if (def_grid.k() < 4 /* cubics */) {
        LOG4CXX_FATAL(log,
            "k >= 4 required to compute two non-trivial spatial derivatives");
        return EXIT_FAILURE;
    }

    LOG4CXX_INFO(log, "Creating new restart file " << restart_file);
    esio_file_create(esioh, restart_file.c_str(), false /* no clobber */);
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing basic scenario parameters");
    {
        esio_line_establish(esioh, 1, 0, 1); // Store as lines to allow comments

        const double Re = def_scenario.Re();
        esio_line_write(esioh, "Re", &Re, 0,
                def_scenario.options().find("Re",false).description().c_str());

        const double Pr = def_scenario.Pr();
        esio_line_write(esioh, "Pr", &Pr, 0,
                def_scenario.options().find("Pr",false).description().c_str());

        esio_line_write(esioh, "Ma", &Ma, 0,
                options.options().find("Ma",false).description().c_str());

        const double gamma = def_scenario.gamma();
        esio_line_write(esioh, "gamma", &gamma, 0,
                def_scenario.options().find("gamma",false).description().c_str());

        esio_line_write(esioh, "M", &M, 0,
                options.options().find("M",false).description().c_str());

        esio_line_write(esioh, "R", &R, 0, "Specific gas constant in J/kg/K");

        const double beta = def_scenario.beta();
        esio_line_write(esioh, "beta", &beta, 0,
                def_scenario.options().find("beta",false).description().c_str());

        const double Lx = def_scenario.Lx();
        esio_line_write(esioh, "Lx", &Lx, 0,
                def_scenario.options().find("Lx",false).description().c_str());

        const double Ly = def_scenario.Ly();
        esio_line_write(esioh, "Ly", &Ly, 0,
                def_scenario.options().find("Ly",false).description().c_str());

        const double Lz = def_scenario.Lz();
        esio_line_write(esioh, "Lz", &Lz, 0,
                def_scenario.options().find("Lz",false).description().c_str());
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing basic grid parameters");
    {
        esio_line_establish(esioh, 1, 0, 1); // Store as lines to allow comments

        const int Nx = numeric_cast<int>(def_grid.Nx());
        esio_line_write(esioh, "Nx", &Nx, 0,
                def_grid.options().find("Nx",false).description().c_str());

        const double DAFx = def_grid.DAFx();
        esio_line_write(esioh, "DAFx", &DAFx, 0,
                def_grid.options().find("DAFx",false).description().c_str());

        const int Ny = numeric_cast<int>(def_grid.Ny());
        esio_line_write(esioh, "Ny", &Ny, 0,
                def_grid.options().find("Ny",false).description().c_str());

        const int k = numeric_cast<int>(def_grid.k());
        esio_line_write(esioh, "k", &k, 0,
                def_grid.options().find("k",false).description().c_str());

        const int Nz = numeric_cast<int>(def_grid.Nz());
        esio_line_write(esioh, "Nz", &Nz, 0,
                def_grid.options().find("Nz",false).description().c_str());

        const double DAFz = def_grid.DAFz();
        esio_line_write(esioh, "DAFz", &DAFz, 0,
                def_grid.options().find("DAFz",false).description().c_str());
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing B-spline basis of uniform order "
                      << (def_grid.k() - 1) << " on [0, Ly] with "
                      << def_grid.Ny() << " DOF");
    {
        scoped_array<real_t> buf(new real_t[def_grid.Ny()]);

        // Compute and store breakpoint locations
        const int nbreak = def_grid.Ny() + 2 - def_grid.k();
        sz::math::linspace(0.0, 1.0, nbreak, buf.get()); // Uniform [0, 1]
        for (int i = 0; i < nbreak; ++i) {               // Stretch 'em out
            buf[i] = def_scenario.Ly()
                   * suzerain_htstretch2(htdelta, 1.0, buf[i]);
        }
        esio_line_establish(esioh, nbreak, 0, nbreak);
        esio_line_write(esioh, "breakpoints", buf.get(), 0,
                "Breakpoint locations used to build B-spline basis");

        // Generate the B-spline workspace based on order and breakpoints
        // Maximum non-trivial derivative operators included
        bspw = make_shared<sz::bspline>(
                def_grid.k(), def_grid.k() - 2, nbreak, buf.get());
        assert(static_cast<unsigned>(bspw->ndof()) == def_grid.Ny());

        // Store collocation points to restart file
        bspw->collocation_points(buf.get(), 1);
        esio_line_establish(esioh, bspw->ndof(), 0, bspw->ndof());
        esio_line_write(esioh, "colpoints", buf.get(), 0,
                "Collocation points used to build discrete operators");
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log,
            "Storing B-spline derivative operators in general band format");
    {
        char name[8];
        char comment[127];
        for (int k = 0; k <= bspw->nderivatives(); ++k) {
            snprintf(name, sizeof(name)/sizeof(name[0]), "D%d", k);
            snprintf(comment, sizeof(comment)/sizeof(comment[0]),
                    "Wall-normal derivative D%d(i,j) = D%d[i+ku+j*(kl+ku)]"
                    " for 0 <= i,j < N when (j-ku <= i && i <= j+kl)", k, k);
            const int aglobal = (bspw->ku(k) + 1 + bspw->kl(k)) * bspw->ndof();
            esio_line_establish(esioh, aglobal, 0, aglobal);
            esio_line_write(esioh, name, bspw->D(k), 0, comment);
            esio_attribute_write(esioh, name, "kl", bspw->kl(k));
            esio_attribute_write(esioh, name, "ku", bspw->ku(k));
            esio_attribute_write(esioh, name, "N",  bspw->ndof());
        }
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Computing derived, dimensional scenario parameters");
    real_t T_wall;
    {
        if (def_scenario.gamma() != 1.4) {
            LOG4CXX_WARN(log, "Using air viscosity vs temperature for non-air");
        }

        tsolver p;
        p.Re     = def_scenario.Re();
        p.L      = def_scenario.Ly();
        p.gamma  = def_scenario.gamma();
        p.R      = R;
        p.Ma     = Ma;
        p.p_wall = p_wall;
        p.s = suzerain_svehla_air_mu_vs_T();
        assert(p.s);

        gsl_function F;
        F.function = &f_tsolver;
        F.params   = &p;

        gsl_root_fsolver * solver
            = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        assert(solver);
        gsl_root_fsolver_set(solver, &F, 100 /* K */, 5000 /* K */);

        real_t low, high;
        int status, iter, maxiter = 100;
        do {
            status = gsl_root_fsolver_iterate(solver);
            T_wall = gsl_root_fsolver_root(solver);
            low    = gsl_root_fsolver_x_lower(solver);
            high   = gsl_root_fsolver_x_upper(solver);
            LOG4CXX_INFO(log, "Looking for desired wall temperature within ["
                    << low << "," << high << "");
            status = gsl_root_test_interval(low, high, 0, 1e-8);
            ++iter;
        } while (status == GSL_CONTINUE && iter < maxiter);

        gsl_root_fsolver_free(solver);
        gsl_spline_free(p.s);
    }

    // Initialize B-spline workspace to find coeffs from collocation points
    bspluzw = make_shared<sz::bspline_luz>(*bspw);
    bspluzw->form_mass(*bspw);

    LOG4CXX_INFO(log, "Closing newly initialized restart file");
    esio_file_close(esioh);
}
