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
static const pb::ScenarioDefinition<real_t> scenario(
        /* default_Re    */ 100,
        /* default_Pr    */ real_t(7)/real_t(10),
        /* default_gamma */ real_t(14)/real_t(10),
        /* default_beta  */ real_t(2)/real_t(3),
        /* default_Lx    */ 4*pi<real_t>(),
        /* default_Ly    */ 2,
        /* default_Lz    */ 4*pi<real_t>()/3);
static const pb::GridDefinition<real_t> grid(
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
    const double const_contrib = 1 / (p->gamma * (p->gamma - 1));
    const double U = f_msolver(y, p);
    const double vel_contrib   = 0.5 * U * U / (p->gamma * p->R * p->T_wall);
    return const_contrib + vel_contrib;
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
    real_t p_wall = GSL_CONST_MKSA_STD_ATMOSPHERE / 100000;
    std::string restart_file;
    real_t htdelta;
    sz::ProgramOptions options(
            "Suzerain-based compressible channel initialization");
    {
        // Cast away const so options processing can modify settings
        options.add_definition(
                const_cast<pb::ScenarioDefinition<real_t>& >(scenario));
        options.add_definition(
                const_cast<pb::GridDefinition<real_t>& >(grid));

        using ::suzerain::validation::ensure_positive;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_positive(ensure_positive<real_t>);
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
                ->notifier(std::bind2nd(ptr_fun_ensure_positive,"htdelta"))
                ->default_value(7),
             "Hyperbolic tangent stretching parameter for breakpoints")
        ;
        options.process(argc, argv);
    }
    const real_t R = GSL_CONST_MKSA_MOLAR_GAS * GSL_CONST_NUM_KILO / M;

    if (grid.k < 4 /* cubics */) {
        LOG4CXX_FATAL(log,
            "k >= 4 required to compute two non-trivial spatial derivatives");
        return EXIT_FAILURE;
    }

    if (grid.Nx != 1) {
        LOG4CXX_FATAL(log, argv[0] << " can only handle Nx == 1");
        return EXIT_FAILURE;
    }

    if (grid.Nz != 1) {
        LOG4CXX_FATAL(log, argv[0] << " can only handle Nz == 1");
        return EXIT_FAILURE;
    }

    LOG4CXX_INFO(log, "Creating new restart file " << restart_file);
    esio_file_create(esioh, restart_file.c_str(),
                     options.variables().count("clobber"));
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing basic scenario parameters");
    {
        esio_line_establish(esioh, 1, 0, 1); // Store as lines to allow comments

        esio_line_write(esioh, "Re", &scenario.Re, 0,
                scenario.options().find("Re",false).description().c_str());

        esio_line_write(esioh, "Pr", &scenario.Pr, 0,
                scenario.options().find("Pr",false).description().c_str());

        esio_line_write(esioh, "Ma", &Ma, 0,
                options.options().find("Ma",false).description().c_str());

        esio_line_write(esioh, "gamma", &scenario.gamma, 0,
                scenario.options().find("gamma",false).description().c_str());

        esio_line_write(esioh, "M", &M, 0,
                options.options().find("M",false).description().c_str());

        esio_line_write(esioh, "R", &R, 0, "Specific gas constant in J/kg/K");

        esio_line_write(esioh, "beta", &scenario.beta, 0,
                scenario.options().find("beta",false).description().c_str());

        esio_line_write(esioh, "Lx", &scenario.Lx, 0,
                scenario.options().find("Lx",false).description().c_str());

        esio_line_write(esioh, "Ly", &scenario.Ly, 0,
                scenario.options().find("Ly",false).description().c_str());

        esio_line_write(esioh, "Lz", &scenario.Lz, 0,
                scenario.options().find("Lz",false).description().c_str());
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing basic grid parameters");
    {
        esio_line_establish(esioh, 1, 0, 1); // Store as lines to allow comments

        const int Nx = numeric_cast<int>(grid.Nx);
        esio_line_write(esioh, "Nx", &Nx, 0,
                grid.options().find("Nx",false).description().c_str());

        esio_line_write(esioh, "DAFx", &grid.DAFx, 0,
                grid.options().find("DAFx",false).description().c_str());

        const int Ny = numeric_cast<int>(grid.Ny);
        esio_line_write(esioh, "Ny", &Ny, 0,
                grid.options().find("Ny",false).description().c_str());

        const int k = numeric_cast<int>(grid.k);
        esio_line_write(esioh, "k", &k, 0,
                grid.options().find("k",false).description().c_str());

        const int Nz = numeric_cast<int>(grid.Nz);
        esio_line_write(esioh, "Nz", &Nz, 0,
                grid.options().find("Nz",false).description().c_str());

        esio_line_write(esioh, "DAFz", &grid.DAFz, 0,
                grid.options().find("DAFz",false).description().c_str());
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing B-spline basis of uniform order "
                      << (grid.k - 1) << " on [0, Ly] with "
                      << grid.Ny << " DOF");
    {
        scoped_array<real_t> buf(new real_t[grid.Ny]);

        // Compute and store breakpoint locations
        const int nbreak = grid.Ny + 2 - grid.k;
        sz::math::linspace(0.0, 1.0, nbreak, buf.get()); // Uniform [0, 1]
        for (int i = 0; i < nbreak; ++i) {               // Stretch 'em out
            buf[i] = scenario.Ly * suzerain_htstretch2(htdelta, 1.0, buf[i]);
        }
        esio_line_establish(esioh, nbreak, 0, nbreak);
        esio_line_write(esioh, "breakpoints", buf.get(), 0,
                "Breakpoint locations used to build B-spline basis");

        // Generate the B-spline workspace based on order and breakpoints
        // Maximum non-trivial derivative operators included
        bspw = make_shared<sz::bspline>(
                grid.k, grid.k - 2, nbreak, buf.get());
        assert(static_cast<unsigned>(bspw->ndof()) == grid.Ny);

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
            snprintf(name, sizeof(name)/sizeof(name[0]), "Dy%d", k);
            snprintf(comment, sizeof(comment)/sizeof(comment[0]),
                    "Wall-normal derivative Dy%d(i,j) = D%d[j,ku+i-j] for"
                    " 0 <= j < n, max(0,j-ku-1) <= i < min(m,j+kl)", k, k);
            const int lda = bspw->ku(k) + 1 + bspw->kl(k);
            esio_plane_establish(
                    esioh, bspw->ndof(), 0, bspw->ndof(), lda, 0, lda);
            esio_plane_write(esioh, name, bspw->D(k), 0, 0, comment);
            esio_attribute_write(esioh, name, "kl", bspw->kl(k));
            esio_attribute_write(esioh, name, "ku", bspw->ku(k));
            esio_attribute_write(esioh, name, "m",  bspw->ndof());
            esio_attribute_write(esioh, name, "n",  bspw->ndof());
        }
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Storing wavenumber vectors for Fourier bases");
    {
        const int Nx = numeric_cast<int>(grid.Nx);
        const int Nz = numeric_cast<int>(grid.Nz);
        const int N  = std::max(Nx, Nz);
        scoped_array<complex_t> buf(new complex_t[N]);

        // Obtain wavenumbers via computing 1*(i*kx)/i
        std::fill_n(buf.get(), N, complex_t(1,0));
        sz::diffwave::apply(1, 0, complex_t(0,-1), buf.get(),
                scenario.Lx, scenario.Lz, 1, Nx, Nx, 0, Nx, 1, 1, 0, 1);
        esio_line_establish(esioh, Nx, 0, Nx);
        esio_line_write(esioh, "kx", reinterpret_cast<real_t *>(buf.get()),
                2, "Wavenumbers in streamwise X direction"); // Re(buf)

        // Obtain wavenumbers via computing 1*(i*kz)/i
        std::fill_n(buf.get(), N, complex_t(1,0));
        sz::diffwave::apply(1, 0, complex_t(0,-1), buf.get(),
                scenario.Lx, scenario.Lz, 1, 1, 1, 0, 1, Nz, Nz, 0, Nz);
        esio_line_establish(esioh, Nz, 0, Nz);
        esio_line_write(esioh, "kz", reinterpret_cast<real_t *>(buf.get()),
                2, "Wavenumbers in spanwise Z direction"); // Re(buf)
    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Computing derived, dimensional reference parameters");
    real_t T_wall, rho_wall;
    {
        if (scenario.gamma != 1.4) {
            LOG4CXX_WARN(log, "Using air viscosity values for non-air!");
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

        LOG4CXX_INFO(log, "Solving " << tsolver_desc << " for T_wall");
        LOG4CXX_INFO(log, "Using Re = " << p.Re
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
                LOG4CXX_FATAL(log, "Likely no solution in interval ["
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
        LOG4CXX_INFO(log, std::setprecision(numeric_limits<real_t>::digits10)
                        << "Wall temperature "
                        << T_wall
                        << " K gives residual "
                        << std::scientific
                        << GSL_FN_EVAL(&F,T_wall));

        rho_wall = p_wall / R / T_wall;
        LOG4CXX_INFO(log, "Wall density is " << rho_wall << " kg/m^3");

        gsl_root_fsolver_free(solver);
        gsl_spline_free(p.s);
    }

    // Initialize B-spline workspace to find coeffs from collocation points
    bspluzw = make_shared<sz::bspline_luz>(*bspw);
    bspluzw->form_mass(*bspw);

    LOG4CXX_INFO(log, "Computing nondimensional mean profiles for restart");
    {
        const int Ny = numeric_cast<int>(grid.Ny);
        esio_field_establish(esioh, 1, 0, 1, 1, 0, 1, Ny, 0, Ny);
        scoped_array<complex_t> buf(new complex_t[Ny]);

        // Nondimensional spanwise and wall-normal velocities are zero
        std::fill_n(buf.get(), Ny, complex_t(0,0));
        esio_field_writev(esioh, "rhov",
                reinterpret_cast<real_t *>(buf.get()), 0, 0, 0, 2,
                "Nondimensional Y momentum coefficients stored row-major ZXY");
        esio_field_writev(esioh, "rhow",
                reinterpret_cast<real_t *>(buf.get()), 0, 0, 0, 2,
                "Nondimensional Z momentum coefficients stored row-major ZXY");

        // Nondimensional density is the constant one
        std::fill_n(buf.get(), Ny, complex_t(1,0));
        esio_field_writev(esioh, "rho",
                reinterpret_cast<real_t *>(buf.get()), 0, 0, 0, 2,
                "Nondimensional density coefficient stored row-major ZXY");

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
        esio_field_writev(esioh, "rhou",
                reinterpret_cast<real_t *>(buf.get()), 0, 0, 0, 2,
                "Nondimensional X momentum coefficients stored row-major ZXY");

        // Find total energy coefficients
        F.function = &f_esolver;
        bspw->find_interpolation_problem_rhs(&F, rhs.get());
        for (int i = 0; i < Ny; ++i) buf[i] = rhs[i];
        bspluzw->solve(1, buf.get(), 1, Ny);
        esio_field_writev(esioh, "rhoe",
                reinterpret_cast<real_t *>(buf.get()), 0, 0, 0, 2,
                "Nondimensional total energy coefficients stored row-major ZXY");

    }
    esio_file_flush(esioh);

    LOG4CXX_INFO(log, "Closing newly initialized restart file");
    esio_file_close(esioh);
}
