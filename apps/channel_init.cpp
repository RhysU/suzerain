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
#include <suzerain/error.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/program_options.hpp>
#include <suzerain/utility.hpp>

#include "logger.hpp"
#include "precision.hpp"
#include "channel.hpp"
#include "nsctpl_rholut.hpp"

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
        /* Ma        */ real_t(115)/real_t(100),
        /* Pr        */ real_t(7)/real_t(10),
        /* bulk_rho  */ 1,
        /* bulk_rhou */ 1,
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
static shared_ptr<nsctpl_rholut::manufactured_solution<real_t> > ms(
            new nsctpl_rholut::manufactured_solution<real_t>);

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

/** Options definitions for tweaking the manufactured solution */
class MSDefinition : public suzerain::problem::IDefinition {

public:

    MSDefinition(nsctpl_rholut::manufactured_solution<real_t> &msoln)
        : IDefinition("Manufactured solution parameters"
                      " (active only when --mms supplied)")
    {
        msoln.rho.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects density field",     _1, _2));
        msoln.u.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects X velocity field",  _1, _2));
        msoln.v.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects Y velocity field",  _1, _2));
        msoln.w.foreach_parameter(boost::bind(option_adder,
                    this->add_options(), "Affects Z velocity field",  _1, _2));
        msoln.T.foreach_parameter(boost::bind(option_adder,
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
    bool   clobber = false;
    real_t mms  = -1;
    {
        suzerain::ProgramOptions options(
                "Suzerain-based compressible channel initialization",
                "RESTART-FILE");

        namespace po = ::boost::program_options;
        using ::suzerain::validation::ensure_nonnegative;
        ::std::pointer_to_binary_function<real_t,const char*,void>
            ptr_fun_ensure_nonnegative(ensure_nonnegative<real_t>);

        nsctpl_rholut::isothermal_channel(*ms);
        MSDefinition msdef(*ms);

        options.add_definition(scenario);
        options.add_definition(grid);
        options.add_definition(msdef);
        options.add_options()
            ("clobber", "Overwrite an existing restart file?")
            ("mms",
             boost::program_options::value<real_t>(&mms)
                ->notifier(std::bind2nd(ptr_fun_ensure_nonnegative, "mms")),
             "If given, prepare a manufactured solution at the specified time.")
        ;
        std::vector<std::string> positional = options.process(argc, argv);

        if (positional.size() != 1) {
            FATAL0("Exactly one restart file name must be specified");
            return EXIT_FAILURE;
        }
        restart_file = positional[0];

        clobber = options.variables().count("clobber");
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
        ms->alpha = scenario.alpha;
        ms->beta  = scenario.beta;
        ms->gamma = scenario.gamma;
        ms->Ma    = scenario.Ma;
        ms->Re    = scenario.Re;
        ms->Pr    = scenario.Pr;
        ms->Lx    = scenario.Lx;
        ms->Ly    = scenario.Ly;
        ms->Lz    = scenario.Lz;
    } else {
        ms.reset();
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
    channel::store(esioh, scenario, ms);
    esio_file_flush(esioh);

    INFO0("Initializing B-spline workspaces");
    bopluz = make_shared<suzerain::bsplineop_luz>(*bop);
    bopluz->form_mass(*bop);

    INFO0("Initializing pencil_grid to obtain parallel decomposition details");
    dgrid = make_shared<suzerain::pencil_grid>(grid.dN, grid.P);
    assert((grid.dN == dgrid->global_physical_extent).all());

    INFO0("Initializing OperatorBase to access decomposition-ready utilities");
    suzerain::OperatorBase<real_t> obase(scenario, grid, *dgrid, *b, *bop);

    INFO0("Allocating and clearing storage for the distributed state fields");
    state_type swave(suzerain::to_yxz(
                channel::field::count, dgrid->local_wave_extent));
    suzerain::multi_array::fill(swave, 0);

    INFO0("Initializing data on collocation points values in physical space");

    // State viewed as a 2D Eigen::Map ordered (F, Y*Z*X).
    channel::physical_view<channel::field::count>::type sphys
        = channel::physical_view<channel::field::count>::create(*dgrid, swave);

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

            const real_t z = obase.z(j);

            for (int i = dgrid->local_physical_start.x();
                i < dgrid->local_physical_end.x();
                ++i, /* NB */ ++offset) {

                const real_t x = obase.z(i);

                // Initialize primitive state for...
                real_t rho, u, v, w, T;
                if (mms < 0) {
                    // ...a very simple parabolic velocity profile.
                    rho = 1;
                    u   = 6*y*(scenario.Ly - y)/(scenario.Ly*scenario.Ly);
                    v   = 0;
                    w   = 0;
                    T   = 1;
                } else {
                    // ...the manufactured solution at t = mms.
                    rho = ms->rho(x, y, z, mms);
                    u   = ms->u  (x, y, z, mms);
                    v   = ms->v  (x, y, z, mms);
                    w   = ms->w  (x, y, z, mms);
                    T   = ms->T  (x, y, z, mms);
                }

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

    INFO0("Converting state to wave space coefficients");
    {
        // Build FFT normalization constant into Y direction's mass matrix
        suzerain::bsplineop_luz massluz(*bop);
        const complex_t scale_factor = grid.dN.x() * grid.dN.z();
        massluz.form(1, &scale_factor, *bop);

        for (std::size_t i = 0; i < channel::field::count; ++i) {
            dgrid->transform_physical_to_wave(&sphys(i, 0));     // X, Z
            obase.bop_solve(massluz, swave, i);                  // Y
        }
    }

    INFO0("Writing state fields to restart file");
    channel::store(esioh, swave, grid, *dgrid);
    esio_file_flush(esioh);

    if (mms < 0) {
        INFO0("Storing new simulation time of zero");
        channel::store_time(esioh, 0);
    } else {
        INFO0("Storing simulation time to match manufactured solution");
        channel::store_time(esioh, mms);
    }
    esio_file_flush(esioh);

    INFO0("Closing newly initialized restart file");
    esio_file_close(esioh);
}
