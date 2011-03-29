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
// channel_mean.cpp: Dump mean statistics for one or more restart files
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <Eigen/Core>
#include <esio/esio.h>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/program_options.hpp>

#include "logger.hpp"
#include "precision.hpp"
#include "channel.hpp"

#pragma warning(disable:383 1572)

// X-macros employed per http://drdobbs.com/blogs/cpp/228700289

#define FORALL_COORDS(apply)     \
    apply(t)                     \
    apply(y)                     \
    apply(yplus)

#define FORALL_STATE_PLUS(apply) \
    apply(uplus)                 \
    apply(vplus)                 \
    apply(wplus)

#define FORALL_STATE_CONS(apply) \
    apply(rho)                   \
    apply(rhou)                  \
    apply(rhov)                  \
    apply(rhow)                  \
    apply(rhoe)

#define FORALL_STATE_PRIM(apply) \
    apply(u)                     \
    apply(v)                     \
    apply(w)                     \
    apply(e)                     \
    apply(p)                     \
    apply(T)                     \
    apply(mu)                    \
    apply(nu)

namespace column {
    enum {
#define X(a) a,
    FORALL_COORDS(X)
    FORALL_STATE_PLUS(X)
    FORALL_STATE_CONS(X)
    FORALL_STATE_PRIM(X)
#undef X
#define X(a) a##_y,
    FORALL_STATE_CONS(X)
    FORALL_STATE_PRIM(X)
#undef X
#define X(a) a##_yy,
    FORALL_STATE_CONS(X)
    FORALL_STATE_PRIM(X)
#undef X
        COUNT
    };

#define STRINGIFY(s) #s
    static const char * name[COUNT] = {
#define X(a) STRINGIFY(a),
    FORALL_COORDS(X)
    FORALL_STATE_PLUS(X)
    FORALL_STATE_CONS(X)
    FORALL_STATE_PRIM(X)
#undef X
#define X(a) STRINGIFY(a##_y),
    FORALL_STATE_CONS(X)
    FORALL_STATE_PRIM(X)
#undef X
#define X(a) STRINGIFY(a##_yy),
    FORALL_STATE_CONS(X)
    FORALL_STATE_PRIM(X)
#undef X
    };
#undef STRINGIFY
}

// Introduce shorthand for common names
using boost::make_shared;
using boost::math::constants::pi;
using boost::numeric_cast;
using boost::shared_ptr;
using std::numeric_limits;

// Used to format output data
static const Eigen::IOFormat iofmt(Eigen::FullPrecision, 0, ", ", "\n");

// Compute quantities based on real-valued mean state coefficients
static
Eigen::Array<real_t, Eigen::Dynamic, column::COUNT>
process(const Eigen::ArrayXXr &s_coeffs,
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::bspline &bspw,
        const real_t time);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit

    // Process incoming arguments
    std::vector<std::string> restart_files;
    {
        suzerain::ProgramOptions options(
                "Suzerain-based channel mean quantity computations",
                "[RESTART-FILE]...");
        restart_files = options.process(argc, argv);
    }

    // Ensure that we're running in a single processor environment
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) > 1) {
        FATAL(argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // Process each command line argument as a file name
    BOOST_FOREACH(const std::string& filename, restart_files) {

        // Create a file-specific ESIO handle using RAII
        shared_ptr<boost::remove_pointer<esio_handle>::type> h(
                esio_handle_initialize(MPI_COMM_WORLD), esio_handle_finalize);

        DEBUG("Loading file " << filename);
        esio_file_open(h.get(), filename.c_str(), 0 /* read-only */);

        // Load time, scenario, grid, and B-spline details from file
        real_t time;
        suzerain::problem::ScenarioDefinition<real_t> scenario;
        suzerain::problem::GridDefinition grid;
        shared_ptr<const suzerain::bspline> bspw;
        channel::load_time(h.get(), time);
        channel::load(h.get(), scenario);
        channel::load(h.get(), grid);
        channel::load(h.get(), bspw);
        assert(bspw->ndof() == grid.N.y());

        // Load zero-zero mode coefficients for all state variables
        Eigen::ArrayXXr s_coeffs(grid.N.y(), field_names.size());
        s_coeffs.setZero();
        {
            Eigen::VectorXc tmp(grid.N.y());
            esio_field_establish(h.get(), grid.N.z(),     0, 1,
                                        grid.N.x()/2+1, 0, 1,
                                        grid.N.y(),     0, grid.N.y());
            for (int i = 0; i < s_coeffs.cols(); ++i) {
                channel::complex_field_read(
                        h.get(), field_names[i], tmp.data());
                assert(tmp.imag().squaredNorm() == 0); // Purely real!
                s_coeffs.col(i) = tmp.real();
            }
        }

        // Close the data file
        esio_file_close(h.get());

        // Compute the quantities of interest
        Eigen::Array<real_t, Eigen::Dynamic, column::COUNT> s
            = process(s_coeffs, scenario, grid, *bspw, time);

        // Save nondimensional quantities to `basename filename .h5`.mean
        std::string outname;
        if (filename.rfind(".h5") == filename.length() - 3) {
            outname = filename.substr(0, filename.length() - 3) + ".mean";
        } else {
            outname = filename + ".mean";
        }
        DEBUG("Saving nondimensional quantities to " << outname);
        {
            std::ofstream outfile(outname.c_str());
            for (size_t i = 0; i < column::COUNT; ++i) {  // Headings
                outfile << std::setw(numeric_limits<real_t>::digits10 + 7)
                        << column::name[i];
                if (i < column::COUNT - 1) outfile << ", ";
            }
            outfile << std::endl;
            outfile << s.format(iofmt) << std::endl;;
            outfile.close();
        }

    }

    return EXIT_SUCCESS;
}

static
Eigen::Array<real_t, Eigen::Dynamic, column::COUNT>
process(const Eigen::ArrayXXr &s_coeffs,
        const suzerain::problem::ScenarioDefinition<real_t> &scenario,
        const suzerain::problem::GridDefinition &grid,
        const suzerain::bspline &bspw,
        const real_t time)
{
    SUZERAIN_UNUSED(grid);

    // Declare storage for all of the quantities of interest
    Eigen::Array<real_t, Eigen::Dynamic, column::COUNT>
        s(s_coeffs.rows(), (int) column::COUNT);
    s.setZero();

    // Populate point-like information of (t,y \in (0, Ly))
    // We'll compute y^{+} information later
    s.col(column::t).setConstant(time);
    bspw.collocation_points(&s.col(column::y)[0], 1);

    // Compute derivatives of conserved state at collocation points.
    bspw.accumulate_operator(0, field_names.static_size,
            1.0, &s_coeffs(0,0), 1, s_coeffs.stride(),
            0.0, &s.col(column::rho)[0], 1, s.stride());
    bspw.accumulate_operator(1, field_names.static_size,
            1.0, &s_coeffs(0,0), 1, s_coeffs.stride(),
            0.0, &s.col(column::rho_y)[0], 1, s.stride());
    bspw.accumulate_operator(2, field_names.static_size,
            1.0, &s_coeffs(0,0), 1, s_coeffs.stride(),
            0.0, &s.col(column::rho_yy)[0], 1, s.stride());

    // Compute specific and primitive state at collocation points
    {
        using namespace column;
        s.col(u)  = s.col(rhou) / s.col(rho);
        s.col(v)  = s.col(rhov) / s.col(rho);
        s.col(w)  = s.col(rhow) / s.col(rho);
        s.col(e)  = s.col(rhoe) / s.col(rho);
        s.col(p)  = (scenario.gamma - 1) * (s.col(rhoe)
                            - s.col(rhou) * s.col(u) / 2
                            - s.col(rhov) * s.col(v) / 2
                            - s.col(rhow) * s.col(w) / 2);
        s.col(T)  = scenario.gamma * s.col(p) / s.col(rho);
        s.col(mu) = s.col(T).pow(scenario.beta);
        s.col(nu) = s.col(mu) / s.col(rho);
    }

    // Compute derivatives of specific and primitive state. Better would be
    // using conserved state derivatives directly, but that's quite a PITA.

    // Copy primitive state collocation values to storage for 1st derivs
#define X(a) s.col(column::a##_y) = s.col(column::a);
    FORALL_STATE_PRIM(X)
#undef X

    // Form mass matrix and obtain coefficients for primitive state
    suzerain::bspline_lu mass(bspw);
    mass.form_mass(bspw);
#define X(a) mass.solve(1, &s.col(column::a##_y)[0], 1, s.stride());
    FORALL_STATE_PRIM(X)
#undef X

    // Copy primitive state coefficients to storage for 2nd derivatives
#define X(a) s.col(column::a##_yy) = s.col(column::a##_y);
    FORALL_STATE_PRIM(X)
#undef X

    // Apply 1st derivative operator to coefficients
#define X(a) bspw.apply_operator(1, 1, 1.0, &s.col(column::a##_y)[0], \
                                 1, s.stride());
    FORALL_STATE_PRIM(X)
#undef X

    // Apply 2nd derivative operator to coefficients
#define X(a) bspw.apply_operator(2, 1, 1.0, &s.col(column::a##_yy)[0], \
                                 1, s.stride());
    FORALL_STATE_PRIM(X)
#undef X

    // Compute wall shear stress, friction velocity, and viscous length scale.
    // These are "almost correct" as they are off by reference factors.
    // We correct for those factors later during {y,u,v,w}^{+} computation.
    const real_t rho_w    = s.col(column::rho)[0];
    const real_t nu_w     = s.col(column::nu)[0];
    const real_t tau_w    = rho_w * nu_w * s.col(column::u_y)[0];
    const real_t u_tau    = std::sqrt(tau_w / rho_w);
    const real_t delta_nu = nu_w / u_tau;

    // Compute {y,u,v,w}^{+}
    {
        using namespace column;
        s.col(yplus) = s.col(y) / delta_nu * std::sqrt(scenario.Re);
        s.col(uplus) = s.col(u) / u_tau    * std::sqrt(scenario.Re);
        s.col(vplus) = s.col(v) / u_tau    * std::sqrt(scenario.Re);
        s.col(wplus) = s.col(w) / u_tau    * std::sqrt(scenario.Re);
    }

    return s;
}
