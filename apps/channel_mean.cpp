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
#include <suzerain/bspline.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/math.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/orthonormal.hpp>
#include <suzerain/scenario_definition.hpp>

#include "logger.hpp"
#include "channel_common.hpp"

#pragma warning(disable:383 1572)

// Introduce shorthand for common names
using boost::make_shared;
using boost::math::constants::pi;
using boost::numeric_cast;
using boost::shared_ptr;
using std::numeric_limits;

static const Eigen::IOFormat iofmt(Eigen::FullPrecision, 0, ", ", "\n");

static bool process(const char * filename);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);                         // Initialize MPI
    atexit((void (*) ()) MPI_Finalize);             // Finalize MPI at exit

    // Ensure that we're running in a single processor environment
    if (suzerain::mpi::comm_size(MPI_COMM_WORLD) > 1) {
        FATAL(argv[0] << " only intended to run on single rank");
        return EXIT_FAILURE;
    }

    // We do not use suzerain::ProgramOptions as it eats the command line
    // Process each command line argument as a file name
    int retval = EXIT_SUCCESS;
    for (int i = 1; i < argc; ++i) {
        if (!process(argv[i])) retval = EXIT_FAILURE;
    }
    return retval;
}

static bool process(const char * filename)
{
    // Create a file-specific ESIO handle using RAII
    shared_ptr<boost::remove_pointer<esio_handle>::type> h(
            esio_handle_initialize(MPI_COMM_WORLD), esio_handle_finalize);

    INFO("Loading file " << filename);
    esio_file_open(h.get(), filename, 0 /* read-only */);

    // Load time, scenario, grid, and B-spline details from file
    real_t t;
    suzerain::problem::ScenarioDefinition<real_t> scenario;
    suzerain::problem::GridDefinition grid;
    shared_ptr<suzerain::bspline> bspw;
    load_time(h.get(), t);
    load(h.get(), scenario);
    load(h.get(), grid);
    load(h.get(), bspw);
    assert(static_cast<unsigned>(bspw->ndof()) == grid.N.y());

    // Load zero-zero mode coefficients for all state variables
    Eigen::ArrayXr s_coeffs(grid.N.y(), field_names.size());
    s_coeffs.setZero();
    {
        Eigen::VectorXc tmp(grid.N.y());
        esio_field_establish(h.get(), grid.N.z(), 0, 1,
                                      grid.N.x(), 0, 1,
                                      grid.N.y(), 0, grid.N.y());
        for (int i = 0; i < s_coeffs.cols(); ++i) {
            complex_field_read(h.get(), field_names[i], tmp.data());
            assert(tmp.imag().squaredNorm() == 0); // Purely real!
            s_coeffs.col(i) = tmp.real();
        }
    }

    // Close the data file
    esio_file_close(h.get());

    // Prepare column names and indices for later use
    std::vector<std::string> column_names;
#define QUANTITY(SYM)                            \
        const int n_##SYM = column_names.size(); \
        column_names.push_back(#SYM)
    QUANTITY(t);    // Point-like information
    QUANTITY(y);
    QUANTITY(rho);  // Conserved state
    QUANTITY(rhou);
    QUANTITY(rhov);
    QUANTITY(rhow);
    QUANTITY(rhoe);
    QUANTITY(u);    // Specific state
    QUANTITY(v);
    QUANTITY(w);
    QUANTITY(e);
    QUANTITY(p);    // Primitive state
    QUANTITY(T);
    QUANTITY(mu);
    QUANTITY(nu);
#undef QUANTITY
#define QUANTITY_DERIVATIVE_1(SYM)                            \
        const int n_##SYM##_y = column_names.size();          \
        column_names.push_back(column_names[n_##SYM] + "_y")
    QUANTITY_DERIVATIVE_1(rho);  // Conserved state
    QUANTITY_DERIVATIVE_1(rhou);
    QUANTITY_DERIVATIVE_1(rhov);
    QUANTITY_DERIVATIVE_1(rhow);
    QUANTITY_DERIVATIVE_1(rhoe);
    QUANTITY_DERIVATIVE_1(u);    // Specific state
    QUANTITY_DERIVATIVE_1(v);
    QUANTITY_DERIVATIVE_1(w);
    QUANTITY_DERIVATIVE_1(e);
    QUANTITY_DERIVATIVE_1(p);    // Primitive state
    QUANTITY_DERIVATIVE_1(T);
    QUANTITY_DERIVATIVE_1(mu);
    QUANTITY_DERIVATIVE_1(nu);
#undef QUANTITY_DERIVATIVE_1
#define QUANTITY_DERIVATIVE_2(SYM)                            \
        const int n_##SYM##_yy = column_names.size();         \
        column_names.push_back(column_names[n_##SYM] + "_yy")
    QUANTITY_DERIVATIVE_2(rho);  // Conserved state
    QUANTITY_DERIVATIVE_2(rhou);
    QUANTITY_DERIVATIVE_2(rhov);
    QUANTITY_DERIVATIVE_2(rhow);
    QUANTITY_DERIVATIVE_2(rhoe);
    QUANTITY_DERIVATIVE_2(u);    // Specific state
    QUANTITY_DERIVATIVE_2(v);
    QUANTITY_DERIVATIVE_2(w);
    QUANTITY_DERIVATIVE_2(e);
    QUANTITY_DERIVATIVE_2(p);    // Primitive state
    QUANTITY_DERIVATIVE_2(T);
    QUANTITY_DERIVATIVE_2(mu);
    QUANTITY_DERIVATIVE_2(nu);
#undef QUANTITY_DERIVATIVE_2

    INFO("Computing " << column_names.size() << " nondimensional quantities");

    // Declare storage for all of the quantities of interest
    Eigen::ArrayXr s(s_coeffs.rows(), column_names.size());
    s.setZero();

    // Populate point-like information of (t,y \in (-L/2, L/2))
    s.col(n_t).setConstant(t);
    bspw->collocation_points(&s.col(n_y)[0], 1);
    s.col(n_y).array() -= scenario.Ly / 2;

    // Compute 0th, 1st, and 2nd derivatives of conserved state
    // at collocation points
    assert(n_rhou == n_rho + 1);
    assert(n_rhov == n_rho + 2);
    assert(n_rhow == n_rho + 3);
    assert(n_rhoe == n_rho + 4);
    bspw->accumulate_operator(0, field_names.static_size,
            1.0, &s_coeffs(0,0),     1, s_coeffs.stride(),
            0.0, &s.col(n_rho)[0],   1, s.stride());
    bspw->accumulate_operator(1, field_names.static_size,
            1.0, &s_coeffs(0,0),     1, s_coeffs.stride(),
            0.0, &s.col(n_rho_y)[0], 1, s.stride());
    (void) n_rhou_y;  // Unused
    (void) n_rhov_y;
    (void) n_rhow_y;
    (void) n_rhoe_y;
    bspw->accumulate_operator(2, field_names.static_size,
            1.0, &s_coeffs(0,0),      1, s_coeffs.stride(),
            0.0, &s.col(n_rho_yy)[0], 1, s.stride());
    (void) n_rhou_yy; // Unused
    (void) n_rhov_yy;
    (void) n_rhow_yy;
    (void) n_rhoe_yy;

    // Compute specific and primitive state at collocation points
    s.col(n_u)  = s.col(n_rhou) / s.col(n_rho);
    s.col(n_v)  = s.col(n_rhov) / s.col(n_rho);
    s.col(n_w)  = s.col(n_rhow) / s.col(n_rho);
    s.col(n_e)  = s.col(n_rhoe) / s.col(n_rho);
    s.col(n_p)  = (scenario.gamma - 1) * (s.col(n_rhoe)
                        - s.col(n_rhou) * s.col(n_u) / 2
                        - s.col(n_rhov) * s.col(n_v) / 2
                        - s.col(n_rhow) * s.col(n_w) / 2);
    s.col(n_T)  = scenario.gamma * s.col(n_p) / s.col(n_rho);
    s.col(n_mu) = s.col(n_T).pow(scenario.beta);
    s.col(n_nu) = s.col(n_mu) / s.col(n_rho);

    // Compute derivatives of specific and primitive state. Better would be
    // using conserved state derivatives directly, but that's quite a PITA.

    s.col(n_u_y)  = s.col(n_u);   // Copy values to storage for 1st derivs
    s.col(n_v_y)  = s.col(n_v);
    s.col(n_w_y)  = s.col(n_w);
    s.col(n_e_y)  = s.col(n_e);
    s.col(n_p_y)  = s.col(n_p);
    s.col(n_T_y)  = s.col(n_T);
    s.col(n_mu_y) = s.col(n_mu);
    s.col(n_nu_y) = s.col(n_nu);
    assert(n_v_y  == n_u_y + 1);
    assert(n_w_y  == n_u_y + 2);
    assert(n_e_y  == n_u_y + 3);
    assert(n_p_y  == n_u_y + 4);
    assert(n_T_y  == n_u_y + 5);
    assert(n_mu_y == n_u_y + 6);
    assert(n_nu_y == n_u_y + 7);

    suzerain::bspline_lu mass(*bspw); // Form mass matrix and obtain coeffs
    mass.form_mass(*bspw);
    mass.solve(n_nu_y - n_u_y + 1, &s.col(n_u_y)[0], 1, s.stride());

    s.col(n_u_yy)  = s.col(n_u_y);   // Copy coeffs to storage for 2nd derivs
    s.col(n_v_yy)  = s.col(n_v_y);
    s.col(n_w_yy)  = s.col(n_w_y);
    s.col(n_e_yy)  = s.col(n_e_y);
    s.col(n_p_yy)  = s.col(n_p_y);
    s.col(n_T_yy)  = s.col(n_T_y);
    s.col(n_mu_yy) = s.col(n_mu_y);
    s.col(n_nu_yy) = s.col(n_nu_y);
    assert(n_v_yy  == n_u_yy + 1);
    assert(n_w_yy  == n_u_yy + 2);
    assert(n_e_yy  == n_u_yy + 3);
    assert(n_p_yy  == n_u_yy + 4);
    assert(n_T_yy  == n_u_yy + 5);
    assert(n_mu_yy == n_u_yy + 6);
    assert(n_nu_yy == n_u_yy + 7);

    // Apply 1st and 2nd derivative operator to coefficients
    bspw->apply_operator(1, n_nu_y - n_u_y + 1,
            1.0, &s.col(n_u_y)[0], 1, s.stride());
    bspw->apply_operator(2, n_nu_yy - n_u_yy + 1,
            1.0, &s.col(n_u_yy)[0], 1, s.stride());

    // Save nondimensional quantities to filename.star
    INFO("Saving nondimensional quantities to " << filename << ".star");
    {
        std::ofstream starfile((std::string(filename) + ".star").c_str());
        for (size_t i = 0; i < column_names.size(); ++i) {  // Headings
            starfile << std::setw(iofmt.precision + 6) << column_names[i];
            if (i < column_names.size() - 1) starfile << ", ";
        }
        starfile << std::endl;
        starfile << s.format(iofmt) << std::endl;;
        starfile.close();
    }

    INFO("Computing quantities in plus units");

    // Re-adjust collocation point offsets so lowest point is at y = 0
    s.col(n_y) -= s.col(n_y)[0];

    // Compute wall shear stress, friction velocity, and viscous length scale.
    // These are "almost correct" as they are off by reference factors.
    // We correct for those factors later.
    const real_t rho_w    = s.col(n_rho)[0];
    const real_t nu_w     = s.col(n_nu)[0];
    const real_t tau_w    = rho_w * nu_w * s.col(n_u_y)[0];
    const real_t u_tau    = std::sqrt(tau_w / rho_w);
    const real_t delta_nu = nu_w / u_tau;

    // We only use data from the lower half of the channel
    const int nplus = (grid.N.y() + 1) / 2;

    // Compute the quantities in plus units
    Eigen::ArrayXr r(nplus, 1 /* t */ + 1 /* y */ + 1 /* y+ */ +  3);
    r.setZero();
    for (int i = 0; i < nplus; ++i) {
        r(i,0) = t;
        r(i,1) = s.col(n_y)[i];
        r(i,2) = s.col(n_y)[i] / delta_nu * std::sqrt(scenario.Re);
        r(i,3) = s.col(n_u)[i] / u_tau    * std::sqrt(scenario.Re);
        r(i,4) = s.col(n_v)[i] / u_tau    * std::sqrt(scenario.Re);
        r(i,5) = s.col(n_w)[i] / u_tau    * std::sqrt(scenario.Re);
        // TODO: \partial{} U^{+} / \partial{} y^{+} and friends
    }

    INFO("Saving plus unit quantities to " << filename << ".plus");
    {
        std::ofstream plusfile((std::string(filename) + ".plus").c_str());
        plusfile << std::setw(iofmt.precision + 7) << "t, ";
        plusfile << std::setw(iofmt.precision + 7) << "y, ";
        plusfile << std::setw(iofmt.precision + 7) << "y+, ";
        plusfile << std::setw(iofmt.precision + 7) << "u+, ";
        plusfile << std::setw(iofmt.precision + 7) << "v+, ";
        plusfile << std::setw(iofmt.precision + 6) << "w+" << std::endl;
        plusfile << r.format(iofmt) << std::endl;
        plusfile.close();
    }

    return true;
}
