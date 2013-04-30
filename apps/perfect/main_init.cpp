//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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
 * Application executing \ref suzerain::perfect::driver_init::run.
 */

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/validation.hpp>

#include "driver.hpp"

#pragma warning(disable:1419)

namespace suzerain { namespace perfect {

/** Application for initializing new restart files. */
struct driver_init : public driver
{
    driver_init(const std::string& revstr)
        : driver("Compressible, perfect gas simulation initialization",
                 "RESTART-FILE",
                 "",
                 revstr)
        , who("init")
    {}

    /** Implementation below in this file */
    int run(int argc, char **argv);

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

} /* namespace perfect */ } /* namespace suzerain */

// Provided by main_init_svnrev.{c,h} so revstr updates are merely relinking
extern "C" const char revstr[];

/** Instantiate and invoke the application */
int main(int argc, char **argv)
{
    suzerain::perfect::driver_init app(revstr);
    return app.run(argc, argv);
}

int
suzerain::perfect::driver_init::run(int argc, char **argv)
{
    using boost::math::constants::pi;
    using boost::math::tgamma;
    using std::pow;
    using std::sin;

    // Establish default grid and domain extents
    grid->L.x()   = 4 * pi<real_t>();
    grid->L.y()   = 2;
    grid->L.z()   = 4 * pi<real_t>() / 3;
    grid->Nx(1);
    grid->DAFx(1.5);
    grid->Ny(32);
    grid->k       = 8;
    grid->htdelta = 3;
    grid->Nz(1);
    grid->DAFz(1.5);

    // Establish default scenario parameters
    scenario->Re         = 100;
    scenario->Ma         = 1.5;
    scenario->Pr         = 0.7;
    scenario->bulk_rho   = 1;
    scenario->bulk_rho_u = 1;
    scenario->alpha      = 0;
    scenario->beta       = real_t(2) / 3;
    scenario->gamma      = 1.4;

    // Establish default isothermal boundary conditions
    isothermal->lower_T = isothermal->upper_T = 1;
    isothermal->lower_u = isothermal->upper_u = 0;
    isothermal->lower_v = isothermal->upper_v = 0;
    isothermal->lower_w = isothermal->upper_w = 0;

    // Establish default time step aggressiveness
    timedef = make_shared<support::time_definition>(/* per Venugopal */ 0.72);

    // Establish default MMS parameters and plug into program options
    msoln = make_shared<manufactured_solution>(
            manufactured_solution::default_caption
            + " (active only when --mms supplied)");
    options.add_definition(msoln->isothermal_channel());

    // Establish binary-specific options
    std::pointer_to_binary_function<real_t,const char*,void>
        ensure_real_tnonnegative(validation::ensure_nonnegative<real_t>);
    real_t mms    = -1;
    real_t npower =  1;
    options.add_options()
        ("clobber", "Overwrite an existing restart file?")
        ("npower",
            boost::program_options::value(&npower)->default_value(npower),
            "Power n in (0, 1] used to control the flatness of the"
            " \"parabolic\" streamwise velocity profile (y*(L-y))^n.")
        ("mms",
            boost::program_options::value(&mms)
            ->notifier(std::bind2nd(ensure_real_tnonnegative, "mms")),
            "If given, prepare a manufactured solution at the specified time.")
    ;

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros becomes usable)
    std::vector<std::string> positional = initialize(argc, argv);
    if (positional.size() != 1) {
        FATAL0("Exactly one restart file name must be specified");
        return EXIT_FAILURE;
    }
    const std::string restart_file = positional[0];
    const bool clobber = options.variables().count("clobber");
    if (npower < 0 || npower > 1) {
        FATAL0("npower in [0,1] required");
        return EXIT_FAILURE;
    }
    if (grid->k < 4) {
        FATAL0("k >= 4 required for two non-trivial wall-normal derivatives");
        return EXIT_FAILURE;
    }

    DEBUG0(who, "Establishing runtime parallel infrastructure and resources");
    establish_ieee_mode();
    load_grid_and_operators(NULL);
    log_discretization_quality();
    establish_decomposition();
    establish_state_storage(fields.size(), fields.size());

    DEBUG0(who, "Initializing state_linear to contain the data in wave space");
    if (mms >= 0) {

        INFO0(who, "Manufactured solution will be initialized at t = " << mms);
        INFO0(who, "Disabling bulk_rho and bulk_rho_u constraints"
                   " due to manufactured solution use");
        scenario->bulk_rho   = std::numeric_limits<real_t>::quiet_NaN();
        scenario->bulk_rho_u = std::numeric_limits<real_t>::quiet_NaN();

        accumulate_manufactured_solution(
                1, *msoln, 0, *state_nonlinear, *grid, *dgrid, *cop, *b, mms);
        state_linear->assign_from(*state_nonlinear);

    } else if (!dgrid->has_zero_zero_modes()) {

        DEBUG("Saving no non-zero information on non-zero-zero rank");
        msoln.reset();
        fill(*state_linear, 0);

    } else {

        INFO("Parabolic profile will be initialized with npower = " << npower);
        msoln.reset();
        fill(*state_linear, 0);

        INFO("Base field uses constant rho = " << scenario->bulk_rho);
        ArrayXr rho = ArrayXr::Constant(grid->N.y(), scenario->bulk_rho);
        ArrayXr v   = ArrayXr::LinSpaced(grid->N.y(),
                                         isothermal->lower_v,
                                         isothermal->upper_v);
        ArrayXr w   = ArrayXr::LinSpaced(grid->N.y(),
                                         isothermal->lower_w,
                                         isothermal->upper_w);
        ArrayXr T   = ArrayXr::LinSpaced(grid->N.y(),
                                         isothermal->lower_T,
                                         isothermal->upper_T);
        INFO("Base field uses v from " << v[0] << " to " << v.tail<1>()[0]);
        INFO("Base field uses w from " << w[0] << " to " << w.tail<1>()[0]);
        INFO("Base field uses T from " << T[0] << " to " << T.tail<1>()[0]);

        INFO("Finding normalization so u = (y*(L-y))^npower integrates to 1");
        real_t normalization = std::numeric_limits<real_t>::quiet_NaN();
        if (npower == 1) {
            // Mathematica: (Integrate[(x (L-x)),{x,0,L}]/L)^(-1)
            normalization = 6 / pow(grid->L.y(), 2);
        } else if (npower > 0) {
            // Mathematica: (Integrate[(x (L - x))^n, {x, 0, L}]/L)^(-1)
            //      -  (Gamma[-n] Gamma[3/2+n])
            //       / (2^(-1-2 n) L^(1+2 n) \[Pi]^(3/2) Csc[n \[Pi]])
            const real_t num1   = sin(npower * pi<real_t>());
            const real_t num2   = tgamma(-npower);
            const real_t num3   = tgamma(real_t(3)/2 + npower);
            const real_t denom1 = pow(           2, -1-2*npower);
            const real_t denom2 = pow( grid->L.y(),  1+2*npower);
            const real_t denom3 = pow(pi<real_t>(),  real_t(3)/2);
            normalization = - (num1 * num2 * num3 * grid->L.y())
                          /   (denom1 * denom2 * denom3);
        } else {
            // Degenerate case
            normalization = 1;
        }

        INFO("Preparing the wall-normal streamwise velocity profile");
        ArrayXr u(grid->N.y());
        for (int j = 0; j < u.size(); ++j) {
            const real_t y_j = b->collocation_point(j);
            u[j] = scenario->bulk_rho_u
                 * normalization
                 * pow(y_j * (grid->L.y() - y_j), npower);
        }

        INFO("Preparing specific internal energy using the equation of state");
        ArrayXr E = T / (scenario->gamma*(scenario->gamma - 1))
                  + (scenario->Ma*scenario->Ma/2) * (u*u + v*v + w*w);

        INFO("Converting to non-specific, conserved state");
        E *= rho; // Notice "E" now contains rho_E
        u *= rho; // Notice "u" now contains rho_u
        v *= rho; // Notice "v" now contains rho_v
        w *= rho; // Notice "w" now contains rho_w

        INFO("Converting conserved state to B-spline coefficients");
        suzerain::bsplineop_lu masslu(*cop);
        masslu.factor_mass(*cop);
        masslu.solve(1, E  .data(), 1, E  .size());
        masslu.solve(1, u  .data(), 1, u  .size());
        masslu.solve(1, v  .data(), 1, v  .size());
        masslu.solve(1, w  .data(), 1, w  .size());
        masslu.solve(1, rho.data(), 1, rho.size());

        INFO("Copying the coefficients directly into the zero-zero modes");
        Map<VectorXc>((*state_linear)[ndx::e  ].origin(), grid->N.y())
                = (E  ).cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::mx ].origin(), grid->N.y())
                = (u  ).cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::my ].origin(), grid->N.y())
                = (v  ).cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::mz ].origin(), grid->N.y())
                = (w  ).cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::rho].origin(), grid->N.y())
                = (rho).cast<complex_t>();

    }

    INFO0(who, "Saving the newly initialized state to disk");
    if (msoln) {
        save_restart(mms, restart_file, clobber);
    } else {
        save_restart(0,   restart_file, clobber);
    }

    return EXIT_SUCCESS;
}
