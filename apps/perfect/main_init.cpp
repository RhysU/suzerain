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
#include <suzerain/exprparse.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
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

/** Provides parsing and validation of several options. */
static void
parse_nonnegative(const std::string& s, suzerain::real_t *t, const char *n)
{
    const suzerain::real_t v = suzerain::exprparse<suzerain::real_t>(s, n);
    suzerain::validation::ensure_nonnegative(v, n);
    *t = v;
}

/** Instantiate and invoke the application. */
int main(int argc, char **argv)
{
    suzerain::perfect::driver_init app(revstr);
    return app.run(argc, argv);
}

/** The logic executed by main(). */
int
suzerain::perfect::driver_init::run(int argc, char **argv)
{
    using boost::math::constants::pi;
    using boost::math::tgamma;
    using std::numeric_limits;
    using std::pow;
    using std::sin;
    using std::string;

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
    real_t mms    = numeric_limits<real_t>::quiet_NaN();
    real_t npower = numeric_limits<real_t>::quiet_NaN();
    options.add_options()
        ("clobber",
         boost::program_options::bool_switch(),
         "Overwrite an existing restart file?")
        ("npower",
         boost::program_options::value<string>()
         ->default_value("1")
         ->notifier(boost::bind(&parse_nonnegative, _1, &npower, "npower")),
         "Power n in [0, 1] used to control the flatness of the"
         " \"parabolic\" streamwise velocity profile (y*(L-y))^n.")
        ("mms",
         boost::program_options::value<string>()
         ->notifier(boost::bind(&parse_nonnegative, _1, &mms, "mms")),
         "If given, prepare a manufactured solution at the specified time.")
    ;

    real_t acoustic_strength  = numeric_limits<real_t>::quiet_NaN();
    real_t acoustic_stiffness = numeric_limits<real_t>::quiet_NaN();
    real_t acoustic_support   = numeric_limits<real_t>::quiet_NaN();
    {
        boost::program_options::options_description pulse_acoustic(
                "Add an acoustic pulse similarly to Baum et al. JCP 1994");
        pulse_acoustic.add_options()
            ("acoustic_strength",
             boost::program_options::value<string>()
             ->notifier(boost::bind(&parse_nonnegative, _1,
                                    &acoustic_strength, "acoustic_strength")),
             "Strength of the acoustic pulse measured in velocity units")
            ("acoustic_stiffness",
             // See notebooks/Bump_Function_Power.nb for why 4 is chosen
             boost::program_options::value<string>()
             ->default_value("4")
             ->notifier(boost::bind(&parse_nonnegative, _1,
                                    &acoustic_stiffness, "acoustic_stiffness")),
             "Stiffness of the acoustic pulse measured by bump function power")
            ("acoustic_support",
             boost::program_options::value<string>()
             ->default_value("1/4")
             ->notifier(boost::bind(&parse_nonnegative, _1,
                                    &acoustic_support, "acoustic_support")),
             "Fraction of the domain over which the pulse is supported")
        ;
        options.options().add(pulse_acoustic);
    }
    real_t entropy_strength  = numeric_limits<real_t>::quiet_NaN();
    real_t entropy_stiffness = numeric_limits<real_t>::quiet_NaN();
    real_t entropy_support   = numeric_limits<real_t>::quiet_NaN();
    {
        boost::program_options::options_description pulse_entropy(
                "Add an entropy pulse similarly to Baum et al. JCP 1994");
        pulse_entropy.add_options()
            ("entropy_strength",
             boost::program_options::value<string>()
             ->notifier(boost::bind(&parse_nonnegative, _1,
                                    &entropy_strength, "entropy_strength")),
             "Strength of the entropy pulse measured in temperature units")
            ("entropy_stiffness",
             // See notebooks/Bump_Function_Power.nb for why 4 is chosen
             boost::program_options::value<string>()
             ->default_value("4")
             ->notifier(boost::bind(&parse_nonnegative, _1,
                                    &entropy_stiffness, "entropy_stiffness")),
             "Stiffness of the entropy pulse measured by bump function power")
            ("entropy_support",
             boost::program_options::value<string>()
             ->default_value("1/4")
             ->notifier(boost::bind(&parse_nonnegative, _1,
                                    &entropy_support, "entropy_support")),
             "Fraction of the domain over which the pulse is supported")
        ;
        options.options().add(pulse_entropy);
    }

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros becomes usable)
    std::vector<std::string> positional = initialize(argc, argv);
    if (positional.size() != 1) {
        FATAL0("Exactly one restart file name must be specified");
        return EXIT_FAILURE;
    }
    const std::string restart_file = positional[0];
    const bool clobber = options.variables()["clobber"].as<bool>();
    if (npower < 0 || npower > 1) {
        FATAL0("npower in [0,1] required");
        return EXIT_FAILURE;
    }
    if (grid->k < 4) {
        FATAL0("k >= 4 required for two non-trivial wall-normal derivatives");
        return EXIT_FAILURE;
    }
    options.conflicting_options("acoustic_strength", "mms");
    options.conflicting_options("entropy_strength",  "mms");
    options.conflicting_options("acoustic_strength", "entropy_strength");

    DEBUG0(who, "Establishing runtime parallel infrastructure and resources");
    establish_ieee_mode();
    load_grid_and_operators(NULL);
    log_discretization_quality();
    establish_decomposition();
    establish_state_storage(fields.size(), fields.size());

    // Shorthand
    const real_t Ly = grid->L.y();
    const int    Ny = grid->N.y();

    DEBUG0(who, "Initializing state_linear to contain the data in wave space");
    fill(*state_linear, 0);
    if (options.variables().count("mms")) {

        INFO0(who, "Manufactured solution will be initialized at t = " << mms);
        INFO0(who, "Disabling bulk_rho and bulk_rho_u constraints"
                   " due to manufactured solution use");
        scenario->bulk_rho   = numeric_limits<real_t>::quiet_NaN();
        scenario->bulk_rho_u = numeric_limits<real_t>::quiet_NaN();

        accumulate_manufactured_solution(
                1, *msoln, 0, *state_nonlinear, *grid, *dgrid, *cop, *b, mms);
        state_linear->assign_from(*state_nonlinear);

    } else if (!dgrid->has_zero_zero_modes()) {

        DEBUG("Saving no non-zero information on non-zero-zero rank");
        msoln.reset(); // No manufactured solution in use

    } else {

        msoln.reset(); // No manufactured solution in use

        INFO("Computing mean base field for use in zero-zero modes");
        // Computation of pressure per suzerain::rholut::p(...)
        ArrayXr u = ArrayXr::LinSpaced(Ny, isothermal->lower_u,
                                           isothermal->upper_u);
        ArrayXr v = ArrayXr::LinSpaced(Ny, isothermal->lower_v,
                                           isothermal->upper_v);
        ArrayXr w = ArrayXr::LinSpaced(Ny, isothermal->lower_w,
                                           isothermal->upper_w);
        ArrayXr T = ArrayXr::LinSpaced(Ny, isothermal->lower_T,
                                           isothermal->upper_T);
        ArrayXr p = scenario->bulk_rho * T / scenario->gamma;
        INFO("Base field uses constant rho = " << scenario->bulk_rho);
        INFO("Base field uses u from " << u(0) << " to " << u(Ny - 1));
        INFO("Base field uses v from " << v(0) << " to " << v(Ny - 1));
        INFO("Base field uses w from " << w(0) << " to " << w(Ny - 1));
        INFO("Base field uses T from " << T(0) << " to " << T(Ny - 1));
        INFO("Base field uses p from " << p(0) << " to " << p(Ny - 1));

        INFO("Parabolic profile will be added with npower = " << npower);
        INFO("Finding normalization so u = (y*(L-y))^npower integrates to 1");
        real_t normalization = numeric_limits<real_t>::quiet_NaN();
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
            const real_t denom1 = pow(2, -1-2*npower);
            const real_t denom2 = pow(Ly,  1+2*npower);
            const real_t denom3 = pow(pi<real_t>(),  real_t(3)/2);
            normalization = - (num1 * num2 * num3 * Ly)
                          /   (denom1 * denom2 * denom3);
        } else {
            // Degenerate npower == 0 case to avoid tgamma(...) domain issues
            normalization = 1;
        }

        INFO("Adding the requested parabolic streamwise velocity profile");
        for (int j = 0; j < Ny; ++j) {
            const real_t y_j = b->collocation_point(j);
            u(j) += (   scenario->bulk_rho_u
                      - (isothermal->lower_u + isothermal->upper_u)/2)
                  * normalization
                  * pow(y_j * (Ly - y_j), npower);
        }

        if (options.variables().count("acoustic_strength")) {
            INFO("Adding acoustic pulse with strength " << acoustic_strength);
            // TODO
        }
        if (options.variables().count("entropy_strength")) {
            INFO("Adding entropy pulse with strength " << entropy_strength);
            // TODO
        }

        INFO("Computing density from pressure and temperature");
        ArrayXr rho = scenario->gamma * p / T;
        T.resize(0); // Mark irrelevant for further use

        INFO("Computing total energy and momentum from primitive state");
        ArrayXr rho_E(Ny);
        for (int j = 0; j < Ny; ++j) {
            const Vector3r m(rho(j) * u(j),
                             rho(j) * v(j),
                             rho(j) * w(j));
            rho_E(j) = rholut::energy_internal(scenario->gamma, p(j))
                     + rholut::energy_kinetic (scenario->Ma, rho(j), m);
        }
        ArrayXr rho_u = rho * u;
        ArrayXr rho_v = rho * v;
        ArrayXr rho_w = rho * w;
        u.resize(0); // Mark irrelevant for further use
        v.resize(0); // Mark irrelevant for further use
        w.resize(0); // Mark irrelevant for further use
        p.resize(0); // Mark irrelevant for further use

        INFO("Converting conserved state to B-spline coefficients");
        suzerain::bsplineop_lu masslu(*cop);
        masslu.factor_mass(*cop);
        masslu.solve(1, rho_E.data(), 1, rho_E.size());
        masslu.solve(1, rho_u.data(), 1, rho_u.size());
        masslu.solve(1, rho_v.data(), 1, rho_v.size());
        masslu.solve(1, rho_w.data(), 1, rho_w.size());
        masslu.solve(1, rho  .data(), 1, rho  .size());

        INFO("Copying the coefficients directly into the zero-zero modes");
        Map<VectorXc>((*state_linear)[ndx::e  ].origin(), Ny)
                = rho_E.cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::mx ].origin(), Ny)
                = rho_u.cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::my ].origin(), Ny)
                = rho_v.cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::mz ].origin(), Ny)
                = rho_w.cast<complex_t>();
        Map<VectorXc>((*state_linear)[ndx::rho].origin(), Ny)
                = rho  .cast<complex_t>();

    }

    INFO0(who, "Saving the newly initialized state to disk");
    if (msoln) {
        save_restart(mms, restart_file, clobber);
    } else {
        save_restart(0,   restart_file, clobber);
    }

    return EXIT_SUCCESS;
}
