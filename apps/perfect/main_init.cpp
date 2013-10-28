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
#include <suzerain/blasius.h>
#include <suzerain/exprparse.hpp>
#include <suzerain/math.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/physical_view.hpp>
#include <suzerain/rholut.hpp>
#include <suzerain/support/grid_definition.hpp>
#include <suzerain/support/isothermal_definition.hpp>
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
"Initializes channel flow profile per --npower option unless --Re_x is\n"
"supplied.  In that case, a Blasius-like boundary layer profile is prepared.\n"
"The Blasius solution will uniquely fix both --upper_u and --upper_v.\n"
"Boundary layers homogenized by Largo employ a linear ramp in v.\n"
"When in doubt, please read through the source code.\n",
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
    using std::sqrt;
    using std::string;

    // Establish default grid and domain extents
    grid.reset(new support::grid_definition( 4 * pi<real_t>()     // Lx
                                           , 1                    // Nx
                                           , real_t(3) / 2        // DAFx
                                           , 2                    // Ly
                                           , 48                   // Ny
                                           , 8                    // k
                                           , 3                    // htdelta
                                           , 4 * pi<real_t>() / 3 // Lz
                                           , 1                    // Nz
                                           , real_t(3) / 2        // DAFz
            ));

    // Establish default scenario parameters
    scenario.reset(new scenario_definition(
                  100                                   // Re
                , real_t(3) / 2                         // Ma
                , real_t(7) / 10                        // Pr
                , 1                                     // bulk_rho
                , 1                                     // bulk_rho_u
                , numeric_limits<real_t>::quiet_NaN()   // bulk_rho_E, disabled
                , 0                                     // alpha
                , real_t(2) / 3                         // beta
                , real_t(14) / 10                       // gamma
            ));

    // Establish default isothermal boundary conditions
    isothermal.reset(new support::isothermal_definition(/* wall_T */ 1));

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
    real_t Re_x   = numeric_limits<real_t>::quiet_NaN();
    options.add_options()
        ("clobber",
         boost::program_options::bool_switch(),
         "Overwrite an existing restart file?")
        ("npower",
         boost::program_options::value<string>()
         ->default_value("1")
         ->notifier(boost::bind(&parse_nonnegative, _1, &npower, "npower")),
         "Power n in [0, 1] used to control the flatness of the channel"
         " \"parabolic\" streamwise velocity profile (y*(L-y))^n."
         " Using zero disables adding an initial streamwise profile.")
        ("Re_x",
         boost::program_options::value<string>()
         ->notifier(boost::bind(&parse_nonnegative, _1, &Re_x, "Re_x")),
         "Local Reynolds number Re_x = u_\\infty / \\nu / x controlling"
         " the Blasius profile used to initialize the boundary layer.")
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
             ->default_value("1/8")
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
             ->default_value("1/8")
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
    if (grid->two_sided() && options.variables().count("Re_x")) {
        WARN0("Conflicting two-sided grid and Blasius profile requested?");
    }
    options.conflicting_options("Re_x", "upper_u"); // Blasius sets freestream
    options.conflicting_options("Re_x", "upper_v"); // Blasius sets freestream

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
    suzerain::multi_array::fill(*state_linear, 0);
    if (options.variables().count("mms")) {

        INFO0(who, "Manufactured solution will be initialized at t = " << mms);
        INFO0(who, "Disabling bulk rho, rho_u, and rho_E constraints"
                   " due to manufactured solution use");
        scenario->bulk_rho   = numeric_limits<real_t>::quiet_NaN();
        scenario->bulk_rho_u = numeric_limits<real_t>::quiet_NaN();
        scenario->bulk_rho_E = numeric_limits<real_t>::quiet_NaN();

        accumulate_manufactured_solution(
                1, *msoln, 0, *state_nonlinear, *grid, *dgrid, *cop, *b, mms);
        state_linear->assign_from(*state_nonlinear);

    } else if (!dgrid->has_zero_zero_modes()) {

        DEBUG("Saving no non-zero information on non-zero-zero rank");
        msoln.reset(); // No manufactured solution in use

    } else {

        msoln.reset(); // No manufactured solution in use

        // Primitive state to be initialized depending on scenario
        // Pressure and temperature convenient for later adding pulses
        ArrayXr u = ArrayXr::Constant(Ny, numeric_limits<real_t>::quiet_NaN());
        ArrayXr v = ArrayXr::Constant(Ny, numeric_limits<real_t>::quiet_NaN());
        ArrayXr w = ArrayXr::Constant(Ny, numeric_limits<real_t>::quiet_NaN());
        ArrayXr T = ArrayXr::Constant(Ny, numeric_limits<real_t>::quiet_NaN());
        ArrayXr p = ArrayXr::Constant(Ny, numeric_limits<real_t>::quiet_NaN());

        if (    options.variables().count("Re_x")
             && !sg->formulation.enabled()) {
            INFO(who, "Initializing a classical boundary layer profile");

            INFO("Blasius base profile for u, v, and T uses Re_x = "
                 << Re_x << ", Pr = " << scenario->Pr);
            INFO("Linear ramp uses for w from "
                 << isothermal->lower_w << " to " << isothermal->upper_w);
            shared_ptr<gsl_spline> fit_u(suzerain_blasius_u(Re_x),
                                         gsl_spline_free);
            shared_ptr<gsl_spline> fit_v(suzerain_blasius_v(Re_x),
                                         gsl_spline_free);
            shared_ptr<gsl_spline> fit_T(suzerain_blasius_T(Re_x, scenario->Pr),
                                         gsl_spline_free);
            shared_ptr<gsl_interp_accel> accel(gsl_interp_accel_alloc(),
                                               gsl_interp_accel_free);

            INFO0(who, "Computing mean freestream behavior per plate scenario");
            // See plate_treatment.tex for origin of T-dependent scaling
            // Blasius profile fixes upper_{u,v} for the freestream...
            isothermal->upper_u   = sqrt(isothermal->upper_T); // Eqn (6)
            isothermal->upper_v   = gsl_spline_eval(fit_v.get(),
                                                    grid->L.y(), accel.get())
                                  * isothermal->upper_u;
            isothermal->upper_rho = pow(isothermal->upper_T,
                                        scenario->beta - 0.5); // Eqn (8)
            // ...where we've blissfully ignored lower_{u,v} (tra la la)

            // Denormalize calls honor lower, upper values for u, w, T.
            // Wall-normal v does not use one as fit_v does not achieve 1
            for (int j = 0; j < Ny; ++j) {
                const double y_j = b->collocation_point(j);
                u[j] = isothermal->denormalize_u(
                        gsl_spline_eval(fit_u.get(), y_j, accel.get()));
                v[j] = gsl_spline_eval(fit_v.get(), y_j, accel.get())
                     * isothermal->upper_u;
                w[j] = isothermal->denormalize_w(y_j / grid->L.y());
                T[j] = isothermal->denormalize_T(
                        gsl_spline_eval(fit_T.get(), y_j, accel.get()));
            }

            // Find pressure using constant density taken from freestream
            p = isothermal->upper_rho * T / scenario->gamma;
            INFO("Base field uses p from " << p(0) << " to " << p(Ny - 1));

        } else if (    options.variables().count("Re_x")
                    && sg->formulation.enabled()) {
            INFO(who, "Initializing approximate homogenized boundary layer");

            INFO("Blasius base profile for u and T uses Re_x = "
                 << Re_x << ", Pr = " << scenario->Pr);
            INFO("Linear ramp uses for v from "
                 << isothermal->lower_v << " to " << isothermal->upper_v);
            INFO("Linear ramp uses for w from "
                 << isothermal->lower_w << " to " << isothermal->upper_w);
            shared_ptr<gsl_spline> fit_u(suzerain_blasius_u(Re_x),
                                         gsl_spline_free);
            shared_ptr<gsl_spline> fit_T(suzerain_blasius_T(Re_x, scenario->Pr),
                                         gsl_spline_free);
            shared_ptr<gsl_interp_accel> accel(gsl_interp_accel_alloc(),
                                               gsl_interp_accel_free);

            INFO(who, "Computing mean freestream behavior per plate scenario");
            // See plate_treatment.tex for origin of T-dependent scaling
            // Blasius profile fixes upper_u for the freestream...
            isothermal->upper_u   = sqrt(isothermal->upper_T); // Eqn (6)
            isothermal->upper_rho = pow(isothermal->upper_T,
                                        scenario->beta - 0.5); // Eqn (8)
            // ...where we've blissfully ignored lower_u (tra la la)

            // Denormalize calls honor lower, upper values for u, v, w, T.
            for (int j = 0; j < Ny; ++j) {
                const double y_j = b->collocation_point(j);
                u[j] = isothermal->denormalize_u(
                        gsl_spline_eval(fit_u.get(), y_j, accel.get()));
                v[j] = isothermal->denormalize_v(y_j / grid->L.y());
                w[j] = isothermal->denormalize_w(y_j / grid->L.y());
                T[j] = isothermal->denormalize_T(
                        gsl_spline_eval(fit_T.get(), y_j, accel.get()));
            }

            // Find pressure using constant density taken from freestream
            p = isothermal->upper_rho * T / scenario->gamma;
            INFO("Base field uses p from " << p(0) << " to " << p(Ny - 1));

        } else {
            INFO(who, "Initializing channel flow profile");

            INFO("Base field starts from linear ramps for u, v, w, and T");
            for (int j = 0; j < Ny; ++j) {
                const double norm_y_j = b->collocation_point(j) / grid->L.y();
                u[j] = isothermal->denormalize_u(norm_y_j);
                v[j] = isothermal->denormalize_v(norm_y_j);
                w[j] = isothermal->denormalize_w(norm_y_j);
                T[j] = isothermal->denormalize_T(norm_y_j);
            }
            INFO("Base field uses u from " << u(0) << " to " << u(Ny - 1));
            INFO("Base field uses v from " << v(0) << " to " << v(Ny - 1));
            INFO("Base field uses w from " << w(0) << " to " << w(Ny - 1));
            INFO("Base field uses T from " << T(0) << " to " << T(Ny - 1));

            INFO("Parabolic profile will be added with npower = " << npower);
            INFO("Finding normalization so u = (y*(L-y))^npower integral is 1");
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
                // Degenerate npower == 0 case avoids tgamma(...) domain issues
                normalization = 0;
            }

            INFO("Adding the requested parabolic streamwise velocity profile");
            for (int j = 0; j < Ny; ++j) {
                const real_t y_j = b->collocation_point(j);
                u(j) += (   scenario->bulk_rho_u
                          - (isothermal->lower_u + isothermal->upper_u)/2)
                      * normalization
                      * pow(y_j * (Ly - y_j), npower);
            }

            // Find pressure from temperature and constant density assumption
            p = scenario->bulk_rho * T / scenario->gamma;
            INFO("Base field uses p from " << p(0) << " to " << p(Ny - 1));

        }

        if (options.variables().count("acoustic_strength")) {
            const real_t left  = (Ly / 2) * (1 - acoustic_support);
            const real_t right = (Ly / 2) * (1 + acoustic_support);
            INFO("Adding wall-normal acoustic pulse with strength "
                 << acoustic_strength
                 << " on (" << left << ", " << right << ")");
            for (int j = 0; j < Ny; ++j) {// Baum et al JCP 1994 pp 254-5
                const real_t v0   = v(j);
                const real_t rho0 = scenario->gamma * p(j) / T(j);
                v(j) += math::bump::shifted(b->collocation_point(j),
                                            left, right, acoustic_strength,
                                            real_t(0), acoustic_stiffness);
                p(j) += rho0 * std::sqrt(T(j)) * scenario->Ma * (v(j) - v0);
                T(j)  = (scenario->gamma * p(j))
                      / (rho0 + scenario->Ma * (v(j) - v0));
            }
        }

        if (options.variables().count("entropy_strength")) {
            const real_t left  = (Ly / 2) * (1 - entropy_support);
            const real_t right = (Ly / 2) * (1 + entropy_support);
            INFO("Adding wall-normal entropy pulse with strength "
                 << entropy_strength
                 << " on (" << left << ", " << right << ")");
            for (int j = 0; j < Ny; ++j) {  // Baum et al JCP 1994 pp 254-5
                T(j) += math::bump::shifted(b->collocation_point(j),
                                            left, right, entropy_strength,
                                            real_t(0), entropy_stiffness);
            }
        }

        INFO("Computing density from pressure and temperature profiles");
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
