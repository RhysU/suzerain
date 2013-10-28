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
 * Application executing \ref suzerain::perfect::driver_advance::run.
 */

#include <esio/esio.h>
#include <largo/largo.h>

#include <suzerain/common.hpp>
#include <suzerain/constraint.hpp>
#include <suzerain/constraint_treatment.hpp>
#include <suzerain/isothermal_specification.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/noise_definition.hpp>
#include <suzerain/zgbsv_specification.hpp>

#include "driver.hpp"
#include "hybrid_operator.hpp"
#include "isothermal_mass_operator.hpp"
#include "nonlinear_operator.hpp"
#include "nonreflecting_treatment.hpp"
#include "perfect.hpp"

#pragma warning(disable:1419)

namespace suzerain {

namespace perfect {

/** Application for initializing new restart files. */
struct driver_advance : public driver
{
    driver_advance(const std::string& revstr)
        : driver("Compressible, perfect gas simulation time advancement",
                 "RESTART-FILE",
                 "",
                 revstr)
        , who("advance")
    {}

    /** Implementation below in this file */
    int run(int argc, char **argv);

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

} // namespace perfect

} // namespace suzerain

// Provided by main_advance_svnrev.{c,h} so revstr updates are merely relinking
extern "C" const char revstr[];

/** Instantiate and invoke the application */
int main(int argc, char **argv)
{
    suzerain::perfect::driver_advance app(revstr);
    return app.run(argc, argv);
}

int
suzerain::perfect::driver_advance::run(int argc, char **argv)
{
    using namespace std;
    using boost::math::isnan;

    // Storage for binary-specific options
    const support::noise_definition noisedef;
    string solver_spec(static_cast<string>(suzerain::zgbsv_specification()));
    string implicit("rhome_xyz");
    string undriven;

    // Register binary-specific options
    options.add_definition(const_cast<support::noise_definition&>(noisedef));
    options.add_options()
        ("explicit", boost::program_options::bool_switch(),
                     "Use purely explicit operators")
        ("implicit", boost::program_options::value(&implicit)
                         ->implicit_value(implicit),
                     "Use hybrid implicit/explicit operators, optionally"
                     " choosing 'rhome_xyz' or 'rhome_y' linearized treatment")
        ("solver",   boost::program_options::value(&solver_spec)
                         ->default_value(solver_spec),
                     "Use the specified algorithm for any --implicit solves")
        ("undriven", boost::program_options::value(&undriven)
                         ->implicit_value("all"),
                     "Disable all or some physics-related driving forces"
                     " by specifying 'all' or 'rho', 'rho_u', etc.")
    ;

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros becomes usable)
    vector<string> positional = initialize(argc, argv);

    // Select type of timestepping operators to use (default implicit)
    options.conflicting_options("explicit", "implicit");
    options.conflicting_options("explicit", "solver");
    const bool use_explicit =  options.variables()["explicit"].as<bool>();
    const bool use_implicit =  options.variables().count("implicit")
                            || !use_explicit;

    // Validate implicit option; implicit_value "captures" following argument
    // Establish common_block.linearization based on the result
    if (use_implicit) {
        boost::algorithm::trim(implicit);
        if (implicit == "rhome_xyz") {
            common_block.linearization = linearize::rhome_xyz;
        } else if (implicit == "rhome_y") {
            common_block.linearization = linearize::rhome_y;
        } else {
            FATAL0("Unknown --implicit argument:  " << implicit);
            return EXIT_FAILURE;
        }
    } else {
        common_block.linearization = linearize::none;
        implicit = "none";
    }

    if (positional.size() != 1) {
        FATAL0("Exactly one restart file name must be specified");
        return EXIT_FAILURE;
    }
    const string restart_file = positional[0];

    INFO0(who, "Loading restart file: " << restart_file);
    real_t initial_t = numeric_limits<real_t>::quiet_NaN();
    {
        // Preserve exact restart file details via Push/Pop/Merge below
        shared_ptr<scenario_definition> restart_scenario
                = make_shared<scenario_definition>();

        // Load the restart details with state going into state_linear
        shared_ptr<boost::remove_pointer<esio_handle>::type> h( // RAII
                esio_handle_initialize(MPI_COMM_WORLD),
                esio_handle_finalize);
        esio_file_open(h.get(), restart_file.c_str(), 0);
        restart_scenario.swap(scenario);                        // Push
        load_restart(h.get(), initial_t);
        restart_scenario.swap(scenario);                        // Pop
        scenario->populate(*restart_scenario, true);            // Merge

        // Adjust total energy as necessary to account for any scenario change
        state_nonlinear->assign_from(*state_linear);
        adjust_scenario(*state_nonlinear, *scenario, *grid, *dgrid, *cop,
                        restart_scenario->Ma, restart_scenario->gamma);
        state_linear->assign_from(*state_nonlinear);
    }

    if (msoln) {
        INFO0(who, "Restart file prescribes a manufactured solution");
        if (!(isnan)(scenario->bulk_rho)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho = "
                  << scenario->bulk_rho);
        }
        if (!(isnan)(scenario->bulk_rho_u)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho_u = "
                  << scenario->bulk_rho_u);
        }
        if (!(isnan)(scenario->bulk_rho_E)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho_E = "
                  << scenario->bulk_rho_E);
        }
    }

    // If requested, add noise to the momentum fields at startup (expensive).
    if (noisedef.percent > 0) {
        state_nonlinear->assign_from(*state_linear);
        add_noise(*state_nonlinear, noisedef, *scenario,
                  *grid, *dgrid, *cop, *b);
        state_linear->assign_from(*state_nonlinear);
    }

    // Initialize any requested slow growth forcing workspace
    sg->workspace = NULL; // Defensive
    common_block.slow_treatment = slowgrowth::none;
    if (sg->formulation.enabled()) {
        common_block.slow_treatment = slowgrowth::largo;
        const std::string& model_name = sg->formulation.name();
        INFO0("Allocating Largo model \"" << model_name
              << "\" for 5 state variables");
        largo_allocate(&sg->workspace, model_name.c_str(), 5, 0, 0, "dns");
        if ((isnan)(sg->grdelta)) {
            WARN0("Slow growth rate grdelta is NaN");
        }
    }

    // Prepare any necessary, problem-specific constraints
    shared_ptr<constraint::treatment<operator_common_block> > constrainer(
            new constraint::treatment<operator_common_block>(
                    scenario->Ma, *dgrid, *b, common_block));
    if        (grid->two_sided()) { // Channel per channel_treatment.tex

        INFO0(who, "Establishing driving, channel-like state constraints");
        constrainer->physical[ndx::rho].reset(
                new constraint::reference_bulk(scenario->bulk_rho  , *b));
        constrainer->physical[ndx::mx ].reset(
                new constraint::reference_bulk(scenario->bulk_rho_u, *b));
        constrainer->physical[ndx::e  ].reset(
                new constraint::reference_bulk(scenario->bulk_rho_E, *b));

    } else if (grid->one_sided()) { // Flat plate per plate_treatment.tex

        INFO0(who, "Computing mean freestream behavior per plate scenario");
        const real_t T_inf   = isothermal->upper_T;                 // Brevity
        const real_t u_inf   = sqrt(T_inf);                         // Eqn ( 6)
        const real_t mx_inf  = pow(T_inf, scenario->beta);          // Eqn ( 7)
        const real_t rho_inf = pow(T_inf, scenario->beta - 0.5);    // Eqn ( 8)
        const real_t e_inf   = pow(T_inf, scenario->beta - 0.5) * ( // Eqn (11)
                                   T_inf/(scenario->gamma*(scenario->gamma-1))
                                 + pow(scenario->Ma, 2) / 2 * (
                                       T_inf
                                     + pow(isothermal->upper_v, 2)
                                     + pow(isothermal->upper_w, 2)
                                   )
                               );

        INFO0(who, "Setting freestream reference state on upper boundary");
        isothermal->upper_rho = rho_inf;
        isothermal->upper_u   = u_inf;
        if (sg->formulation.enabled()) {
            // Homogenization modifies the inviscid characteristics.
            // See discussion (or placeholder thereof) at Redmine #2982
            isothermal->upper_v = isothermal->lower_v
                                - grid->L.y()*sg->grdelta;
            INFO0(who, "Matching upper_v with prescribed"
                       " grdelta, lower_v, and Ly: " << isothermal->upper_v);
        } else {
            DEBUG0(who, "Keeping reference upper_v = " << isothermal->upper_v);
        }
        DEBUG0(who, "Keeping reference upper_w = " << isothermal->upper_w);
        DEBUG0(who, "Keeping reference upper_T = " << isothermal->upper_T);

        if (scenario->bulk_rho) {
            WARN0(who, "Removing channel-like bulk_rho setting");
            scenario->bulk_rho   = numeric_limits<real_t>::quiet_NaN();
        }
        if (scenario->bulk_rho_u) {
            WARN0(who, "Removing channel-like bulk_rho_u setting");
            scenario->bulk_rho_u = numeric_limits<real_t>::quiet_NaN();
        }
        if (scenario->bulk_rho_E) {
            WARN0(who, "Removing channel-like bulk_rho_E setting");
            scenario->bulk_rho_E = numeric_limits<real_t>::quiet_NaN();
        }

        INFO0(who, "Establishing driving, freestream-like state constraints");
        constrainer->physical[ndx::e  ].reset(
                new constraint::constant_upper(e_inf,   *cop, 0));
        constrainer->physical[ndx::mx ].reset(
                new constraint::constant_upper(mx_inf,  *cop, 0));
        constrainer->physical[ndx::rho].reset(
                new constraint::constant_upper(rho_inf, *cop, 0));

        if (isothermal->upper_v > 0 && (boost::math::isinf)(scenario->Re)) {
            WARN0(who, "Nonreflecting viscous outflow boundary problematic"
                       << " (Redmine #2983)");
        }

    } else {

        FATAL0(who, "Sanity error in constraint selection");
        return EXIT_FAILURE;

    }

    // Prepare spatial operators depending on requested advancement type.
    // Notice constrainer always wraps the implicit operator
    // and that the same nonlinear_operator is used pervasively.
    L = constrainer;
    N.reset(new nonlinear_operator(
                *scenario, *grid, *dgrid, *cop, *b, common_block, *sg, msoln));
    if (use_explicit) {

        INFO0(who, "Initializing fully explicit spatial operators");
        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            shared_ptr<nonreflecting_treatment> nonreflecting(
                    new nonreflecting_treatment(
                        *scenario, *isothermal, *grid, *dgrid, *cop, *b));
            nonreflecting->N = N;
            N = nonreflecting;
        }
        constrainer->L.reset(new isothermal_mass_operator(
                    *scenario, *isothermal,
                    *grid, *dgrid, *cop, *b, common_block));

    } else if (use_implicit) {

        INFO0(who, "Initializing hybrid implicit/explicit spatial operators");
        INFO0(who, "Implicit linearization employed: " << implicit);
        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            FATAL0(who, "Nonreflecting upper boundary treatment"
                        " not currently usable with implicit advance");
            return EXIT_FAILURE;
        }
        constrainer->L.reset(new isothermal_hybrid_linear_operator(
                    solver_spec, *scenario, *isothermal,
                    *grid, *dgrid, *cop, *b, common_block));

    } else {

        FATAL0(who, "Sanity error in operator selection");
        return EXIT_FAILURE;

    }

    // Use --undriven as a testing- and debugging-related tool.
    // For example, to investigate nonreflecting boundary condition behavior.
    if (options.variables().count("undriven")) {
        boost::algorithm::trim(undriven);
        const size_t j = distance(
                ndx::identifier.begin(),
                find(ndx::identifier.begin(), ndx::identifier.end(), undriven));
        if (j < ndx::identifier.static_size) {
            INFO0(who, "Disabling any " << ndx::description[j]
                       << " physical constraint per --undriven=" << undriven);
            constrainer->physical[j] = constrainer->none;
        } else if (undriven == "all") {
            INFO0(who, "Disabling all physical constraints per --undriven="
                       << undriven);
            fill(constrainer->physical.begin(), constrainer->physical.end(),
                 constrainer->none);
        } else {
            FATAL0("Unknown --undriven argument:  " << undriven);
            return EXIT_FAILURE;
        }
    }

    // Perform final housekeeping and then advance time as requested
    establish_ieee_mode();
    log_discretization_quality();
    prepare_controller(initial_t, dgrid->chi());
    save_metadata();
    const real_t elapsed_wall_time = advance_controller(); // Negative on error

    // If we advanced by any time steps, log the observed linearization error
    if (   !use_explicit
        && elapsed_wall_time >= 0
        && controller->current_nt() > 0) {
        log_linearization_error(build_timeprefix(controller->current_t(),
                                                 controller->current_nt()),
                                controller->current_t(),
                                controller->current_nt());
    }

    // Deallocate any slow growth forcing workspace
    if (sg->formulation.enabled()) {
        if (sg->workspace) {
            largo_deallocate(&sg->workspace);
        }
        sg->workspace = NULL; // Defensive
    }

    // Report error to the OS iff advance_control reported an error
    return elapsed_wall_time < 0;
}
