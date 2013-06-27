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

#include <suzerain/common.hpp>
#include <suzerain/constraint.hpp>
#include <suzerain/constraint_treatment.hpp>
#include <suzerain/error.h>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/noise_definition.hpp>
#include <suzerain/zgbsv_specification.hpp>
#include <suzerain/hybrid_residual_operator.hpp>

#include "driver.hpp"
//#include "channel_treatment.hpp"
#include "explicit_operator.hpp"

#include "hybrid_operator.hpp"
#include "nonreflecting_treatment.hpp"

#ifdef SUZERAIN_HAVE_LARGO
#include "largo/largo.h" 
#endif

#pragma warning(disable:1419)

namespace suzerain { namespace reacting {

/** Application for initializing new restart files. */
struct driver_advance : public driver
{
    driver_advance(const std::string& revstr)
        : driver("Compressible, reacting flow simulation time advancement",
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

} /* namespace reacting */ } /* namespace suzerain */

// Provided by main_advance_svnrev.{c,h} so revstr updates are merely relinking
extern "C" const char revstr[];

/** Instantiate and invoke the application */
int main(int argc, char **argv)
{
    suzerain::reacting::driver_advance app(revstr);
    return app.run(argc, argv);
}

int
suzerain::reacting::driver_advance::run(int argc, char **argv)
{
    using boost::math::isnan;
    using std::numeric_limits;
    using std::string;
    using std::vector;

    // Storage for binary-specific options
    const support::noise_definition noisedef;
    string solver_spec(static_cast<string>(suzerain::zgbsv_specification()));
    string filter_spec("none");
    string undriven;

    // Register binary-specific options
    options.add_definition(const_cast<support::noise_definition&>(noisedef));
    options.add_options()
        ("explicit", boost::program_options::bool_switch(),
                     "Use purely explicit operators")
        ("implicit", boost::program_options::bool_switch(),
                     "Use hybrid implicit/explicit operators")
        ("solver",   boost::program_options::value(&solver_spec)
                         ->default_value(solver_spec),
                     "Use the specified algorithm for any implicit solves")
        ("filter",   boost::program_options::value(&filter_spec)
                         ->default_value(filter_spec),
                     "Use the specified type to construct filter source")
        ("undriven", boost::program_options::value(&undriven)
                         ->implicit_value("all"),
                         "Disable all or some constraint-based driving forces"
                         " by specifying 'all', 'rho', 'rho_u', or 'rho_E'");
    ;

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros becomes usable)
    vector<string> positional = initialize(argc, argv);

    // Select type of timestepping operators to use (default implicit)
    options.conflicting_options("explicit", "implicit");
    const bool use_explicit =  options.variables()["explicit"].as<bool>();
    const bool use_implicit =  options.variables()["implicit"].as<bool>()
                            || !use_explicit;

    // Only one implicit option now, so just use if we got '--implicit'
    if (use_implicit) {
        common_block.linearization = linearize::rhome_y;
    } else {
        common_block.linearization = linearize::none;
    }

    // Select type of filtering to use (default none)
    const bool use_filter =  options.variables().count("filter");

    if (use_filter) {
        boost::algorithm::trim(filter_spec);
        if (filter_spec == "none") {
            common_block.filter_treatment = filter::none;
        } else if (filter_spec == "cook") {
            common_block.filter_treatment = filter::cook;
        } else if (filter_spec == "viscous") {
            common_block.filter_treatment = filter::viscous;
        } else {
            FATAL0("Unknown --filter argument:  " << filter_spec);
            return EXIT_FAILURE;
        }
    } else {
        common_block.filter_treatment = filter::none;
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
        shared_ptr<channel_definition> restart_chdef
                = make_shared<channel_definition>();

        // Load the restart details with state going into state_linear
        shared_ptr<boost::remove_pointer<esio_handle>::type> h( // RAII
                esio_handle_initialize(MPI_COMM_WORLD),
                esio_handle_finalize);
        esio_file_open(h.get(), restart_file.c_str(), 0);
        restart_chdef.swap(chdef);
        load_restart(h.get(), initial_t);
        restart_chdef.swap(chdef);
        chdef->populate(*restart_chdef, true);

        // FIXME: filter source
        //        initialize with info from restart here (?)

        // Adjust total energy as necessary to account for any scenario change
        state_nonlinear->assign_from(*state_linear);
        // FIXME: Currently no functionality to adjust scenario
        //adjust_scenario(*state_nonlinear, *scenario, *grid, *dgrid, *cop, restart_scenario->Ma, restart_scenario->gamma);
        state_linear->assign_from(*state_nonlinear);
    }

    INFO0(who, "Initializing antioch_constitutive");
    cmods->init_antioch();

    common_block.Ns = cmods->Ns();

    // However, if msoln was provided, match its contents to other members
    // Do here b/c cmods doesn't make sense until after init_antioch
    if (msoln) {
        if (cmods) msoln->match(*cmods);
        if (grid)  msoln->match(*grid);
    }

    if (msoln) {
        INFO0(who, "Restart file prescribes a manufactured solution");
        if (!(isnan)(chdef->bulk_rho)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho = "
                  << chdef->bulk_rho);
        }
        if (!(isnan)(chdef->bulk_rho_u)) {
            WARN0(who, "Manufactured solution incompatible with bulk_rho_u = "
                  << chdef->bulk_rho_u);
        }
    }

    // If requested, add noise to the momentum fields at startup (expensive).
    if (noisedef.percent > 0) {
        state_nonlinear->assign_from(*state_linear);
        add_noise(*state_nonlinear, noisedef,
                  *grid, *dgrid, *cop, *b);
        state_linear->assign_from(*state_nonlinear);
    }


    // Prepare any necessary, problem-specific constraints
    shared_ptr<constraint::treatment<operator_common_block> > constrainer(
            new constraint::treatment<operator_common_block>(
                    1.0, *dgrid, *b, common_block));
    if        (grid->two_sided()) { // Channel per channel_treatment.tex

        INFO0(who, "Establishing driving, channel-like state constraints");
        (*constrainer)[ndx::rho].reset(
                new constraint::reference_bulk(chdef->bulk_rho  , *b));
        (*constrainer)[ndx::mx ].reset(
                new constraint::reference_bulk(chdef->bulk_rho_u, *b));
        (*constrainer)[ndx::e  ].reset(
                new constraint::reference_bulk(chdef->bulk_rho_E, *b));

    } else if (grid->one_sided()) { // Flat plate

        INFO0(who, "Computing mean freestream behavior per plate scenario");
        using namespace std;
        const real_t T_inf   = isothermal->upper_T;  
        const real_t u_inf   = isothermal->upper_u;
        const real_t rho_inf = isothermal->upper_rho;
        const real_t mx_inf  = u_inf * rho_inf;

        // Compute total energy at infinity based on input upper values
        // Neglect contribution of wall-normal velocity to kinetic energy,
        // and assume that species concentrations are at equilibrium
        // FIXME: adjust the value of target rho_E dynamically
        //        to account for wall-normal velocity and variations in
        //        species mass fractions
        const size_t Ns = cmods->Ns();
        std::vector<real_t> mass_fractions(Ns);
        for (unsigned int s=0; s<Ns; ++s) {
            mass_fractions[s]      = isothermal->upper_cs[s];
        }
        const real_t e_inf   = (cmods->e_from_T(T_inf, mass_fractions)
                                + 0.5 * (u_inf * u_inf)) * rho_inf;

        INFO0(who, "Setting constraints using freestream reference state"
                   " on upper boundary");

        if (chdef->bulk_rho) {
            WARN0(who, "Removing channel-like bulk_rho setting");
            chdef->bulk_rho   = numeric_limits<real_t>::quiet_NaN();
        }
        if (chdef->bulk_rho_u) {
            WARN0(who, "Removing channel-like bulk_rho_u setting");
            chdef->bulk_rho_u = numeric_limits<real_t>::quiet_NaN();
        }
        if (chdef->bulk_rho_E) {
            WARN0(who, "Removing channel-like bulk_rho_E setting");
            chdef->bulk_rho_E = numeric_limits<real_t>::quiet_NaN();
        }

        INFO0(who, "Establishing driving, freestream-like state constraints");
        (*constrainer)[ndx::rho].reset(
                new constraint::constant_upper(rho_inf, *b));
        (*constrainer)[ndx::mx ].reset(
                new constraint::constant_upper(mx_inf,  *b));
        (*constrainer)[ndx::e  ].reset(
                new constraint::constant_upper(e_inf,   *b));

    } else {

        FATAL0(who, "Sanity error in constraint selection");
        return EXIT_FAILURE;

    }


    // Nonreflecting must mutate chdef/isothermal before instantiating
    // linear operator
    if (grid->one_sided()) {
        // FIXME: Plate implementation in progress.
        // Need to set reference value for rho (rho_ref)
        
        INFO0(who, "Preparing nonreflecting upper boundary treatment");

        WARN0(who, "Non-reflecting boundary treatment enabling boundary" 
                   "layer simulations is still under development.      "
                   "Proceed at your own risk.");

        
        // Assign some reference values
        const size_t Ns = cmods->Ns();
        
        // SUZERAIN_ENSURE(!(isnan)(isothermal->upper_rho));
        if ((isnan)(isothermal->upper_rho)) {
            FATAL0(who, "upper_rho for nonreflecting boundary"
                        " not specified and not in restart file");
            return EXIT_FAILURE;
        }
        
        common_block.rho_ref = isothermal->upper_rho;
        common_block.T_ref   = isothermal->upper_T;
        common_block.u_ref   = isothermal->upper_u;
        common_block.v_ref   = isothermal->upper_v;
        common_block.w_ref   = isothermal->upper_w;
        
        common_block.cs_ref.resize(Ns);
        std::vector<real_t> mass_fractions(Ns);
        for (unsigned int s=0; s<Ns; ++s) {
            common_block.cs_ref(s) = isothermal->upper_cs[s];
            mass_fractions[s]      = isothermal->upper_cs[s];
        }
        
        common_block.etots_ref.resize(Ns);
        cmods->evaluate_for_nonreflecting(common_block.T_ref,
                                          common_block.cs_ref,
                                          common_block.a_ref,
                                          common_block.gamma_ref,
                                          common_block.R_ref,
                                          common_block.etots_ref);
        
        common_block.Cv_ref =  common_block.R_ref
                            / (common_block.gamma_ref - 1);
        
        // FIXME: Remove this info or make it a debug output option
        INFO0(who, "rho_ref   = " << common_block.rho_ref  );
        INFO0(who, "T_ref     = " << common_block.T_ref    );
        INFO0(who, "a_ref     = " << common_block.a_ref    );
        INFO0(who, "gamma_ref = " << common_block.gamma_ref);
        INFO0(who, "R_ref     = " << common_block.R_ref    );
        INFO0(who, "E_ref     = " << cmods->e_from_T(common_block.T_ref,
                                                     mass_fractions));
    }

    // Prepare spatial operators depending on requested advance type
    // Notice constrainer always wraps the implicit operator
    // and that the same nonlinear_operator is used pervasively.
    L = constrainer;
    if (use_explicit) {
        INFO0(who, "Initializing fully explicit spatial operators");

        N.reset(new explicit_nonlinear_operator(
                    *cmods, *grid, *dgrid, *cop, *b, common_block, *fsdef,
                    *sgdef, msoln));

        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            shared_ptr<nonreflecting_treatment> nonreflecting(
                    new nonreflecting_treatment(
                        *grid, *dgrid, *cop, *b, common_block));
            nonreflecting->N = N;
            N = nonreflecting;
        }

        constrainer->L.reset(new isothermal_mass_operator(
                    *cmods, *isothermal, 
                    *chdef, *grid, *dgrid, *cop, *b, common_block));
                   
    } else if (use_implicit) {
        INFO0(who, "Initializing hybrid implicit/explicit spatial operators");

        constrainer->L.reset(new isothermal_hybrid_linear_operator(
                    solver_spec, *cmods, *isothermal, *chdef, *grid, *dgrid,
                    *cop, *b, common_block));

        shared_ptr<suzerain::hybrid_residual_operator>
            tmp_hybrid( new hybrid_residual_operator(dgrid->chi()) );

        tmp_hybrid->R.reset(new explicit_nonlinear_operator(
                                *cmods, *grid, *dgrid, *cop, *b,
                                common_block, *fsdef, *sgdef, msoln));

        tmp_hybrid->L = this->L;

        this->N = tmp_hybrid;

        if (grid->one_sided()) {
            WARN0(who, "Non-reflecting boundary treatment with hybrid       " 
                       "implicit/explict time marching is very experimental."
                       "Don't say we didn't warn you.");

            shared_ptr<nonreflecting_treatment> nonreflecting(
                new nonreflecting_treatment(
                    *grid, *dgrid, *cop, *b, common_block));
            nonreflecting->N = this->N;
            this->N = nonreflecting;
        }


    } else {
        FATAL0(who, "Sanity error in operator selection");
        return EXIT_FAILURE;
    }

#ifdef SUZERAIN_HAVE_LARGO
    // Allocate slow growth workspace
    {
        int state_count = (int) cmods->Ns()+4;
        int Ns = (int) cmods->Ns()-1; 
        largo_bl_temporal_allocate(&sgdef->workspace, state_count, Ns);
    }
#endif

    // Use --undriven as a testing- and debugging-related tool.
    // For example, to investigate nonreflecting boundary condition behavior.
    if (options.variables().count("undriven")) {
        boost::algorithm::trim(undriven);
        shared_ptr<constraint::base> disabled(new constraint::disabled(*b));
        if (undriven == "all") {
            INFO0(who, "Disabling all established constraints per --undriven");
            (*constrainer)[ndx::e  ] = disabled;
            (*constrainer)[ndx::mx ] = disabled;
            (*constrainer)[ndx::my ] = disabled;
            (*constrainer)[ndx::mz ] = disabled;
            (*constrainer)[ndx::rho] = disabled;
        } else if (undriven == "rho") {
            INFO0(who, "Disabling any density constraint per --undriven");
            (*constrainer)[ndx::rho] = disabled;
        } else if (undriven == "rho_u") {
            INFO0(who, "Disabling any momentum constraint per --undriven");
            (*constrainer)[ndx::mx ] = disabled;
        } else if (undriven == "rho_E") {
            INFO0(who, "Disabling any total energy constraint per --undriven");
            (*constrainer)[ndx::e  ] = disabled;
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

    // Report error to the OS iff advance_control reported an error
    return elapsed_wall_time < 0;
}
