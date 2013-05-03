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
#include <suzerain/support/logging.hpp>
#include <suzerain/support/noise_definition.hpp>
#include <suzerain/zgbsv_specification.hpp>
#include <suzerain/hybrid_residual_operator.hpp>

#include "driver.hpp"
#include "channel_treatment.hpp"
#include "explicit_operator.hpp"

#include "hybrid_operator.hpp"
#include "nonreflecting_treatment.hpp"


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
                     "Use the specified type to construct filter source");
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

    // Prepare spatial operators depending on requested advance type
    if (use_explicit) {
        INFO0(who, "Initializing fully explicit spatial operators");
        L.reset(new channel_treatment<isothermal_mass_operator>(
                    *cmods, *isothermal, *chdef, *grid, *dgrid, *cop, *b, common_block));
        N.reset(new explicit_nonlinear_operator(
                    *cmods, *grid, *dgrid, *cop, *b, common_block, *fsdef, msoln));

        // nonreflecting
        if (grid->one_sided()) {
            INFO0(who, "Preparing nonreflecting upper boundary treatment");
            shared_ptr<nonreflecting_treatment> nonreflecting(
                    new nonreflecting_treatment(
                        *grid, *dgrid, *cop, *b, common_block));
            nonreflecting->N = N;
            N = nonreflecting;
            chdef->bulk_rho      = numeric_limits<real_t>::quiet_NaN();
            chdef->bulk_rho_u    = numeric_limits<real_t>::quiet_NaN();
            isothermal->upper_T  = numeric_limits<real_t>::quiet_NaN();
            isothermal->upper_u  = numeric_limits<real_t>::quiet_NaN();
            isothermal->upper_v  = numeric_limits<real_t>::quiet_NaN();
            isothermal->upper_w  = numeric_limits<real_t>::quiet_NaN();
	    // TODO Consider appropriate BC for species in this context
        }

    } else if (use_implicit) {
        INFO0(who, "Initializing hybrid implicit/explicit spatial operators");

        // nonreflecting
        if (grid->one_sided()) {
            FATAL0(who, "Nonreflecting upper boundary treatment"
                        " not usable with implicit advance");
            return EXIT_FAILURE;
        }


        L.reset(new channel_treatment<isothermal_hybrid_linear_operator>(
                    solver_spec, *cmods, *isothermal, *chdef, *grid, *dgrid,
                    *cop, *b, common_block));

        shared_ptr<suzerain::hybrid_residual_operator>
            tmp_hybrid( new hybrid_residual_operator(dgrid->chi()) );

        tmp_hybrid->R.reset(new explicit_nonlinear_operator(
                                *cmods, *grid, *dgrid, *cop, *b,
                                common_block, *fsdef, msoln));

        tmp_hybrid->L = this->L;

        this->N = tmp_hybrid;
    } else {
        FATAL0(who, "Sanity error in operator selection");
        return EXIT_FAILURE;
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
