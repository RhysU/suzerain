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

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <esio/esio.h>

#include <suzerain/common.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/noise_definition.hpp>
#include <suzerain/zgbsv_specification.hpp>

#include "driver.hpp"
#include "channel_treatment.hpp"
#include "explicit_operator.hpp"
#include "hybrid_operator.hpp"

#pragma warning(disable:1419)

namespace suzerain { namespace perfect {

/** Application for initializing new restart files. */
struct driver_advance : public driver
{
    driver_advance(const std::string& revstr)
        : driver("Compressible, perfect gas simulation time advancement",
                 "RESTART-FILE",
                 "",
                 revstr)
    {}

    /** Implementation below in this file */
    int run(int argc, char **argv);
};

} /* namespace perfect */ } /* namespace suzerain */

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
    using boost::math::isnan;
    using std::numeric_limits;
    using std::string;
    using std::vector;

    // Storage for binary-specific options
    const support::noise_definition noisedef;
    bool use_explicit  = false;
    bool use_implicit  = false;
    string solver_spec(static_cast<string>(suzerain::zgbsv_specification()));

    // Register binary-specific options
    options.add_definition(const_cast<support::noise_definition&>(noisedef));
    options.add_options()
        ("explicit", "Use purely explicit operators")
        ("implicit", "Use hybrid implicit/explicit operators")
        ("solver",   boost::program_options::value<string>(&solver_spec)
                         ->default_value(solver_spec),
                     "Use the specified algorithm for any implicit solves")
    ;

    // Initialize application and then process binary-specific options
    // (henceforth suzerain::support::logging macros becomes usable)
    vector<string> positional = initialize(argc, argv);

    // Select type of timestepping operators to use (default implicit)
    options.conflicting_options("implicit", "explicit");
    if (options.variables().count("explicit")) {
        use_explicit = true;
    } else {
        use_implicit = true;
    }

    if (positional.size() != 1) {
        FATAL0("Exactly one restart file name must be specified");
        return EXIT_FAILURE;
    }
    const std::string restart_file = positional[0];

    INFO0("Loading restart file: " << restart_file);
    real_t initial_t = numeric_limits<real_t>::quiet_NaN();
    {
        // Obtain exact restart scenario via "push/pop" pair below
        shared_ptr<scenario_definition> restart_scenario
                = make_shared<scenario_definition>();

        // Load the restart details with state going into state_linear
        esio_handle esioh = esio_handle_initialize(MPI_COMM_WORLD);
        esio_file_open(esioh, restart_file.c_str(), 0 /* read-only */);
        restart_scenario.swap(scenario);  // "push"
        load_restart(esioh, initial_t);
        restart_scenario.swap(scenario);  // "pop"
        esio_file_close(esioh);
        esio_handle_finalize(esioh);

        // If necessary, adjust total energy to account for scenario changes
        bool necessary = false;
        if (    (isnan)(scenario->Ma)
             || scenario->Ma != restart_scenario->Ma) {
            necessary = true;
            scenario->Ma = restart_scenario->Ma;
        }
        if (    (isnan)(scenario->gamma)
            || scenario->gamma != restart_scenario->gamma) {
            necessary = true;
            scenario->gamma = restart_scenario->gamma;
        }
        if (necessary) {
            state_nonlinear->assign(*state_linear);
            adjust_scenario(*state_nonlinear, *scenario, *grid,
                            *dgrid, *cop, restart_scenario->Ma,
                            restart_scenario->gamma);
            state_linear->assign(*state_nonlinear);
        }
    }

    if (msoln) {
        INFO0("Restart file prescribes a manufactured solution");
        if (!(isnan)(scenario->bulk_rho)) {
            WARN0("Manufactured solution incompatible with bulk_rho = "
                  << scenario->bulk_rho);
        }
        if (!(isnan)(scenario->bulk_rho_u)) {
            WARN0("Manufactured solution incompatible with bulk_rho_u = "
                  << scenario->bulk_rho_u);
        }
    }

    // If requested, add noise to the momentum fields at startup (expensive).
    if (noisedef.percent > 0) {
        state_nonlinear->assign(*state_linear);
        add_noise(*state_nonlinear, noisedef, *scenario,
                  *grid, *dgrid, *cop, *b);
        state_linear->assign(*state_nonlinear);
    }

    // Prepare spatial operators depending on request advance type
    if (use_explicit) {
        INFO0("Initializing explicit timestepping operators");
        L.reset(new channel_treatment<isothermal_mass_operator>(
                    *scenario, *grid, *dgrid, *cop, *b, common_block));
        N.reset(new explicit_nonlinear_operator(
                    *scenario, *grid, *dgrid, *cop, *b, common_block, msoln));
    } else if (use_implicit) {
        INFO0("Initializing hybrid implicit/explicit timestepping operators");
        L.reset(new channel_treatment<isothermal_hybrid_linear_operator>(
                    solver_spec, *scenario, *grid, *dgrid,
                    *cop, *b, common_block));
        N.reset(new hybrid_nonlinear_operator(
                    *scenario, *grid, *dgrid, *cop, *b, common_block, msoln));
    } else {
        FATAL0("Sanity error in operator selection");
        return EXIT_FAILURE;
    }

    // Perform final housekeeping and then advance time as requested
    establish_ieee_mode();
    log_discretization_quality();
    prepare_controller(initial_t);
    save_metadata();
    const real_t elapsed_wall_time = advance_controller(); // Negative on error

    // If we advanced by any time steps, log the observed linearization error
    if (elapsed_wall_time >= 0 && controller->current_nt() > 0) {
        log_linearization_error(build_timeprefix(controller->current_t(),
                                                 controller->current_nt()),
                                controller->current_t(),
                                controller->current_nt());
    }

    return elapsed_wall_time >= 0;
}
