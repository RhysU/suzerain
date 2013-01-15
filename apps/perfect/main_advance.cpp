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
    using std::string;
    using std::vector;

    // Storage for binary-specific options
    const support::noise_definition noisedef;
    bool use_explicit  = false;
    bool use_implicit  = false;
    bool use_yang11    = false;
    bool use_smr91     = false;
    string solver_spec(static_cast<string>(suzerain::zgbsv_specification()));

    // Register binary-specific options
    options.add_definition(const_cast<support::noise_definition&>(noisedef));
    options.add_options()
        ("explicit", "Use purely explicit operators")
        ("implicit", "Use hybrid implicit/explicit operators")
        ("smr91",    "Advance time per Spalart, Moser, and Rogers 1991")
        ("yang11",   "Advance time per Shan Yang's 2011 thesis")
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

    // Select type of timestepping schemes to use (default smr91)
    options.conflicting_options("smr91", "yang11");
    if (options.variables().count("yang11")) {
        use_yang11 = true;
    } else {
        use_smr91 = true;
    }

    if (positional.size() != 1) {
        FATAL0("Exactly one restart file name must be specified");
        return EXIT_FAILURE;
    }
    const std::string restart_file = positional[0];

    // FIXME STARTHERE

    DEBUG0("Establishing runtime parallel infrastructure and resources");
    establish_ieee_mode();
    load_grid_and_operators(NULL);
    establish_decomposition();
    establish_state_storage(fields.size(), fields.size());

    return EXIT_SUCCESS;
}
