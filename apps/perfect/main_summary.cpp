//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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
 * Application executing \ref suzerain::driver_base::summary_run
 */

#include "driver.hpp"

#include <suzerain/bspline.hpp>
#include <suzerain/profile.hpp>
#include <suzerain/summary.hpp>
#include <suzerain/support/driver_base.hpp>
#include <suzerain/support/support.hpp>

#include "instantaneous.hpp"

#pragma warning(disable:1419)

namespace suzerain {

namespace perfect {

/** Summarize mean statistics for one or more restart files. */
struct driver_summary : public driver
{
    driver_summary()
        : driver("Compressible, perfect gas simulation summarization",
                 driver_base::summary_argument_synopsis,
                 driver_base::summary_description,
                 REVISIONSTR)
    {}

    /** Logging requirements are simpler than what superclass provides. */
    virtual std::string log4cxx_config() { return support::log4cxx_config; }

    /** Invoked by \c main just below. */
    int run(int argc, char **argv)
    {
        // Invoke driver-agnostic summarization logic
        summary_pool_type pool;
        summary final;
        support::shared_esio_handle esioh;
        const int status = summary_run(argc, argv, pool, final, esioh);

        // Perform driver-specific post-processing on nontrivial results
        if (pool.size()) {

            // Convert final from collocation point values to coefficients
            bsplineop_lu masslu(*cop);
            masslu.factor_mass(*cop);
            masslu.solve(final.nongrid().cols(),
                         final.nongrid().data(),
                         final.nongrid().innerStride(),
                         final.nongrid().outerStride());

            // Compute QoI given final profile information
            profile prof;
            prof = final;
            instantaneous inst;
            inst = final;
            log_quantities_of_interest(
                    build_timeprefix_description("summary", "summary"),
                    prof, inst, esioh);
        }

        return status;
    }

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;
};

} // namespace perfect

} // namespace suzerain

/** Instantiate and invoke the application. */
int main(int argc, char **argv)
{
    suzerain::perfect::driver_summary app;
    return app.run(argc, argv);
}
