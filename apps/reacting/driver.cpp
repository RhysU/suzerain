//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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
 * @copydoc driver.hpp
 */

#include "driver.hpp"

#include <suzerain/diffwave.hpp>
#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

#include "reacting.hpp"

namespace suzerain {

namespace reacting {

driver::driver(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &revstr)
    : driver_base(application_synopsis,
                  argument_synopsis,
                  description,
                  revstr)
    , isothermal(make_shared<support::definition_isothermal>())
    , chdef(make_shared<definition_channel>())
    , cmods(make_shared<antioch_constitutive>())
    , fsdef(make_shared<definition_filter>())
    , sgdef(make_shared<support::definition_largo>())
    , who("reacting")
{
    // Sets up usual 5 fields.  If necessary, species are added later.
    this->fields = default_fields();
}

std::vector<std::string>
driver::initialize(
        int argc,
        char **argv)
{
    // msoln is not used by all binaries and is therefore not added below
    options.add_definition(*chdef);
    options.add_definition(*cmods);
    options.add_definition(*fsdef);
    options.add_definition(*isothermal);
    options.add_definition(*sgdef);

    // Delegate to superclass initialization
    std::vector<std::string> positional = super::initialize(argc, argv);

    return positional;
}

void
driver::log_manufactured_solution_absolute_error(
        const std::string& timeprefix,
        const time_type t,
        const step_type nt)
{
    SUZERAIN_UNUSED(nt);

    // NOP whenever a manufactured_solution is not in use
    if (!msoln) return;

    // Avoid computational cost when logging is disabled
    support::logging::logger_type mms_abserr
            = support::logging::get_logger("mms.abserr");
    if (!INFO0_ENABLED(mms_abserr)) return;

    // Compute $L^2_{xyz}$ of error of state against manufactured solution
    state_nonlinear->assign_from(*state_linear);
    accumulate_manufactured_solution(
            1, *msoln, -1, *state_nonlinear,
            *grid, *dgrid, *cop, *b, t);
    const std::vector<field_L2xyz> L2
        = compute_field_L2xyz(*state_nonlinear, *grid, *dgrid, *gop);

    // Output absolute global errors for each field
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(L2[k].total());
    }
    INFO0(mms_abserr, msg.str());
}

void
driver::compute_statistics(
        driver::time_type t)
{
    SUZERAIN_TIMER_SCOPED("driver::compute_statistics");

    // Obtain mean samples from instantaneous fields stored in state_linear
    state_nonlinear->assign_from(*state_linear);
    mean = reacting::sample_quantities(*cmods, *grid, *dgrid,
                                       *cop, *state_nonlinear, t);

    // Obtain mean quantities computed via implicit forcing (when possible)
    if (common_block.means.rows() == mean.storage.rows()) {
        mean.f().col(0)     = common_block.fx();
        mean.f().col(1)     = common_block.fy();
        mean.f().col(2)     = common_block.fz();
        mean.f_dot_u()      = common_block.f_dot_u();
        mean.qb()           = common_block.qb();
        mean.CrhoE()        = common_block.CrhoE();
        mean.Crhou().col(0) = common_block.Crhou();
        mean.Crhou().col(1) = common_block.Crhov();
        mean.Crhou().col(2) = common_block.Crhow();
        mean.Crho()         = common_block.Crho();
        mean.Crhou_dot_u()  = common_block.Crhou_dot_u();
    } else {
        WARN0(who, "Could not obtain mean quantities set by implicit forcing");
        mean.storage.middleCols<quantities::nscalars::implicit>(
                quantities::start::implicit).setConstant(
                    std::numeric_limits<real_t>::quiet_NaN());
    }

    // Convert implicit values at collocation points to expansion coefficients
    bsplineop_lu mass(*cop);
    mass.opform_mass(*cop);
    mass.factor();
    mass.solve(quantities::nscalars::implicit,
               mean.storage.middleCols<quantities::nscalars::implicit>(
                   quantities::start::implicit).data(),
               mean.storage.innerStride(), mean.storage.outerStride());
}

bool
driver::log_status_hook(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    const bool retval = super::log_status_hook(timeprefix, t, nt);
    log_manufactured_solution_absolute_error(timeprefix, t, nt);
    return retval;
}

void
driver::save_metadata_hook(
        const esio_handle esioh)
{
    super::save_metadata_hook(esioh);
    chdef->save(esioh);
    cmods->save(esioh);
    fsdef->save(esioh);
    fsdef->save_filteropz(esioh, dgrid->global_wave_extent.y());
    isothermal->save(esioh);
    sgdef->save(esioh);
    save(esioh, msoln, *cmods, *grid);
    return;
}

void
driver::load_metadata_hook(
        const esio_handle esioh)
{
    super::load_metadata_hook(esioh);
    chdef->load(esioh);
    cmods->load(esioh);
    fsdef->load(esioh);
    isothermal->load(esioh);
    sgdef->load(esioh);

    // After cmods->load, have valid species names, so we can update
    // fields to include species
    add_species_fields(cmods->species_names, this->fields);

    load(esioh, msoln, *cmods, *grid);
    return;
}

bool
driver::save_statistics_hook(
        const esio_handle esioh,
        const driver::time_type t)
{
    DEBUG0(who, "In driver::save_statistics_hook");

    // Should we compute fresh statistics or re-use cached values?
#pragma warning(push,disable:1572)
    const bool use_cached = controller && (t == mean.t);
#pragma warning(pop)

    // Compute a prefix for logging purposes, including a trailing space
    std::string prefix;
    if (controller) {
        prefix = build_timeprefix(t, controller->current_nt());
        prefix.append(1, ' ');
    }

    // Compute statistics whenever necessary
    if (use_cached) {
        DEBUG0(who, prefix << "Cowardly refusing to re-compute statistics");
    } else {
        const double starttime = MPI_Wtime();
        compute_statistics(t);
        const double elapsed = MPI_Wtime() - starttime;
        INFO0(who, prefix << "Computed statistics in "
              << elapsed << " seconds");
    }

    // Save statistics and invoke superclass hook
    if (!mean.save(esioh)) {
        WARN0(who, "Incomplete statistics saved at time " << t);
    }
    return super::save_statistics_hook(esioh, t);
}

void
driver::load_statistics_hook(
        const esio_handle esioh)
{
    support::load_time(esioh, mean.t);
    if (!mean.load(esioh)) {
        // Notice this forces use_cached == false in save_statistics_hook
        WARN0(who, "Incomplete statistics loaded from time " << mean.t);
        mean.t = std::numeric_limits<real_t>::quiet_NaN();
    }
    return super::load_statistics_hook(esioh);
}

// Please review and possibly synchronize apps/perfect/driver.cpp when changing
void
driver::default_restart_interval(
        time_type& t,
        step_type&)
{
    using std::max;
    using std::abs;

    // Look for the largest magnitude, problem-dependent, macro velocity scale
    real_t velocity = std::numeric_limits<real_t>::quiet_NaN();

    // In a channel...
    if (grid && grid->two_sided()) {

        // ...any approximate bulk velocity from a driving force...
        if (chdef) {

            if ((boost::math::isnan)(chdef->bulk_rho)) {
                TRACE0(who, "No bulk density scale available so assuming 1");
                velocity = max(velocity,
                               abs(chdef->bulk_rho_u) / /*rho*/ 1);
            } else {
                velocity = max(velocity,
                               abs(chdef->bulk_rho_u) / chdef->bulk_rho);
            }

        }

        // ...may be trumped by driving the upper and lower walls.
        if (isothermal) {
            velocity = max(velocity, abs(isothermal->upper_u));
            velocity = max(velocity, abs(isothermal->lower_u));
        }

    }

    // On a plate...
    if (grid && grid->one_sided()) {

        if (sgdef->formulation.enabled() && sgdef->baseflow && fields.size()) {

            // ...prefer wall velocity from inviscid baseflow as freestream
            // surrogate as it is independent of the wall-normal domain size...
            const size_t nvalues = fields.size() + /*pressure*/1;
            std::vector<real_t> buf(3 * nvalues);
            sgdef->baseflow->conserved(0.0,               // wall
                                       &buf[0*nvalues],   // base
                                       &buf[1*nvalues],   // dybase
                                       &buf[2*nvalues]);  // dxbase
            const real_t base_rho   = buf[0];
            const real_t base_rho_u = buf[1];
            if (base_rho > 0) {
                velocity = max(velocity, abs(base_rho_u) / base_rho);
            } else {
                TRACE0(who, "Trivial baseflow detected; using upper velocity");
                velocity = max(velocity, abs(isothermal->upper_u));
            }

        } else if (isothermal) {
            // ...taking freestream reference if-and-only-if no baseflow...
            velocity = max(velocity, abs(isothermal->upper_u));
        }

        // ...and permit a driven lower wall velocity to trump.
        if (isothermal) {
            velocity = max(velocity, abs(isothermal->lower_u));
        }

    }

    // If we have a domain size and a velocity scale, compute a timescale
    // and default to saving eight restarts across that timescale.
    // (e.g. in a channel resulting in 80 restarts in 10 flow throughs)
    if (boost::math::isnormal(velocity)) {
        if (grid) {
            t = (grid->L.x() / velocity) / 8;
        } else {
            DEBUG0(who, "No grid details for default_restart_interval");
        }
    } else {
        DEBUG0(who, "No macro scale velocity for default_restart_interval");
    }
}

} // end namespace reacting

} // end namespace suzerain
