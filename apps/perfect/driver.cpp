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

#include <suzerain/bl.h>
#include <suzerain/channel.h>
#include <suzerain/diffwave.hpp>
#include <suzerain/error.h>
#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_tools.hpp>
#include <suzerain/profile.hpp>
#include <suzerain/samples.hpp>
#include <suzerain/state.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

#include "perfect.hpp"

namespace suzerain {

namespace perfect {

driver::driver(
        const std::string &application_synopsis,
        const std::string &argument_synopsis,
        const std::string &description,
        const std::string &revstr)
    : driver_base(application_synopsis,
                  argument_synopsis,
                  description,
                  revstr)
    , scenario(make_shared<definition_scenario>())
    , isothermal(make_shared<support::definition_isothermal>())
    , sg(make_shared<support::definition_largo>())
    , rad(make_shared<support::definition_radialflow>(/*deltae*/ 1.0))
    , who("perfect")
{
    this->fields = default_fields();
}

std::vector<std::string>
driver::initialize(
        int argc,
        char **argv)
{
    // Only add groups of options when non-trivial at initialization
    if (scenario  ) options.add_definition(*scenario  ) ;
    if (isothermal) options.add_definition(*isothermal) ;
    if (sg        ) options.add_definition(*sg        ) ;
    if (rad       ) options.add_definition(*rad       ) ;
    if (msoln     ) options.add_definition(*msoln     ) ;

    // Delegate to superclass initialization
    std::vector<std::string> positional = super::initialize(argc, argv);

    // If msoln was provided, match its contents to other members
    if (msoln) {
        if (scenario) msoln->match(*scenario);
        if (grid)     msoln->match(*grid);
    }

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
driver::log_linearization_error(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    SUZERAIN_UNUSED(nt);

    // Avoid computational cost when logging is disabled
    support::logging::logger_type lin_abserr
            = support::logging::get_logger("lin.abserr");
    if (!INFO0_ENABLED(lin_abserr)) return;

    SUZERAIN_TIMER_SCOPED("driver::log_linearization_error");

    // Algorithm proceeds directly from necessary interface guarantees.
    // Rewriting it symbolically from the code and grokking it is a good
    // exercise towards understanding the low storage time interfaces.
    //
    // Success requires that the following preconditions all hold:
    //    i) At substep zero, this->N computes N(u) - Lu including
    //       the necessary reference quantities.
    //   ii) The reference quantities are stored in common_block and may be
    //       expressly zeroed.  Zeroing them nukes linearized contributions
    //       from I{Nonlinear,Linear}Operator but not linear contributions
    //       from this->L.
    //  iii) On non-zero substeps for zero reference values, this->N computes
    //       only N(u) and this->L computes only
    //       \left(M+\varphi\mathscr{L}\right)u.
    //   iv) At the beginning of this process, state_linear contains
    //       valid information and adheres to boundary conditions.
    //    v) The explicit_operator_nonlinear application requires an auxiliary
    //       scaling factor "chi" to account for Fourier transform
    //       normalization needs.
    //   vi) Using diffwave::apply zeros wavenumbers used only for dealiasing
    //       purposes-- this data is unimportant from the perspective of
    //       measuring actual linearization error.
    state_nonlinear->assign_from(*state_linear);
    common_block.set_zero(grid->dN.y());  // Defensive
    N->apply_operator(t, *state_nonlinear, *method, /*substep*/0);
    L->accumulate_mass_plus_scaled_operator(
            1, *state_linear, dgrid->chi(), *state_nonlinear, /*substep*/0);
    common_block.set_zero(grid->dN.y());  // Zero reference quantities
    L->accumulate_mass_plus_scaled_operator(
            1, *state_linear, -1, *state_nonlinear, /*substep*/0);
    for (size_t k = 0; k < fields.size(); ++k) {
        diffwave::apply(0, 0, 1, (*state_nonlinear)[k].origin(),
            grid->L.x(), grid->L.z(), dgrid->global_wave_extent.y(),
            grid->N.x(), grid->dN.x(),
            dgrid->local_wave_start.x(), dgrid->local_wave_end.x(),
            grid->N.z(), grid->dN.z(),
            dgrid->local_wave_start.z(), dgrid->local_wave_end.z());
    }
    state_nonlinear->exchange(*state_linear);
    N->apply_operator(t, *state_nonlinear, *method, /*substep*/1);
    for (size_t k = 0; k < fields.size(); ++k) {
        diffwave::apply(0, 0, 1, (*state_nonlinear)[k].origin(),
            grid->L.x(), grid->L.z(), dgrid->global_wave_extent.y(),
            grid->N.x(), grid->dN.x(),
            dgrid->local_wave_start.x(), dgrid->local_wave_end.x(),
            grid->N.z(), grid->dN.z(),
            dgrid->local_wave_start.z(), dgrid->local_wave_end.z());
    }
    state_nonlinear->add_scaled(1/dgrid->chi(), *state_linear);

    // When this->L and this->N match perfectly, state_nonlinear should be
    // identically zero.  Seeing more than acceptable (e.g. 1e-10) floating
    // point error indicates something is amiss.
    const std::vector<field_L2xyz> L2
        = compute_field_L2xyz(*state_nonlinear, *grid, *dgrid, *gop);
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(L2[k].total());
    }
    INFO0(lin_abserr, msg.str());
}

void driver::log_quantities_of_interest(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    SUZERAIN_UNUSED(t);
    SUZERAIN_UNUSED(nt);

    if (!grid) return;

    SUZERAIN_TIMER_SCOPED("driver::log_quantities_of_interest");

    // When possible, any operator_tools superclass of N is reused so that
    // take_samples may benefit from any cached factorizations.
    shared_ptr<operator_tools> otool
            = dynamic_pointer_cast<operator_tools>(N);
    if (!otool) {
        otool = make_shared<operator_tools>(*grid, *dgrid, *cop);
    }

    // If possible, use existing information from mean quantities
    // Otherwise compute from instantaneous fields stored in state_linear
    profile prof;
    if (controller && mean && mean->t == t) {
        prof = *mean;
    } else {
        state_nonlinear->assign_from(*state_linear);
        prof = *take_profile(*scenario, *otool, *state_nonlinear);
    }

    if (grid->one_sided()) {

        // Compute many details about the boundary layer for logging
        suzerain_bl_local       wall;
        suzerain_bl_viscous     viscous;
        suzerain_bl_local       edge;
        suzerain_bl_local       edge99;
        suzerain_bl_thicknesses thick;
        suzerain_bl_reynolds    reynolds;
        suzerain_bl_qoi         qoi;
        suzerain_bl_pg          pg;
        summarize_boundary_layer_nature(prof, *scenario, sg,
                                        *otool->masslu(), *b,
                                        wall, viscous, thick,
                                        edge, edge99,
                                        reynolds, qoi, pg);

        // Log messages using application-agnostic superclass functionality
        this->log_quantities_boundary_layer(timeprefix,
                                            &wall, &viscous, &thick,
                                            &edge, &edge99,
                                            &reynolds, &qoi, &pg);

    } else if (grid->two_sided()) {

        // Compute many details about the boundary layer for logging
        suzerain_channel_local   wall;
        suzerain_channel_viscous viscous;
        suzerain_channel_local   center;
        suzerain_channel_qoi     qoi;
        summarize_channel_nature(prof, *scenario, *b,
                                 wall, viscous, center, qoi);

        // Log messages using application-agnostic superclass functionality
        this->log_quantities_channel(timeprefix,
                                     &wall, &viscous, &center, &qoi);

    } else {

        SUZERAIN_ERROR_REPORT_UNIMPLEMENTED();

    }
}

void
driver::compute_statistics(
        driver::time_type t)
{
    SUZERAIN_TIMER_SCOPED("driver::compute_statistics");

    // Obtain mean samples from instantaneous fields stored in state_linear
    state_nonlinear->assign_from(*state_linear);
    // When possible, any operator_tools superclass of N is reused so that
    // take_samples may benefit from any cached factorizations.
    shared_ptr<operator_tools> otool
            = dynamic_pointer_cast<operator_tools>(N);
    if (!otool) {
        otool = make_shared<operator_tools>(*grid, *dgrid, *cop);
    }
    mean = take_samples(*scenario, *otool, *state_nonlinear, t);

    // Obtain mean quantities computed via implicit forcing (when possible)
    if (common_block.implicits.rows() == mean->storage.rows()) {
        mean->SrhoE()        = common_block.SrhoE();
        mean->Srhou().col(0) = common_block.Srhou();
        mean->Srhou().col(1) = common_block.Srhov();
        mean->Srhou().col(2) = common_block.Srhow();
        mean->Srho()         = common_block.Srho();
        mean->Srhou_dot_u()  = common_block.Srhou_dot_u();
        mean->f().col(0)     = common_block.fx();
        mean->f().col(1)     = common_block.fy();
        mean->f().col(2)     = common_block.fz();
        mean->f_dot_u()      = common_block.f_dot_u();
        mean->qb()           = common_block.qb();
        mean->CrhoE()        = common_block.CrhoE();
        mean->Crhou().col(0) = common_block.Crhou();
        mean->Crhou().col(1) = common_block.Crhov();
        mean->Crhou().col(2) = common_block.Crhow();
        mean->Crho()         = common_block.Crho();
        mean->Crhou_dot_u()  = common_block.Crhou_dot_u();
    } else {
        WARN0(who, "Could not obtain mean quantities set by implicit forcing");
        common_block.implicits.setConstant(
                std::numeric_limits<real_t>::quiet_NaN());
    }

    // Convert implicit values at collocation points to expansion coefficients
    bsplineop_lu mass(*cop);
    mass.opform_mass(*cop);
    mass.factor();
    mass.solve(samples::nscalars::implicit,
               mean->storage.middleCols<samples::nscalars::implicit>(
                   samples::start::implicit).data(),
               mean->storage.innerStride(), mean->storage.outerStride());
}

bool
driver::log_status_hook(
        const std::string& timeprefix,
        const real_t t,
        const std::size_t nt)
{
    const bool retval = super::log_status_hook(timeprefix, t, nt);
    log_manufactured_solution_absolute_error(timeprefix, t, nt);
    log_quantities_of_interest(timeprefix, t, nt);
    return retval;
}

void
driver::save_metadata_hook(
        const esio_handle esioh)
{
    super::save_metadata_hook(esioh);
    if (scenario)   scenario->save(esioh);
    if (isothermal) isothermal->save(esioh);
    if (sg)         sg->save(esioh);
    if (rad)        rad->save(esioh);
    if (msoln)      save(esioh, msoln, *scenario, *grid);
    return;
}

void
driver::load_metadata_hook(
        const esio_handle esioh)
{
    super::load_metadata_hook(esioh);

    if (!scenario) {
        scenario = make_shared<definition_scenario>();
    }
    scenario->load(esioh);

    if (!isothermal) {
        isothermal = make_shared<support::definition_isothermal>();
    }
    isothermal->load(esioh);

    if (!sg) {
        sg = make_shared<support::definition_largo>();
    }
    sg->load(esioh);

    if (!rad) {
        rad = make_shared<support::definition_radialflow>();
    }
    rad->load(esioh);

    load(esioh, msoln, *scenario, *grid);
    return;
}

bool
driver::save_statistics_hook(
        const esio_handle esioh,
        const driver::time_type t)
{
    // Should we compute fresh statistics or re-use cached values?
#pragma warning(push,disable:1572)
    const bool use_cached = controller && mean && (mean->t == t);
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
    // FIXME Also save boundary layer quantities of interest when one-sided?
    if (!support::save_samples(esioh, *mean)) {
        WARN0(who, "Incomplete statistics saved at time " << t);
    }
    return super::save_statistics_hook(esioh, t);
}

void
driver::load_statistics_hook(
        const esio_handle esioh)
{
    if (!mean) {
        mean.reset(new samples());
    }
    support::load_time(esioh, mean->t);

    if (!support::load_samples(esioh, *mean)) {
        // Notice this forces use_cached == false in save_statistics_hook
        WARN0(who, "Incomplete statistics loaded from time " << mean->t);
        mean->t = std::numeric_limits<real_t>::quiet_NaN();
    }
    return super::load_statistics_hook(esioh);
}

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
        if (scenario) {

            if ((boost::math::isnan)(scenario->bulk_rho)) {
                TRACE0(who, "No bulk density scale available so assuming 1");
                velocity = max(velocity,
                               abs(scenario->bulk_rho_u) / /*rho*/ 1);
            } else {
                velocity = max(velocity,
                               abs(scenario->bulk_rho_u) / scenario->bulk_rho);
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

        // ...prefer wall velocity from the inviscid baseflow as a surrogate of
        // freestream as it is independent of the wall-normal domain size...
        if (sg->formulation.enabled() && sg->baseflow) {
            largo_state freestream, dontcare;
            sg->baseflow->conserved(0.0, freestream.as_is(),
                                    dontcare.as_is(), dontcare.as_is());
            velocity = max(velocity, abs(freestream.u()));
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

} // end namespace perfect

} // end namespace suzerain
