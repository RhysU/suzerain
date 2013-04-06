//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Rhys Ulerich
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
    , chdef(make_shared<channel_definition>())
    , cmods(make_shared<antioch_constitutive>())
    , fsdef(make_shared<filter_definition>())
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

    // Compute L2 of error of state against manufactured solution
    state_nonlinear->assign_from(*state_linear);
    accumulate_manufactured_solution(
            1, *msoln, -1, *state_nonlinear,
            *grid, *dgrid, *cop, *b, t);
    const std::vector<field_L2> L2
        = compute_field_L2(*state_nonlinear, *grid, *dgrid, *gop);

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
    //    v) The explicit_nonlinear_operator application requires an auxiliary
    //       scaling factor "chi" to account for Fourier transform
    //       normalization needs.
    //   vi) Using diffwave::apply zeros wavenumbers used only for dealiasing
    //       purposes-- this data is unimportant from the perspective of
    //       measuring actual linearization error.
    const real_t chi = real_t(1) / (grid->dN.x() * grid->dN.z());
    state_nonlinear->assign_from(*state_linear);
    common_block.set_zero(grid->dN.y());  // Defensive
    N->apply_operator(t, *state_nonlinear,
            method->evmaxmag_real(), method->evmaxmag_imag(), /*substep*/0);
    L->accumulate_mass_plus_scaled_operator(
            1, *state_linear, chi,  *state_nonlinear, /*substep*/0);
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
    N->apply_operator(t, *state_nonlinear,
            method->evmaxmag_real(), method->evmaxmag_imag(), /*substep*/1);
    for (size_t k = 0; k < fields.size(); ++k) {
        diffwave::apply(0, 0, 1, (*state_nonlinear)[k].origin(),
            grid->L.x(), grid->L.z(), dgrid->global_wave_extent.y(),
            grid->N.x(), grid->dN.x(),
            dgrid->local_wave_start.x(), dgrid->local_wave_end.x(),
            grid->N.z(), grid->dN.z(),
            dgrid->local_wave_start.z(), dgrid->local_wave_end.z());
    }
    state_nonlinear->add_scaled(1/chi, *state_linear);

    // When this->L and this->N match perfectly, state_nonlinear should be
    // identically zero.  Seeing more than acceptable (e.g. 1e-10) floating
    // point error indicates something is amiss.
    const std::vector<field_L2> L2
        = compute_field_L2(*state_nonlinear, *grid, *dgrid, *gop);
    std::ostringstream msg;
    msg << timeprefix;
    for (size_t k = 0; k < L2.size(); ++k) {
        msg << ' ' << fullprec<>(L2[k].total());
    }
    INFO0(lin_abserr, msg.str());
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
        mean.f().col(0) = common_block.f();  // Only streamwise momentum...
        mean.f().rightCols<2>().setZero();   // ...not wall-normal, spanwise
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
        using std::numeric_limits;
        mean.f          ().setConstant(numeric_limits<real_t>::quiet_NaN());
        mean.f_dot_u    ().setConstant(numeric_limits<real_t>::quiet_NaN());
        mean.qb         ().setConstant(numeric_limits<real_t>::quiet_NaN());
        mean.CrhoE      ().setConstant(numeric_limits<real_t>::quiet_NaN());
        mean.Crhou      ().setConstant(numeric_limits<real_t>::quiet_NaN());
        mean.Crho       ().setConstant(numeric_limits<real_t>::quiet_NaN());
        mean.Crhou_dot_u().setConstant(numeric_limits<real_t>::quiet_NaN());
    }

    // Convert implicit values at collocation points to expansion coefficients
    bsplineop_lu mass(*cop);
    mass.opform_mass(*cop);
    mass.factor();
    mass.solve(quantities::nscalars::implicit,
               mean.storage.middleCols<quantities::nscalars::implicit>(
                   quantities::nscalars::physical).data(),
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
    mean.save(esioh);
    return super::save_statistics_hook(esioh, t);
}

void
driver::load_statistics_hook(
        const esio_handle esioh)
{
    mean.load(esioh);
    support::load_time(esioh, mean.t);
    return super::load_statistics_hook(esioh);
}

void
driver::default_restart_interval(
        time_type& t,
        step_type&)
{
    real_t flowthrough_time = std::numeric_limits<real_t>::quiet_NaN();
    if (grid && chdef) {
        flowthrough_time =  grid->L.x()
                         / (chdef->bulk_rho_u / chdef->bulk_rho);
    }
    if (boost::math::isnormal(flowthrough_time)) {
        t = flowthrough_time;
    }
}

} // end namespace reacting

} // end namespace suzerain
