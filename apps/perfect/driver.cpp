//--------------------------------------------------------------------------
//
// Copyright (C) 2008, 2009, 2010, 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// driver.cpp: A driver for nondimensional perfect gas simulations
// $Id$

#include "driver.hpp"

#include <suzerain/format.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/support/logging.hpp>
#include <suzerain/support/support.hpp>

#include "perfect.hpp"

namespace suzerain {

namespace perfect {

driver::driver(
        const std::string &application_synopsis,
        const std::string &description,
        const std::string &revstr)
    : driver_base(application_synopsis,
                  description,
                  revstr)
{
    this->fields = default_fields();
}

std::vector<std::string>
driver::initialize(
        int argc,
        char **argv)
{
    // msoln is not used by all binaries and is therefore not added below
    options.add_definition(*scenario);

    return super::initialize(argc, argv);
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
    state_nonlinear->assign(*state_linear);
    accumulate_manufactured_solution(
            1, *msoln, -1, *state_nonlinear,
            *grid, *dgrid, *b, *cop, t);
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
driver::compute_statistics(
        driver::time_type t)
{
    SUZERAIN_TIMER_SCOPED("driver::compute_statistics");

    // Obtain mean samples from instantaneous fields stored in state_linear
    state_nonlinear->assign(*state_linear);
    mean = perfect::sample_mean_quantities(
            *scenario, *grid, *dgrid, *b, *cop, *state_nonlinear, t);

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
        WARN0("Could not obtain mean quantities set by implicit forcing");
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
    mass.solve(mean_quantities::nscalars::implicit,
               mean.storage.middleCols<mean_quantities::nscalars::implicit>(
                   mean_quantities::nscalars::physical).data(),
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
    save(esioh, *scenario);
    save(esioh, msoln, *scenario, *grid);
    return;
}

void
driver::load_metadata_hook(
        const esio_handle esioh)
{
    super::load_metadata_hook(esioh);
    load(esioh, *scenario);
    load(esioh, msoln, *scenario, *grid);
    return;
}

bool
driver::save_statistics_hook(
        const esio_handle esioh)
{
    // Either compute mean quantities or re-use existing data in this->mean
    const time_type t  = controller->current_t();
#pragma warning(push,disable:1572)
    if (t == mean.t) {
#pragma warning(pop)
        DEBUG0("Cowardly refusing to re-sample statistics at t = " << t);
    } else {
        const double starttime = MPI_Wtime();
        compute_statistics(t);
        const double elapsed = MPI_Wtime() - starttime;
        const step_type nt = controller->current_nt();
        const std::string timeprefix(build_timeprefix(t, nt));
        INFO0(timeprefix << " Computed statistics in "
              << elapsed << " seconds");
    }

    save(esioh, mean);
    return super::save_statistics_hook(esioh);
}

void
driver::load_statistics_hook(
        const esio_handle esioh)
{
    load(esioh, mean);
    support::load_time(esioh, mean.t);
    return super::load_statistics_hook(esioh);
}

} // end namespace perfect

} // end namespace suzerain
