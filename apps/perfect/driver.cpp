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

#include <suzerain/support/logging.hpp>

#include "nsctpl_rholut.hpp"
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
driver::compute_statistics()
{
    const real_t NaN = std::numeric_limits<real_t>::quiet_NaN();
    const real_t t   = controller ? controller->current_t() : NaN;

    // Defensively avoid multiple invocations with no intervening changes
#pragma warning(push,disable:1572)
    if (controller && t == mean.t) {
#pragma warning(pop)
        DEBUG0("Cowardly refusing to re-sample statistics at t = " << t);
        return;
    }

    SUZERAIN_TIMER_SCOPED("driver::compute_statistics");

    // Obtain mean samples from instantaneous fields stored in state_linear
    const double starttime = MPI_Wtime();
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
        mean.f          ().setConstant(NaN);
        mean.f_dot_u    ().setConstant(NaN);
        mean.qb         ().setConstant(NaN);
        mean.CrhoE      ().setConstant(NaN);
        mean.Crhou      ().setConstant(NaN);
        mean.Crho       ().setConstant(NaN);
        mean.Crhou_dot_u().setConstant(NaN);
    }

    // Convert implicit values at collocation points to expansion coefficients
    bsplineop_lu mass(*cop);
    mass.opform_mass(*cop);
    mass.factor();
    mass.solve(mean_quantities::nscalars::implicit,
               mean.storage.middleCols<mean_quantities::nscalars::implicit>(
                   mean_quantities::nscalars::physical).data(),
               mean.storage.innerStride(), mean.storage.outerStride());

    const double elapsed = MPI_Wtime() - starttime;
    if (controller) {
        std::string timeprefix(build_timeprefix(t, controller->current_nt()));
        INFO0(timeprefix << " Computed statistics in " << elapsed << " seconds");
    }
}

// FIXME Implement

} // end namespace perfect

} // end namespace suzerain
