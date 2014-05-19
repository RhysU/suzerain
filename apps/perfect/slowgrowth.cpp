//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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
 * @copydoc slowgrowth.hpp
 */

#include "slowgrowth.hpp"

#include <largo/largo.h>

#include <suzerain/common.hpp>
#include <suzerain/error.h>
#include <suzerain/ndx.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timers.h>

#include "instantaneous.hpp"

namespace suzerain {

namespace perfect {

slowgrowth::slowgrowth()
    : meanrms(0)
    , meanrms_y(0)
{
}

void
slowgrowth::initialize(
        const type slow_treatment,
        const specification_largo &sg,
        const real_t code_Ma,
        const std::size_t substep_index)
{
    switch (slow_treatment) {
    default:
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();

    case slowgrowth::none:
        break;

    case slowgrowth::largo:
        // If necessary, perform initialization calls for Largo
        if (substep_index == 0) {

            SUZERAIN_TIMER_SCOPED("slowgrowth::initialize");

            // Initialize the slow growth workspace
            // Avoids debugging-related memory allocation if at all possible
            SUZERAIN_ENSURE(sg.gramp_mean.size());
            SUZERAIN_ENSURE(sg.gramp_mean.size() == sg.gramp_rms.size());
            if (SUZERAIN_UNLIKELY(sg.ignore_gramp_mean||sg.ignore_gramp_rms)) {
                std::vector<real_t> z(sg.gramp_mean.size(), 0);
                largo_init(sg.workspace,
                           sg.grdelta,
                           sg.ignore_gramp_mean ? &z[0] : &sg.gramp_mean[0],
                           sg.ignore_gramp_rms  ? &z[0] : &sg.gramp_rms [0]);
            } else {
                largo_init(sg.workspace,
                           sg.grdelta,
                           &sg.gramp_mean[0],
                           &sg.gramp_rms[0]);
            }

            // Prepare any necessary baseflow information at the wall
            // assuming that baseflow state already matched code_Ma.
            largo_state dy, dx;
            if (sg.baseflow) {
                sg.baseflow->conserved(
                        0.0, basewall.as_is(), dy.as_is(), dx.as_is());
                sg.baseflow->pressure (
                        0.0, basewall.p, dy.p, dx.p); // as_is()
            }

            // TODO Implement, though currently always trivially zero
            largo_state dt, src;

            // Present the baseflow information to Largo
            const real_t inv_Ma2 = 1 / (code_Ma * code_Ma);
            largo_init_wall_baseflow(sg.workspace,
                                     basewall.rescale(inv_Ma2),
                                     dy      .rescale(inv_Ma2),
                                     dt      .rescale(inv_Ma2),
                                     dx      .rescale(inv_Ma2),
                                     src     .rescale(inv_Ma2));
        }
        break;

    }
}

void
slowgrowth::gather_wavexz(
        const slowgrowth::type slow_treatment,
        const operator_base &o,
        const contiguous_state<4,complex_t> &swave)
{
    switch (slow_treatment) {
    default:
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();

    case slowgrowth::none:
        break;

    case slowgrowth::largo:
        {
            SUZERAIN_TIMER_SCOPED("slowgrowth::gather_wavexz");

            // With state being Fourier in X and Z but collocation in Y...
            // ...compute L^2_{xz} of state at each collocation point
            meanrms = compute_field_L2xz(swave, o.grid, o.dgrid);

            // ...and rescale to convert to root-mean-square (RMS) fluctuations
            // (mean L2 values are unused so also defensively NaN that storage)
            const real_t rms_adjust = 1 / sqrt(o.grid.L.x() * o.grid.L.z());
            for (size_t i = 0; i < meanrms.size(); ++i) {
                meanrms[i].mean.setConstant(
                        std::numeric_limits<real_t>::quiet_NaN());
                meanrms[i].fluctuating *= rms_adjust;
            }
        }
        break;
    }
}

void
slowgrowth::gather_physical_cons(
        const slowgrowth::type slow_treatment,
        const operator_base &o,
        const instantaneous &inst)
{
    switch (slow_treatment) {
    default:
        SUZERAIN_ERROR_VOID_UNIMPLEMENTED();

    case slowgrowth::none:
        break;

    case slowgrowth::largo:
        {
            SUZERAIN_TIMER_SCOPED("slowgrowth::gather_physical_cons");

            // Slow growth requires mean conserved state at collocation points.
            // Abuse unused pieces within 'meanrms' to avoid more allocations.
            SUZERAIN_ENSURE(meanrms.size() == 5);
            meanrms[ndx::e  ].mean = inst.rhoE();
            meanrms[ndx::mx ].mean = inst.rhou();
            meanrms[ndx::my ].mean = inst.rhov();
            meanrms[ndx::mz ].mean = inst.rhow();
            meanrms[ndx::rho].mean = inst.rho();

            // Slow growth requires mean pressure and its fluctuation profiles.
            // Abuse 'meanrms' to store the information after conserved state.
            // The RMS of fluctuating pressure computation can be numerically
            // noisy hence require a non-negative result prior to the sqrt.
            meanrms.resize(meanrms.size() + 1);
            meanrms.back().mean        = inst.p();
            meanrms.back().fluctuating = (inst.p2() - inst.p().square())
                                         .max(0).sqrt();

            // Wall-normal derivatives of every mean RMS quantity are required
            meanrms_y.resize(meanrms.size());
            ArrayX2r tmp(o.dgrid.global_physical_extent.y(), 2);
            std::vector<field_L2xz>::const_iterator src = meanrms.begin();
            std::vector<field_L2xz>::const_iterator end = meanrms.end();
            std::vector<field_L2xz>::      iterator dst = meanrms_y.begin();
            for (/*just above*/; src != end; ++src, ++dst) {
                tmp.col(0) = (*src).mean;
                tmp.col(1) = (*src).fluctuating;
                o.masslu()->solve(tmp.cols(), tmp.data(),
                                  tmp.innerStride(), tmp.outerStride());
                (*dst).mean       .resizeLike(tmp.col(0));
                (*dst).fluctuating.resizeLike(tmp.col(1));
                o.cop.accumulate(1, 1.0, tmp.col(0).data(), tmp.innerStride(),
                                    0.0, (*dst).mean.data(), 1);
                o.cop.accumulate(1, 1.0, tmp.col(1).data(), tmp.innerStride(),
                                    0.0, (*dst).fluctuating.data(), 1);
            }
        }
        break;
    }
}


} // namespace perfect

} // namespace suzerain
