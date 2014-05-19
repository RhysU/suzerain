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

#include <largo/largo.h>

#include <suzerain/common.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/specification_largo.hpp>
#include <suzerain/state.hpp>
#include <suzerain/timers.h>

#include "slowgrowth.hpp"

namespace suzerain {

namespace perfect {

slowgrowth::slowgrowth()
    : meanrms(0)
{
}

void
slowgrowth::gatherwave_rms(const slowgrowth::type slow_treatment,
                           const operator_base &o,
                           const contiguous_state<4,complex_t> &swave)
{
    // Quick return sans timers whenever no work required
    if (slow_treatment == slowgrowth::none) return;

    SUZERAIN_TIMER_SCOPED("slowgrowth::gatherwave_rms");

    // With conserved state being Fourier in X and Z but collocation in Y...
    // ...collectively compute L^2_{xz} of state at each collocation point
    meanrms = compute_field_L2xz(swave, o.grid, o.dgrid);

    // ...and rescale to convert to root-mean-square (RMS) fluctuations
    // (mean L2 values are uninteresting so also defensively NaN storage).
    const real_t rms_adjust = 1 / sqrt(o.grid.L.x() * o.grid.L.z());
    for (size_t i = 0; i < meanrms.size(); ++i) {
        meanrms[i].mean.setConstant(std::numeric_limits<real_t>::quiet_NaN());
        meanrms[i].fluctuating *= rms_adjust;
    }
}

} // namespace perfect

} // namespace suzerain
