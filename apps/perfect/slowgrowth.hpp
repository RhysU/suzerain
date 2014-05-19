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

#ifndef SUZERAIN_PERFECT_SLOWGROWTH_HPP
#define SUZERAIN_PERFECT_SLOWGROWTH_HPP

/** @file
 * Declarations related to slow growth formulation selection.
 */

#include <suzerain/largo_state.hpp>
#include <suzerain/l2.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class operator_base;
class specification_largo;

namespace perfect {

/** Provides scoping semantics for slowgrowth::type */
class slowgrowth {

public:

    /**
     * What slow growth sources are employed?
     * All valid types evaluate to true in a boolean context.
     */
    enum type {
        none  = 1,  ///< No slow growth sources
        largo       ///< Slow growth sources computed by Largo library
    };

    slowgrowth();

    void
    initialize(const type slow_treatment,
               const specification_largo &sg,
               const real_t code_Ma,
               const std::size_t substep_index = 0);

    void
    gather_wavexz_rms(const type slow_treatment,
                      const operator_base &o,
                      const contiguous_state<4,complex_t> &swave);

    //gatherphys_mean
    //differentiate_meanrms
    //inner_y
    //inner_xz
    //forcing_xyz

private:

    largo_state basewall;

    std::vector<field_L2xz> meanrms;

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_SLOWGROWTH_HPP */
