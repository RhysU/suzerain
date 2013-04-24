//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

#ifndef SUZERAIN_NONREFLECTING_TREATMENT_HPP
#define SUZERAIN_NONREFLECTING_TREATMENT_HPP

/** @file
 * A physics-agnostic, low-storage B-spline mass operator
 */

#include <suzerain/common.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;
class grid_specification;
class pencil_grid;

namespace perfect {

// Forward declarations
class operator_common_block;
class scenario_definition;

/**
 * Provides Giles-like nonreflecting boundary conditions at the upper boundary.
 * Background and notation set within <tt>writeups/perfect_gas.tex</tt> section
 * entitled "Nonreflecting freestream boundary conditions".
 */
class nonreflecting_treatment
    : public operator_base
    , public timestepper::nonlinear_operator< contiguous_state<4,complex_t> >
{
public:

    /**
     * Constructor.
     * After construction, #N must be provided.
     */
    nonreflecting_treatment(
            const scenario_definition& scenario,
            const grid_specification& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b,
            operator_common_block& common);

    /** Applies Giles conditions delegating most processing to #N */
    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const;

    /** The operator whose behavior is modified by this instance. */
    shared_ptr<timestepper::nonlinear_operator<
                contiguous_state<4,complex_t>
            > > N;

protected:

    /** The scenario in which the operator is used. */
    const scenario_definition &scenario;

    /** Provides reference values used when computing the conditions. */
    operator_common_block &common;

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    // Using boost::noncopyable trips Intel non-virtual base destructor warnings.
    nonreflecting_treatment(const nonreflecting_treatment&);
    nonreflecting_treatment& operator=(const nonreflecting_treatment&);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_NONREFLECTING_TREATMENT_HPP */
