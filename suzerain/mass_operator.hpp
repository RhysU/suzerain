//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
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
// mass_operator.hpp: A physics-agnostic, low-storage B-spline mass operator
// $Id$

#ifndef SUZERAIN_MASS_OPERATOR_HPP
#define SUZERAIN_MASS_OPERATOR_HPP

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/grid_specification.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/state_fwd.hpp>
#include <suzerain/timestepper.hpp>

namespace suzerain {

/** A linear operator which applies or inverts a B-spline mass matrix. */
class mass_operator
  : public operator_base,
    public timestepper::lowstorage::linear_operator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    >
{
public:

    mass_operator(
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b);

    virtual ~mass_operator();

    virtual void apply_mass_plus_scaled_operator(
             const complex_t &phi,
             multi_array::ref<complex_t,4> &state,
             const timestepper::lowstorage::method_interface<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void accumulate_mass_plus_scaled_operator(
             const complex_t &phi,
             const multi_array::ref<complex_t,4> &input,
             const complex_t &beta,
             contiguous_state<4,complex_t> &output,
             const timestepper::lowstorage::method_interface<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index) const;

     virtual void invert_mass_plus_scaled_operator(
             const complex_t &phi,
             multi_array::ref<complex_t,4> &state,
             const timestepper::lowstorage::method_interface<complex_t> &method,
             const component delta_t,
             const std::size_t substep_index,
             multi_array::ref<complex_t,4> *ic0 = NULL) const;

private:

     /** Precomputed mass matrix factorization */
    bsplineop_luz massluz;

};

} // namespace suzerain

#endif  /* SUZERAIN_MASS_OPERATOR_HPP */
