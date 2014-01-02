//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
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

#ifndef SUZERAIN_REACTING_HYBRID_OPERATOR_HPP
#define SUZERAIN_REACTING_HYBRID_OPERATOR_HPP

/** @file
 * Hybrid implicit/explicit Navier--Stokes operators.
 */

#include <suzerain/lowstorage.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/reacting_imexop.h>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bsmbsm_solver;
class grid_specification;
class pencil_grid;
class isothermal_specification;
class zgbsv_specification;

namespace reacting {

// Forward declarations
class operator_common_block;
class antioch_constitutive;
class channel_definition;

/**
 * A hybrid implicit operator that provides no slip, isothermal walls.  It
 * requires interoperation with nonlinear_operator via operator_common_block.
 */
class isothermal_hybrid_linear_operator
  : public operator_base,
    public lowstorage::linear_operator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    >
{
public:

    isothermal_hybrid_linear_operator(
            const zgbsv_specification& spec,
            const antioch_constitutive &cmods,
            const isothermal_specification &isospec,
            const channel_definition &chdef,
            const grid_specification &grid,
            const pencil_grid &dgrid,
            const bsplineop &cop,
            bspline &b,
            operator_common_block &common);

    ~isothermal_hybrid_linear_operator();

    virtual void apply_mass_plus_scaled_operator(
             const complex_t &phi,
             multi_array::ref<complex_t,4> &state,
             const std::size_t substep_index) const;

    virtual void accumulate_mass_plus_scaled_operator(
            const complex_t &phi,
            const multi_array::ref<complex_t,4> &input,
            const complex_t &beta,
            contiguous_state<4,complex_t> &output,
            const std::size_t substep_index) const;

    /**
     * Sets no-slip, isothermal boundary conditions at first and last
     * wall-normal collocation point using that
     * rho_wall = e_wall * gamma * (gamma - 1).
     */
    virtual void invert_mass_plus_scaled_operator(
            const complex_t &phi,
            multi_array::ref<complex_t,4> &state,
            const lowstorage::method_interface<complex_t> &method,
            const component delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4> *ic0 = NULL) const;

protected:

    /** Controls the solves performed during invert_mass_plus_scaled_operator */
    shared_ptr<bsmbsm_solver> flow_solver;

    std::vector<shared_ptr<bsmbsm_solver> > species_solver;

    /** Provides info about isothermal specifications */
    const isothermal_specification &isospec;

    /** Provides info about channel case, e.g., T_wall */
    const channel_definition &chdef;

    /** Access to constitutive laws */
    const antioch_constitutive &cmods;

    /** Houses data required for operator application and inversion */
    operator_common_block &common;

private:

    /** Pack scenario parameters for reacting_imexop.h usage */
    suzerain_reacting_imexop_scenario imexop_s() const;

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

} // namespace reacting

} // namespace suzerain

#endif  /* SUZERAIN_REACTING_HYBRID_OPERATOR_HPP */
