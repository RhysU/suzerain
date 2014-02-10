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

#ifndef SUZERAIN_PERFECT_OPERATOR_HYBRID_HPP
#define SUZERAIN_PERFECT_OPERATOR_HYBRID_HPP

/** @file
 * A hybrid implicit/explicit Navier--Stokes operator for isothermal wall(s).
 */

#include <suzerain/common.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/rholut_imexop.h>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bsmbsm_solver;
class specification_grid;
class specification_isothermal;
class pencil_grid;
class specification_zgbsv;

namespace perfect {

// Forward declarations
class operator_common_block;
class definition_scenario;

/**
 * A hybrid implicit/explicit operator that providing isothermal wall conditions
 * and possibly a nonreflecting freestream. It requires interoperation with
 * operator_nonlinear, set on member #N, via operator_common_block.
 */
class operator_hybrid_isothermal
  : public operator_base,
    public lowstorage::linear_operator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    >
  , public lowstorage::operator_nonlinear<
        contiguous_state<4,complex_t>
    >
{
public:

    /** The superclass specifying the linear operator interface. */
    typedef lowstorage::linear_operator<
        multi_array::ref<complex_t,4>,
        contiguous_state<4,complex_t>
    > linear;

    /** The superclass specifying the nonlinear operator interface. */
    typedef lowstorage::operator_nonlinear<
        contiguous_state<4,complex_t>
    > nonlinear;

    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Construct either a one-sided or two-sided instance.
     *
     * When specification_grid::two_sided(), both the lower and upper boundaries
     * are constrained to be isothermal.  When specification_grid::one_sided(),
     * the lower boundary is constrained in this manner and a nonreflecting
     * boundary condition is applied at the upper boundary.
     *
     * After construction, member #N must be populated by a
     * \ref operator_nonlinear instance.
     */
    operator_hybrid_isothermal(
            const specification_zgbsv& spec,
            const definition_scenario& scenario,
            const specification_isothermal& isothermal,
            const specification_grid& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b,
            operator_common_block& common);

    ~operator_hybrid_isothermal();

    virtual void apply_mass_plus_scaled_operator(
             const complex_t& phi,
             multi_array::ref<complex_t,4> &state,
             const std::size_t substep_index) const;

    virtual void accumulate_mass_plus_scaled_operator(
            const complex_t& phi,
            const multi_array::ref<complex_t,4> &input,
            const complex_t& beta,
            contiguous_state<4,complex_t> &output,
            const std::size_t substep_index) const;

    /**
     *
     * Sets no-slip, isothermal boundary conditions at first and possibly the
     * last wall-normal collocation point using rho_wall = e_wall * gamma *
     * (gamma - 1).
     */
    virtual void invert_mass_plus_scaled_operator(
            const complex_t& phi,
            multi_array::ref<complex_t,4> &state,
            const lowstorage::method_interface<complex_t> &method,
            const linear::component delta_t,
            const std::size_t substep_index,
            multi_array::ref<complex_t,4>* ic0 = NULL) const;

    /**
     * Delegates to #N, adjusting the result to achieve
     * a nonreflecting upper boundary when grid.one_sided().
     */
    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const;

    /** The operator whose behavior is modified by this instance. */
    shared_ptr<lowstorage::operator_nonlinear<
                contiguous_state<4,complex_t>
            > > N;

protected:

    /** Controls the solves performed during invert_mass_plus_scaled_operator */
    shared_ptr<bsmbsm_solver> solver;

    /** The scenario in which the operator is used */
    const definition_scenario& scenario;

    /** The isothermal wall boundary details */
    const specification_isothermal& isothermal;

    /** Houses data required for operator application and inversion */
    operator_common_block& common;

    /**
     * When <code>grid.one_sided()</code>, set to the 5x5 column major matrix
     * \f${R^Y}^{-1} \left[V^L S\right]^{-1} \left[P^G C^G\right] \left[V^L
     * S\right] {R^Y}\f$ for upper NRBC on substep 0 during \ref
     * apply_operator().  Otherwise, ignored.
     *
     * \see suzerain_rholut_imexop_accumulate for details behind \c a suffix.
     */
    Matrix5r upper_nrbc_a;

    /**
     * When <code>grid.one_sided()</code>, set to the 5x5 column major matrix
     * \f${R^Y}^{-1} \left[V^L S\right]^{-1} \left[P^G B^G\right] \left[V^L
     * S\right] {R^Y}\f$ for upper NRBC on substep 0 during \ref apply_operator.
     * Otherwise, ignored.
     *
     * \see suzerain_rholut_imexop_accumulate for details behind \c b suffix.
     */
    Matrix5r upper_nrbc_b;

    /**
     * When <code>grid.one_sided()</code>, set to the 5x5 column major matrix
     * \f${R^Y}^{-1} \left[V^L S\right]^{-1} \left[P^G    \right] \left[V^L
     * S\right] {R^Y}\f$ for upper NRBC on substep 0 during \ref apply_operator.
     * Otherwide, ignored.
     *
     * \see suzerain_rholut_imexop_accumulate for details behind \c c suffix.
     */
    Matrix5r upper_nrbc_c;

    /**
     * When <code>grid.one_sided()</code>, set to the 5x5 column major matrix
     * \f${R^Y}^{-1} \left[V^L S\right]^{-1} \left[I - P^G\right] \left[V^L
     * S\right] {R^Y}\f$ for upper NRBC on substep 0 during \ref apply_operator.
     * Otherwide, ignored.
     */
    Matrix5r upper_nrbc_n;

    /**
     * Invokes compute_giles_matrices() for upper boundary.
     * Broken out separately to help document reference state choices.
     */
    void compute_giles_matrices_upper();

private:

    /** Pack scenario parameters for rholut_imexop.h usage */
    suzerain_rholut_imexop_scenario imexop_s() const;

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_OPERATOR_HYBRID_HPP */
