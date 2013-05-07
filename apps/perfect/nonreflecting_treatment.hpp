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
 * Provides Giles-like nonreflecting boundary conditions at the upper \f$y =
 * L_y\f$ boundary.  Background and notation set within
 * <tt>writeups/perfect_gas.tex</tt> section entitled "Nonreflecting freestream
 * boundary conditions".
 */
class nonreflecting_treatment
    : public operator_base
    , public timestepper::nonlinear_operator< contiguous_state<4,complex_t> >
{
public:
    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
            const operator_common_block& common);

    /**
     * Applies Giles conditions delegating most processing to #N.
     * That is, compute
     * \f[
     * \underbrace{
     *   {R^Y}^{-1}
     *   \left[V^L S\right]^{-1}
     *   \left(
     *     \chi^{-1}
     *     \left( \ii k_x \left[C^G\right] + \ii k_z \left[B^G\right] \right)
     *     \left[V^L S\right] R^Y M
     *     \hat{V}
     *     +
     *     \left(I - P^G\right)
     *     \left[V^L S\right] R^Y
     *     N(\hat{V})
     *   \right)
     * }_{N_E^G\left(\hat{V}\right)}
     * \f]
     * assuming that #N implements \f$N(\hat{V})\f$.
     *
     * Uses, and on <tt>substep_index == 0</tt> updates, all of
     * <ol>
     *   <li>#VL_S_RY</li>
     *   <li>#BG_VL_S_RY_by_chi</li>
     *   <li>#CG_VL_S_RY_by_chi</li>
     *   <li>#ImPG_CG_VL_S_RY</li>
     *   <li>#invVL_S_RY</li>
     * </ol>
     * during invocation.
     */
    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const real_t evmaxmag_real,
            const real_t evmaxmag_imag,
            const std::size_t substep_index) const;

    /**
     * Compute subsonic nonreflecting boundary condition matrices given
     * reference state information and the \ref scenario_definition provided at
     * construction time.  Updates all of
     * <ol>
     *   <li>#VL_S_RY</li>
     *   <li>#BG_VL_S_RY_by_chi</li>
     *   <li>#CG_VL_S_RY_by_chi</li>
     *   <li>#ImPG_CG_VL_S_RY</li>
     *   <li>#invVL_S_RY</li>
     * </ol>
     * during invocation.
     *
     * @param ref_rho Nondimensional reference density \f$rho\f$.
     * @param ref_u   Nondimensional reference streamwise velocity \f$u\f$.
     * @param ref_v   Nondimensional reference boundary-normal velocity \f$v\f$.
     *                Inflow conditions are used when this velocity is strictly
     *                negative.  Otherwise, outflow conditions are used.
     * @param ref_w   Nondimensional reference spanwise velocity \f$w\f$.
     * @param ref_a   Nondimensional reference sound speed \f$a\f$.
     */
    void compute_subsonic_matrices(
            const real_t ref_rho,
            const real_t ref_u,
            const real_t ref_v,
            const real_t ref_w,
            const real_t ref_a);

    /** The operator whose behavior is modified by this instance. */
    shared_ptr<timestepper::nonlinear_operator<
                contiguous_state<4,complex_t>
            > > N;

protected:

    /** The scenario in which the operator is used. */
    const scenario_definition &scenario;

    /** Provides reference values used when computing the conditions. */
    const operator_common_block &common;

    /**
     * Wavenumber-independent matrix
     * \f$ \left[V^L S\right] R^Y \f
     * used by apply_operator().
     */
    Matrix5r VL_S_RY;

    /**
     * Wavenumber-independent matrix
     * \f$ \chi^{-1} \left[B^G\right] \left[V^L S\right] R^Y \f$
     * used by apply_operator().
     */
    Matrix5r BG_VL_S_RY_by_chi;

    /**
     * Wavenumber-independent matrix
     * \f$ \chi^{-1} \left[C^G\right] \left[V^L S\right] R^Y \f$
     * used by apply_operator().
     */
    Matrix5r CG_VL_S_RY_by_chi;

    /**
     * Wavenumber-independent matrix
     * \f$ \left(I - P^G\right) \left[V^L S\right] R^Y \f$
     * used by apply_operator().
     */
    Matrix5r ImPG_VL_S_RY;

    /**
     * Wavenumber-independent matrix
     * \f$ {R^Y}^{-1} \left[V^L S\right]^{-1} \f$
     * used by apply_operator().
     */
    Matrix5r inv_VL_S_RY;

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
