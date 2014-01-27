//--------------------------------------------------------------------------
//
// Copyright (C) 2013-2014 Rhys Ulerich
// Copyright (C) 2013-2014 The PECOS Development Team
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

#ifndef SUZERAIN_PERFECT_TREATMENT_NONREFLECTING_HPP
#define SUZERAIN_PERFECT_TREATMENT_NONREFLECTING_HPP

/** @file
 * A physics-agnostic, low-storage B-spline mass operator
 */

#include <suzerain/common.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/specification_isothermal.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;
class specification_grid;
class pencil_grid;

namespace perfect {

// Forward declarations
class definition_scenario;

/**
 * Provides Giles-like nonreflecting boundary conditions at the upper \f$y =
 * L_y\f$ boundary.  Background and notation set within
 * <tt>writeups/perfect_gas.tex</tt> section entitled "Nonreflecting freestream
 * boundary conditions".
 */
class treatment_nonreflecting
    : public operator_base
    , public lowstorage::operator_nonlinear< contiguous_state<4,complex_t> >
{
public:
    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Constructor.
     * After construction, #N must be provided.
     */
    treatment_nonreflecting(
            const definition_scenario& scenario,
            const specification_isothermal& isothermal,
            const specification_grid& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b);

    /**
     * Applies Giles conditions delegating most processing to #N.
     * That is, compute
     * \f[
     * \underbrace{
     *   {R^Y}^{-1}
     *   \left[V^L S\right]^{-1}
     *   \left(
     *     \chi^{-1}
     *     \left(   \ii k_x \left[P^G C^G\right]
     *            + \ii k_z \left[P^G B^G\right] \right)
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
     *   <li>#PG_BG_VL_S_RY_by_chi</li>
     *   <li>#PG_CG_VL_S_RY_by_chi</li>
     *   <li>#ImPG_CG_VL_S_RY</li>
     *   <li>#invVL_S_RY</li>
     * </ol>
     * during invocation.
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

    /** The scenario in which the operator is used. */
    const definition_scenario &scenario;

    /**
     * The lower and upper isothermal boundary specification.
     * Provides the reference values used for nonreflecting treatment.
     */
    const specification_isothermal &isothermal;

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
    Matrix5r PG_BG_VL_S_RY_by_chi;

    /**
     * Wavenumber-independent matrix
     * \f$ \chi^{-1} \left[C^G\right] \left[V^L S\right] R^Y \f$
     * used by apply_operator().
     */
    Matrix5r PG_CG_VL_S_RY_by_chi;

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

    /**
     * Compute Medida's Giles-like nonreflecting boundary condition matrices
     * given reference state information and the \ref definition_scenario
     * provided at construction time.  Inflow vs outflow and subsonic vs
     * supersonic conditions are determined using \c ref_v and \c normal_sign.
     * Updates all of
     * <ol>
     *   <li>#VL_S_RY</li>
     *   <li>#PG_BG_VL_S_RY_by_chi</li>
     *   <li>#PG_CG_VL_S_RY_by_chi</li>
     *   <li>#ImPG_CG_VL_S_RY</li>
     *   <li>#inv_VL_S_RY</li>
     * </ol>
     * during invocation.
     *
     * @param ref_rho     Nondimensional reference density \f$rho\f$.
     * @param ref_u       Nondimensional reference streamwise velocity \f$u\f$.
     * @param ref_v       Nondimensional reference velocity \f$v\f$.
     * @param ref_w       Nondimensional reference spanwise velocity \f$w\f$.
     * @param ref_a       Nondimensional reference sound speed \f$a\f$.
     * @param normal_sign Sign of a boundary-normal vector.
     *                    For the boundary at \f$y=0\f$ this should be negative.
     *                    For the boundary at \f$y=L_y\f$, it must be positive.
     */
    void compute_giles_matrices(
            const real_t ref_rho,
            const real_t ref_u,
            const real_t ref_v,
            const real_t ref_w,
            const real_t ref_a,
            const real_t normal_sign);

    /**
     * Invokes compute_giles_matrices() for lower boundary.
     * Broken out for brevity and to document reference state choices.
     */
    void compute_giles_matrices_lower()
    {
        return compute_giles_matrices(isothermal.lower_rho,
                                      isothermal.lower_u,
                                      isothermal.lower_v,
                                      isothermal.lower_w,
                                      std::sqrt(isothermal.lower_T),
                                      -1);
    }

    /**
     * Invokes compute_giles_matrices() for upper boundary.
     * Broken out for brevity and to document reference state choices.
     */
    void compute_giles_matrices_upper()
    {
        return compute_giles_matrices(isothermal.upper_rho,
                                      isothermal.upper_u,
                                      isothermal.upper_v,
                                      isothermal.upper_w,
                                      std::sqrt(isothermal.upper_T),
                                      +1);
    }

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    treatment_nonreflecting(const treatment_nonreflecting&);
    treatment_nonreflecting& operator=(const treatment_nonreflecting&);

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_TREATMENT_NONREFLECTING_HPP */
