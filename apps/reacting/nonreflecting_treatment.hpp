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

#ifndef SUZERAIN_NONREFLECTING_TREATMENT_HPP
#define SUZERAIN_NONREFLECTING_TREATMENT_HPP

/** @file
 * A physics-agnostic, low-storage B-spline mass operator
 */

#include <suzerain/common.hpp>
#include <suzerain/lowstorage.hpp>
#include <suzerain/operator_base.hpp>
#include <suzerain/state_fwd.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;
class grid_specification;
class pencil_grid;

namespace reacting {

// Forward declarations
class operator_common_block;
// class scenario_definition;

/**
 * Provides Giles-like nonreflecting boundary conditions at the upper boundary.
 * Background and notation set within <tt>writeups/reacting_gas.tex</tt> section
 * entitled "Nonreflecting freestream boundary conditions".
 */
class nonreflecting_treatment
    : public operator_base
    , public lowstorage::nonlinear_operator< contiguous_state<4,complex_t> >
{
public:

    /**
     * Constructor.
     * After construction, #N must be provided.
     */
    nonreflecting_treatment(
//             const scenario_definition& scenario,
            const grid_specification& grid,
            const pencil_grid& dgrid,
            const bsplineop& cop,
            bspline& b,
            operator_common_block& common
            );

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
     */
    virtual std::vector<real_t> apply_operator(
            const real_t time,
            contiguous_state<4,complex_t> &swave,
            const lowstorage::method_interface<complex_t> &method,
            const std::size_t substep_index) const;

    /** The operator whose behavior is modified by this instance. */
    shared_ptr<lowstorage::nonlinear_operator<
                contiguous_state<4,complex_t>
            > > N;

protected:

    /** Provides reference values used when computing the conditions. */
    operator_common_block &common;

private:

    /** Helps to identify from whom logging messages are being emitted. */
    std::string who;

    // boost::noncopyable trips Intel non-virtual base destructor warnings.
    nonreflecting_treatment(const nonreflecting_treatment&);
    nonreflecting_treatment& operator=(const nonreflecting_treatment&);

    /** FIXME: Numbrer of species. */
    size_t Ns;

    /** Mean species densities at upper boundary. */
    mutable VectorXr rhos;
};

} // namespace reacting

} // namespace suzerain

#endif  /* SUZERAIN_NONREFLECTING_TREATMENT_HPP */
