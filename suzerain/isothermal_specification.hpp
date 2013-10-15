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

#ifndef SUZERAIN_ISOTHERMAL_SPECIFICATION_HPP
#define SUZERAIN_ISOTHERMAL_SPECIFICATION_HPP

/** @file
 * Provides \ref isothermal_specification.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Holds parameters for specifying simple isothermal boundary conditions on the
 * \f$y=0\f$ and \f$y=L_y\f$ planes.  These conditions include no-slip walls,
 * transpiring walls, and constant freestream boundaries.  Code consuming
 * these values may treat <tt>NaN</tt> as a do-not-enforce-constraint flag.
 */
class isothermal_specification
{
public:

    /**
     * Construct an instance with all parameters set to NaN and zero-length
     * #lower_cs and #upper_cs.  Clients can use NaN as a not-yet-specified or
     * use-the-default value.
     */
    isothermal_specification();

    /**
     * Specify no-slip walls with given \f$T\f$ at both boundaries
     * for a single "dilluter" species.
     *
     * @param wall_T \f$T\f$ for both walls.
     */
    isothermal_specification(real_t wall_T);

    /**
     * Specify no-slip walls with given \f$T\f$ at both boundaries.
     * @copydetails isothermal_specification(real_t)
     * @param wall_cs Mass fractions, one per species, at both walls.
     */
    isothermal_specification(real_t wall_T,
                             const std::vector<real_t>& wall_cs);

    /**
     * Specify a transpiring wall with given inflow velocity \f$v\f$ and
     * temperature \f$T\f$ at both boundaries for a single "diluter" species.
     *
     * @param wall_T          \f$T\f$ for both walls.
     * @param inflow_velocity \f$v\f$ for both walls with positive values
     *                        oriented blowing into the domain interior.
     */
    isothermal_specification(real_t wall_T,
                             real_t inflow_velocity);

    /**
     * Specify a transpiring wall with given inflow velocity \f$v\f$ and
     * temperature \f$T\f$ at both boundaries.
     * @copydetails isothermal_specification(real_t,real_t)
     * @param wall_cs Mass fractions, one per species, at both boundaries.
     */
    isothermal_specification(real_t wall_T,
                             real_t inflow_velocity,
                             const std::vector<real_t>& wall_cs);

    /**
     * Specify two different temperatures and wall blowing velocities
     * for a single "diluter" species.
     *
     * @param lower_T \f$T\f$ at \f$y=0\f$
     * @param lower_v \f$v\f$ at \f$y=0\f$
     * @param upper_T \f$T\f$ at \f$y=L_y\f$
     * @param upper_v \f$v\f$ at \f$y=L_y\f$
     */
    isothermal_specification(real_t lower_T,
                             real_t lower_v,
                             real_t upper_T,
                             real_t upper_v);

    /**
     * Specify two different temperatures and wall blowing velocities.
     * @copydetails isothermal_specification(real_t,real_t,real_t,real_t)
     * @param lower_cs Mass fractions, one per species, at \f$y=0\f$ boundary.
     * @param upper_cs Mass fractions, one per species, at \f$y=L_y\f$ boundary.
     * @throw invalid_argument if <tt>lower_cs.size() != upper_cs.size()</tt>
     */
    isothermal_specification(real_t lower_T,
                             real_t lower_v,
                             const std::vector<real_t>& lower_cs,
                             real_t upper_T,
                             real_t upper_v,
                             const std::vector<real_t>& upper_cs);

    /**
     * Specify two different temperatures, wall blowing velocities,
     * and associated streamwise velocities for a single "diluter" species.
     *
     * <tt>NaN</tt> may be specified for \c lower_v or \c upper_v to indicate
     * that the boundary-normal velocity should be determined by some other
     * means.  Such usage requires auxiliary information (e.g. an instantaneous
     * mean velocity at the boundary) in conjunction with a boundary condition
     * implementation (e.g. \ref isothermal_mass_operator).
     *
     * @param lower_T \f$T\f$ at \f$y=0\f$
     * @param lower_u \f$u\f$ at \f$y=0\f$
     * @param lower_v \f$v\f$ at \f$y=0\f$
     * @param upper_T \f$T\f$ at \f$y=L_y\f$
     * @param upper_u \f$u\f$ at \f$y=L_y\f$
     * @param upper_v \f$v\f$ at \f$y=L_y\f$
     */
    isothermal_specification(real_t lower_T,
                             real_t lower_u,
                             real_t lower_v,
                             real_t upper_T,
                             real_t upper_u,
                             real_t upper_v);

    /**
     * Specify two different temperatures, wall blowing velocities,
     * and associated streamwise velocities.
     * @copydetails isothermal_specification(real_t,real_t,real_t,real_t,real_t,real_t)
     * @param lower_cs Mass fractions, one per species, at \f$y=0\f$ boundary.
     * @param upper_cs Mass fractions, one per species, at \f$y=L_y\f$ boundary.
     * @throw invalid_argument if <tt>lower_cs.size() != upper_cs.size()</tt>
     */
    isothermal_specification(real_t lower_T,
                             real_t lower_u,
                             real_t lower_v,
                             const std::vector<real_t>& lower_cs,
                             real_t upper_T,
                             real_t upper_u,
                             real_t upper_v,
                             const std::vector<real_t>& upper_cs);

    /**
     * Specify two different temperatures, wall blowing velocities,
     * associated streamwise velocities, and reference densities.
     * @copydetails isothermal_specification(real_t,real_t,real_t,real_t,real_t,real_t,real_t,real_t)
     * @param lower_cs Mass fractions, one per species, at \f$y=0\f$ boundary.
     * @param upper_cs Mass fractions, one per species, at \f$y=L_y\f$ boundary.
     * @throw invalid_argument if <tt>lower_cs.size() != upper_cs.size()</tt>
     */
    isothermal_specification(real_t lower_T,
                             real_t lower_u,
                             real_t lower_v,
                             real_t lower_rho,
                             const std::vector<real_t>& lower_cs,
                             real_t upper_T,
                             real_t upper_u,
                             real_t upper_v,
                             real_t upper_rho,
                             const std::vector<real_t>& upper_cs);

    /**
     * Conditions on the \f$y=0\f$ boundary.
     * @{
     */

    /**
     * Temperature \f$T\f$ at \f$y=0\f$.
     * Though not enforced, this should be strictly positive for realizability.
     */
    real_t lower_T;

    /**
     * Streamwise velocity \f$u\f$ at \f$y=0\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * This would permit, for example, slip boundaries.
     */
    real_t lower_u;

    /**
     * Wall-normal velocity \f$v\f$ at \f$y=0\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * This would permit, for example, evolving the reference
     * velocity of an inflow or outflow boundary.
     */
    real_t lower_v;

    /**
     * Spanwise velocity \f$w\f$ at \f$y=0\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * This would permit, for example, slip boundaries.
     */
    real_t lower_w;

    /**
     * Density \f$rho\f$ at \f$y=0\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * An assigned value may serve as a reference value for
     * a nonreflecting boundary.
     */
    real_t lower_rho;

    /**
     * Species mass fractions \f$c_s\f$ at \f$y=0\f$ including the "diluter"
     * species as index zero.  Though not enforced, these should sum to one for
     * consistency.
     */
    std::vector<real_t> lower_cs;

    /**@}*/

    /**
     * Conditions on the \f$y=L_y\f$ boundary.
     * @{
     */

    /**
     * Temperature \f$T\f$ at \f$y=L_y\f$.
     * Though not enforced, this should be strictly positive for realizability.
     */
    real_t upper_T;

    /**
     * Streamwise velocity \f$u\f$ at \f$y=L_y\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * This would permit, for example, slip boundaries.
     */
    real_t upper_u;

    /**
     * Wall-normal velocity \f$v\f$ at \f$y=L_y\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * This would permit, for example, evolving the reference
     * velocity of an inflow or outflow boundary.
     */
    real_t upper_v;

    /**
     * Spanwise velocity \f$w\f$ at \f$y=L_y\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * This would permit, for example, slip boundaries.
     */
    real_t upper_w;

    /**
     * Density \f$rho\f$ at \f$y=L_y\f$.
     * Code may treat <tt>NaN</tt> as a do-not-enforce flag.
     * An assigned value may serve as a reference value for
     * a nonreflecting boundary.
     */
    real_t upper_rho;

    /**
     * Species mass fractions \f$c_s\f$ at \f$y=L_y\f$ including the "diluter"
     * species as index zero.  Though not enforced, these should sum to one for
     * consistency.
     */
    std::vector<real_t> upper_cs;

    /**@}*/

    /**
     * Shorthand for computing the increment from \f$y=0\f$ to \f$y=L_y\f$.
     * @{
     */
    real_t delta_T  () const { return upper_T   - lower_T  ; }
    real_t delta_u  () const { return upper_u   - lower_u  ; }
    real_t delta_v  () const { return upper_v   - lower_v  ; }
    real_t delta_w  () const { return upper_w   - lower_w  ; }
    real_t delta_rho() const { return upper_rho - lower_rho; }
    /**@}*/

    /**
     * Shorthand for computing the ratio of \f$y=L_y\f$ to \f$y=0\f$.
     * @{
     */
    real_t ratio_T  () const { return upper_T   / lower_T  ; }
    real_t ratio_u  () const { return upper_u   / lower_u  ; }
    real_t ratio_v  () const { return upper_v   / lower_v  ; }
    real_t ratio_w  () const { return upper_w   / lower_w  ; }
    real_t ratio_rho() const { return upper_rho / lower_rho; }
    /**@}*/

};

} // namespace suzerain

#endif // SUZERAIN_ISOTHERMAL_SPECIFICATION_HPP
