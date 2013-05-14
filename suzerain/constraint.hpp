//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_PERFECT_CONSTRAINT_HPP
#define SUZERAIN_PERFECT_CONSTRAINT_HPP

/** @file
 * Provides \ref constraint.
 */

#include <suzerain/common.hpp>

namespace suzerain {

// Forward declarations
class bspline;

/**
 * Encapsulates the ways \ref constraint_treatment can constrain a value.
 */
class constraint
{
public:

    /**
     * Tag types marking in which ways a particular value can be constrained.
     * @{
     */

    /** Enforce mean collocation value at \f$y=0\f$ */
    class value_lower {};

    /** Enforce mean collocation value at \f$y=L_y\f$ */
    class value_upper {};

    /** Enforce bulk value across \f$y=\left[0,L_y\right]\f$ */
    class value_bulk  {};

    /**@}*/

    /** Default constructor creates an NOP constraint. */
    constraint();

    /** Construct a \ref value_lower constraint achieving \ref target. */
    constraint(const value_lower&, const real_t target, bspline &b);

    /** Construct a \ref value_upper constraint achieving \ref target. */
    constraint(const value_upper&, const real_t target, bspline &b);

    /** Construct a \ref value_bulk constraint achieving \ref target. */
    constraint(const value_bulk&, const real_t target, bspline &b);

    /** What scalar value should be targeted by the constraint? */
    real_t target;

    /**
     * What functional, when dotted with B-spline coefficients, computes the
     * quantity of interest?
     * */
    VectorXr coeff;

    /** Is this constraint computable as specified? */
    bool enabled() const
    { return !(boost::math::isnan)(target); }

};

} // namespace suzerain

#endif /* SUZERAIN_PERFECT_CONSTRAINT_HPP */
