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

namespace constraint {

// TODO The base::coeff API would be better done using Eigen 3.2 Refs

/** Abstract base class for constraints expressible as a functional. */
class base
{
public:

    /**
     * What functional, when dotted with B-spline coefficients, computes the
     * quantity of interest?
     * */
    VectorXr coeff;

    /**
     * What scalar value should be targeted by the constraint?  This value may
     * be nontrivial to compute and should be cached when used repeated.
     */
    virtual real_t target() const = 0;

    /** Is this constraint computable as specified? */
    virtual bool enabled() const = 0;

    /** Virtual destructor to permit use as a base class. */
    virtual ~base();
};

/** Constrain the collocation value at \f$y=0\f$. */
class lower : public virtual base
{
public:

    explicit lower(bspline &b);

};

/** Constrain the collocation value at \f$y=L_y\f$. */
class upper : public virtual base
{
public:

    explicit upper(bspline &b);

};

/** Constrain the bulk value across \f$y=\left[0,L_y\right]\f$. */
class bulk : public virtual base
{
public:

    explicit bulk(bspline &b);

};

/** Target some constant value. */
class constant : public virtual base
{
public:

    /** Specify the constant value to target. */
    explicit constant(const real_t target);

    virtual real_t target() const;

    virtual bool enabled() const;

private:

    /** Stores the constant target value. */
    real_t t;

};

/**
 * Constrain the collocation value at \f$y=0\f$ to be some given constant.
 */
class constant_lower : public virtual constant, public virtual lower
{
public:

    constant_lower(const real_t target, bspline &b);

};

/**
 * Constrain the collocation value at \f$y=L_y\f$ to be some given constant.
 */
class constant_upper : public virtual constant, public virtual upper
{
public:

    constant_upper(const real_t target, bspline &b);

};

/**
 * Constrain the bulk value across \f$y=\left[0,L_y\right]\f$ to be some
 * constant.
 */
class constant_bulk : public virtual constant, public virtual bulk
{
public:

    constant_bulk(const real_t target, bspline &b);

};

} // namespace constraint

} // namespace suzerain

#endif /* SUZERAIN_PERFECT_CONSTRAINT_HPP */
