//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#ifndef SUZERAIN_PERFECT_CONSTRAINT_HPP
#define SUZERAIN_PERFECT_CONSTRAINT_HPP

/** @file
 * Provides \ref constraint.
 */

#include <suzerain/common.hpp>

namespace suzerain {

// Forward declarations
class bspline;
class bsplineop;

/** Provide logic for constraining mean state in some fashion. */
namespace constraint {

// TODO The base::coeff API would be better done using Eigen 3.2 Refs

/** Abstract base class for constraints expressible as a functional. */
class base
{
public:

    /**
     * What functional, when dotted with B-spline coefficients, computes the
     * quantity of interest?  This member should only be used when enabled()
     * returns \c true.
     * */
    VectorXr coeff;

    /**
     * What shape, stored as relative magnitudes at collocation points,
     * expresses the wall-normal forcing profile to use when enforcing
     * constraints?  For example, a uniform profile would be expressed using
     * all ones.  A point source would be all zeros with the exception of a
     * single one.  Unlike #coeff, this member must be usable even when
     * enabled() returns \c false.
     */
    ArrayXr shape;

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

/**
 * A building block preparing a uniform forcing profile in #shape.
 */
struct uniform : virtual base
{
    explicit uniform(const bspline &b);

    explicit uniform(const bsplineop &bop);
};

/**
 * A building block preparing a linear forcing profile in #shape.
 *
 * The profile takes value \c lower at the 0th collocation point
 * and value \c upper at the last collocation point.
 */
struct linear : virtual base
{
    linear(const bspline &b, const real_t lower, const real_t upper);
};

/**
 * A building block preparing a forcing profile in #shape targeting the
 * <tt>i</tt>-th B-spline coefficient.
 */
struct coefficient : virtual base
{
    explicit coefficient(const bsplineop &bop, const int i);
};

/**
 * A building block to constrain the <tt>nderiv</tt>-th function derivative at
 * location \f$y=0\f$.
 */
struct lower : virtual base
{
    explicit lower(const bsplineop& bop, const int nderiv);
};

/**
 * A building block to constrain the <tt>nderiv</tt>-th function derivative at
 * location \f$y=L_y\f$.
 */
struct upper : virtual base
{
    explicit upper(const bsplineop& bop, const int nderiv);
};

/**
 * A building block to constrain the bulk value across
 * \f$y=\left[0,L_y\right]\f$.
 */
struct bulk : virtual base
{
    explicit bulk(bspline &b);
};

/** A building block to target some constant value. */
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

/** A building block to target some externally maintained reference. */
class reference : public virtual base
{
public:

    /**
     * Specify the reference to target.
     * The reference must have lifetime longer than this instance.
     */
    explicit reference(const real_t& target);

    virtual real_t target() const;

    virtual bool enabled() const;

private:

    /** Stores a pointer to the referenced value. */
    const real_t *t;

};

/**
 * A constraint reporting that it is always disabled.
 */
struct disabled : uniform
{
    explicit disabled(bspline &b);

    /** Always returns \c NaN. */
    virtual real_t target() const;

    /** Always returns \c false. */
    virtual bool enabled() const;
};

/**
 * Constrain the <tt>nderiv</tt>-th derivative at \f$y=0\f$ to be some given
 * constant using a uniform forcing profile.
 */
struct constant_lower : constant, lower, uniform
{
    constant_lower(const real_t target,
                   const bsplineop &bop,
                   const int nderiv);
};

/**
 * Constrain the <tt>nderiv</tt>-th derivative at \f$y=L_y\f$ to be some given
 * constant using a uniform forcing profile.
 */
struct constant_upper : constant, upper, uniform
{
    constant_upper(const real_t target,
                   const bsplineop &bop,
                   const int nderiv);
};

/**
 * Constrain the bulk value across \f$y=\left[0,L_y\right]\f$ to be some
 * constant using a uniform forcing profile.
 */
struct constant_bulk : constant, bulk, uniform
{
    constant_bulk(const real_t target, bspline &b);
};

/**
 * Constrain the <tt>nderiv</tt>-th derivative at \f$y=0\f$ to be some given
 * reference using a uniform forcing profile.
 */
struct reference_lower : reference, lower, uniform
{
    reference_lower(const real_t target,
                    const bsplineop &bop,
                    const int nderiv = 0);
};

/**
 * Constrain the <tt>nderiv</tt>-th derivative at \f$y=L_y\f$ to track some
 * referenced value using a uniform forcing profile.
 */
struct reference_upper : reference, upper, uniform
{
    reference_upper(const real_t target,
                    const bsplineop &bop,
                    const int nderiv = 0);
};

/**
 * Constrain the bulk value across \f$y=\left[0,L_y\right]\f$ to track some
 * referenced value using a uniform forcing profile.
 */
struct reference_bulk : reference, bulk, uniform
{
    reference_bulk(const real_t& target, bspline &b);
};

/**
 * Constrain the <tt>nderiv</tt>-th derivative at \f$y=0\f$ to be some target
 * value by manipulating the <tt>(nderiv+1)</tt>-th B-spline basis coefficient.
 * This approach modifies the <tt>nderiv</tt>-th and higher derivatives but
 * does not change the behavior of lower ones.
 */
struct constant_lower_derivative : constant, lower, coefficient
{
    constant_lower_derivative(
            const real_t target,
            const bsplineop &bop,
            const int nderiv);
};

/**
 * Constrain the <tt>nderiv</tt>-th derivative at \f$y=L_y\f$ to be some target
 * value by manipulating the <tt>(nderiv+1)</tt>-th B-spline basis coefficient.
 * \copydetails constant_lower_derivative
 */
struct constant_upper_derivative : constant, upper, coefficient
{
    constant_upper_derivative(
            const real_t target,
            const bsplineop &bop,
            const int nderiv);
};

} // namespace constraint

} // namespace suzerain

#endif /* SUZERAIN_PERFECT_CONSTRAINT_HPP */
