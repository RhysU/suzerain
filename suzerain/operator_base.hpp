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

#ifndef SUZERAIN_OPERATOR_BASE_HPP
#define SUZERAIN_OPERATOR_BASE_HPP

/** @file
 * Useful base classes for operator implementations.
 */

#include <suzerain/common.hpp>
#include <suzerain/operator_tools.hpp>

namespace suzerain {

/**
 * Provides potentially heavyweight B-spline and parallel FFT infrastructure.
 * Intended as a base class for nonlinear_operator implementations and other
 * logic requiring physical space coordinate information.
 */
class operator_base : public operator_tools
{
public:

    /**
     * @copydoc operator_tools::operator_tools
     *
     * @param b B-spline workspace for obtaining necessary details,
     *          e.g. integration coefficients.
     */
    operator_base(const specification_grid &grid,
                  const pencil_grid &dgrid,
                  const bsplineop &cop,
                  bspline &b);

    /** Virtual destructor to permit use as a base class */
    virtual ~operator_base();

    /**
     * Return the <tt>i</tt>th \c globally-indexed x grid point.
     */
    real_t x(std::size_t i) const
    {
        return i * grid.L.x() / grid.dN.x() - grid.L.x() / 2;
    }

    /**
     * Return the <tt>j</tt>th \c globally-indexed y grid point.
     * Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y()
     */
    real_t y(std::size_t j) const
    {
        return y_[j];
    }

    /**
     * Return the <tt>k</tt>th \c globally-indexed z grid point.
     */
    real_t z(std::size_t k) const
    {
        return k * grid.L.z() / grid.dN.z() - grid.L.z() / 2;
    }

    /**
     * Return the <tt>j</tt>th \c globally-indexed y grid spacing.
     * Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y()
     */
    real_t one_over_delta_y(std::size_t j) const
    {
        return one_over_delta_y_[j];
    }

    /**
     * Return the globally-indexed maximum pure imaginary eigenvalue estimate
     * for the first derivative operator in y at the <tt>j</tt>th collocation
     * point.  Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y().
     */
    real_t lambda1_y(std::size_t j) const
    {
        return lambda1_y_[j];
    }

    /**
     * Return the globally-indexed maximum pure real eigenvalue estimate
     * for the second derivative operator in y at the <tt>j</tt>th collocation
     * point.  Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y().
     */
    real_t lambda2_y(std::size_t j) const
    {
        return lambda2_y_[j];
    }

    /** Uniform grid spacing in x */
    const real_t one_over_delta_x;

    /** Maximum pure imaginary eigenvalue magnitude for first x derivative */
    const real_t lambda1_x;

    /** Maximum pure real eigenvalue magnitude for second x derivative */
    const real_t lambda2_x;

    /** Uniform grid spacing in z */
    const real_t one_over_delta_z;

    /** Maximum pure imaginary eigenvalue magnitude for first z derivative */
    const real_t lambda1_z;

    /** Maximum pure real eigenvalue magnitude for second z derivative */
    const real_t lambda2_z;

private:

    /** Stores y grid points on this rank in physical space */
    boost::multi_array<real_t,1> y_;

    /** Stores y grid spacing on this rank in physical space */
    boost::multi_array<real_t,1> one_over_delta_y_;

    /** Stores pure imaginary eigenvalue magnitudes for first y derivative */
    boost::multi_array<real_t,1> lambda1_y_;

    /** Stores pure real eigenvalue magnitudes for second y derivative */
    boost::multi_array<real_t,1> lambda2_y_;
};

} // namespace suzerain

#endif  /* SUZERAIN_OPERATOR_BASE_HPP */
