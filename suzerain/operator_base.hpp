//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
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
// operator_base.hpp: useful base class for operator implementation work
// $Id$

#ifndef SUZERAIN_OPERATOR_BASE_HPP
#define SUZERAIN_OPERATOR_BASE_HPP

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/grid_definition.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/timers.h>

namespace suzerain {

/**
 * Provides common double-precision B-spline and parallel FFT infrastructure
 * useful across many ILinearOperator and INonlinearOperator implementations.
 */
class OperatorBase
{
public:

    /**
     * Construct an instance based upon the provided numeric utilities.  The
     * instance makes the constant arguments available by reference throughout
     * its lifetime.  Accordingly, the constant reference arguments should have
     * a lifetime longer than that of this instance.
     *
     * @param grid     Grid definition to store.
     * @param dgrid    Decomposition providing parallel grid details.
     * @param b        B-spline workspace for obtaining necessary details,
     *                 e.g. integration coefficients.
     * @param bop      B-spline operators to use.
     */
    OperatorBase(const suzerain::problem::GridDefinition &grid,
                 const suzerain::pencil_grid &dgrid,
                 suzerain::bspline &b,
                 const suzerain::bsplineop &bop);

    /** Virtual destructor to permit use as a base class */
    virtual ~OperatorBase();

    /**
     * Perform scaled operator accumulation on two state fields.
     *
     * @see suzerain::bsplineop::accumulate
     */
    template<typename AlphaType, typename MultiArrayX,
             typename BetaType,  typename MultiArrayY>
    int bop_accumulate(
            int nderiv,
            const AlphaType& alpha, const MultiArrayX &x, int ndx_x,
            const BetaType& beta,         MultiArrayY &y, int ndx_y) const
    {
        struct Guard {
            Guard()  { SUZERAIN_TIMER_BEGIN("OperatorBase::bop_accumulate"); }
            ~Guard() { SUZERAIN_TIMER_END  ("OperatorBase::bop_accumulate"); }
        } scope_guard;

        assert(x.shape()[1] == (unsigned) bop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2] );
        assert((unsigned) y.strides()[3] == y.shape()[2] * y.strides()[2] );
        assert(std::equal(x.shape() + 1, x.shape() + 4, y.shape() + 1));

        return bop.accumulate(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx_x].origin(), x.strides()[1], x.strides()[2],
                beta,   y[ndx_y].origin(), y.strides()[1], y.strides()[2]);
    }

    /**
     * Perform scaled operator application on one state field.
     *
     * @see suzerain::bsplineop::apply
     */
    template<typename AlphaType, typename MultiArray>
    int bop_apply(
            int nderiv, const AlphaType& alpha, MultiArray &x, int ndx) const
    {
        struct Guard {
            Guard()  { SUZERAIN_TIMER_BEGIN("OperatorBase::bop_apply"); }
            ~Guard() { SUZERAIN_TIMER_END  ("OperatorBase::bop_apply"); }
        } scope_guard;

        assert(x.shape()[1] == (unsigned) bop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return bop.apply(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx].origin(), x.strides()[1], x.strides()[2]);
    }

    /**
     * Perform real-valued B-spline operator inversion on one state field.
     *
     * @see suzerain::bsplineop_lu::solve
     */
    template<typename MultiArray>
    int bop_solve(
            const suzerain::bsplineop_lu &lu, MultiArray &x, int ndx) const
    {
        struct Guard {
            Guard()  { SUZERAIN_TIMER_BEGIN("OperatorBase::bop_solve(real)"); }
            ~Guard() { SUZERAIN_TIMER_END  ("OperatorBase::bop_solve(real)"); }
        } scope_guard;

        assert(x.shape()[1] == (unsigned) lu.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return lu.solve(x.shape()[2]*x.shape()[3], x[ndx].origin(),
                        x.strides()[1], x.strides()[2]);
    }

    /**
     * Perform complex-valued B-spline operator inversion on one state field.
     *
     * @see suzerain::bsplineop_luz::solve
     */
    template<typename MultiArray>
    int bop_solve(
            const suzerain::bsplineop_luz &luz, MultiArray &x, int ndx) const
    {
        struct Guard {
            Guard()  { SUZERAIN_TIMER_BEGIN("OperatorBase::bop_solve(cmplx)"); }
            ~Guard() { SUZERAIN_TIMER_END  ("OperatorBase::bop_solve(cmplx)"); }
        } scope_guard;

        assert(x.shape()[1] == (unsigned) luz.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return luz.solve(x.shape()[2]*x.shape()[3], x[ndx].origin(),
                         x.strides()[1], x.strides()[2]);
    }

    /**
     * Perform wave space-based differentiation using accumulation.
     *
     * @see suzerain::diffwave::accumulate
     */
    template<typename MultiArrayX, typename MultiArrayY>
    void diffwave_accumulate(int dxcnt,
                             int dzcnt,
                             const typename MultiArrayX::element& alpha,
                             const MultiArrayX &x,
                             int ndx_x,
                             const typename MultiArrayY::element& beta,
                             MultiArrayY &y,
                             int ndx_y) const
    {
        struct Guard {
            Guard()
            {SUZERAIN_TIMER_BEGIN( "OperatorBase::diffwave_accumulate"); }

            ~Guard()
            { SUZERAIN_TIMER_END  ("OperatorBase::diffwave_accumulate"); }
        } scope_guard;

        assert(std::equal(x.shape()   + 1, x.shape()   + 4, y.shape()   + 1));
        assert(std::equal(x.strides() + 1, x.strides() + 4, y.strides() + 1));

        return suzerain::diffwave::accumulate(
                dxcnt, dzcnt,
                alpha, x[ndx_x].origin(),
                beta,  y[ndx_y].origin(),
                grid.L.x(), grid.L.z(),
                dgrid.global_wave_extent.y(),
                grid.N.x(),
                grid.dN.x(),
                dgrid.local_wave_start.x(),
                dgrid.local_wave_end.x(),
                grid.N.z(),
                grid.dN.z(),
                dgrid.local_wave_start.z(),
                dgrid.local_wave_end.z());

    }

    /**
     * Perform wave space-based differentiation using application.
     *
     * @see suzerain::diffwave::apply
     */
    template<typename MultiArray>
    void diffwave_apply(int dxcnt,
                        int dzcnt,
                        const typename MultiArray::element& alpha,
                        MultiArray &x,
                        int ndx_x) const
    {
        struct Guard {
            Guard()  { SUZERAIN_TIMER_BEGIN("OperatorBase::diffwave_apply"); }
            ~Guard() { SUZERAIN_TIMER_END  ("OperatorBase::diffwave_apply"); }
        } scope_guard;

        return suzerain::diffwave::apply(
                dxcnt, dzcnt,
                alpha, x[ndx_x].origin(),
                grid.L.x(), grid.L.z(),
                dgrid.global_wave_extent.y(),
                grid.N.x(),
                grid.dN.x(),
                dgrid.local_wave_start.x(),
                dgrid.local_wave_end.x(),
                grid.N.z(),
                grid.dN.z(),
                dgrid.local_wave_start.z(),
                dgrid.local_wave_end.z());
    }

    /** Zero wave-space modes present only for dealiasing purposes */
    template<typename MultiArray>
    void zero_dealiasing_modes(MultiArray &x,
                               int ndx_x) const
    {
        struct Guard {
            Guard()
            { SUZERAIN_TIMER_BEGIN("OperatorBase::zero_dealiasing_modes"); }
            ~Guard()
            { SUZERAIN_TIMER_END  ("OperatorBase::zero_dealiasing_modes"); }
        } scope_guard;

        // Applying the dx^0 dz^0 operator with a scale factor of 1
        // zeros the dealiasing-only wavenumbers within the field
        return this->diffwave_apply(0, 0, 1.0, x, ndx_x);
    }

    /**
     * Return the <tt>i</tt>th \c globally-indexed x grid point.
     */
    double x(std::size_t i) const
    {
        return i * grid.L.x() / grid.dN.x() - grid.L.x() / 2;
    }

    /**
     * Return the <tt>j</tt>th \c globally-indexed y grid point.
     * Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y()
     */
    double y(std::size_t j) const
    {
        return y_[j];
    }

    /**
     * Return the <tt>k</tt>th \c globally-indexed z grid point.
     */
    double z(std::size_t k) const
    {
        return k * grid.L.z() / grid.dN.z() - grid.L.z() / 2;
    }

    /**
     * Return the <tt>j</tt>th \c globally-indexed y grid spacing.
     * Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y()
     */
    double one_over_delta_y(std::size_t j) const
    {
        return one_over_delta_y_[j];
    }

    /**
     * Return the globally-indexed maximum pure imaginary eigenvalue estimate
     * for the first derivative operator in y at the <tt>j</tt>th collocation
     * point.  Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y().
     */
    double lambda1_y(std::size_t j) const
    {
        return lambda1_y_[j];
    }

    /**
     * Return the globally-indexed maximum pure real eigenvalue estimate
     * for the second derivative operator in y at the <tt>j</tt>th collocation
     * point.  Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y().
     */
    double lambda2_y(std::size_t j) const
    {
        return lambda2_y_[j];
    }

    /** Uniform grid spacing in x */
    const double one_over_delta_x;

    /** Maximum pure imaginary eigenvalue magnitude for first derivative in x */
    const double lambda1_x;

    /** Maximum pure real eigenvalue magnitude for second derivatives in x */
    const double lambda2_x;

    /** Uniform grid spacing in z */
    const double one_over_delta_z;

    /** Maximum pure imaginary eigenvalue magnitude for first derivative in z */
    const double lambda1_z;

    /** Maximum pure real eigenvalue magnitude for second derivatives in z */
    const double lambda2_z;

    /** The grid in which the operator is used */
    const suzerain::problem::GridDefinition &grid;

    /** The parallel decomposition grid in which the operator is used */
    const suzerain::pencil_grid &dgrid;

    /** The B-spline operators with which the operator is used */
    const suzerain::bsplineop &bop;

private:

    /** Stores y grid points on this rank in wave space */
    boost::multi_array<double,1> y_;

    /** Stores y grid spacing on this rank in wave space */
    boost::multi_array<double,1> one_over_delta_y_;

    /** Stores pure imaginary eigenvalue magnitudes for y first derivatives */
    boost::multi_array<double,1> lambda1_y_;

    /** Stores pure real eigenvalue magnitudes for y second derivatives */
    boost::multi_array<double,1> lambda2_y_;

    // Noncopyable
    OperatorBase(const OperatorBase&);
    OperatorBase& operator=(const OperatorBase&);
};

} // namespace suzerain

#endif  /* SUZERAIN_OPERATOR_BASE_HPP */
