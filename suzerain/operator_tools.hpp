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

#ifndef SUZERAIN_OPERATOR_TOOLS_HPP
#define SUZERAIN_OPERATOR_TOOLS_HPP

/** @file
 * Provides combined, timed B-spline and parallel FFT infrastructure routines.
 */

#include <suzerain/common.hpp>
#include <suzerain/bspline.hpp>
#include <suzerain/diffwave.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/specification_grid.hpp>
#include <suzerain/timers.h>

namespace suzerain {

/**
 * Provides combined, timed B-spline and parallel FFT infrastructure routines.
 * Intended as lightweight, syntactic sugar atop coordinated use of several
 * other classes.
 */
class operator_tools
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
     * @param cop      B-spline operators to use.
     */
    operator_tools(const grid_specification &grid,
                   const pencil_grid &dgrid,
                   const bsplineop &cop)
        : grid(grid)
        , dgrid(dgrid)
        , cop(cop)
    {}

    /** Gain the benefits and curses of virtual behavior. */
    virtual ~operator_tools() {}

    /**
     * Perform scaled operator accumulation on two state fields.
     *
     * @see bsplineop::accumulate
     */
    template<typename AlphaType, typename MultiArrayX,
             typename BetaType,  typename MultiArrayY>
    int bop_accumulate(
            int nderiv,
            const AlphaType& alpha, const MultiArrayX &x, int ndx_x,
            const BetaType& beta,         MultiArrayY &y, int ndx_y) const
    {
        SUZERAIN_TIMER_SCOPED("operator_tools::bop_accumulate");

        assert(x.shape()[1] == (unsigned) cop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2] );
        assert((unsigned) y.strides()[3] == y.shape()[2] * y.strides()[2] );
        assert(std::equal(x.shape() + 1, x.shape() + 4, y.shape() + 1));

        return cop.accumulate(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx_x].origin(), x.strides()[1], x.strides()[2],
                beta,   y[ndx_y].origin(), y.strides()[1], y.strides()[2]);
    }

    /**
     * Perform scaled operator application on one state field.
     *
     * @see bsplineop::apply
     */
    template<typename AlphaType, typename MultiArray>
    int bop_apply(
            int nderiv, const AlphaType& alpha, MultiArray &x, int ndx) const
    {
        SUZERAIN_TIMER_SCOPED("operator_tools::bop_apply");

        assert(x.shape()[1] == (unsigned) cop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return cop.apply(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx].origin(), x.strides()[1], x.strides()[2]);
    }

    /**
     * Perform real-valued B-spline operator inversion on one state field.
     *
     * @see bsplineop_lu::solve
     */
    template<typename MultiArray>
    int bop_solve(
            const bsplineop_lu &lu, MultiArray &x, int ndx) const
    {
        SUZERAIN_TIMER_SCOPED("operator_tools::bop_solve(real)");

        assert(x.shape()[1] == (unsigned) lu.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return lu.solve(x.shape()[2]*x.shape()[3], x[ndx].origin(),
                        x.strides()[1], x.strides()[2]);
    }

    /**
     * Perform complex-valued B-spline operator inversion on one state field.
     *
     * @see bsplineop_luz::solve
     */
    template<typename MultiArray>
    int bop_solve(
            const bsplineop_luz &luz, MultiArray &x, int ndx) const
    {
        SUZERAIN_TIMER_SCOPED("operator_tools::bop_solve(cmplx)");

        assert(x.shape()[1] == (unsigned) luz.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return luz.solve(x.shape()[2]*x.shape()[3], x[ndx].origin(),
                         x.strides()[1], x.strides()[2]);
    }

    /**
     * Perform wave space-based differentiation using accumulation.
     *
     * @see diffwave::accumulate
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
        SUZERAIN_TIMER_SCOPED("operator_tools::diffwave_accumulate");

        assert(std::equal(x.shape()   + 1, x.shape()   + 4, y.shape()   + 1));
        assert(std::equal(x.strides() + 1, x.strides() + 4, y.strides() + 1));

        return diffwave::accumulate(
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
     * @see diffwave::apply
     */
    template<typename MultiArray>
    void diffwave_apply(int dxcnt,
                        int dzcnt,
                        const typename MultiArray::element& alpha,
                        MultiArray &x,
                        int ndx_x) const
    {
        SUZERAIN_TIMER_SCOPED("operator_tools::diffwave_apply");

        return diffwave::apply(
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
        SUZERAIN_TIMER_SCOPED("operator_tools::zero_dealiasing_modes");

        // Applying the dx^0 dz^0 operator with a scale factor of 1
        // zeros the dealiasing-only wavenumbers within the field
        return this->diffwave_apply(0, 0, 1.0, x, ndx_x);
    }

    /** The grid in which the operator is used */
    const grid_specification &grid;

    /** The parallel decomposition grid in which the operator is used */
    const pencil_grid &dgrid;

    /** The B-spline operators with which the operator is used */
    const bsplineop &cop;

    /** Retrieve a, possibly cached, factorized, real-valued mass matrix. */
    shared_ptr<const bsplineop_lu> masslu() const
    {
        if (SUZERAIN_UNLIKELY(!masslu_)) {
            bsplineop_lu *p = new bsplineop_lu(cop);
            SUZERAIN_ENSURE(0 == p->factor_mass(cop));
            masslu_.reset(p);
        }
        return masslu_;
    }

    /** Retrieve a, possibly cached, factorized, complex-valued mass matrix. */
    shared_ptr<const bsplineop_luz> massluz() const
    {
        if (SUZERAIN_UNLIKELY(!massluz_)) {
            bsplineop_luz *p = new bsplineop_luz(cop);
            SUZERAIN_ENSURE(0 == p->factor_mass(cop));
            massluz_.reset(p);
        }
        return massluz_;
    }

private:

    mutable shared_ptr<bsplineop_lu>  masslu_;

    mutable shared_ptr<bsplineop_luz> massluz_;

};

} // namespace suzerain

#endif  /* SUZERAIN_OPERATOR_TOOLS_HPP */
