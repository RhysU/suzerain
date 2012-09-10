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
#include <suzerain/scenario_definition.hpp>

// TODO OperatorBase should be templated on, e.g. ContiguousState

namespace suzerain {

/**
 * Provides common B-spline and parallel FFT infrastructure useful across
 * many ILinearOperator and INonlinearOperator implementations.
 *
 * @tparam FPT Floating point type to employ.
 */
template<typename FPT>
class OperatorBase
{
public:

    /**
     * Construct an instance based upon the provided scenario and numerics.
     * The instance makes the constant arguments available by reference
     * throughout its lifetime.  Accordingly, the constant reference arguments
     * should have a lifetime longer than that of this instance.
     *
     * @param scenario Scenario definition to store.
     * @param grid     Grid definition to store.
     * @param dgrid    Decomposition providing parallel grid details.
     * @param b        B-spline workspace for obtaining necessary details,
     *                 e.g. integration coefficients.
     * @param bop      B-spline operators to use.
     */
    OperatorBase(
            const typename suzerain::problem::ScenarioDefinition<FPT> &scenario,
            const suzerain::problem::GridDefinition &grid,
            const suzerain::pencil_grid &dgrid,
            suzerain::bspline &b,
            const suzerain::bsplineop &bop)
        : one_over_delta_x(grid.N.x() /* !dN.x() */ / scenario.Lx),
          lambda1_x(boost::math::constants::pi<FPT>() * one_over_delta_x),
          lambda2_x(lambda1_x * lambda1_x),
          one_over_delta_z(grid.N.z() /* !dN.z() */ / scenario.Lz),
          lambda1_z(boost::math::constants::pi<FPT>() * one_over_delta_z),
          lambda2_z(lambda1_z * lambda1_z),
          scenario(scenario),
          grid(grid),
          dgrid(dgrid),
          bop(bop),
          y_(boost::extents[boost::multi_array_types::extent_range(
                  dgrid.local_physical_start.y(),
                  dgrid.local_physical_end.y())]),
          one_over_delta_y_(
                  boost::extents[boost::multi_array_types::extent_range(
                        dgrid.local_physical_start.y(),
                        dgrid.local_physical_end.y())]),
          lambda1_y_(
                  boost::extents[boost::multi_array_types::extent_range(
                        dgrid.local_physical_start.y(),
                        dgrid.local_physical_end.y())]),
          lambda2_y_(
                  boost::extents[boost::multi_array_types::extent_range(
                        dgrid.local_physical_start.y(),
                        dgrid.local_physical_end.y())])
    {
        const FPT pi = boost::math::constants::pi<FPT>();

        // Compute the B-spline-dependent correction factor to obtain
        // good maximum eigenvalue estimates given wall-normal inhomogeneity.
        // See model document for definitions of C^{(1)} and C^{(2)}.
        double C1, Clow1, Chigh1;
        suzerain_bspline_htstretch2_evdeltascale(
                1, b.k(), grid.htdelta, b.n(), &C1, &Clow1, &Chigh1);
        double C2, Clow2, Chigh2;
        suzerain_bspline_htstretch2_evdeltascale(
                2, b.k(), grid.htdelta, b.n(), &C2, &Clow2, &Chigh2);

        // Compute collocation point-based information local to this rank
        for (int j = dgrid.local_physical_start.y();
             j < dgrid.local_physical_end.y();
             ++j) {

            // Collocation point locations
            y_[j] = b.collocation_point(j);

            // Inverse spacing to next adjacent collocation point
            one_over_delta_y_[j] = 1 / b.spacing_collocation_point(j);

            // Eigenvalue estimates using collocation point spacing
            const double inverse_spacing = one_over_delta_y_[j];

            // Estimating wall-normal first derivative eigenvalue magnitudes
            // TODO: Use estimates of C^{(1)} after working out safety factor
            lambda1_y_[j]  = pi * inverse_spacing / 4;

            // Estimating wall-normal second derivative eigenvalue magnitudes
            // TODO: Use estimates of C^{(2)} after working out safety factor
            lambda2_y_[j]  = pi * inverse_spacing / 1;
            lambda2_y_[j] *= lambda2_y_[j];
        }
    }

    /** Virtual destructor for peace of mind */
    virtual ~OperatorBase() {}

    /** Shorthand for scaled operator accumulation */
    template<typename AlphaType, typename MultiArrayX,
             typename BetaType,  typename MultiArrayY>
    int bop_accumulate(
            int nderiv,
            const AlphaType& alpha, const MultiArrayX &x, int ndx_x,
            const BetaType& beta,         MultiArrayY &y, int ndx_y) const
    {
#if defined(SUZERAIN_HAVE_GRVY) && defined(GRVY_LIB_VERSION)
        struct Guard {
            Guard()  { grvy_timer_begin("OperatorBase::bop_accumulate"); }
            ~Guard() { grvy_timer_end(  "OperatorBase::bop_accumulate"); }
        } scope_guard;
#endif
        assert(x.shape()[1] == (unsigned) bop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2] );
        assert((unsigned) y.strides()[3] == y.shape()[2] * y.strides()[2] );
        assert(std::equal(x.shape() + 1, x.shape() + 4, y.shape() + 1));

        return bop.accumulate(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx_x].origin(), x.strides()[1], x.strides()[2],
                beta,   y[ndx_y].origin(), y.strides()[1], y.strides()[2]);
    }

    /** Shorthand for scaled operator application */
    template<typename AlphaType, typename MultiArray>
    int bop_apply(
            int nderiv, const AlphaType& alpha, MultiArray &x, int ndx) const
    {
#if defined(SUZERAIN_HAVE_GRVY) && defined(GRVY_LIB_VERSION)
        struct Guard {
            Guard()  { grvy_timer_begin("OperatorBase::bop_apply"); }
            ~Guard() { grvy_timer_end(  "OperatorBase::bop_apply"); }
        } scope_guard;
#endif
        assert(x.shape()[1] == (unsigned) bop.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return bop.apply(
                nderiv, x.shape()[2] * x.shape()[3],
                alpha,  x[ndx].origin(), x.strides()[1], x.strides()[2]);
    }

    /** Shorthand for real-valued operator inversion */
    template<typename MultiArray>
    int bop_solve(
            const suzerain::bsplineop_lu &lu, MultiArray &x, int ndx) const
    {
#if defined(SUZERAIN_HAVE_GRVY) && defined(GRVY_LIB_VERSION)
        struct Guard {
            Guard()  { grvy_timer_begin("OperatorBase::bop_solve(real)"); }
            ~Guard() { grvy_timer_end(  "OperatorBase::bop_solve(real)"); }
        } scope_guard;
#endif
        assert(x.shape()[1] == (unsigned) lu.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return lu.solve(x.shape()[2]*x.shape()[3], x[ndx].origin(),
                        x.strides()[1], x.strides()[2]);
    }

    /** Shorthand for complex-valued operator inversion */
    template<typename MultiArray>
    int bop_solve(
            const suzerain::bsplineop_luz &luz, MultiArray &x, int ndx) const
    {
#if defined(SUZERAIN_HAVE_GRVY) && defined(GRVY_LIB_VERSION)
        struct Guard {
            Guard()  { grvy_timer_begin("OperatorBase::bop_solve(cmplx)"); }
            ~Guard() { grvy_timer_end(  "OperatorBase::bop_solve(cmplx)"); }
        } scope_guard;
#endif
        assert(x.shape()[1] == (unsigned) luz.n());
        assert((unsigned) x.strides()[3] == x.shape()[2] * x.strides()[2]);

        return luz.solve(x.shape()[2]*x.shape()[3], x[ndx].origin(),
                         x.strides()[1], x.strides()[2]);
    }

    /** Shorthand for wave space-based differentiation accumulation */
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
#if defined(SUZERAIN_HAVE_GRVY) && defined(GRVY_LIB_VERSION)
        struct Guard {
            Guard()  { grvy_timer_begin("OperatorBase::diffwave_accumulate"); }
            ~Guard() { grvy_timer_end(  "OperatorBase::diffwave_accumulate"); }
        } scope_guard;
#endif
        assert(std::equal(x.shape()   + 1, x.shape()   + 4, y.shape()   + 1));
        assert(std::equal(x.strides() + 1, x.strides() + 4, y.strides() + 1));

        return suzerain::diffwave::accumulate(
                dxcnt, dzcnt,
                alpha, x[ndx_x].origin(),
                beta,  y[ndx_y].origin(),
                scenario.Lx, scenario.Lz,
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

    /** Shorthand for wave space-based differentiation application */
    template<typename MultiArray>
    void diffwave_apply(int dxcnt,
                        int dzcnt,
                        const typename MultiArray::element& alpha,
                        MultiArray &x,
                        int ndx_x) const
    {
#if defined(SUZERAIN_HAVE_GRVY) && defined(GRVY_LIB_VERSION)
        struct Guard {
            Guard()  { grvy_timer_begin("OperatorBase::diffwave_apply"); }
            ~Guard() { grvy_timer_end(  "OperatorBase::diffwave_apply"); }
        } scope_guard;
#endif
        return suzerain::diffwave::apply(
                dxcnt, dzcnt,
                alpha, x[ndx_x].origin(),
                scenario.Lx, scenario.Lz,
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

    /** Shorthand for zeroing modes present only for dealiasing */
    template<typename MultiArray>
    void zero_dealiasing_modes(
                        MultiArray &x,
                        int ndx_x) const
    {
#if defined(SUZERAIN_HAVE_GRVY) && defined(GRVY_LIB_VERSION)
        struct Guard {
            Guard()  { grvy_timer_begin("OperatorBase::zero_dealiasing_modes"); }
            ~Guard() { grvy_timer_end(  "OperatorBase::zero_dealiasing_modes"); }
        } scope_guard;
#endif
        // Applying the dx^0 dz^0 operator with a scale factor of 1
        // zeros the dealiasing-only wavenumbers within the field
        return this->diffwave_apply(0, 0, 1.0, x, ndx_x);
    }

    /**
     * Return the <tt>i</tt>th \c globally-indexed x grid point.
     */
    FPT x(std::size_t i) const {
        return i * scenario.Lx / grid.dN.x() - scenario.Lx / 2;
    }

    /**
     * Return the <tt>j</tt>th \c globally-indexed y grid point.
     * Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y()
     */
    FPT y(std::size_t j) const {
        return y_[j];
    }

    /**
     * Return the <tt>k</tt>th \c globally-indexed z grid point.
     */
    FPT z(std::size_t k) const {
        return k * scenario.Lz / grid.dN.z() - scenario.Lz / 2;
    }

    /**
     * Return the <tt>j</tt>th \c globally-indexed y grid spacing.
     * Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y()
     */
    FPT one_over_delta_y(std::size_t j) const {
        return one_over_delta_y_[j];
    }

    /**
     * Return the globally-indexed maximum pure imaginary eigenvalue estimate
     * for the first derivative operator in y at the <tt>j</tt>th collocation
     * point.  Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y().
     */
    FPT lambda1_y(std::size_t j) const {
        return lambda1_y_[j];
    }

    /**
     * Return the globally-indexed maximum pure real eigenvalue estimate
     * for the second derivative operator in y at the <tt>j</tt>th collocation
     * point.  Only valid for j \f$\in\f$ dgrid.local_physical_{start,end}.y().
     */
    FPT lambda2_y(std::size_t j) const {
        return lambda2_y_[j];
    }

    /** Uniform grid spacing in x */
    const FPT one_over_delta_x;

    /** Maximum pure imaginary eigenvalue magnitude for first derivatives in x */
    const FPT lambda1_x;

    /** Maximum pure real eigenvalue magnitude for second derivatives in x */
    const FPT lambda2_x;

    /** Uniform grid spacing in z */
    const FPT one_over_delta_z;

    /** Maximum pure imaginary eigenvalue magnitude for first derivatives in z */
    const FPT lambda1_z;

    /** Maximum pure real eigenvalue magnitude for second derivatives in z */
    const FPT lambda2_z;

    /** The scenario in which the operator is used */
    const typename suzerain::problem::ScenarioDefinition<FPT> &scenario;

    /** The grid in which the operator is used */
    const suzerain::problem::GridDefinition &grid;

    /** The parallel decomposition grid in which the operator is used */
    const suzerain::pencil_grid &dgrid;

    /** The B-spline operators with which the operator is used */
    const suzerain::bsplineop &bop;

private:

    /** Stores y grid points on this rank in wave space */
    boost::multi_array<FPT,1> y_;

    /** Stores y grid spacing on this rank in wave space */
    boost::multi_array<FPT,1> one_over_delta_y_;

    /** Stores pure imaginary eigenvalue magnitudes for y first derivatives */
    boost::multi_array<FPT,1> lambda1_y_;

    /** Stores pure real eigenvalue magnitudes for y second derivatives */
    boost::multi_array<FPT,1> lambda2_y_;

    // Noncopyable
    OperatorBase(const OperatorBase&);
    OperatorBase& operator=(const OperatorBase&);
};

} // namespace suzerain

#endif  /* SUZERAIN_OPERATOR_BASE_HPP */
