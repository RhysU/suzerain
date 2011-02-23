/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * grid_definition.hpp: classes handling grid definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_GRID_DEFINITION_HPP
#define __SUZERAIN_GRID_DEFINITION_HPP

#include <suzerain/common.hpp>
#include <suzerain/problem.hpp>
#include <suzerain/types.hpp>
#include <suzerain/validation.hpp>

/** @file
 * Provides classes handling problem grid definitions.
 */

namespace suzerain {

namespace problem {

/**
 * Holds basic three dimensional computational grid details for a distributed,
 * mixed Fourier/B-spline method.  The B-spline representation is used in the
 * wall-normal Y direction.
 */
template< typename FPT = double >
class GridDefinition : public IDefinition, public integral_types
{
public:
    /**
     * Construct an instance with the given default values.  Setting default
     * values of zero changes option semantics so the parameters become
     * optional.
     *
     * @param default_Nx   Default logical grid size in the X direction.
     * @param default_DAFx Default dealiasing factor in the X direction.
     * @param default_Ny   Default logical grid size in the Y direction.
     *                     This is the number of B-spline basis functions
     *                     (equivalently, wall-normal degrees of freedom)
     *                     to use.
     * @param default_k    Default uniform B-spline basis order plus one.
     *                     Piecewise cubics correspond to
     *                     <tt>default_k == 4</tt>.
     * @param default_Nz   Default logical grid size in the Z direction.
     * @param default_DAFz Default dealiasing factor in the Z direction.
     */
    explicit GridDefinition(
            size_type default_Nx,
            FPT       default_DAFx,
            size_type default_Ny,
            size_type default_k,
            size_type default_Nz,
            FPT       default_DAFz);

    /**
     * Retrieve global logical computational grid extents.  It does
     * not account for dealiasing factors.
     *
     * @return the global logical grid extents in the X, Y, and Z directions.
     */
    const size_type_3d& global_extents() const { return global_extents_; }

    /**
     * Retrieve global dealiased computational grid extents. These are the
     * global_extents() multiplied by the dealiasing factors DAFx(), one,
     * and DAFz().
     *
     * @return the global dealiased grid extents in the X, Y, and Z directions.
     */
    size_type_3d dealiased_extents() const {
        size_type_3d retval = global_extents_;
        retval[0] *= DAFx_;
        retval[2] *= DAFz_;
        return retval;
    }

    /**
     * Retrieve computational grid size in the X direction.  This is the number
     * of points in the domain without accounting for any dealiasing.
     *
     * @return the logical grid size in the X direction.
     */
    size_type Nx() const { return global_extents_[0]; }

    /** \copydoc Nx() const */
    size_type& Nx() { return global_extents_[0]; }

    /**
     * Retrieve the dealiasing factor for the X direction.  This factor
     * should be multiplied times Nx() to obtain an extent for Fourier
     * transformations.
     *
     * @return the dealiasing factor for the X direction.
     */
    FPT DAFx() const { return DAFx_; }

    /** \copydoc DAFx() const */
    FPT& DAFx() { return DAFx_; }

    /**
     * Retrieve computational grid size in the Y direction.  This is the number
     * of B-spline basis functions (equivalently, wall-normal degrees of
     * freedom) in use.  This direction does not support dealiasing.
     *
     * @return the logical grid size in the Y direction.
     */
    size_type Ny() const { return global_extents_[1]; }

    /** \copydoc Ny() const */
    size_type& Ny() { return global_extents_[1]; }

    /**
     * Retrieve the B-spline basis order plus one.  For example,
     * piecewise cubics have <tt>k() == 4</tt>.
     *
     * @return the B-spline basis order.
     */
    size_type k() const { return k_; }

    /** \copydoc k() const */
    size_type& k() { return k_; }

    /**
     * Retrieve grid size in the Z direction.  This is the number
     * of points in the domain without accounting for any dealiasing.
     *
     * @return the logical grid size in the Z direction.
     */
    size_type Nz() const { return global_extents_[2]; }

    /** \copydoc Nz() const */
    size_type& Nz() { return global_extents_[2]; }

    /**
     * Retrieve the dealiasing factor for the Z direction.  This factor
     * should be multiplied times Nz() to obtain an extent for Fourier
     * transformations.
     *
     * @return the dealiasing factor for the Z direction.
     */
    FPT DAFz() const { return DAFz_; }

    /** \copydoc DAFz() const */
    FPT& DAFz() { return DAFz_; }

    /**
     * Retrieve the two dimensional processor grid extents.
     *
     * In physical space, \f$ P_A \f$ is the grid extent in the Z direction
     * and \f$ P_B \f$ is the grid extent in the Y direction.
     * In wave space, \f$ P_A \f$ is the grid extent in the X direction
     * and \f$ P_B \f$  is the grid extent in the Z direction.
     *
     * @return the processor grid extents in the \f$ P_A \f$ and \f$ P_B \f$
     *         directions.
     */
    const size_type_2d& processor_grid() const { return processor_grid_; }

    /**
     * Retrieve the processor grid extent in the \f$ P_A \f$ direction.
     *
     * @return the processor grid extents in the \f$ P_A \f$ direction.
     * @see processor_grid() for more details.
     */
    size_type Pa() const { return processor_grid_[0]; }

    /** \copydoc Pa() const */
    size_type& Pa() { return processor_grid_[0]; }

    /**
     * Retrieve the processor grid extent in the \f$ P_B \f$ direction.
     *
     * @return the processor grid extents in the \f$ P_B \f$ direction.
     * @see processor_grid() for more details.
     */
    size_type Pb() const { return processor_grid_[1]; }

    /** \copydoc Pb() const */
    size_type& Pb() { return processor_grid_[1]; }

private:

    /** Stores the computational grid extents */
    size_type_3d global_extents_;

    FPT DAFx_;     /**< Stores the X direction dealiasing factor */
    size_type k_;  /**< Stores the B-spline basis order */
    FPT DAFz_;     /**< Stores the Z direction dealiasing factor */

    /** Stores the processor grid extents */
    size_type_2d processor_grid_;
};

template< typename FPT >
GridDefinition<FPT>::GridDefinition(size_type default_Nx,
                                    FPT       default_DAFx,
                                    size_type default_Ny,
                                    size_type default_k,
                                    size_type default_Nz,
                                    FPT       default_DAFz)
    : IDefinition("Mixed Fourier/B-spline computational grid definition")
{
    using ::std::auto_ptr;
    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;
    using ::boost::program_options::typed_value;
    using ::boost::program_options::value;

    // Created to solve ambiguous type issues below
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_positive_FPT(ensure_positive<FPT>);
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_nonnegative_FPT(ensure_nonnegative<FPT>);

    // Complicated add_options() calls done to allow changing the validation
    // routine in use when the default provided value is zero.  Zero is
    // generally used a NOP value by some client code.

    { // Nx
        auto_ptr<typed_value<size_type> > v(value(&global_extents_[0]));
        if (default_Nx) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),   "Nx"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"Nx"));
        }
        v->default_value(default_Nx);
        this->add_options()("Nx", v.release(),
                "Spectral coefficient count in streamwise X direction");
    }

    { // DAFx
        auto_ptr<typed_value<FPT> > v(value(&DAFx_));
        if (default_DAFx) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,   "DAFx"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT,"DAFx"));
        }
        v->default_value(default_DAFx);
        this->add_options()("DAFx", v.release(),
            "Dealiasing factor in streamwise X direction");
    }

    { // Ny
        auto_ptr<typed_value<size_type> > v(value(&global_extents_[1]));
        if (default_Ny) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),   "Ny"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"Ny"));
        }
        v->default_value(default_Ny);
        this->add_options()("Ny", v.release(),
                "Collocation point count in wall-normal Y direction");
    }

    { // k
        auto_ptr<typed_value<size_type> > v(value(&k_));
        if (default_k) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),   "k"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"k"));
        }
        v->default_value(default_k);
        this->add_options()("k", v.release(),
                "B-spline basis order where k = 4 indicates piecewise cubics");
    }

    { // Nz
        auto_ptr<typed_value<size_type> > v(value(&global_extents_[2]));
        if (default_Nz) {
            v->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),   "Nz"));
        } else {
            v->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"Nz"));
        }
        v->default_value(default_Nz);
        this->add_options()("Nz", v.release(),
                "Spectral coefficient count in spanwise Z direction");
    }

    { // DAFz
        auto_ptr<typed_value<FPT> > v(value(&DAFz_));
        if (default_DAFz) {
            v->notifier(bind2nd(ptr_fun_ensure_positive_FPT,   "DAFz"));
        } else {
            v->notifier(bind2nd(ptr_fun_ensure_nonnegative_FPT,"DAFz"));
        }
        v->default_value(default_DAFz);
        this->add_options()("DAFz", v.release(),
                "Dealiasing factor in spanwise Z direction");
    }

    { // Pa
        auto_ptr<typed_value<size_type> > v(value(&processor_grid_[0]));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"Pa"));
        v->default_value(0);
        this->add_options()("Pa", v.release(),
            "Processor count in the P_A decomposition direction");
    }

    { // Pb
        auto_ptr<typed_value<size_type> > v(value(&processor_grid_[1]));
        v->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"Pb"));
        v->default_value(0);
        this->add_options()("Pb", v.release(),
            "Processor count in the P_B decomposition direction");
    }
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_GRID_DEFINITION_HPP
