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
     * Construct an instance with the given default values.
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

    /**
     * Retrieve the dealiasing factor for the X direction.  This factor
     * should be multiplied times Nx() to obtain an extent for Fourier
     * transformations.
     *
     * @return the dealiasing factor for the X direction.
     */
    FPT DAFx() const { return DAFx_; }

    /**
     * Retrieve computational grid size in the Y direction.  This is the number
     * of B-spline basis functions (equivalently, wall-normal degrees of
     * freedom) in use.  This direction does not support dealiasing.
     *
     * @return the logical grid size in the Y direction.
     */
    size_type Ny() const { return global_extents_[1]; }

    /**
     * Retrieve the B-spline basis order plus one.  For example,
     * piecewise cubics have <tt>k() == 4</tt>.
     *
     * @return the B-spline basis order.
     */
    size_type k() const { return k_; }

    /**
     * Retrieve grid size in the Z direction.  This is the number
     * of points in the domain without accounting for any dealiasing.
     *
     * @return the logical grid size in the Z direction.
     */
    size_type Nz() const { return global_extents_[2]; }

    /**
     * Retrieve the dealiasing factor for the Z direction.  This factor
     * should be multiplied times Nz() to obtain an extent for Fourier
     * transformations.
     *
     * @return the dealiasing factor for the Z direction.
     */
    FPT DAFz() const { return DAFz_; }


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

    /**
     * Retrieve the processor grid extent in the \f$ P_B \f$ direction.
     *
     * @return the processor grid extents in the \f$ P_B \f$ direction.
     * @see processor_grid() for more details.
     */
    size_type Pb() const { return processor_grid_[1]; }

    /*! @copydoc IDefinition::options */
    const boost::program_options::options_description& options() {
        return options_;
    }

private:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;

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
    : options_("Mixed Fourier/B-spline grid definition"),
      DAFx_(default_DAFx),
      k_(default_k),
      DAFz_(default_DAFz)
{
    global_extents_[0] = default_Nx;
    global_extents_[1] = default_Ny;
    global_extents_[2] = default_Nz;
    std::fill(processor_grid_.begin(), processor_grid_.end(), 0);

    namespace po = ::boost::program_options;

    using ::std::bind2nd;
    using ::std::ptr_fun;
    using ::suzerain::validation::ensure_nonnegative;
    using ::suzerain::validation::ensure_positive;

    // Created to solve ambiguous type issues below
    ::std::pointer_to_binary_function<FPT,const char*,void>
        ptr_fun_ensure_positive_FPT(ensure_positive<FPT>);

    options_.add_options()
        ("Nx", po::value<size_type>(&global_extents_[0])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Nx"))
            ->default_value(default_Nx),
        "Spectral coefficient count in streamwise X direction")
        ("DAFx", po::value<FPT>(&DAFx_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"DAFx"))
            ->default_value(default_DAFx),
        "Dealiasing factor in streamwise X direction")
        ("Ny", po::value<size_type>(&global_extents_[1])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Ny"))
            ->default_value(default_Ny),
        "Collocation point count in wall-normal Y direction")
        ("k", po::value<size_type>(&k_)
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"k"))
            ->default_value(default_k),
        "B-spline basis order; k = 4 indicates piecewise cubics")
        ("Nz", po::value<size_type>(&global_extents_[2])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Nz"))
            ->default_value(default_Nz),
        "Spectral coefficient count in spanwise Z direction")
        ("DAFz", po::value<FPT>(&DAFz_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"DAFz"))
            ->default_value(default_DAFz),
        "Dealiasing factor in spanwise Z direction")
        ("Pa", po::value<size_type>(&processor_grid_[0])
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"Pa"))
            ->default_value(0),
        "Processor count in the P_A decomposition direction; 0 for automatic")
        ("Pb", po::value<size_type>(&processor_grid_[1])
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"Pb"))
            ->default_value(0),
        "Processor count in the P_B decomposition direction; 0 for automatic")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_GRID_DEFINITION_HPP
