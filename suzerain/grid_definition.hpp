/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
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
 * Provides classes handling grid definitions, which are runtime
 * arguments or scenario parameters used to perform calculations.
 */

namespace suzerain {

namespace problem {

/**
 * Holds basic three dimensional computational grid dimensions,
 * including the two dimensional processor grid definition.
 */
template< typename FPT = double >
class GridDefinition : public IDefinition, public integral_types
{
public:
    typedef FPT floating_point_type; /**< Floating point type in use */

    /**
     * Construct an instance with the given default size in each direction
     * and the given default length.
     *
     * @param default_size Default grid size in the X, Y, and Z directions.
     * @param default_length Default grid length in the X, Y, and Z directions.
     */
    GridDefinition(size_type default_size = 16, FPT default_length = 2.0*M_PI);

    /**
     * Retrieve the domain length in the X direction.
     *
     * @return the domain's X length.
     */
    FPT lx() const { return lx_; }

    /**
     * Retrieve the domain length in the Y direction.
     *
     * @return the domain's Y length.
     */
    FPT ly() const { return ly_; }

    /**
     * Retrieve the domain length in the Z direction.
     *
     * @return the domain's Z length.
     */
    FPT lz() const { return lz_; }

    /**
     * Retrieve global computational grid extents.
     *
     * @return the global grid extents in the X, Y, and Z directions.
     */
    const size_type_3d& global_extents() const { return global_extents_; }

    /**
     * Retrieve computational grid size in the X direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the X direction.
     */
    size_type nx() const { return global_extents_[0]; }

    /**
     * Retrieve computational grid size in the Y direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the Y direction.
     */
    size_type ny() const { return global_extents_[1]; }

    /**
     * Retrieve computational grid size in the Z direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the Z direction.
     */
    size_type nz() const { return global_extents_[2]; }

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
    size_type pa() const { return processor_grid_[0]; }

    /**
     * Retrieve the processor grid extent in the \f$ P_B \f$ direction.
     *
     * @return the processor grid extents in the \f$ P_B \f$ direction.
     * @see processor_grid() for more details.
     */
    size_type pb() const { return processor_grid_[1]; }

    /*! @copydoc IDefinition::options */
    const boost::program_options::options_description& options() {
        return options_;
    }

private:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;

    FPT lx_;  /**< Stores the X direction length */
    FPT ly_;  /**< Stores the Y direction length */
    FPT lz_;  /**< Stores the Z direction length */

    /** Stores the computational grid extents */
    size_type_3d global_extents_;

    /** Stores the processor grid extents */
    size_type_2d processor_grid_;
};

template< typename FPT >
GridDefinition<FPT>::GridDefinition(size_type default_size,
                                    FPT default_length)
    : options_("Grid definition"),
      lx_(default_length),
      ly_(default_length),
      lz_(default_length)
{
    std::fill(global_extents_.begin(), global_extents_.end(), default_size);
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
        ("lx", po::value<FPT>(&lx_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"lx"))
            ->default_value(default_length),
        "Nondimensional grid length in X (streamwise) direction")
        ("ly", po::value<FPT>(&ly_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"ly"))
            ->default_value(default_length),
        "Nondimensional grid length in Y (wall normal) direction")
        ("lz", po::value<FPT>(&lz_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"lz"))
            ->default_value(default_length),
        "Nondimensional grid length in Z (spanwise) direction")
        ("nx", po::value<size_type>(&global_extents_[0])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"nx"))
            ->default_value(default_size),
        "Number of grid points in X (streamwise) direction")
        ("ny", po::value<size_type>(&global_extents_[1])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"ny"))
            ->default_value(default_size),
        "Number of grid points in Y (wall normal) direction")
        ("nz", po::value<size_type>(&global_extents_[2])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"nz"))
            ->default_value(default_size),
        "Number of grid points in Z (spanwise) direction")
        ("pa", po::value<size_type>(&processor_grid_[0])
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"pa"))
            ->default_value(0),
        "Processor count in the P_A decomposition direction; 0 for automatic")
        ("pb", po::value<size_type>(&processor_grid_[1])
            ->notifier(bind2nd(ptr_fun(ensure_nonnegative<size_type>),"pb"))
            ->default_value(0),
        "Processor count in the P_B decomposition direction; 0 for automatic")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_GRID_DEFINITION_HPP
