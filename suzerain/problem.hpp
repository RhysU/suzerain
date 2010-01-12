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
 * problem.hpp: classes handling problem definitions
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_PROBLEM_HPP
#define __SUZERAIN_PROBLEM_HPP

#include <suzerain/common.hpp>
#include <suzerain/types.hpp>

/** @file
 * Provides classes handling problem definitions, which are runtime
 * arguments or scenario parameters used to perform calculations.
 */

namespace suzerain {

/**
 * Provides classes handling problem definitions, which are runtime
 * arguments or scenario parameters used to perform calculations.
 */
namespace problem {

/**
 * An interface describing a problem definition.
 */
class IDefinition
{
public:
    /** Constructor appropriate for an abstract base class */
    IDefinition() {};

    /** Virtual destructor appropriate for an abstract base class */
    virtual ~IDefinition() {};

    /**
     * Obtain a Boost options_description encompassing all information
     * in the definition.
     *
     * @return A reference suitable for <tt>add</tt>-ing to a
     *         <tt>boost::program_options::options_description</tt> instance.
     *
     * @see <a href="http://www.boost.org/doc/libs/release/libs/program_options">
     *      Boost.Program_options</a> for more information.
     */
    virtual const boost::program_options::options_description& options() = 0;
};

namespace detail {

template< typename T >
void ensure_positive(const T &t)
throw(std::invalid_argument)
{
    BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);
    if (boost::is_signed<T>::value && t <= T(0)) {
        std::ostringstream msg;
        msg << "Value " << t << " is non-positive";
        throw std::invalid_argument(msg.str());
    }
}

} // namespace detail

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
     * In physical space, \f$ P_0 \f$ is the grid extent in the Z direction
     * and \f$ P_1 \f$ is the grid extent in the Y direction.
     * In wave space, \f$ P_0 \f$ is the grid extent in the X direction
     * and \f$ P_1 \f$  is the grid extent in the Z direction.
     *
     * @return the processor grid extents in the \f$ P_0 \f$ and \f$ P_1 \f$
     *         directions.
     */
    const size_type_2d& processor_grid() const { return processor_grid_; }

    /**
     * Retrieve the processor grid extent in the \f$ P_0 \f$ direction.
     *
     * @return the processor grid extents in the \f$ P_0 \f$ direction.
     * @see processor_grid() for more details.
     */
    size_type pg0() const { return processor_grid_[0]; }

    /**
     * Retrieve the processor grid extent in the \f$ P_1 \f$ direction.
     *
     * @return the processor grid extents in the \f$ P_1 \f$ direction.
     * @see processor_grid() for more details.
     */
    size_type pg1() const { return processor_grid_[1]; }

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

    options_.add_options()
        ("lx",
         po::value<FPT>(&lx_)
                ->notifier(detail::ensure_positive<FPT>)
                ->default_value(default_length),
        "Nondimensional grid length in X (streamwise) direction")
        ("ly",
         po::value<FPT>(&ly_)
                ->notifier(detail::ensure_positive<FPT>)
                ->default_value(default_length),
        "Nondimensional grid length in Y (wall normal) direction")
        ("lz",
         po::value<FPT>(&lz_)
                ->notifier(detail::ensure_positive<FPT>)
                ->default_value(default_length),
        "Nondimensional grid length in Z (spanwise) direction")
        ("nx",
         po::value<size_type>(&global_extents_[0])
                ->notifier(detail::ensure_positive<size_type>)
                ->default_value(default_size),
        "Number of grid points in X (streamwise) direction")
        ("ny",
         po::value<size_type>(&global_extents_[1])
                ->notifier(detail::ensure_positive<size_type>)
                ->default_value(default_size),
        "Number of grid points in Y (wall normal) direction")
        ("nz",
         po::value<size_type>(&global_extents_[2])
                ->notifier(detail::ensure_positive<size_type>)
                ->default_value(default_size),
        "Number of grid points in Z (spanwise) direction")
        ("pg0",
         po::value<size_type>(&processor_grid_[0])
                ->notifier(detail::ensure_positive<size_type>)
                ->default_value(0),
        "Processor count in the P_0 decomposition direction; 0 for automatic")
        ("pg1",
         po::value<size_type>(&processor_grid_[1])
                ->notifier(detail::ensure_positive<size_type>)
                ->default_value(0),
        "Processor count in the P_1 decomposition direction; 0 for automatic")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_PROBLEM_HPP
