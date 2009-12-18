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
 * Holds basic three dimensional computational grid dimensions.
 */
template< typename FPT           = double,
          typename SizeType      = int >
class GridDefinition : public IDefinition
{
public:
    typedef FPT      floating_point_type; /**< Floating point type in use */
    typedef SizeType size_type;           /**< Integral size type in use */

    /**
     * Construct an instance with the given default size in each direction
     * and the given default length.
     *
     * @param default_size Default grid size in the X, Y, and Z directions.
     * @param default_length Default grid length in the X, Y, and Z directions.
     */
    GridDefinition(SizeType default_size = 16, FPT default_length = 2.0*M_PI);

    /**
     * Retrieve the domain length in the X direction.
     *
     * @return the domain's X length.
     */
    FPT Lx() const { return Lx_; }

    /**
     * Retrieve the domain length in the Y direction.
     *
     * @return the domain's Y length.
     */
    FPT Ly() const { return Ly_; }

    /**
     * Retrieve the domain length in the Z direction.
     *
     * @return the domain's Z length.
     */
    FPT Lz() const { return Lz_; }

    /**
     * Retrieve computational grid size in the X direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the X direction.
     */
    SizeType Nx() const { return Nx_; }

    /**
     * Retrieve computational grid size in the Y direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the Y direction.
     */
    SizeType Ny() const { return Ny_; }

    /**
     * Retrieve computational grid size in the Z direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the Z direction.
     */
    SizeType Nz() const { return Nz_; }

    /*! @copydoc IDefinition::options */
    const boost::program_options::options_description& options() {
        return options_;
    }

private:

    /** Stores the program options processing information */
    boost::program_options::options_description options_;

    FPT Lx_;        /**< Stores the X direction length */
    FPT Ly_;        /**< Stores the Y direction length */
    FPT Lz_;        /**< Stores the Z direction length */

    SizeType Nx_;   /**< Stores the X direction grid size */
    SizeType Ny_;   /**< Stores the Y direction grid size */
    SizeType Nz_;   /**< Stores the Z direction grid size */
};

template< typename FPT, typename SizeType >
GridDefinition<FPT,SizeType>::GridDefinition(SizeType default_size,
                                             FPT default_length)
    : options_("Grid definition"),
      Lx_(default_length),
      Ly_(default_length),
      Lz_(default_length),
      Nx_(default_size),
      Ny_(default_size),
      Nz_(default_size)
{
    namespace po = ::boost::program_options;

    options_.add_options()
        ("Lx",
         po::value<FPT>(&Lx_)
                ->notifier(detail::ensure_positive<FPT>)
                ->default_value(default_length),
        "Nondimensional grid length in X (streamwise) direction")
        ("Ly",
         po::value<FPT>(&Ly_)
                ->notifier(detail::ensure_positive<FPT>)
                ->default_value(default_length),
        "Nondimensional grid length in Y (wall normal) direction")
        ("Lz",
         po::value<FPT>(&Lz_)
                ->notifier(detail::ensure_positive<FPT>)
                ->default_value(default_length),
        "Nondimensional grid length in Z (spanwise) direction")
        ("Nx",
         po::value<SizeType>(&Nx_)
                ->notifier(detail::ensure_positive<SizeType>)
                ->default_value(default_size),
        "Number of grid points in X (streamwise) direction")
        ("Ny",
         po::value<SizeType>(&Ny_)
                ->notifier(detail::ensure_positive<SizeType>)
                ->default_value(default_size),
        "Number of grid points in Y (wall normal) direction")
        ("Nz",
         po::value<SizeType>(&Nz_)
                ->notifier(detail::ensure_positive<SizeType>)
                ->default_value(default_size),
        "Number of grid points in Z (spanwise) direction")
    ;
}

} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_PROBLEM_HPP
