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
 * @todo Document this class
 */
template< typename FPT           = double,
          typename SizeType      = int >
class GridDefinition : public IDefinition
{
public:
    typedef SizeType size_type;
    typedef FPT      floating_point_type;

    GridDefinition(SizeType default_size = 16, FPT default_length = 2.0*M_PI);

    FPT Lx() const { return Lx_; }
    FPT Ly() const { return Ly_; }
    FPT Lz() const { return Lz_; }

    SizeType Nx() const { return Nx_; }
    SizeType Ny() const { return Ny_; }
    SizeType Nz() const { return Nz_; }

    /*! @copydoc IDefinition::options */
    const boost::program_options::options_description& options() {
        return options_;
    }

private:
    boost::program_options::options_description options_;
    FPT Lx_;
    FPT Ly_;
    FPT Lz_;
    SizeType Nx_;
    SizeType Ny_;
    SizeType Nz_;
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
