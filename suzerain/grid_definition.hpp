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
    /**
     * Construct an instance with the given default size in each direction
     * and the given default length.
     *
     * @param default_Nx Default grid size in the X direction.
     * @param default_Ny Default grid size in the Y direction.
     * @param default_Nz Default grid size in the Z direction.
     * @param default_Lx Default domain length in the X direction.
     * @param default_Ly Default domain length in the Y direction.
     * @param default_Lz Default domain length in the Z direction.
     */
    GridDefinition(size_type default_Nx = 16,
                   size_type default_Ny = 16,
                   size_type default_Nz = 16,
                   FPT default_Lx = 2.0*M_PI,
                   FPT default_Ly = 2.0*M_PI,
                   FPT default_Lz = 2.0*M_PI);

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
    size_type Nx() const { return global_extents_[0]; }

    /**
     * Retrieve computational grid size in the Y direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the Y direction.
     */
    size_type Ny() const { return global_extents_[1]; }

    /**
     * Retrieve computational grid size in the Z direction.
     * This is the number of points in the domain.
     *
     * @return the grid size in the Z direction.
     */
    size_type Nz() const { return global_extents_[2]; }

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

    FPT Lx_;  /**< Stores the X direction length */
    FPT Ly_;  /**< Stores the Y direction length */
    FPT Lz_;  /**< Stores the Z direction length */

    /** Stores the computational grid extents */
    size_type_3d global_extents_;

    /** Stores the processor grid extents */
    size_type_2d processor_grid_;
};

template< typename FPT >
GridDefinition<FPT>::GridDefinition(size_type default_Nx,
                                    size_type default_Ny,
                                    size_type default_Nz,
                                    FPT default_Lx,
                                    FPT default_Ly,
                                    FPT default_Lz)
    : options_("Grid definition"),
      Lx_(default_Lx),
      Ly_(default_Ly),
      Lz_(default_Lz)
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
        ("Lx", po::value<FPT>(&Lx_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"Lx"))
            ->default_value(default_Lx),
        "Nondimensional grid length in X (streamwise) direction")
        ("Ly", po::value<FPT>(&Ly_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"Ly"))
            ->default_value(default_Ly),
        "Nondimensional grid length in Y (wall normal) direction")
        ("Lz", po::value<FPT>(&Lz_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"Lz"))
            ->default_value(default_Lz),
        "Nondimensional grid length in Z (spanwise) direction")
        ("Nx", po::value<size_type>(&global_extents_[0])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Nx"))
            ->default_value(default_Nx),
        "Number of grid points in X (streamwise) direction")
        ("Ny", po::value<size_type>(&global_extents_[1])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Ny"))
            ->default_value(default_Ny),
        "Number of grid points in Y (wall normal) direction")
        ("Nz", po::value<size_type>(&global_extents_[2])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Nz"))
            ->default_value(default_Nz),
        "Number of grid points in Z (spanwise) direction")
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

/**
 * Holds three dimensional computational grid dimensions for a channel problem.
 * This is just a subclass of GridDefinition which provides some better default
 * values.
 */
template< typename FPT = double >
class ChannelDefinition : public GridDefinition<FPT>
{
public:
    /**
     * Construct an instance with the given default size in each direction
     * and the given default length.
     *
     * @param default_Nx Default grid size in the X direction.
     * @param default_Ny Default grid size in the Y direction.
     * @param default_Nz Default grid size in the Z direction.
     * @param default_Lx Default domain length in the X direction.
     * @param default_Ly Default domain length in the Y direction.
     * @param default_Lz Default domain length in the Z direction.
     */
    ChannelDefinition(typename GridDefinition<FPT>::size_type default_Nx = 16,
                      typename GridDefinition<FPT>::size_type default_Ny = 16,
                      typename GridDefinition<FPT>::size_type default_Nz = 16,
                      FPT default_Lx = 4.0*M_PI,
                      FPT default_Ly = 2.0,
                      FPT default_Lz = 4.0*M_PI/3.0)
        : GridDefinition<FPT>(default_Nx, default_Ny, default_Nz,
                              default_Lx, default_Ly, default_Lz);
};


} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_GRID_DEFINITION_HPP
