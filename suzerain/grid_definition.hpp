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
 * Provides classes handling grid definitions, which are runtime
 * arguments or scenario parameters used to perform calculations.
 */

namespace suzerain {

namespace problem {

/**
 * Holds basic three dimensional computational grid dimensions, including
 * physical and discrete grid sizes, dealiasing factors, and the two
 * dimensional processor grid definition.
 */
template< typename FPT = double >
class GridDefinition : public IDefinition, public integral_types
{
public:
    /**
     * Construct an instance with the given default logical, size
     * in each direction, the given default length, and the given 
     * multiplicative dealiasing factors.
     *
     * @param default_Nx Default logical grid size in the X direction.
     * @param default_Ny Default logical grid size in the Y direction.
     * @param default_Nz Default logical grid size in the Z direction.
     * @param default_Lx Default domain length in the X direction.
     * @param default_Ly Default domain length in the Y direction.
     * @param default_Lz Default domain length in the Z direction.
     * @param default_DAFx Default dealiasing factor in the X direction.
     * @param default_DAFy Default dealiasing factor in the Y direction.
     * @param default_DAFz Default dealiasing factor in the Z direction.
     */
    explicit GridDefinition(
            size_type default_Nx = 16,
            size_type default_Ny = 16,
            size_type default_Nz = 16,
            FPT default_Lx = 2*boost::math::constants::pi<FPT>(),
            FPT default_Ly = 2*boost::math::constants::pi<FPT>(),
            FPT default_Lz = 2*boost::math::constants::pi<FPT>(),
            FPT default_DAFx = 1,
            FPT default_DAFy = 1,
            FPT default_DAFz = 1);

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
     * Retrieve global logical computational grid extents.  It does
     * not account for dealiasing factors.
     *
     * @return the global logical grid extents in the X, Y, and Z directions.
     */
    const size_type_3d& global_extents() const { return global_extents_; }

    /**
     * Retrieve global dealiased computational grid extents. These are the
     * global_extents() multiplied by the dealiasing factors DAFx(), DAFy(),
     * and DAFz().
     *
     * @return the global dealiased grid extents in the X, Y, and Z directions.
     */
    size_type_3d dealiased_extents() const {
        size_type_3d retval = global_extents_;
        retval[0] *= DAFx_;
        retval[1] *= DAFy_;
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
     * Retrieve computational grid size in the Y direction.  This is the number
     * of points in the domain without accounting for any dealiasing.
     *
     * @return the logical grid size in the Y direction.
     */
    size_type Ny() const { return global_extents_[1]; }

    /**
     * Retrieve grid size in the Z direction.  This is the number
     * of points in the domain without accounting for any dealiasing.
     *
     * @return the logical grid size in the Z direction.
     */
    size_type Nz() const { return global_extents_[2]; }

    /**
     * Retrieve the dealiasing factor for the X direction.  This factor
     * should be multiplied times Nx() to obtain an extent for Fourier
     * transformations.
     *
     * @return the dealiasing factor for the X direction.
     */
    FPT DAFx() const { return DAFx_; }

    /**
     * Retrieve the dealiasing factor for the Y direction.  This factor
     * should be multiplied times Ny() to obtain an extent for Fourier
     * transformations.
     *
     * @return the dealiasing factor for the Y direction.
     */
    FPT DAFy() const { return DAFy_; }

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

    FPT Lx_;  /**< Stores the X direction length */
    FPT Ly_;  /**< Stores the Y direction length */
    FPT Lz_;  /**< Stores the Z direction length */

    /** Stores the computational grid extents */
    size_type_3d global_extents_;

    FPT DAFx_;  /**< Stores the X direction dealiasing factor */
    FPT DAFy_;  /**< Stores the Y direction dealiasing factor */
    FPT DAFz_;  /**< Stores the Z direction dealiasing factor */

    /** Stores the processor grid extents */
    size_type_2d processor_grid_;
};

template< typename FPT >
GridDefinition<FPT>::GridDefinition(size_type default_Nx,
                                    size_type default_Ny,
                                    size_type default_Nz,
                                    FPT default_Lx,
                                    FPT default_Ly,
                                    FPT default_Lz,
                                    FPT default_DAFx,
                                    FPT default_DAFy,
                                    FPT default_DAFz)
    : options_("Grid definition"),
      Lx_(default_Lx),
      Ly_(default_Ly),
      Lz_(default_Lz),
      DAFx_(default_DAFx),
      DAFy_(default_DAFy),
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
        "Grid point count in X (streamwise) direction")
        ("Ny", po::value<size_type>(&global_extents_[1])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Ny"))
            ->default_value(default_Ny),
        "Grid point count in Y (wall normal) direction")
        ("Nz", po::value<size_type>(&global_extents_[2])
            ->notifier(bind2nd(ptr_fun(ensure_positive<size_type>),"Nz"))
            ->default_value(default_Nz),
        "Grid point count in Z (spanwise) direction")
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
        ("DAFx", po::value<FPT>(&DAFx_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"DAFx"))
            ->default_value(default_DAFx),
        "Dealiasing factor in X (streamwise) direction")
        ("DAFy", po::value<FPT>(&DAFy_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"DAFy"))
            ->default_value(default_DAFy),
        "Dealiasing factor in Y (wall normal) direction")
        ("DAFz", po::value<FPT>(&DAFz_)
            ->notifier(bind2nd(ptr_fun_ensure_positive_FPT,"DAFz"))
            ->default_value(default_DAFz),
        "Dealiasing factor in Z (spanwise) direction")
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
                      FPT default_Lx = 4*boost::math::constants::pi<FPT>(),
                      FPT default_Ly = 2,
                      FPT default_Lz = 4*boost::math::constants::pi<FPT>()/3,
                      FPT default_DAFx = 3/FPT(2),
                      FPT default_DAFy = 1,
                      FPT default_DAFz = 3/FPT(2))
        : GridDefinition<FPT>(default_Nx,   default_Ny,   default_Nz,
                              default_Lx,   default_Ly,   default_Lz,
                              default_DAFx, default_DAFy, default_DAFz) {}
};


} // namespace problem

} // namespace suzerain

#endif // __SUZERAIN_GRID_DEFINITION_HPP
