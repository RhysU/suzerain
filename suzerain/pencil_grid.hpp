/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
 * pencil_grid.hpp: Class to manage data layout concerns for P3DFFT usage
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_PENCIL_GRID_H
#define __SUZERAIN_PENCIL_GRID_H

#include <suzerain/common.hpp>
#include <suzerain/mpi.hpp>

namespace suzerain
{

// Forward declarations for different implementations
class pencil_grid_p3dfft;
class pencil_grid_underling;

// Choose pencil grid implementation based on compile-time availability
// Uses public #defines present from suzerain/suzerain-config.h
#if defined(SUZERAIN_HAVE_P3DFFT)
/** Use P3DFFT-based \c pencil_grid implementation */
typedef pencil_grid_p3dfft pencil_grid;
#elif defined(SUZERAIN_HAVE_UNDERLING)
/** Use underling-based \c pencil_grid implementation */
typedef pencil_grid_underling pencil_grid;
#else
# error "Found neither suzerain-p3dfft nor underling libraries"
#endif

/**
 * Encapsulates P3DFFT %pencil grid details, including the global grid extent
 * and the processor grid decomposition parameters.  Appropriately handles
 * munging this information to obtain P3DFFT calls which are stride one
 * in Y in wave space.  Unless otherwise noted, all indices start from
 * zero with X, Y, and Z having indices 0, 1, and 2, respectively.
 */
class pencil_grid_p3dfft
{
public:
    // See http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /**
     * Constructs an instance for the given global physical grid extents.
     * Under the covers, P3DFFT's <tt>p3dfft_setup</tt> is invoked to determine
     * pencil decomposition parameters for the local process.
     *
     * @param global_physical_extent Global physical grid extents in the
     *        streamwise (X), wall-normal (Y), and spanwise (Z)
     *        directions.
     * @param processor_grid The processor grid decomposition to use in the
     *        \f$ P_0 \f$ and \f$ P_1 \f$ directions.  Providing a zero
     *        for either value causes the value to be determined automatically.
     */
    template < typename RandomAccessContainer1,
               typename RandomAccessContainer2 >
    pencil_grid_p3dfft(const RandomAccessContainer1 &global_physical_extent,
                       const RandomAccessContainer2 &processor_grid)
    {
        using boost::numeric_cast;
        this->construct_(numeric_cast<int>(global_physical_extent[0]),
                         numeric_cast<int>(global_physical_extent[1]),
                         numeric_cast<int>(global_physical_extent[2]),
                         numeric_cast<int>(processor_grid[0]),
                         numeric_cast<int>(processor_grid[1]));
    }

    /**
     * Tears down an instance.
     * Under the covers, P3DFFT's <tt>p3dfft_clean</tt> is invoked.
     */
    ~pencil_grid_p3dfft();

    /** Global grid extents using physical space sizes. */
    Eigen::Array3i global_physical_extent;

    /** Global grid extents using wave space sizes. */
    Eigen::Array3i global_wave_extent;

    /**
     * The processor grid extents.
     *
     * In physical space, \f$ P_0 \f$ is the grid extent in the Z direction
     * and \f$ P_1 \f$ is the grid extent in the Y direction.
     * In wave space, \f$ P_0 \f$ is the grid extent in the X direction
     * and \f$ P_1 \f$  is the grid extent in the Z direction.
     *
     * @return the processor grid extents in the \f$ P_0 \f$
     *         and \f$ P_1 \f$ directions as indices 0 and 1, respectively.
     */
    Eigen::Array2i processor_grid;

    /**
     * Local pencil physical space starting indices (inclusive) within
     * the global extents.
     */
    Eigen::Array3i local_physical_start;

    /**
     * Local pencil physical space ending indices (exclusive)
     * within the global extents.
     */
    Eigen::Array3i local_physical_end;

    /**
     * Local pencil physical space extents.
     */
    Eigen::Array3i local_physical_extent;

    /**
     * Retrieve the number of contiguous real scalars required to store
     * a pencil's worth of data.  This accounts for any padding required
     * due to differences in the local physical and wave space extents
     * and the fact that wave storage requires complex scalars.
     *
     * @return Number of real-valued scalars (i.e. <tt>double</tt>s)
     *         required to store one pencil's contiguous data.
     */
    std::size_t local_physical_storage() const;

    /**
     * Local pencil wave space starting indices (inclusive) within
     * the global extents.
     */
    Eigen::Array3i local_wave_start;

    /**
     * Local pencil wave space ending indices (exclusive)
     * within the global extents.
     */
    Eigen::Array3i local_wave_end;

    /**
     * Local pencil wave space extents.
     */
    Eigen::Array3i local_wave_extent;

    /**
     * Retrieve the number of contiguous complex scalars required to store
     * a pencil's worth of data.  This accounts for any padding required
     * due to differences in the local physical and wave space extents
     * and the fact that physical space storage requires real-valued scalars.
     *
     * @return Number of complex-valued scalars (i.e. <tt>double[2]</tt>s)
     *         required to store one pencil's contiguous data.
     */
    std::size_t local_wave_storage() const;

    /**
     * Collectively transform a field from wave space to physical space.  The
     * input will be complex-valued and stored according to local_wave_start()
     * and local_wave_end().  The output will be real-valued and stored
     * according to local_physical_state() and local_physical_end().
     *
     * @param inout Field to transform in place.
     */
    void transform_wave_to_physical(double * inout) const;

    /**
     * Collectively transform a field from physical space to wave space.  The
     * input will be complex-valued and stored according to
     * local_physical_start() and local_physical_end().  The output will be
     * real-valued and stored according to local_wave_state() and
     * local_wave_end().
     *
     * @param inout Field to transform in place.
     */
    void transform_physical_to_wave(double * inout) const;

private:

    /** Was p3dfft_setup successfully called during construction? */
    bool p3dfft_setup_called_;

    /**
     * Internal routine performing construction-like tasks
     */
    void construct_(int Nx, int Ny, int Nz, int Pa, int Pb);
};

} // namespace suzerain

#endif // __SUZERAIN_PENCIL_GRID_H
