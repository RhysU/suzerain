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
#include <suzerain/types.hpp>
#include <suzerain/mpi.hpp>
#include <p3dfft_d.h>

namespace suzerain
{

/**
 * Encapsulates P3DFFT %pencil grid details, including the global grid extent
 * and the processor grid decomposition parameters.  Appropriately handles
 * munging this information to obtain P3DFFT calls which are stride one
 * in Y in wave space.  Unless otherwise noted, all indices start from
 * zero with X, Y, and Z having indices 0, 1, and 2, respectively.
 */
class pencil_grid : public integral_types
{

public:
    /**
     * Constructs an instance for the given global grid extents.
     * Under the covers, P3DFFT's <tt>p3dfft_setup</tt> is
     * invoked to determine pencil decomposition parameters for
     * the local process.
     *
     * @param global_extents Global grid extents in the
     *        streamwise (X), wall-normal (Y), and spanwise (Z)
     *        directions.
     * @param processor_grid The processor grid decomposition to use in the
     *        \f$ P_0 \f$ and \f$ P_1 \f$ directions.  Providing a zero
     *        for either value causes that value to be determined automatically.
     */
    pencil_grid(const size_type_3d &global_extents,
                const size_type_2d &processor_grid);

    /**
     * Tears down an instance.  Under the covers, P3DFFT's <tt>p3dfft_clean</tt>
     * is invoked.
     */
    ~pencil_grid();

    /**
     * Retrieve global computational grid extents.
     *
     * @return the global grid extents in the X, Y, and Z directions.
     */
    const size_type_3d& global_extents() const { return global_extents_; }

    /**
     * Retrieve the processor grid extents.
     *
     * In physical space, \f$ P_0 \f$ is the grid extent in the Z direction
     * and \f$ P_1 \f$ is the grid extent in the Y direction.
     * In wave space, \f$ P_0 \f$ is the grid extent in the X direction
     * and \f$ P_1 \f$  is the grid extent in the Z direction.
     *
     * @return the processor grid extents in the \f$ P_0 \f$
     *         and \f$ P_1 \f$ directions as indices 0 and 1, respectively.
     */
    const size_type_2d& processor_grid() const { return processor_grid_; }

    /**
     * Retrieve local pencil physical space starting indices (inclusive) within
     * the global extents.
     *
     * @return local pencil physical space starting indices in the X, Y, and Z
     * directions.
     */
    const index_3d& local_physical_start() const { return pstart_; }

    /**
     * Retrieve local pencil physical space ending indices (exclusive)
     * within the global extents.
     *
     * @return local pencil physical space starting indices in the X, Y, and Z
     * directions.
     */
    const index_3d& local_physical_end() const { return pend_; }

    /**
     * Retrieve local pencil physical space extents.
     *
     * @return local pencil physical space extents in the X, Y, and Z
     * directions.
     */
    const index_3d& local_physical_extent() const { return pextent_; }

    /**
     * Retrieve the number of contiguous real scalars required to store
     * a pencil's worth of data.  This accounts for any padding required
     * due to differences in the local physical and wave space extents
     * and the fact that wave storage requires complex scalars.
     *
     * @return Number of real-valued scalars (i.e. <tt>double</tt>s)
     *         required to store one pencil's contiguous data.
     */
    size_type local_physical_storage() const;

    /**
     * Retrieve local pencil wave space starting indices (inclusive) within
     * the global extents.
     *
     * @return local pencil wave space starting indices in the X, Y, and Z
     * directions.
     */
    const index_3d& local_wave_start() const { return wstart_; }

    /**
     * Retrieve local pencil wave space ending indices (exclusive)
     * within the global extents.
     *
     * @return local pencil wave space starting indices in the X, Y, and Z
     * directions.
     */
    const index_3d& local_wave_end() const { return wend_; }

    /**
     * Retrieve local pencil wave space extents.
     *
     * @return local pencil wave space extents in the X, Y, and Z
     * directions.
     */
    const index_3d& local_wave_extent() const { return wextent_; }

    /**
     * Retrieve the number of contiguous complex scalars required to store
     * a pencil's worth of data.  This accounts for any padding required
     * due to differences in the local physical and wave space extents
     * and the fact that physical space storage requires real-valued scalars.
     *
     * @return Number of complex-valued scalars (i.e. <tt>double[2]</tt>s)
     *         required to store one pencil's contiguous data.
     */
    size_type local_wave_storage() const;

    /**
     * Collectively transform a field from wave space to physical space.  The
     * input will be complex-valued and stored according to local_wave_start()
     * and local_wave_end().  The output will be real-valued and stored
     * according to local_physical_state() and local_physical_end().
     *
     * @param inout Field to transform in place.
     */
    void transform_wave_to_physical(double * inout) const {
        p3dfft_btran_c2r(inout, inout);
    }

    /**
     * Collectively transform a field from wave space to physical space.  The
     * input will be complex-valued and stored according to local_wave_start()
     * and local_wave_end().  The output will be real-valued and stored
     * according to local_physical_state() and local_physical_end().
     *
     * @param in  Input complex-valued field.
     * @param out Output real-valued field.
     */
    void transform_wave_to_physical(double * in, double * out) const {
        p3dfft_btran_c2r(in, out);
    }

    /**
     * Collectively transform a field from physical space to wave space.  The
     * input will be complex-valued and stored according to
     * local_physical_start() and local_physical_end().  The output will be
     * real-valued and stored according to local_wave_state() and
     * local_wave_end().
     *
     * @param inout Field to transform in place.
     */
    void transform_physical_to_wave(double * inout) const {
        p3dfft_ftran_r2c(inout, inout);
    }

    /**
     * Collectively transform a field from physical space to wave space.  The
     * input will be complex-valued and stored according to
     * local_physical_start() and local_physical_end().  The output will be
     * real-valued and stored according to local_wave_state() and
     * local_wave_end().
     *
     * @param in  Input real-valued field.
     * @param out Output complex-valued field.
     */
    void transform_physical_to_wave(double * in, double * out) const {
        p3dfft_ftran_r2c(in, out);
    }

private:
    /**
     * Global grid extent in the streamwise, wall-normal,
     * and spanwise directions.
     **/
    size_type_3d global_extents_;

    /** Processor grid extent in the \f$ P_0 \f$ and \f$ P_1 \f$ directions. */
    size_type_2d processor_grid_;

    /** Physical space starting indices for local storage within global grid */
    index_3d pstart_;

    /** Physical space ending indices for local storage within global grid */
    index_3d pend_;

    /** Physical space dimensions for local storage */
    index_3d pextent_;

    /** Wave space starting indices for local storage within global grid */
    index_3d wstart_;

    /** Wave space ending indices for local storage within global grid */
    index_3d wend_;

    /** Wave space dimensions for local storage */
    index_3d wextent_;

    /** Was p3dfft_setup successfully called during construction? */
    bool p3dfft_setup_called_;
};

} // namespace suzerain

#endif // __SUZERAIN_PENCIL_GRID_H
