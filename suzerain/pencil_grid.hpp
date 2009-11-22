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
#ifndef __SUZERAIN_PENCILGRID_H
#define __SUZERAIN_PENCILGRID_H

#include <suzerain/common.hpp>
#include <suzerain/exceptions.hpp>

namespace suzerain
{

/** Encapsulates P3DFFT %pencil grid details, including the global grid size
 * and the processor grid decomposition parameters.  These are the same details
 * provided to P3DFFT's \c p3dfft_setup method.
 */
template < typename I = int >
class pencil_grid
{

public:
    typedef I dim_type; /**< Dimension type used to specify the grid */

    /**
     * Constructs an instance per the \c p3dff_setup method.
     *
     * @param proc_dims Processor grid sizes in \f$ P_1, P_2 \f$ directions
     * @param nx Global grid size in the streamwise direction
     * @param ny Global grid size in the wall-normal direction
     * @param nz Global grid size in the spanwise direction
     * @exception suzerain::invalid_argument on negative input
     */
    pencil_grid(const I proc_dims[2], const I nx, const I ny, const I nz)
    throw(suzerain::invalid_argument);

    const dim_type pg1; /**< Processor grid size in \f$ P_1 \f$ direction */
    const dim_type pg2; /**< Processor grid size in \f$ P_2 \f$ direction */

    const dim_type nx;  /**< Global grid size in the streamwise direction */
    const dim_type ny;  /**< Global grid size in the wall-normal direction */
    const dim_type nz;  /**< Global grid size in the spanwise direction */

};

template<typename I>
pencil_grid<I>::pencil_grid(
    const I proc_dims[2],
    const I nx, const I ny, const I nz)
throw(suzerain::invalid_argument)
        : pg1(proc_dims[0]), pg2(proc_dims[1]),
        nx(nx), ny(ny), nz(nz)
{
    using suzerain::invalid_argument;

    if (pg1 < 0) throw invalid_argument("pg1 must be nonnegative");

    if (pg2 < 0) throw invalid_argument("pg2 must be nonnegative");

    if (nx  < 0) throw invalid_argument("nx must be nonnegative");

    if (ny  < 0) throw invalid_argument("ny must be nonnegative");

    if (nz  < 0) throw invalid_argument("nz must be nonnegative");
}

} // namespace suzerain

#endif // __SUZERAIN_PENCILGRID_H
