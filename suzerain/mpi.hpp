/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2009 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * mpi.hpp: miscellaneous utilities for working with MPI's C interface
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_MPI_HPP
#define __SUZERAIN_MPI_HPP

#include <suzerain/common.hpp>
#include <mpi.h>

namespace suzerain {

/**
 * Provides miscellaneous utilities for working with MPI's C interface.
 */
namespace mpi {

/**
 * Return an error string for the given MPI error code.
 *
 * @param errorcode to retrieve information about via
 *                  <tt>MPI_Error_string</tt>.
 *
 * @return either the MPI error string or a string describing why no
 * information was available.
 * @see MPI_Error_string for more details.
 */
std::string error_string(const int errorcode);

/**
 * Determine the size of the group associated with a communicator.
 * On error, throws an appropriate exception.
 *
 * @param comm The communicator to query, e.g. MPI_COMM_WORLD
 *
 * @return the size of the group.
 * @see MPI_Comm_size for more details.
 */
int comm_size(MPI_Comm comm);

/**
 * Determine the rank of the calling process in the communicator.
 * On error, throws an appropriate exception.
 *
 * @param comm The communicator to query, e.g. MPI_COMM_WORLD
 *
 * @return the rank of the calling process in the group of \c comm.
 * @see MPI_Comm_rank for more details.
 */
int comm_rank(MPI_Comm comm);

/**
 * Create a division of processors in a Cartesian grid.
 * The dimensionality is chosen via the templated second parameter type.
 * On error, throws an appropriate exception.
 *
 * @param nnodes The number of nodes in a grid
 * @param dims   Specifies the number of nodes in each dimension.  A
 *               value of zero indicates that the method should fill in
 *               a suitable value.
 *
 * @see MPI_Dims_create for more details.
 */
template<std::size_t N, class Integer1, class Integer2>
void dims_create(const Integer1 nnodes, boost::array<Integer2,N> &dims) {
    using boost::numeric::converter;

    const int int_N      = boost::numeric_cast<int>(N);
    const int int_nnodes = boost::numeric_cast<int>(nnodes);
    boost::array<int,N> int_dims;
    std::transform(dims.begin(), dims.end(),
                   int_dims.begin(), converter<int,Integer2>());

    const int status = MPI_Dims_create(int_nnodes, int_N, int_dims.data());
    if (status != MPI_SUCCESS) throw std::runtime_error(error_string(status));

    std::transform(int_dims.begin(), int_dims.end(),
                    dims.begin(), converter<Integer2,int>());
}

} // namespace mpi

} // namespace suzerain

#endif // __SUZERAIN_MPI_HPP
