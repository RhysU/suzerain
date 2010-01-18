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
 * The dimensionality is chosen via the templated second parameter type.  On
 * error, throws an appropriate exception.  Any zero value contained in the
 * sequence (<tt>dimsBegin</tt>,<tt>dimsEnd</tt>] indicates that the method
 * should choose a suitable value.  The dimensionality of the problem is taken
 * from <tt>std::distance(dimsBegin,dimsEnd)</tt>.
 *
 * @param nnodes The number of nodes in a grid
 * @param dimsBegin The beginning iterator for the dimension lengths to choose.
 * @param dimsEnd   The ending iterator for the dimension lengths to choose.
 *
 * @see MPI_Dims_create for more details.
 */
template<class Integer, class ForwardIterator>
void dims_create(const Integer nnodes,
                 const ForwardIterator dimsBegin,
                 const ForwardIterator dimsEnd) {

    typedef typename
        std::iterator_traits<ForwardIterator>::value_type value_type;

    // Cast/copy input types as type 'int' to match MPI API
    using boost::numeric_cast;
    const int int_N      = numeric_cast<int>(std::distance(dimsBegin,dimsEnd));
    const int int_nnodes = numeric_cast<int>(nnodes);
    std::vector<int> int_dims(int_nnodes);
    std::transform(dimsBegin, dimsEnd, int_dims.begin(),
            boost::numeric::converter<int,value_type>());

    // Invoke the MPI API
    const int status = MPI_Dims_create(int_nnodes, int_N, int_dims.data());
    if (status != MPI_SUCCESS) throw std::runtime_error(error_string(status));

    // Copy the output back to the caller
    std::transform(int_dims.begin(), int_dims.end(), dimsBegin,
            boost::numeric::converter<value_type,int>());
}

} // namespace mpi

} // namespace suzerain

#endif // __SUZERAIN_MPI_HPP
