//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_MPI_HPP
#define SUZERAIN_MPI_HPP

/** @file
 * C++ utilities for working with MPI's C-based interface.
 */

#include <suzerain/common.hpp>
#include <suzerain/mpi.h>

namespace suzerain {

/**
 * C++ utilities for working with MPI's C-based interface.
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
 * Ensure that MPI_Init has already been called.
 * In debug mode, entering this method when <tt>!MPI_Initialized</tt>
 * results in a failed assertion.  Otherwise, entering this
 * method when <tt>!MPI_Initialized</tt> causes a std::logic_error.
 *
 * @throw std::logic_error if MPI_Init has not been called.
 */
void ensure_mpi_initialized() throw (std::logic_error);

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
 * Create a reasonable identifier for the current process within the given MPI
 * communicator.  The result is based on the communicator's name, if present,
 * and the rank of the current process.
 *
 * @param comm Communicator to use.
 *
 * @return an identifier munging the communicator and rank into a
 * human-readable string.
 */
std::string comm_rank_identifier(MPI_Comm comm);

/**
 * Create a division of processors in a Cartesian grid.
 * The dimensionality is chosen via the distance between \c dimsBegin and \c
 * dimsEnd.  On error, throws an appropriate exception.  Any zero value
 * contained in the sequence (<tt>dimsBegin</tt>,<tt>dimsEnd</tt>] indicates
 * that the method should choose a suitable value.  The dimensionality of the
 * problem is taken from <tt>std::distance(dimsBegin,dimsEnd)</tt>.
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

    // numeric_cast and/or copy input types to type int to match MPI API
    using boost::numeric_cast;
    const int int_nnodes = numeric_cast<int>(nnodes);
    const int int_ndims  = numeric_cast<int>(std::distance(dimsBegin,dimsEnd));
    scoped_array<int> int_dims(new int[int_ndims]);
    std::transform(dimsBegin, dimsEnd, int_dims.get(),
            boost::numeric::converter<int,value_type>());

    // Invoke the MPI API
    const int status = MPI_Dims_create(int_nnodes, int_ndims, int_dims.get());
    if (status != MPI_SUCCESS) throw std::runtime_error(error_string(status));

    // Copy the output back to the caller
    std::transform(int_dims.get(), int_dims.get() + int_ndims, dimsBegin,
            boost::numeric::converter<value_type,int>());
}

/**
 * Create two MPI_Barriers across the communicator and sleep indefinitely
 * between them on the given rank.  Useful when one needs to attach a debugger
 * to a process.
 *
 * @param comm The MPI Communicator used in the two \c MPI_Barrier calls.
 * @param sleep_rank The rank of the process to enter an indefinite sleep loop.
 * @param os   The output stream on which to write the ready message.
 *
 * @return The modified output stream.
 * @see The <a href="http://www.open-mpi.org/faq/">
 *      Open MPI FAQ</a> for more information.
 */
template< typename charT, typename traits >
std::basic_ostream<charT,traits>& sleep_barrier(
        MPI_Comm comm = MPI_COMM_SELF,
        const int sleep_rank = 0,
        std::basic_ostream<charT,traits> &os = std::cerr)
{
    if (sleep_rank < 0) {
        throw (std::runtime_error("sleep_rank < 0"));
    }
    if (sleep_rank >= comm_size(comm)) {
        throw (std::runtime_error("sleep_rank > comm_size(comm)"));
    }

    // Retrieve the current hostname in a beyond-paranoid manner, continuing
    // the sleep_barrier iff the hostname retrieval succeeds.
    char hostname[256];
    {
        const std::string id = comm_rank_identifier(comm);
        // Ensure null termination in event of truncation
        hostname[sizeof(hostname)-1] = 0;
        errno = 0;
        int local_gethostname_err = gethostname(hostname, sizeof(hostname)-1);
        if (local_gethostname_err) {
            os << "Error from gethostname on " << id
               << ":" << strerror(errno)
               << "; skipping sleep_barrier."  << std::endl;
        }
        int global_gethostname_err = 0;
        const int err = MPI_Allreduce(&local_gethostname_err,
                                      &global_gethostname_err,
                                      1, MPI_INT, MPI_LOR, comm);
        if (err) throw std::runtime_error(error_string(err));
        if (global_gethostname_err) return os;
    }

    // Force all processes to sync up prior to sleeping
    {
        const int err = MPI_Barrier(comm);
        if (err) throw std::runtime_error(error_string(err));
    }

    // Sleep indefinitely if we're the lucky process.
    // Set i nonzero using a debugger to exit infinite loop below.
    if (comm_rank(comm) == sleep_rank) {
        os << "sleep_barrier: pid " << getpid() << " on host "
           << hostname  << " ready for attach." << std::endl;
        int i = 0;
        while (0 == i) {
            sleep(5);
        }
    }

    // Force all processes to sync up after sleeping
    {
        const int err = MPI_Barrier(comm);
        if (err) throw std::runtime_error(error_string(err));
    }

    return os;
}

} // namespace mpi

} // namespace suzerain

#endif // SUZERAIN_MPI_HPP
