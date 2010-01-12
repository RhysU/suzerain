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
 * mpi.cpp: miscellaneous utilities for working with MPI's C interface
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#include <suzerain/mpi.hpp>
#include <mpi.h>

namespace suzerain {

namespace mpi {

std::string error_string(const int errorcode) {
    std::ostringstream retval;

    char *errorstring = NULL;
    int resultlen;
    const int status = MPI_Error_string(errorcode, errorstring, &resultlen);
    switch (status) {
        case MPI_SUCCESS:
            retval << errorstring;
            break;
        case MPI_ERR_ARG:
            retval << "No error string known for code " << errorcode;
            break;
        default:
            retval << "Unknown error " << status
                   << " when retrieving string for code " << errorcode;
    }

    return retval.str();
}

int comm_size(MPI_Comm comm) {
    int size;
    const int status = MPI_Comm_size(comm, &size);
    if (status != MPI_SUCCESS) throw std::runtime_error(error_string(status));
    return size;
}

int comm_rank(MPI_Comm comm) {
    int rank;
    const int status = MPI_Comm_rank(comm, &rank);
    if (status != MPI_SUCCESS) throw std::runtime_error(error_string(status));
    return rank;
}

} // namespace mpi

} // namespace suzerain
