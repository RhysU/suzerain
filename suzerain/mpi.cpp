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

/** @file
 * @copydoc mpi.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/mpi.hpp>

namespace suzerain {

namespace mpi {

std::string error_string(const int errorcode)
{
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

void ensure_mpi_initialized() throw (std::logic_error)
{
    int flag;
    MPI_Initialized(&flag);
    assert(flag);
    if (!flag) throw std::logic_error("MPI_Init has not been invoked");
}

int comm_size(MPI_Comm comm)
{
    int size;
    const int error = MPI_Comm_size(comm, &size);
    if (error) throw std::runtime_error(error_string(error));
    return size;
}

int comm_rank(MPI_Comm comm)
{
    int rank;
    const int error = MPI_Comm_rank(comm, &rank);
    if (error) throw std::runtime_error(error_string(error));
    return rank;
}

std::string comm_rank_identifier(MPI_Comm comm)
{
    // All MPI calls below fail for MPI_COMM_NULL,
    // so handle it as a special case.
    if (comm == MPI_COMM_NULL) {
        return "NULL";
    }

    // We want the name for comm as a C-style string at p_buffer
    // with resultlen = strlen(p_buffer)
    char buffer[MPI_MAX_OBJECT_NAME];
    char * p_buffer = buffer;
    int resultlen;
    const int error = MPI_Comm_get_name(comm, buffer, &resultlen);
    if (error) throw std::runtime_error(error_string(error));
    if (resultlen == 0) {
        // Providing a default value if MPI had none.
        resultlen = snprintf(buffer, sizeof(buffer), "c0x%X", comm);
        assert(resultlen >= 0);
    } else {
        // Brevity, trim any leading MPI_COMM_ business by modifying p_buffer
        static const char prefix[]   = "MPI_COMM_";
        static const int  prefix_len = sizeof(prefix) - 1;
        if (0 == strncasecmp(buffer, prefix, prefix_len)) {
            p_buffer += prefix_len;
            resultlen -= prefix_len;
        }
    }

    // Look up comm's size and this process' rank
    int size, rank;
    char sep;
    if (comm == MPI_COMM_SELF) {
        // size and rank are useless for MPI_COMM_SELF,
        // so provide details from MPI_COMM_WORLD instead
        sep  = 'w';
        size = comm_size(MPI_COMM_WORLD);
        rank = comm_rank(MPI_COMM_WORLD);
    } else {
        sep  = 'r';
        size = comm_size(comm);
        rank = comm_rank(comm);
    }

    // Pack up and return a std::string containing name and rank
    // Take care so that all results from the same comm are the same length
    std::ostringstream oss;
    oss.write(p_buffer, resultlen);
    oss << sep << std::setfill('0')
#pragma warning(push,disable:810 2259)
        << std::setw(ceil(log10(size)))
#pragma warning(pop)
        << rank;
    return oss.str();
}

} // namespace mpi

} // namespace suzerain
