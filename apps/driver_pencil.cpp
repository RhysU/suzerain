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
 * driver_pencil.cc: A parallel MPI transpose proving ground
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <log4cxx/logger.h>
#include <suzerain/program_options.hpp>
#include <suzerain/mpi.hpp>

#define ONLYPROC0(expr) if (!procid) { expr ; } else

#pragma warning(push,disable:1418)
template< typename Integer >
std::pair<Integer,Integer> calculate_processor_extents(
        const Integer nproc,
        const Integer procid,
        const Integer npoints)
{
    const Integer remainder = npoints - nproc*(npoints/nproc);
    Integer start = 0;
    for (Integer i = 0; i < procid; ++i){
        start += npoints/nproc;
        if (i >= nproc - remainder) start += 1;
    }

    Integer end = start + npoints/nproc;
    if (procid >= nproc - remainder) end += 1;

    return std::make_pair(start, end); // [start, end)
}
#pragma warning(pop)

int main(int argc, char **argv)
{
    const int NDIM = 2;                      // Dimensionality of the problem
    boost::array<int,NDIM> N = { 16, 16 };   // Global extents in X, Y direction

    MPI_Init(&argc, &argv);                   // Initialize MPI on startup
    atexit((void (*) ()) MPI_Finalize);       // Finalize MPI at exit

    const int nproc  = suzerain::mpi::comm_size(MPI_COMM_WORLD);
    const int procid = suzerain::mpi::comm_rank(MPI_COMM_WORLD);
    log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger(
            suzerain::mpi::comm_rank_identifier(MPI_COMM_WORLD));

    // Process command line arguments in process 0 and broadcast them
    suzerain::ProgramOptions options;
    namespace po = boost::program_options;
    options.add_options()
        ("nx", po::value<int>(&N[0])->default_value(N[0]),
        "Processor grid size in X direction.")
        ("ny", po::value<int>(&N[1])->default_value(N[1]),
        "Processor grid size in Y direction.")
    ;
    if (!procid) {
        options.process(argc, argv);
    }
    MPI_Bcast(&N, 2, MPI_INTEGER, 0, MPI_COMM_WORLD);

    ONLYPROC0(LOG4CXX_INFO(logger, "Number of processors: " << nproc));
    ONLYPROC0(LOG4CXX_INFO(logger, "Global grid dimensions: "
                           << boost::format("(% 4d, % 4d)") % N[0] % N[1]));

    // Allocate the storage we need to house either X or Y pencils
    // First compute the processor extents when long in either direction
    boost::multi_array<int,NDIM> procextents(boost::extents[NDIM][2]);
    for (int i = 0; i < procextents.shape()[0]; ++i) {
        BOOST_AUTO(extents, calculate_processor_extents(nproc, procid, N[i]));
        procextents[i][0] = extents.first;
        procextents[i][1] = extents.second;
    }
    // Second compute the storage required when long in either direction
    int max_extent_size = boost::integer_traits<int>::const_min;
    int max_storage     = boost::integer_traits<int>::const_min;
    for (int i = 0; i < procextents.shape()[0]; ++i) {
        int i_storage = procextents[i][1] - procextents[i][0];
        max_extent_size = std::max(max_extent_size, i_storage);
        for (int j = 0; j < procextents.shape()[0]; ++j) {
            if (j != i) i_storage *= N[j];
        }
        max_storage = std::max(max_storage, i_storage);
    }
    // Third allocate and create views of the required raw storage
    boost::scoped_array<double> raw_storage(new double[max_storage]);
    typedef boost::multi_array_ref<double,NDIM> array_type;
    boost::array<bool,
                 array_type::dimensionality> ascending = {{ true, true }};
    boost::array<array_type::size_type,
                 array_type::dimensionality> ordering_X = {{ 0, 1 }};
    boost::array<array_type::size_type,
                 array_type::dimensionality> ordering_Y = {{ 1, 0 }};
    array_type pencil_X(
            raw_storage.get(),
            boost::extents[N[0]][procextents[0][1] - procextents[0][0]],
            boost::general_storage_order<array_type::dimensionality>(
                ordering_X.begin(), ascending.begin()));
    array_type pencil_Y(
            raw_storage.get(),
            boost::extents[procextents[1][1] - procextents[1][0]][N[1]],
            boost::general_storage_order<array_type::dimensionality>(
                ordering_Y.begin(), ascending.begin()));

    LOG4CXX_INFO(logger, "storage required: " << max_storage);
    LOG4CXX_INFO(logger, "pencil_X extents: " << boost::format("(% 4d, % 4d)")
                 % pencil_X.shape()[0] % pencil_X.shape()[1]);
    LOG4CXX_INFO(logger, "pencil_Y extents: " << boost::format("(% 4d, % 4d)")
                 % pencil_Y.shape()[0] % pencil_Y.shape()[1]);

}
