//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//
// Copyright (C) 2010 The PECOS Development Team
// Based heavily on the Boost Test predicate utilities
//
// Please see http://pecos.ices.utexas.edu for more information.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Suzerain is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------
//
// test_underling_tools.hpp: one-off test tools for underling
//
// $Id$
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef PECOS_SUZERAIN_TEST_UNDERLING_TOOLS_HPP
#define PECOS_SUZERAIN_TEST_UNDERLING_TOOLS_HPP

#include <suzerain/common.hpp>
#include <suzerain/mpi.hpp>
#include <suzerain/underling.hpp>
#include <fftw3-mpi.h>

/** A test fixture to setup and teardown MPI and FFTW MPI */
struct FFTWMPIFixture {

    FFTWMPIFixture() {
        MPI_Init(NULL, NULL); // NULL valid per MPI standard section 8.7
#ifdef HAVE_FFTW3_THREADS
        assert(fftw_init_threads());
        int nthreads = 1;
#if defined HAVE_OPENMP
        const char * const envstr = getenv("OMP_NUM_THREADS");
        if (envstr) {
            const int envnum = atoi(envstr);
            if (envnum > 0) nthreads = envnum;
        }
#elif defined HAVE_PTHREAD
    // TODO Provide sane nthreads default for FFTW pthread environment
#else
#error "Sanity check failed; unknown FFTW threading model in use."
#endif
        BOOST_TEST_MESSAGE("Using FFTW with " << nthreads << " threads.");
        fftw_plan_with_nthreads(nthreads);
#endif
        fftw_mpi_init();
    }

    ~FFTWMPIFixture() {
        fftw_mpi_cleanup();
#ifdef HAVE_FFTW3_THREADS
        fftw_cleanup_threads();
#endif
        MPI_Finalize();
    }
};

/* A test fixture to setup and teardown an underling use case */
struct UnderlingFixture {

    UnderlingFixture(MPI_Comm comm,
                     const int n0, const int n1, const int n2,
                     const int howmany,
                     const unsigned transposed_flags = 0)
        : grid(comm, n0, n1, n2),
          problem(grid, howmany, transposed_flags),
          data((underling_real *) fftw_malloc(
                    problem.local_memory()*sizeof(underling_real)),
                &fftw_free),
          plan(problem, data.get(),
               suzerain::underling::transpose::all, FFTW_ESTIMATE)
    {
        BOOST_REQUIRE(grid);
        BOOST_REQUIRE(problem);
        BOOST_REQUIRE(data);
        BOOST_REQUIRE(plan);
    }

    suzerain::underling::grid grid;
    suzerain::underling::problem problem;
    boost::shared_array<underling_real> data;
    suzerain::underling::plan plan;

};

#endif // PECOS_SUZERAIN_TEST_UNDERLING_TOOLS_HPP
