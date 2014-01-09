//-----------------------------------------------------------------------bl-
// underling 0.3.1: an FFTW MPI-based library for 3D pencil decompositions
// http://red.ices.utexas.edu/projects/underling
//
// Copyright (C) 2010, 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
//
// This file is part of underling.
//
// underling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// underling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with underling.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------el-
// $Id$

#ifndef UNDERLING_FIXTURES_HPP
#define UNDERLING_FIXTURES_HPP

#include <mpi.h>
#include <limits.h>
#include <underling/error.h>
#include <underling/underling.hpp>

/* A test fixture to setup and teardown an underling-based test case */
struct UnderlingFixture {

    UnderlingFixture(MPI_Comm comm,
                     const int n0, const int n1, const int n2,
                     const int howmany,
                     const unsigned transposed_flags = 0,
                     const bool in_place = true)
        : in_place(in_place),
          grid(comm, n0, n1, n2),
          problem(grid, howmany, transposed_flags),
          in(allocate_(problem.local_memory())),
          out(in_place ? in : allocate_(problem.local_memory())),
          plan(problem, in, out, underling::transpose::all, FFTW_ESTIMATE)
    {
        BOOST_REQUIRE(grid);
        BOOST_REQUIRE(problem);
        BOOST_REQUIRE(in);
        BOOST_REQUIRE(out);
        BOOST_REQUIRE(plan);
    }

    ~UnderlingFixture()
    {
        if (in)               fftw_free(in);
        if (!in_place && out) fftw_free(out);
    }

    const bool in_place;
    underling::grid grid;
    underling::problem problem;
    underling_real * const in;
    underling_real * const out;
    underling::plan plan;

    void fill_in_with_NaNs() {
        std::fill_n(in, problem.local_memory(),
                    std::numeric_limits<underling_real>::quiet_NaN());
    }

    void fill_out_with_NaNs() {
        std::fill_n(out, problem.local_memory(),
                    std::numeric_limits<underling_real>::quiet_NaN());
    }

private:
    static underling_real * allocate_(size_t count) {
        return (underling_real *) fftw_malloc(sizeof(underling_real) * count);
    }
};

/** A fixture for the Boost.Test that replaces underling_error */
#pragma warning(push,disable:2017 2021)
class BoostFailErrorHandlerFixture {
public:
    /** A underling_error_handler_t that invokes BOOST_FAIL */
    static void boost_fail_error_handler(
            const char *reason, const char *file, int line, int underling_errno)
    {
        std::ostringstream oss;
        oss << "Encountered '"
            << underling_strerror(underling_errno)
            << "' at "
            << file
            << ':'
            << line
            << " with reason '"
            << reason
            << "'";
        BOOST_FAIL(oss.str());
    }

    BoostFailErrorHandlerFixture()
        : previous_(underling_set_error_handler(&boost_fail_error_handler)) {}

    ~BoostFailErrorHandlerFixture() {
        underling_set_error_handler(previous_);
    }
private:
    underling_error_handler_t * previous_;
};
#pragma warning(pop)

#endif // UNDERLING_FIXTURES_HPP
