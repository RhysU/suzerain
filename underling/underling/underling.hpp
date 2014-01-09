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

#ifndef UNDERLING_HPP
#define UNDERLING_HPP

#include <ostream>
#include <underling/underling.h>

/** @file
 * Provides C++ wrappers for the C-based API in underling.h.  In
 * particular, provides RAII semantics for underling's opaque types and some
 * std::basic_ostream helpers.
 */

/** Provides C++ wrappers for underling's C-based API. */
namespace underling {

/** @see underling_real */
typedef underling_real real;

/** @see underling_extents */
typedef underling_extents extents;

/** Wraps transpose direction flags */
namespace transpose {

    enum { // Anonymous to avoid introducing an unnecessary type

        /** @see UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1 */
        long_n2_to_long_n1 = UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1,

        /** @see UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0 */
        long_n1_to_long_n0 = UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0,

        /** @see UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1 */
        long_n0_to_long_n1 = UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1,

        /** @see UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2 */
        long_n1_to_long_n2 = UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2,

        /** @see UNDERLING_TRANSPOSE_ALL */
        all                =  UNDERLING_TRANSPOSE_ALL,

        /** @see UNDERLING_TRANSPOSE_NONE */
        none               =  UNDERLING_TRANSPOSE_NONE

    };

}

/** Wraps transposed storage flags */
namespace transposed {

    enum { // Anonymous to avoid introducing an unnecessary type

        /** @see UNDERLING_TRANSPOSED_LONG_N2 */
        long_n2 = UNDERLING_TRANSPOSED_LONG_N2,

        /** @see UNDERLING_TRANSPOSED_LONG_N0 */
        long_n0 = UNDERLING_TRANSPOSED_LONG_N0

    };

}

/** A mixin for making a subclass non-copyable via private inheritance. */
class noncopyable
{
protected:
    noncopyable () {}
    ~noncopyable () {}

private:
    noncopyable (const noncopyable &);
    noncopyable & operator = (const noncopyable &);
};

/**
 * Provides a thin RAII wrapper for underling_grid.
 * @see underling_grid.
 */
class grid : private noncopyable {
public:

    /** @see underling_grid_create */
    grid(MPI_Comm comm, int n0, int n1, int n2, int pA = 0, int pB = 0)
        : grid_(underling_grid_create(comm, n0, n1, n2, pA, pB))
    {};

    /** @see underling_grid_create */
    template< class InputIterator1, class InputIterator2 >
    grid(MPI_Comm comm, InputIterator1 n, InputIterator2 p)
        : grid_(NULL)
    {
        const int n0 = *n++; const int n1 = *n++; const int n2 = *n;
        const int pA = *p++; const int pB = *p;

        grid_ = underling_grid_create(comm, n0, n1, n2, pA, pB);
    }

    /** @see underling_grid_create */
    template< class InputIterator >
    grid(MPI_Comm comm, InputIterator n)
        : grid_(NULL)
    {
        const int n0 = *n++; const int n1 = *n++; const int n2 = *n;

        grid_ = underling_grid_create(comm, n0, n1, n2, 0, 0);
    }

    /** @see underling_grid_destroy */
    ~grid() { underling_grid_destroy(grid_); }

    /** @see underling_grid_pA_size */
    int pA_size() const { return underling_grid_pA_size(grid_); }

    /** @see underling_grid_pB_size */
    int pB_size() const { return underling_grid_pB_size(grid_); }

    /** @return The wrapped underling_grid instance. */
    underling_grid get() const { return grid_; }

    /** @return True if the wrapped underling_grid instance is non-NULL. */
    operator bool () const { return grid_ != NULL; };

private:
    underling_grid grid_; /**< The wrapped underling_grid instance */
};

/**
 * Provides a thin RAII wrapper for underling_problem.
 * @see underling_problem.
 */
class problem : private noncopyable {
public:

    /** @see underling_problem_create */
    problem(const grid &g, int howmany, unsigned transposed_flags = 0)
        : problem_(underling_problem_create(
                g.get(), howmany, transposed_flags))
    {};

    /** @see underling_problem_destroy */
    ~problem() { underling_problem_destroy(problem_); }

    /** @return The wrapped underling_problem instance. */
    underling_problem get() const { return problem_; }

    /** @see underling_local_extents */
    underling_extents local_extents(int i) const
    {
        return underling_local_extents(problem_, i);
    }

    /** @see underling_local */
    size_t local(int i,
                 int *start  = NULL,
                 int *size   = NULL,
                 int *stride = NULL,
                 int *order  = NULL) const
    {
        return underling_local(problem_, i, start, size, stride, order);
    }

    /** @see underling_local_memory */
    size_t local_memory() const
    {
        return underling_local_memory(problem_);
    }

    /** @see underling_local_memory_optimum */
    size_t local_memory_optimum() const
    {
        return underling_local_memory_optimum(problem_);
    }

    /** @return True if the wrapped underling_problem instance is non-NULL. */
    operator bool () const { return problem_ != NULL; };

private:
    underling_problem problem_; /**< The wrapped underling_problem instance */
};

/** @see underling_local_memory_maximum */
inline
size_t local_memory_maximum(const grid &g, const problem &p)
{
    return underling_local_memory_maximum(g.get(), p.get());
}

/** @see underling_local_memory_minimum */
inline
size_t local_memory_minimum(const grid &g, const problem &p)
{
    return underling_local_memory_minimum(g.get(), p.get());
}

/** @see underling_global_memory */
inline
size_t global_memory(const grid &g, const problem &p)
{
    return underling_global_memory(g.get(), p.get());
}

/** @see underling_global_memory_optimum */
inline
size_t global_memory_optimum(const grid &g, const problem &p)
{
    return underling_global_memory_optimum(g.get(), p.get());
}

/**
 * Provides a thin RAII wrapper for underling_plan.
 * @see underling_plan.
 */
class plan : private noncopyable {
public:

    /** @see underling_plan_create */
    plan(const problem &p,
         underling_real * in,
         underling_real * out,
         unsigned transpose_flags  = 0,
         unsigned fftw_rigor_flags = 0)
        : plan_(underling_plan_create(p.get(),
                                      in,
                                      out,
                                      transpose_flags,
                                      fftw_rigor_flags))
    {}

    /** @see underling_plan_destroy */
    ~plan() { underling_plan_destroy(plan_); }

    /** @return The wrapped underling_plan instance. */
    underling_plan get() const { return plan_; }

    /** @see underling_execute_long_n2_to_long_n1 */
    int execute_long_n2_to_long_n1(underling_real * in,
                                   underling_real * out) const
    {
        return underling_execute_long_n2_to_long_n1(plan_, in, out);
    }

    /** @see underling_execute_long_n1_to_long_n0 */
    int execute_long_n1_to_long_n0(underling_real * in,
                                   underling_real * out) const
    {
        return underling_execute_long_n1_to_long_n0(plan_, in, out);
    }

    /** @see underling_execute_long_n0_to_long_n1 */
    int execute_long_n0_to_long_n1(underling_real * in,
                                   underling_real * out) const
    {
        return underling_execute_long_n0_to_long_n1(plan_, in, out);
    }

    /** @see underling_execute_long_n1_to_long_n2 */
    int execute_long_n1_to_long_n2(underling_real * in,
                                   underling_real * out) const
    {
        return underling_execute_long_n1_to_long_n2(plan_, in, out);
    }

    /** @return True if the wrapped underling_problem instance is non-NULL. */
    operator bool () const { return plan_ != NULL; };

private:
    underling_plan plan_;  /**< The wrapped underling_plan instance */
};

} // namespace underling

/**
 * Outputs an underling_extents or underling::extents instance
 * as a human-readable string on any std::basic_ostream.
 *
 * @param os On which to output \c e.
 * @param e  To be output.
 *
 * @return The modified \c os.
 */
template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os,
        const underling::extents &e)
{
    return os << '['
              << e.start[0] << ',' << (e.start[0] + e.size[0])
              << ")x["
              << e.start[1] << ',' << (e.start[1] + e.size[1])
              << ")x["
              << e.start[2] << ',' << (e.start[2] + e.size[2])
              << ")x["
              << e.start[3] << ',' << (e.start[3] + e.size[3])
              << ')';
}

/** @see underling_extents_cmp */
inline
bool operator==(const underling::extents &e1,
                const underling::extents &e2) {
    return !underling_extents_cmp(&e1, &e2);
}

#endif // UNDERLING_HPP
