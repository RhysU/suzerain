/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * underling.hpp: C++ wrappers for the C-based underling API
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_UNDERLING_HPP
#define __SUZERAIN_UNDERLING_HPP

#include <ostream>
#include <suzerain/underling.h>

// TODO Document the C++ API
// TODO Add basic_ostream details
// TODO Overload operator<< for underling_extents

namespace suzerain {

/**
 * Provides C++ wrappers for the C-based API in underling.h.  In particular,
 * provides RAII semantics for underling's opaque types and some
 * std::basic_ostream helpers.
 */
namespace underling {

typedef underling_real real;

typedef underling_extents extents;

namespace transpose {
    const unsigned long_n2_to_long_n1 = UNDERLING_TRANSPOSE_LONG_N2_TO_LONG_N1;

    const unsigned long_n1_to_long_n0 = UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N0;

    const unsigned long_n0_to_long_n1 = UNDERLING_TRANSPOSE_LONG_N0_TO_LONG_N1;

    const unsigned long_n1_to_long_n2 = UNDERLING_TRANSPOSE_LONG_N1_TO_LONG_N2;

    const unsigned all                =  UNDERLING_TRANSPOSE_ALL;
};

class grid {
public:

    grid(MPI_Comm comm, int n0, int n1, int n2, int pA = 0, int pB = 0)
        : grid_(underling_grid_create(comm, n0, n1, n2, pA, pB)) {};

    template< class InputIterator1, class InputIterator2 >
    grid(MPI_Comm comm, InputIterator1 n, InputIterator2 p)
        : grid_(NULL)
    {
        const int n0 = *n++; const int n1 = *n++; const int n2 = *n;
        const int pA = *p++; const int pB = *p;

        grid_ = underling_grid_create(comm, n0, n1, n2, pA, pB);
    }

    template< class InputIterator >
    grid(MPI_Comm comm, InputIterator n)
        : grid_(NULL)
    {
        const int n0 = *n++; const int n1 = *n++; const int n2 = *n;

        grid_ = underling_grid_create(comm, n0, n1, n2, 0, 0);
    }

    ~grid() { underling_grid_destroy(grid_); }

    const underling_grid get() const { return grid_; }

    operator bool () const { return grid_ != NULL; };

private:
    underling_grid grid_;
};

class problem {
public:

    problem(const grid &g, int howmany)
        : problem_(underling_problem_create(g.get(), howmany)) {};

    ~problem() { underling_problem_destroy(problem_); }

    const underling_problem get() const { return problem_; }

    underling_extents local_extents(int i) const {
        return underling_local_extents(problem_, i);
    }

    size_t local(int i,
                 int *start,
                 int *size,
                 int *stride,
                 int *strideorder) const {
        return underling_local(problem_, i, start, size, stride, strideorder);
    }

    size_t local_memory() const { return underling_local_memory(problem_); }

    operator bool () const { return problem_ != NULL; };

private:
    underling_problem problem_;
};

inline
size_t local_memory(const problem &p) {
    return p.local_memory();
}

inline
size_t local_memory_optimum(const grid &g, const problem &p) {
    return underling_local_memory_optimum(g.get(), p.get());
}

inline
size_t local_memory_maximum(const grid &g, const problem &p) {
    return underling_local_memory_maximum(g.get(), p.get());
}

inline
size_t local_memory_minimum(const grid &g, const problem &p) {
    return underling_local_memory_minimum(g.get(), p.get());
}

inline
size_t global_memory(const grid &g, const problem &p) {
    return underling_global_memory(g.get(), p.get());
}

inline
size_t global_memory_optimum(const grid &g, const problem &p) {
    return underling_global_memory_optimum(g.get(), p.get());
}

class plan {
public:

    plan(const problem &p,
         underling_real * data,
         unsigned transform_flags  = 0,
         unsigned fftw_rigor_flags = 0)
        : plan_(underling_plan_create(p.get(),
                                      data,
                                      transform_flags,
                                      fftw_rigor_flags)) {}

    ~plan() { underling_plan_destroy(plan_); }

    const underling_plan get() const { return plan_; }

    int execute_long_n2_to_long_n1() const {
        return underling_execute_long_n2_to_long_n1(plan_);
    }

    int execute_long_n1_to_long_n0() const {
        return underling_execute_long_n1_to_long_n0(plan_);
    }

    int execute_long_n0_to_long_n1() const {
        return underling_execute_long_n0_to_long_n1(plan_);
    }

    int execute_long_n1_to_long_n2() const {
        return underling_execute_long_n1_to_long_n2(plan_);
    }

    operator bool () const { return plan_ != NULL; };

private:
    underling_plan plan_;
};

} // namespace underling

} // namespace suzerain

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os,
        const suzerain::underling::extents &e)
{
    return os << '['
              << e.start[0] << ',' << (e.start[0] + e.size[0])
              << ")x["
              << e.start[1] << ',' << (e.start[1] + e.size[1])
              << ")x["
              << e.start[2] << ',' << (e.start[2] + e.size[2])
              << ')';
}

#endif // __SUZERAIN_UNDERLING_HPP
