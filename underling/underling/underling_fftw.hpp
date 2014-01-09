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

#ifndef UNDERLING_FFTW_HPP
#define UNDERLING_FFTW_HPP

#include <underling/underling_fftw.h>
#include <underling/underling.hpp>

/** @file
 * Provides C++ wrappers for the C-based API in underling_fftw.h.  In
 * particular, provides RAII semantics for opaque types.
 */

namespace underling {

namespace fftw {

/** @see underling_fftw_extents */
typedef underling_fftw_extents extents;

/** Wraps packed storage flags */
namespace packed {

    enum { // Anonymous to avoid introducing an unnecessary type

        /** @see UNDERLING_FFTW_PACKED_LONG_N2 */
        long_n2 = UNDERLING_FFTW_PACKED_LONG_N2,

        /** @see UNDERLING_FFTW_PACKED_LONG_N0 */
        long_n0 = UNDERLING_FFTW_PACKED_LONG_N0,

        /** @see UNDERLING_FFTW_PACKED_ALL */
        all     =  UNDERLING_FFTW_PACKED_ALL,

        /** @see UNDERLING_FFTW_PACKED_NONE */
        none    =  UNDERLING_FFTW_PACKED_NONE

    };

}

/**
 * Provides a thin RAII wrapper for underling_fftw_plan.
 * @see underling_fftw_plan.
 */
class plan : public noncopyable {
public:

    /** A tag type used to indicate a complex-to-complex forward transform */
    struct c2c_forward  {};

    /** A tag type used to indicate a complex-to-complex backward transform */
    struct c2c_backward {};

    /** A tag type used to indicate a real-to-complex forward transform */
    struct r2c_forward  {};

    /** A tag type used to indicate a complex-to-real backward transform */
    struct c2r_backward {};

    /** @see underling_fftw_plan_create_c2c_forward */
    plan(const c2c_forward tag,
         const problem &p,
         int long_ni,
         underling_real *in,
         underling_real *out,
         unsigned fftw_rigor_flags = 0,
         unsigned packed_flags     = 0)
        : plan_(underling_fftw_plan_create_c2c_forward(
                    p.get(), long_ni, in, out, fftw_rigor_flags, packed_flags))
    {
        (void) tag; // unused
    }

    /** @see underling_fftw_plan_create_c2c_backward */
    plan(const c2c_backward tag,
         const problem &p,
         int long_ni,
         underling_real *in,
         underling_real *out,
         unsigned fftw_rigor_flags = 0,
         unsigned packed_flags     = 0)
        : plan_(underling_fftw_plan_create_c2c_backward(
                    p.get(), long_ni, in, out, fftw_rigor_flags, packed_flags))
    {
        (void) tag; // unused
    }

    /** @see underling_fftw_plan_create_r2c_forward */
    plan(const r2c_forward tag,
         const problem &p,
         int long_ni,
         underling_real *in,
         underling_real *out,
         unsigned fftw_rigor_flags = 0,
         unsigned packed_flags     = 0)
        : plan_(underling_fftw_plan_create_r2c_forward(
                    p.get(), long_ni, in, out, fftw_rigor_flags, packed_flags))
    {
        (void) tag; // unused
    }

    /** @see underling_fftw_plan_create_c2r_backward */
    plan(const c2r_backward tag,
         const problem &p,
         int long_ni,
         underling_real *in,
         underling_real *out,
         unsigned fftw_rigor_flags = 0,
         unsigned packed_flags     = 0)
        : plan_(underling_fftw_plan_create_c2r_backward(
                    p.get(), long_ni, in, out, fftw_rigor_flags, packed_flags))
    {
        (void) tag; // unused
    }

    /** @see underling_fftw_plan_create_inverse */
    plan(const plan& plan_to_invert,
         underling_real * in,
         underling_real * out,
         unsigned fftw_rigor_flags = 0)
        : plan_(underling_fftw_plan_create_inverse(
                    plan_to_invert.get(), in, out, fftw_rigor_flags))
    {};

    /** @see underling_fftw_plan_destroy */
    ~plan() { underling_fftw_plan_destroy(plan_); }

    /** @return The wrapped underling_fftw_plan instance. */
    underling_fftw_plan get() const { return plan_; }

    /** @see underling_fftw_local_extents_input */
    underling_fftw_extents local_extents_input() const
    {
        return underling_fftw_local_extents_input(plan_);
    }

    /** @see underling_fftw_local_extents_output */
    underling_fftw_extents local_extents_output() const
    {
        return underling_fftw_local_extents_output(plan_);
    }

    /** @see underling_fftw_local_input */
    int local_input(int *start  = NULL,
                    int *size   = NULL,
                    int *stride = NULL,
                    int *order  = NULL) const
    {
        return underling_fftw_local_input(plan_, start, size, stride, order);
    }

    /** @see underling_fftw_local_output */
    int local_output(int *start  = NULL,
                     int *size   = NULL,
                     int *stride = NULL,
                     int *order  = NULL) const {
        return underling_fftw_local_output(plan_, start, size, stride, order);
    }

    /** @return True if wrapped underling_fftw_plan instance is non-NULL. */
    operator bool () const { return plan_ != NULL; };

    /** @see underling_fftw_plan_execute */
    int execute(underling_real *in,
                underling_real *out) const
    {
        return underling_fftw_plan_execute(plan_, in, out);
    }

private:
    underling_fftw_plan plan_; /**< The wrapped underling_fftw_plan instance */
};

} // namespace fftw

} // namespace underling

/**
 * Outputs an underling_extents or underling::fftw::extents instance
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
        const underling::fftw::extents &e)
{
    return os << '['
              << e.start[0] << ',' << (e.start[0] + e.size[0])
              << ")x["
              << e.start[1] << ',' << (e.start[1] + e.size[1])
              << ")x["
              << e.start[2] << ',' << (e.start[2] + e.size[2])
              << ")x["
              << e.start[3] << ',' << (e.start[3] + e.size[3])
              << ")x["
              << e.start[4] << ',' << (e.start[4] + e.size[4])
              << ')';
}

/** @see underling_fftw_extents_cmp */
inline
bool operator==(const underling::fftw::extents &e1,
                const underling::fftw::extents &e2) {
    return !underling_fftw_extents_cmp(&e1, &e2);
}

#endif // UNDERLING_FFTW_HPP
