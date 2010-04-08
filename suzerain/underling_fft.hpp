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
 * underling_fft.hpp: C++ wrappers for the C-based underling_fft API
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_UNDERLING_FFT_HPP
#define __SUZERAIN_UNDERLING_FFT_HPP

#include <suzerain/underling_fft.h>
#include <suzerain/underling.hpp>

/** @file
 * Provides C++ wrappers for the C-based API in underling_fft.h.  In
 * particular, provides RAII semantics for opaque types.
 */

namespace suzerain {

namespace underling {

namespace fft {

/** @see underling_fft_extents */
typedef underling_fft_extents extents;

/**
 * Provides a thin RAII wrapper for underling_fft_plan.
 * @see underling_fft_plan.
 */
class plan : public boost::noncopyable {
public:

    /** A tag type used to indicate a complex-to-complex forward transform */
    struct c2c_forward  {};

    /** A tag type used to indicate a complex-to-complex backward transform */
    struct c2c_backward {};

    /** A tag type used to indicate a real-to-complex forward transform */
    struct r2c_forward  {};

    /** A tag type used to indicate a complex-to-real backward transform */
    struct c2r_backward {};

    /** @see underling_fft_plan_create_c2c_forward */
    plan(const c2c_forward tag,
         const problem &p,
         int long_ni,
         underling_real *data,
         unsigned fftw_rigor_flags)
        : plan_(underling_fft_plan_create_c2c_forward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fft_plan_create_c2c_backward */
    plan(const c2c_backward tag,
         const problem &p,
         int long_ni,
         underling_real *data,
         unsigned fftw_rigor_flags)
        : plan_(underling_fft_plan_create_c2c_backward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fft_plan_create_r2c_forward */
    plan(const r2c_forward tag,
         const problem &p,
         int long_ni,
         underling_real *data,
         unsigned fftw_rigor_flags)
        : plan_(underling_fft_plan_create_r2c_forward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fft_plan_create_c2r_backward */
    plan(const c2r_backward tag,
         const problem &p,
         int long_ni,
         underling_real *data,
         unsigned fftw_rigor_flags)
        : plan_(underling_fft_plan_create_c2r_backward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fft_plan_create_inverse */
    plan(const plan& plan_to_invert,
         underling_real * data,
         unsigned fftw_rigor_flags)
        : plan_(underling_fft_plan_create_inverse(
                    plan_to_invert.get(), data, fftw_rigor_flags)) {};

    /** @see underling_fft_plan_destroy */
    ~plan() { underling_fft_plan_destroy(plan_); }

    /** @return The wrapped underling_fft_plan instance. */
    const underling_fft_plan get() const { return plan_; }

    /** @see underling_fft_local_extents_input */
    underling_fft_extents local_extents_input() const {
        return underling_fft_local_extents_input(plan_);
    }

    /** @see underling_fft_local_extents_output */
    underling_fft_extents local_extents_output() const {
        return underling_fft_local_extents_output(plan_);
    }

    /** @see underling_fft_local_input */
    int local_input(int *start  = NULL,
                    int *size   = NULL,
                    int *stride = NULL,
                    int *order  = NULL) const {
        return underling_fft_local_input(plan_, start, size, stride, order);
    }

    /** @see underling_fft_local_output */
    int local_output(int *start  = NULL,
                     int *size   = NULL,
                     int *stride = NULL,
                     int *order  = NULL) const {
        return underling_fft_local_output(plan_, start, size, stride, order);
    }

    /** @return True if the wrapped underling_fft_plan instance is non-NULL. */
    operator bool () const { return plan_ != NULL; };

    /** @see underling_fft_plan_execute */
    int execute() const {
        return underling_fft_plan_execute(plan_);
    }

private:
    underling_fft_plan plan_; /**< The wrapped underling_fft_plan instance */
};

} // namespace fft

} // namespace underling

} // namespace suzerain

/**
 * Outputs an underling_extents or suzerain::underling::fft::extents instance
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
        const suzerain::underling::fft::extents &e)
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

/** @see underling_fft_extents_cmp */
bool operator==(const suzerain::underling::fft::extents &e1,
                const suzerain::underling::fft::extents &e2) {
    return !underling_fft_extents_cmp(&e1, &e2);
}

#endif // __SUZERAIN_UNDERLING_FFT_HPP
