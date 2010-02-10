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

/**
 * @file Provides C++ wrappers for the C-based API in underling_fft.h.  In
 * particular, provides RAII semantics for opaque types.
 */

namespace suzerain {

namespace underling {

/**
 * Provides a thin RAII wrapper for underling_fftplan.
 * @see underling_fftplan.
 */
class fftplan {
public:

    /** A tag type used to indicate a complex-to-complex forward transform */
    struct c2c_forward  {};

    /** A tag type used to indicate a complex-to-complex backward transform */
    struct c2c_backward {};

    /** A tag type used to indicate a real-to-complex forward transform */
    struct r2c_forward  {};

    /** A tag type used to indicate a complex-to-real backward transform */
    struct c2r_backward {};

    /** @see underling_fftplan_create_c2c_forward */
    fftplan(const c2c_forward tag,
            const problem &p,
            int long_ni,
            underling_real *data,
            unsigned fftw_rigor_flags)
        : fftplan_(underling_fftplan_create_c2c_forward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fftplan_create_c2c_backward */
    fftplan(const c2c_backward tag,
            const problem &p,
            int long_ni,
            underling_real *data,
            unsigned fftw_rigor_flags)
        : fftplan_(underling_fftplan_create_c2c_backward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fftplan_create_r2c_forward */
    fftplan(const r2c_forward tag,
            const problem &p,
            int long_ni,
            underling_real *data,
            unsigned fftw_rigor_flags)
        : fftplan_(underling_fftplan_create_r2c_forward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fftplan_create_c2r_backward */
    fftplan(const c2r_backward tag,
            const problem &p,
            int long_ni,
            underling_real *data,
            unsigned fftw_rigor_flags)
        : fftplan_(underling_fftplan_create_c2r_backward(
                    p.get(), long_ni, data, fftw_rigor_flags)) {};

    /** @see underling_fftplan_destroy */
    ~fftplan() { underling_fftplan_destroy(fftplan_); }

    /** @return The wrapped underling_fftplan instance. */
    const underling_fftplan get() const { return fftplan_; }

    /** @return True if the wrapped underling_fftplan instance is non-NULL. */
    operator bool () const { return fftplan_ != NULL; };

    /** @see underling_fftplan_execute */
    int execute() const {
        return underling_fftplan_execute(fftplan_);
    }

private:
    underling_fftplan fftplan_; /**< The wrapped underling_fftplan instance */
};

} // namespace underling

} // namespace suzerain

#endif // __SUZERAIN_UNDERLING_FFT_HPP
