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
 * fftw.hpp: miscellaneous utilities for working with FFTW's C interface
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_FFTW_HPP
#define __SUZERAIN_FFTW_HPP

#include <suzerain/common.hpp>
#include <fftw3.h>

namespace suzerain {

/**
 * Provides miscellaneous utilities for working with FFTW's C interface.
 */
namespace fftw {

enum rigor {
    estimate    = FFTW_ESTIMATE,
    measure     = FFTW_MEASURE,         // Default
    patient     = FFTW_PATIENT,
    exhaustive  = FFTW_EXHAUSTIVE
};

rigor rigor_from(const char *name);

rigor rigor_from(unsigned flags);

const char * c_str(rigor r);

template< typename charT, typename traits >
::std::basic_ostream<charT,traits>& operator<<(
        ::std::basic_ostream<charT,traits> &os,
        const rigor &r)
{
    return os << c_str(r);
}

} // namespace fftw

} // namespace suzerain

#endif // __SUZERAIN_FFTW_HPP
