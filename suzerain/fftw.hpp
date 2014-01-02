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

#ifndef SUZERAIN_FFTW_HPP
#define SUZERAIN_FFTW_HPP

/** @file
 * C++ utilities for working with FFTW's C-based API.
 */

#include <ostream>

#include <fftw3.h>

namespace suzerain {

/**
 * Provides miscellaneous utilities for working with FFTW's C interface.
 */
namespace fftw {

/**
 * A type safe enum for working with FFTW planning flags.  Meant as a
 * thin veneer to allow conversion to/from human-readable strings.
 *
 * @see See the <a href="http://www.fftw.org/fftw3_doc/Planner-Flags.html">
 *      FFTW manual section on planner flags</a> for more details.
 */
enum rigor {
    /** Corresponds to \c FFTW_ESTIMATE */
    estimate    = FFTW_ESTIMATE,

    /** Corresponds to \c FFTW_MEASURE and is the default value. */
    measure     = FFTW_MEASURE,

    /** Corresponds to \c FFTW_PATIENT */
    patient     = FFTW_PATIENT,

    /** Corresponds to \c FFTW_EXHAUSTIVE */
    exhaustive  = FFTW_EXHAUSTIVE,

    /** Corresponds to \c FFTW_WISDOM_ONLY */
    wisdom_only  = FFTW_WISDOM_ONLY
};

/**
 * Convert from a human-readable name into a rigor flag.  Comparisons are case
 * insensitive and at most only the first two characters are examined.  This
 * causes, for example, 'p', 'P', 'patIENT', and 'pAtRiArch' to all be treated
 * as patient.  Ambiguous values are treated towards less costly planning, e.g.
 * 'e' is treated as estimate and not exhaustive.  Unknown strings are
 * converted to measure.
 *
 * @param name to convert.
 *
 * @return The corresponding rigor value.
 */
rigor rigor_from(const char *name);

/**
 * Extract a rigor instance from FFTW's <tt>define</tt>d flags.
 *
 * @param flags to convert.
 *
 * @return The corresponding rigor value.
 */
rigor rigor_from(const unsigned flags);

/**
 * Return a human-readable string representing the given rigor.  Returned
 * values should be "round trip" capable, meaning they can be fed to rigor_from
 * to determine \c r again.
 *
 * @param r for which to provide a string.
 *
 * @return a human-readable string.
 */
const char * c_str(const rigor r);

/**
 * Retrieve a sensible default number of threads for threaded FFTW planning.
 * Takes into account whether or not threading is available and the desired
 * threading method.  For example, it uses <tt>OMP_NUM_THREADS</tt> when OpenMP
 * threading is used.
 *
 * @return a sensible default number of threads.
 */
int default_nthreads();

/**
 * Enable output of a rigor on a <tt>basic_ostream</tt>.
 *
 * @param os on which to output \c r.
 * @param r value to output.
 *
 * @return the modified output stream.
 */
template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os,
        const rigor &r)
{
    return os << c_str(r);
}

} // namespace fftw

} // namespace suzerain

#endif // SUZERAIN_FFTW_HPP
