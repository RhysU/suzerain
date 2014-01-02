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

/** @file
 * @copydoc fftw.hpp
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/fftw.hpp>

#include <cctype>

namespace suzerain {

namespace fftw {

rigor rigor_from(const char *name)
{
    using std::tolower;
    if (name) {
        switch (tolower(name[0])) {
        default:  break;
        case 'm': return measure;
        case 'p': return patient;
        case 'e':
            switch (tolower(name[1])) {
            default:  // Assume estimate
            case 's': return estimate;
            case 'x': return exhaustive;
            }
        case 'w': return wisdom_only;
        }
    }

    return measure;
}

rigor rigor_from(const unsigned flags)
{
    if (flags & FFTW_WISDOM_ONLY) {
        return wisdom_only;
    } else if (flags & FFTW_EXHAUSTIVE) {
        return exhaustive;
    } else if (flags & FFTW_PATIENT) {
        return patient;
    } else if (flags & FFTW_ESTIMATE) {
        return estimate;
    } else {
        return measure; // Default
    }
}

const char * c_str(const rigor r)
{
    switch (r) {
    case estimate:    return "estimate";
    case patient:     return "patient";
    case exhaustive:  return "exhaustive";
    case wisdom_only: return "wisdom_only";
    default:
    case measure:     return "measure";    // Default
    }
}

} // namespace fftw

} // namespace suzerain
