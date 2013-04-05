//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#ifndef SUZERAIN_ZGBSV_SPECIFICATION_HPP
#define SUZERAIN_ZGBSV_SPECIFICATION_HPP

/** @file
 * Encapsulates parsing and storing complex-valued banded solve options.
 */

#include <cassert>
#include <iosfwd>
#include <string>

namespace suzerain {

/** Encapsulates parsing and storing complex-valued banded solve options. */
class zgbsv_specification
{
public:

    /** Which general banded solver is specified? */
    enum method_type { zgbsvx = 1, zgbsv, zcgbsvx };

    /** Construct a default solver specification */
    zgbsv_specification();

    /** Construct a solver specification by parsing the given \c spec. */
    zgbsv_specification(const std::string &spec);

    method_type  method()   const { return method_; }
    bool         equil()    const { return equil_;  }
    bool         reuse()    const { return reuse_;  }
    int          aiter()    const { return aiter_;  }
    int          siter()    const { return siter_;  }
    int          diter()    const { return diter_;  }
    double       tolsc()    const { return tolsc_;  }
    bool         in_place() const;
    const char * mname()    const;

    operator std::string () const;
    bool operator==(const zgbsv_specification &that) const;

private:

    method_type method_;
    bool        equil_;
    bool        reuse_;
    int         aiter_;
    int         siter_;
    int         diter_;
    double      tolsc_;

};

template< typename CharT, typename Traits >
std::basic_ostream<CharT,Traits>& operator<<(
        std::basic_ostream<CharT,Traits> &os,
        const zgbsv_specification& s)
{
    return os << static_cast<std::string>(s);
}

} // namespace suzerain

#endif // SUZERAIN_ZGBSV_SPECIFICATION_HPP
