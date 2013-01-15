//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

    enum method_type { zgbsvx = 1, zgbsv, zcgbsvx };

    zgbsv_specification();
    zgbsv_specification(const std::string &spec);

    method_type method() const { return method_; }
    bool        equil()  const { return equil_;  }
    bool        reuse()  const { return reuse_;  }
    int         aiter()  const { return aiter_;  }
    int         siter()  const { return siter_;  }
    int         diter()  const { return diter_;  }
    double      tolsc()  const { return tolsc_;  }

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
