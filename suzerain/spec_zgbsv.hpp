//--------------------------------------------------------------------------
//
// Copyright (C) 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// spec_zgbsv.hpp: Encapsulates complex-valued banded solve options
// $Id$

#ifndef SUZERAIN_SPEC_ZGBSV_HPP
#define SUZERAIN_SPEC_ZGBSV_HPP

#include <cassert>
#include <iosfwd>
#include <string>

namespace suzerain {

class spec_zgbsv
{
public:

    enum method_type { zgbsvx = 1, zgbsv, zcgbsvx };

    spec_zgbsv();
    spec_zgbsv(const std::string &spec);

    method_type method() const { return method_; }
    bool        equil()  const { return equil_;  }
    bool        reuse()  const { return reuse_;  }
    int         aiter()  const { return aiter_;  }
    int         siter()  const { return siter_;  }
    int         diter()  const { return diter_;  }
    double      tolsc()  const { return tolsc_;  }

    operator std::string () const;
    bool operator==(const spec_zgbsv &that) const;

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
        const spec_zgbsv& s)
{
    return os << static_cast<std::string>(s);
}

} // namespace suzerain

#endif // SUZERAIN_SPEC_ZGBSV_HPP
