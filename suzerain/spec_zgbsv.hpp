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

#include <iosfwd>
#include <string>

namespace suzerain {

class spec_zgbsv
{
public:

    enum method_type { zgbsvx = 1, zgbsv, zcgbsvx } method;
    bool   reuse;
    bool   equil;
    int    aiter;
    int    siter;
    int    diter;
    double tolsc;

    spec_zgbsv();
    spec_zgbsv(const std::string &spec);
    operator std::string () const;
    bool operator==(const spec_zgbsv &that) const;

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
