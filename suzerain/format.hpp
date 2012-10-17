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
// format.hpp: Provides std::ostream-related formatting-related utilities
// $Id$

#ifndef SUZERAIN_FORMAT_HPP
#define SUZERAIN_FORMAT_HPP

#include <suzerain/common.hpp>

namespace suzerain {

template <
    typename FPT           = real_t,
    std::streamsize Digits = std::numeric_limits<FPT>::digits10,
    std::streamsize Width  = Digits + 5
>
class append_real_traits
{
public:

    typedef FPT type;

    static const std::streamsize digits = Digits;

    static const std::streamsize width = Width;

    static const FPT fixedmax;

    static const FPT fixedmin;

private:

    append_real_traits();

    append_real_traits(const append_real_traits&);

    append_real_traits& operator=(const append_real_traits&);

};

// Magic "2" is the width of a sign and a decimal point
template <typename FPT, std::streamsize Digits, std::streamsize Width>
const FPT append_real_traits<FPT,Digits,Width>::fixedmax
        = std::pow(FPT(10), FPT(Width) - Digits - 2);

// Magic 3 is the width of a sign, leading zero, and decimal point
template <typename FPT, std::streamsize Digits, std::streamsize Width>
const FPT append_real_traits<FPT,Digits,Width>::fixedmin
        = std::pow(FPT(10), -(FPT(Width) - Digits - 3));

//// template <typename Number, std::streamsize Digits, std::streamsize Width>
//// template <class CharT, class Traits>
//// std::basic_ostream<CharT,Traits>& append_real<Number,Digits,Width>::operator()(
////         std::basic_ostream<CharT,Traits>& os,
////         Number value)
//// {
////     // Format in fixed or scientific form as appropriate in given width
////     // Care taken to not perturb observable ostream state after function call
////     std::ios::fmtflags savedflags;
////     std::streamsize savedprec;
////     if (value >= fixedmin && value <= fixedmax) {
////         savedflags = os.setf(std::ios::fixed      | std::ios::right,
////                              std::ios::floatfield | std::ios::adjustfield);
////         savedprec = os.precision(append_real_prec);
////     } else {
////         savedflags = os.setf(std::ios::scientific | std::ios::right,
////                              std::ios::floatfield | std::ios::adjustfield);
////         savedprec = os.precision(append_real_width - 9);
////     }
////     os << std::setw(append_real_width) << static_cast<real_t>(value);
////     os.precision(savedprec);
////     os.setf(savedflags);
////
////     return os;
//// }

} // namespace suzerain

#endif // SUZERAIN_FORMAT_HPP
