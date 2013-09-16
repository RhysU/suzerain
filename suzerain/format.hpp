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

#ifndef SUZERAIN_FORMAT_HPP
#define SUZERAIN_FORMAT_HPP

/** @file
 * Provides std::ostream-related formatting utilities
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * A stream insertion wrapper to provide full precision floating point
 * output.  The wrapper does not perturb observable stream state.
 * Scientific notation is avoided for small-magnitude quantities.
 * Padding is added to facilitate table-like output.
 *
 * A simple use case with the generated output follows:
 * \code
 * cout << fullprec<>(1.23456789012345678e-5) << '\n'; // "   1.23456789012e-05"
 * cout << fullprec<>(1.23456789012345678e-4) << '\n'; // "   1.23456789012e-04"
 * cout << fullprec<>(1.23456789012345678e-3) << '\n'; // "   1.23456789012e-03"
 * cout << fullprec<>(1.23456789012345678e-2) << '\n'; // "   0.012345678901235"
 * cout << fullprec<>(1.23456789012345678e-1) << '\n'; // "   0.123456789012346"
 * cout << fullprec<>(1.23456789012345678e+0) << '\n'; // "   1.234567890123457"
 * cout << fullprec<>(1.23456789012345678e+1) << '\n'; // "  12.345678901234567"
 * cout << fullprec<>(1.23456789012345678e+2) << '\n'; // " 123.456789012345681"
 * cout << fullprec<>(1.23456789012345678e+3) << '\n'; // "   1.23456789012e+03"
 * cout << fullprec<>(1.23456789012345678e+4) << '\n'; // "   1.23456789012e+04"
 * cout << fullprec<>(1.23456789012345678e+5) << '\n'; // "   1.23456789012e+05"
 * \endcode
 *
 * @tparam FPT    Target floating point type.
 * @tparam Digits Number of digits to consider as "full precision".
 *                Default uses <code>std::numeric_limits</code>.
 * @tparam Width  Right-justified output width to use.
 *                Default permits leading sign and exponential notation.
 */
template <
    typename FPT           = real_t,
    std::streamsize Digits = std::numeric_limits<FPT>::digits10,
    std::streamsize Width  = Digits + 5
>
class fullprec
{
public:

    /** Number of digits to output. */
    static const std::streamsize digits = Digits;

    /** Right-justified output width to use. */
    static const std::streamsize width = Width;

    /** Maximum value that will be output in fixed decimal format. */
    static const FPT fixedmax;

    /** Minimum value that will be output in fixed decimal format. */
    static const FPT fixedmin;

    /**
     * Construct an instant outputting \c val.
     *
     * The lifetime of \c val must be longer than the lifetime
     * of the constructed instance.
     *
     * @param val Floating point value to output.
     */
    explicit fullprec(const FPT& val) : val(val) {}

    /**
     * Output the value provided at construction on the given stream.
     * See \ref fullprec for example usage.
     *
     * @tparam CharT  Per <code>std::basic_ostream</code>
     * @tparam Traits Per <code>std::basic_ostream</code>
     * @param  os     Stream on which to output the value.
     * @param  fp     Instance housing the data to be output.
     *
     * @return \c os
     */
    template <class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>& operator<<(
        std::basic_ostream<CharT,Traits>& os,
        const fullprec<FPT,Digits,Width>& fp)
    {
        // Format in fixed or scientific form as appropriate in given width
        // Care taken to not perturb observable state of os after return
        std::ios::fmtflags savedflags;
        std::streamsize savedprec;
        if (fp.val >= fp.fixedmin && fp.val <= fp.fixedmax) {
            savedflags = os.setf(std::ios::fixed      | std::ios::right,
                                 std::ios::floatfield | std::ios::adjustfield);
            savedprec  = os.precision(fp.digits);
        } else {
            savedflags = os.setf(std::ios::scientific | std::ios::right,
                                 std::ios::floatfield | std::ios::adjustfield);
            savedprec  = os.precision(fp.width - 9);
        }
        os << std::setw(fp.width) << fp.val;
        os.precision(savedprec);
        os.setf(savedflags);

        return os;
    }

private:

    /** A reference to the value provided at construction time. */
    const FPT& val;

    /** Private to prevents generated assignment operator. */
    fullprec& operator=(const fullprec&);
};

// Magic "2" is the width of a sign and a decimal point
template <typename FPT, std::streamsize Digits, std::streamsize Width>
const FPT fullprec<FPT,Digits,Width>::fixedmax
        = std::pow(FPT(10), FPT(Width) - Digits - 2);

// Magic 3 is the width of a sign, leading zero, and decimal point
template <typename FPT, std::streamsize Digits, std::streamsize Width>
const FPT fullprec<FPT,Digits,Width>::fixedmin
        = std::pow(FPT(10), -(FPT(Width) - Digits - 3));

} // namespace suzerain

#endif // SUZERAIN_FORMAT_HPP
