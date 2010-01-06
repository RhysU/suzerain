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
 * math.hpp: provides general mathematics routines
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_MATH_HPP
#define __SUZERAIN_MATH_HPP

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Provides general mathematics routines.
 */
namespace math {

/**
 * Computes \f$x^n\f$ efficiently for small integer \f$n\f$, including
 * \f$n <= 0\f$.  No overflow checking is performed.  Algorithm taken
 * from the GNU Scientific Library's \c gsl_pow_int.
 *
 * @param x \f$x\f$
 * @param n \f$n\f$
 *
 * @return \f$x^n\f$
 */
template<typename FPT, typename Integral>
FPT integer_power(FPT x, Integral n)
{
    using std::numeric_limits;

    // Avoid shooting ourselves by accidentally requesting a negative power for
    // an integer input.  Long lines to ensure messages appear in compiler
    // error output when used improperly.
    BOOST_STATIC_ASSERT(numeric_limits<Integral>::is_integer);
    BOOST_STATIC_ASSERT(!numeric_limits<Integral>::is_signed || !numeric_limits<FPT>::is_integer);

    FPT retval = 1;
    // Convert all requests into one involving a positive power
    if (n < 0) {
        x = ((FPT) 1)/x;
        n = -n;
    }
    // Repeated squaring method.  Returns 0.0^0 = 1; continuous in x
    do {
        if (n & 1) retval *= x;  /* for n odd */
        n >>= 1;
        x *= x;
    } while (n);

    return retval;
}

} // namespace math

} // namespace suzerain

#endif // __SUZERAIN_MATH_HPP
