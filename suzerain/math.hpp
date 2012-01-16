/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010, 2011, 2012 The PECOS Development Team
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
 * Return the lesser of \c a and \c b favoring <tt>NaN</tt>s over numbers.
 * See <a href="http://en.wikipedia.org/wiki/IEEE_754_revision#min_and_max">
 * IEEE 754 revision on Wikipedia</a> for why <tt>std::min</tt> is not
 * usable when <tt>NaN</tt>s are preferred.
 */
template<class T>
inline
const T& minnan(const T& a, const T& b)
{
    return (a < b) || SUZERAIN_UNLIKELY((boost::math::isnan)(a)) ? a : b;
}

/**
 * Return the greater of \c a and \c b favoring <tt>NaN</tt>s over numbers.
 * See <a href="http://en.wikipedia.org/wiki/IEEE_754_revision#min_and_max">
 * IEEE 754 revision on Wikipedia</a> for why <tt>std::max</tt> is not
 * usable when <tt>NaN</tt>s are preferred.
 */
template<class T>
inline
const T& maxnan(const T& a, const T& b)
{
    return (a > b) || SUZERAIN_UNLIKELY((boost::math::isnan)(a)) ? a : b;
}

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
template<typename FPT, typename Integer>
typename boost::enable_if< boost::is_signed<Integer>, FPT >::type
integer_power(FPT x, Integer n)
{
    using std::numeric_limits;
    BOOST_STATIC_ASSERT(numeric_limits<Integer>::is_integer);

    // Prevent accidentally requesting a negative power for integer x.
    BOOST_STATIC_ASSERT(!numeric_limits<FPT>::is_integer);

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

/**
 * Computes \f$x^n\f$ efficiently for small integer \f$n\f$.
 * No overflow checking is performed.  Algorithm taken
 * from the GNU Scientific Library's \c gsl_pow_int.
 *
 * @param x \f$x\f$
 * @param n \f$n\f$
 *
 * @return \f$x^n\f$
 */
template<typename FPT, typename Integer>
typename boost::disable_if< boost::is_signed<Integer>, FPT >::type
integer_power(FPT x, Integer n)
{
    using std::numeric_limits;
    BOOST_STATIC_ASSERT(numeric_limits<Integer>::is_integer);

    FPT retval = 1;
    // Repeated squaring method.  Returns 0.0^0 = 1; continuous in x
    do {
        if (n & 1) retval *= x;  /* for n odd */
        n >>= 1;
        x *= x;
    } while (n);

    return retval;
}


namespace {

template<typename T, int N> struct impl_fixed_integer_power;

template<typename T> struct impl_fixed_integer_power<T,0> {
    SUZERAIN_FORCEINLINE static T fixed_integer_power(const T t) {
        SUZERAIN_UNUSED(t);
        return 1;
    }
};

template<typename T> struct impl_fixed_integer_power<T,1> {
    SUZERAIN_FORCEINLINE static T fixed_integer_power(const T t) {
        return t;
    }
};

template<typename T> struct impl_fixed_integer_power<T,2> {
    SUZERAIN_FORCEINLINE static T fixed_integer_power(const T t) {
        return t*t;
    }
};

template<typename T> struct impl_fixed_integer_power<T,3> {
    SUZERAIN_FORCEINLINE static T fixed_integer_power(const T t) {
        return t*t*t;
    }
};

template<typename T> struct impl_fixed_integer_power<T,4> {
    SUZERAIN_FORCEINLINE static T fixed_integer_power(const T t) {
        T tmp = t*t;
        return tmp*tmp;
    }
};

} // anonymous

/**
 * Compute <tt>t^N</tt> for small, fixed integer powers.
 *
 * @param t to process.
 *
 * @return <tt>t^N</tt>.
 */
template<int N, typename T>
T fixed_integer_power(const T t)
{
    return impl_fixed_integer_power<T,N>::fixed_integer_power(t);
}

/**
 * Output \n linearly spaced values spanning the range <tt>[xbegin, xend]</tt>
 * (inclusive).
 *
 * @param xbegin Beginning value
 * @param xend   Ending value
 * @param n      Number of linearly spaced values to use.  Must be
 *               nonnegative.
 * @param x      Output locations
 * @return One plus that last output location.
 */
template<typename FPT, typename Integer, typename OutputIterator>
OutputIterator linspace(const FPT xbegin,
                        const FPT xend,
                        const Integer n,
                        OutputIterator x)
{
    BOOST_STATIC_ASSERT(std::numeric_limits<Integer>::is_integer);
    if (SUZERAIN_UNLIKELY(n <= 0)) throw std::invalid_argument("n <= 0");

    if (SUZERAIN_UNLIKELY(n == 1)) {
        if (SUZERAIN_UNLIKELY(xbegin != xend))
            throw std::invalid_argument("n == 1 && xbegin != xend");
        *x++ = xbegin;
    } else {
        const FPT xh = (xend - xbegin)/(n-1);
        *x++ = xbegin;
        for (Integer i = 1; i < n-1; ++i) {
            *x++ = xh*i + xbegin;
        }
        *x++ = xend;
    }

    return x;
}

/**
 * Output \n logarithmically spaced values spanning the range
 * <tt>[pow(base,xbegin), pow(base,xend)]</tt> (inclusive).
 *
 * @param xbegin Beginning value
 * @param xend   Ending value
 * @param n      Number of linearly spaced values to use.  Must be
 *               nonnegative.
 * @param x      Output locations
 * @param base   Exponent base
 * @return One plus that last output location.
 */
template<typename FPT, typename Integer, typename OutputIterator>
OutputIterator logspace(const FPT xbegin,
                        const FPT xend,
                        const Integer n,
                        OutputIterator x,
                        const FPT base = 10)
{
    BOOST_STATIC_ASSERT(std::numeric_limits<Integer>::is_integer);
    if (SUZERAIN_UNLIKELY(n <= 0)) throw std::invalid_argument("n <= 0");

    if (SUZERAIN_UNLIKELY(n == 1)) {
        if (SUZERAIN_UNLIKELY(xbegin != xend))
            throw std::invalid_argument("n == 1 && xbegin != xend");
        *x++ = std::pow(base, xbegin);
    } else {
        const FPT xh = (xend - xbegin)/(n-1);
        *x++ = pow(base, xbegin);
        for (Integer i = 1; i < n-1; ++i) {
            *x++ = std::pow(base, xh*i + xbegin);
        }
        *x++ = std::pow(base, xend);
    }

    return x;
}

/**
 * Output \n values spanning the range <tt>[xbegin, xend]</tt> (inclusive)
 * stretched linearly according to \f$\Delta{}x_\text{end} = \alpha
 * \Delta{}x_\text{begin}\f$.
 *
 * @param xbegin Beginning value
 * @param xend   Ending value
 * @param n      Number of linearly spaced values to use.  Must be
 *               greater than two.
 * @param alpha  The stretching factor relating the ratio of the last
 *               interval width to the first interval width.  Must be
 *               greater than zero.
 * @param x      Output locations
 * @return One plus that last output location.
 */
template<typename FPT, typename Integer, typename OutputIterator>
OutputIterator stretchspace(const FPT xbegin,
                            const FPT xend,
                            const Integer n,
                            const FPT alpha,
                            OutputIterator x)
{
    BOOST_STATIC_ASSERT(std::numeric_limits<Integer>::is_integer);
    if (SUZERAIN_UNLIKELY(n < 3))
        throw std::invalid_argument("n < 3");
    if (SUZERAIN_UNLIKELY(alpha <= 0))
        throw std::invalid_argument("alpha <= 0");

    // Compute parameters; see notes dated 6 April 2010.  Expressions arise
    // from setting dx_{n-1} = alpha * dx_{0} = beta*(n-1) + dx_{0}, expressing
    // beta in terms of dx_{0}, and solving b-a = \sum_{i=0}^{ninterval} dx_{i}
    // for dx0.
    const Integer ninterval = n - 1;
    const FPT dx0           = (2*(xend - xbegin)) / (ninterval*(1 + alpha));
    const FPT beta          = (alpha - 1) / (ninterval - 1) * dx0;

    // Use parameters to compute the output sequence
    FPT xlast = xbegin;
    for (Integer i = 0; i < ninterval; ++i) {
        *x++ = xlast;
        xlast += beta*i + dx0;
    }
    *x++ = xend;

    return x;
}

} // namespace math

} // namespace suzerain

#endif // __SUZERAIN_MATH_HPP
