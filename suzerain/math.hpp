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

#ifndef SUZERAIN_MATH_HPP
#define SUZERAIN_MATH_HPP

/** @file
 * Provides general mathematical routines.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Provides general mathematical routines.
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

/**
 * Variants of the bump function \f$ \exp\left(-\frac{1}{1-x^2}\right) \f$.
 * @see For example, http://en.wikipedia.org/wiki/Bump_function.
 */
namespace bump {

/**
 * Evaluate the classic \f$ \exp\left(-\frac{1}{1-x^2}\right) \f$
 * with support defined to be only \f$ x \in \left(-1, 1\right) \f$.
 *
 * @param x Location at which to evaluate the bump function.
 */
template<typename FPT>
FPT classic(const FPT x)
{
    if (x <= -1 || 1 <= x) {
        return 0;
    } else {
        using std::exp;
        return exp(-1 / (1 - x*x));
    }
}

/**
 * Evaluate a bump function on \f$ x \in \left(-1, 1\right) \f$ rescaled
 * to give \f$1\f$ when \f$ x = 0 \f$.
 *
 * @param x Location at which to evaluate the bump function.
 */
template<typename FPT>
FPT scaled(const FPT x)
{
    if (x <= -1 || 1 <= x) {
        return 0;
    } else {
        using std::exp;
        const FPT x2 = x*x;
        return exp(x2 / (x2 - 1));
    }
}

/**
 * Evaluate a bump function on \f$ x \in \left(-1, 1\right) \f$ rescaled
 * to give \f$p\f$ when \f$ x = 0 \f$.
 *
 * @param x Location at which to evaluate the bump function.
 * @param p Value to return when \f$ x = 0 \f$.
 */
template<typename FPT>
FPT scaled(const FPT x, const FPT p)
{
    if (x <= -1 || 1 <= x) {
        return 0;
    } else {
        using std::exp;
        const FPT x2 = x*x;
        return p * exp(x2 / (x2 - 1));
    }
}

/**
 * @copydoc scaled(const FPT,const FPT)
 *
 * @param q Value to return outside \f$ \left(-1, 1\right) \f$.
 */
template<typename FPT>
FPT scaled(const FPT x, const FPT p, const FPT q)
{
    if (x <= -1 || 1 <= x) {
        return q;
    } else {
        using std::exp;
        const FPT x2 = x*x;
        return (p - q) * exp(x2 / (x2 - 1)) + q;
    }
}

/**
 * @copydoc scaled(const FPT,const FPT,const FPT)
 *
 * @param n Evaulate the \f$n\f$th power of the bump function.
 */
template<typename FPT>
FPT scaled(const FPT x, const FPT p, const FPT q, const FPT n)
{
    if (x <= -1 || 1 <= x) {
        return q;
    } else {
        using std::exp;
        const FPT x2 = x*x;
        return (p - q) * exp(n * (x2 / (x2 - 1))) + q;
    }
}

/**
 * Evaluate a bump function on \f$ x \in \left(l, r\right) \f$ rescaled
 * to give \f$1\f$ when \f$ x = \frac{l+r}{2} \f$.
 *
 * @param x Location at which to evaluate the bump function.
 * @param l Left boundary, exclusive, of the function's support.
 * @param r Right boundary, exclusive, of the function's support.
 */
template<typename FPT>
FPT shifted(const FPT x, const FPT l, const FPT r)
{
    if (x <= l || r <= x) {
        return 0;
    } else {
        using std::exp;
        FPT t = l + r - 2 * x;
        t *= t;
        t /= 4 * (x - l) * (x - r);
        return exp(t);
    }
}

/**
 *
 * Evaluate a bump function on \f$ x \in \left(l, r\right) \f$ rescaled
 * to give \f$p\f$ when \f$ x = \frac{l+r}{2} \f$.
 *
 * @param p Value to return when \f$ x = \frac{l+r}{2} \f$.
 */
template<typename FPT>
FPT shifted(const FPT x, const FPT l, const FPT r, const FPT p)
{
    if (x <= l || r <= x) {
        return 0;
    } else {
        using std::exp;
        FPT t = l + r - 2 * x;
        t *= t;
        t /= 4 * (x - l) * (x - r);
        return p * exp(t);
    }
}

/**
 * @copydoc shifted(const FPT,const FPT,const FPT,const FPT)
 *
 * @param q Value to return outside \f$ \left(l, r\right) \f$.
 */
template<typename FPT>
FPT shifted(const FPT x, const FPT l, const FPT r,
            const FPT p, const FPT q)
{
    if (x <= l || r <= x) {
        return q;
    } else {
        using std::exp;
        FPT t = l + r - 2 * x;
        t *= t;
        t /= 4 * (x - l) * (x - r);
        return (p - q) * exp(t) + q;
    }
}

/**
 * @copydoc shifted(const FPT,const FPT,const FPT,const FPT,const FPT)
 *
 * @param n Evaulate the \f$n\f$th power of the bump function.
 */
template<typename FPT>
FPT shifted(const FPT x, const FPT l, const FPT r,
            const FPT p, const FPT q, const FPT n)
{
    if (x <= l || r <= x) {
        return q;
    } else {
        using std::exp;
        FPT t = l + r - 2 * x;
        t *= t;
        t /= 4 * (x - l) * (x - r);
        return (p - q) * exp(n * t) + q;
    }
}

} // namespace bump

/**
 * Routines for evaluating polynomial-based interpolants.
 *
 * @see Mathematica's <a
 * href="http://reference.wolfram.com/mathematica/ref/InterpolatingPolynomial.html">
 * InterpolatingPolynomial</a> for an easy way to generate these and other
 * similar results.
 */
namespace interpolate {

/**
 * Interpolate a function given values at two distinct points.
 * This is quite obviously a line.  It has been coded up to be consistent
 * with \ref value_deriv.  The interpolated results are exact for
 * constant and linear functions.
 *
 * @param[in ]  x1 The first distinct point.
 * @param[in ]  y1 The function value at \c x1.
 * @param[in ]  x2 The second distinct point.
 * @param[in ]  y2 The function value at \c x2.
 * @param[in ]  x3 The point of interest.
 * @param[out]  y3 The function value at \c x3.
 */
template<typename FPT>
void value(const FPT&  x1,
           const FPT&  y1,
           const FPT&  x2,
           const FPT&  y2,
           const FPT&  x3,
           FPT&        y3)
{
    y3 = y1 + (y2 - y1) / (x2 - x1) * (x3 - x1);
}

/**
 * Interpolate a function and its derivative given values and derivatives at
 * two distinct points.  The expressions used are all derived from
 * interrogating <code>InterpolatingPolynomial[{ x1, y1, yp1}, {x2, y2, yp2}},
 * x3] </code> with Mathematica.  The interpolated results are exact
 * for constants through cubic polynomials, inclusive.
 *
 * @param[in ]  x1 The first distinct point.
 * @param[in ]  y1 The function value at \c x1.
 * @param[in ] yp1 The function derivative at \c x1.
 * @param[in ]  x2 The second distinct point.
 * @param[in ]  y2 The function value at \c x2.
 * @param[in ] yp2 The function derivative at \c x2.
 * @param[in ]  x3 The point of interest.
 * @param[out]  y3 The function value at \c x3.
 * @param[out] yp3 The function derivative at \c x3.
 */
template<typename FPT>
void value_deriv(const FPT&  x1,
                 const FPT&  y1,
                 const FPT& yp1,
                 const FPT&  x2,
                 const FPT&  y2,
                 const FPT& yp2,
                 const FPT&  x3,
                 FPT&        y3,
                 FPT&       yp3)
{
    const FPT x2m1      = x2 - x1;
    const FPT x2m1p2    = x2m1   * x2m1;
    const FPT x2m1p3    = x2m1p2 * x2m1;
    const FPT invx2m1p3 = 1 / x2m1p3;
    const FPT x3m1      = x3 - x1;
    const FPT x3m1p2    = x3m1 * x3m1;
    const FPT x3m2      = x3 - x2;
    const FPT y2m1      = y2 - y1;

    y3  = invx2m1p3*(  x2m1p3*(y1 + x3m1*yp1)
                     + x3m1p2*(- 2*x3m2*y2m1 - x2m1p2*yp1
                               + x2m1*(y2m1 + x3m2*(yp1 + yp2))));
    yp3 = invx2m1p3*(  x2m1p3*yp1
                     - 2*x3m1*(x3m1*y2m1 + 2*x3m2*y2m1 + x2m1p2*yp1)
                     + x2m1*x3m1*(2*y2m1 + (x3m1 + 2*x3m2)*(yp1 + yp2))) ;
}

} // namespace interpolate

} // namespace math

} // namespace suzerain

#endif // SUZERAIN_MATH_HPP
