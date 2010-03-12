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
 * complex.hpp: utilities for manipulating various complex number formats
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_COMPLEX_HPP
#define __SUZERAIN_COMPLEX_HPP

#include <suzerain/common.hpp>
#include <suzerain/functional.hpp>

// TODO Broken details::assign_* if FFTW3 discovers the C99 _Complex type
// This header does not currently require fftw3.h to function, but fixing
// the above TODO may change that fact.

namespace suzerain {

/**
 * Provides utilities for manipulating various complex number types.
 * Includes capabilities for manipulating C++03 and FFTW3 complex types.
 */
namespace complex {

/**
 * Provides complex type traits suitable for manipulating C++03 and FFTW3
 * complex types.  Many of these build atop the Boost.TypeTraits framework.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/type_traits/">
 * Boost.TypeTraits</a> for more details.
 */
namespace traits {

/**
 * An <tt>is_complex</tt> type trait similar to <tt>boost::is_complex</tt>.  If
 * the template parameter \c T is complex, then <tt>is_complex</tt> will
 * inherit from <tt>boost::true_type</tt>.  Otherwise it inherits from
 * <tt>boost::false_type</tt>.
 */
template<class T, class Enable = void>
struct is_complex
        : public boost::is_complex<T> {};

/**
 * A specialization of <tt>is_complex</tt> extended to recognize FFTW-like
 * complex types.
 */
template<class T>
struct is_complex<T[2],
                  typename boost::enable_if<boost::is_arithmetic<T> >::type>
        : public boost::true_type {};

/**
 * A type trait that, given a valid complex-value type, returns the type of the
 * real scalar component as <tt>type</tt>.
 */
template<class T> struct real {};

/** A specialization to handle <tt>std::complex<T></tt>. */
template<class T> struct real<std::complex<T> > {
    typedef typename std::complex<T>::value_type type;
};

/** A specialization to handle FFTW-like types. */
template<class T> struct real<T[2]> {
    typedef T type;
};

/**
 * An type trait similar to <tt>is_complex</tt> which additionally checks that
 * the underlying real type is <tt>float</tt>.  Inherits from
 * <tt>boost::true_type</tt> if that is the case.  Otherwise it inherits from
 * <tt>boost::false_type</tt>.
 */
template<class T, class Enable = void>
struct is_complex_float : public boost::false_type {};

/** A specialization to handle recognized complex types */
template<class T>
struct is_complex_float<T, typename boost::enable_if<is_complex<T> >::type>
    : public boost::is_same<float,typename real<T>::type> {};

/**
 * An type trait similar to <tt>is_complex</tt> which additionally checks that
 * the underlying real type is <tt>double</tt>.  Inherits from
 * <tt>boost::true_type</tt> if that is the case.  Otherwise it inherits from
 * <tt>boost::false_type</tt>.
 */
template<class T, class Enable = void>
struct is_complex_double : public boost::false_type {};

/** A specialization to handle recognized complex types */
template<class T>
struct is_complex_double<T, typename boost::enable_if<is_complex<T> >::type>
    : public boost::is_same<double,typename real<T>::type> {};

/**
 * An type trait similar to <tt>is_complex</tt> which additionally checks that
 * the underlying real type is <tt>long double</tt>.  Inherits from
 * <tt>boost::true_type</tt> if that is the case.  Otherwise it inherits from
 * <tt>boost::false_type</tt>.
 */
template<class T, class Enable = void>
struct is_complex_long_double : public boost::false_type {};

/** A specialization to handle recognized complex types */
template<class T>
struct is_complex_long_double<T, typename boost::enable_if<is_complex<T> >::type>
    : public boost::is_same<long double,typename real<T>::type> {};

} // namespace traits

/**
 * Returns a <tt>std::complex</tt> instance with matching floating point type
 * whose real and imaginary components are equal to NaN.
 *
 * @return <tt>std::complex<FPT>(quiet_NaN,quiet_NaN)</tt>
 */
template<typename FPT>
typename boost::enable_if<
    boost::is_floating_point<FPT>,
    std::complex<FPT>
>::type NaN()
{
    BOOST_STATIC_ASSERT(std::numeric_limits<FPT>::has_quiet_NaN);
    const FPT quiet_NaN = std::numeric_limits<FPT>::quiet_NaN();
    return std::complex<FPT>(quiet_NaN,quiet_NaN);
}

/**
 * Returns a <tt>std::complex</tt> instance with matching floating point type
 * whose real and imaginary components are equal to NaN.
 *
 * @return <tt>std::complex<FPT>(quiet_NaN,quiet_NaN)</tt>
 */
template<class T>
std::complex<typename traits::real<T>::type>
NaN()
{
    typedef typename traits::real<T>::type real_type;
    BOOST_STATIC_ASSERT(std::numeric_limits<real_type>::has_quiet_NaN);
    const real_type quiet_NaN = std::numeric_limits<real_type>::quiet_NaN();
    return std::complex<real_type>(quiet_NaN,quiet_NaN);
}

/** Import <tt>std::real(std::complex<T>)</tt> */
using std::real;

/** Import <tt>std::imag(std::complex<T>)</tt> */
using std::imag;

/**
 * Obtain the real part of a complex number.
 *
 * @param z complex number stored as a two-element array.
 *
 * @return <tt>Re(z)</tt>
 */
template<typename FPT>
SUZERAIN_FORCEINLINE
FPT& real(FPT (&z)[2]) {
    return z[0];
}

/**
 * Obtain the real part of a complex number.
 *
 * @param z complex number stored as a two-element array.
 *
 * @return <tt>Re(z)</tt>
 */
template<typename FPT>
SUZERAIN_FORCEINLINE
const FPT& real(const FPT (&z)[2]) {
    return z[0];
}

/**
 * Obtain the imaginary part of a complex number.
 *
 * @param z complex number stored as a two-element array.
 *
 * @return <tt>Im(z)</tt>
 */
template<typename FPT>
SUZERAIN_FORCEINLINE
FPT& imag(FPT (&z)[2]) {
    return z[1];
}

/**
 * Obtain the imaginary part of a complex number.
 *
 * @param z complex number stored as a two-element array.
 *
 * @return <tt>Im(z)</tt>
 */
template<typename FPT>
SUZERAIN_FORCEINLINE
const FPT& imag(const FPT (&z)[2]) {
    return z[1];
}

/**
 * Overwrite \c dest with \c src.
 *
 * @param dest destination
 * @param src source
 */
template<class Complex, class Source>
SUZERAIN_FORCEINLINE
typename boost::enable_if<
    traits::is_complex<Source>, void
>::type assign_complex(Complex &dest, const Source &src)
{
    real(dest) = real(src);
    imag(dest) = imag(src);
}

/**
 * Overwrite \c dest with \c src.  The imaginary
 * part of \c dest is set to zero.
 *
 * @param dest destination
 * @param src source
 */
template<class Complex, class Source>
SUZERAIN_FORCEINLINE
typename boost::enable_if<
    boost::is_arithmetic<Source>, void
>::type assign_complex(Complex &dest, const Source &src)
{
    real(dest) = src;
    imag(dest) = 0;
}

/**
 * Overwrite \c dest with <tt>src_real + I*src_imag</tt> where \c I is
 * the imaginary unit.
 *
 * @param dest destination
 * @param src_real real part of the source
 * @param src_imag imag part of the source
 */
template<class Complex, typename FPT1, typename FPT2>
SUZERAIN_FORCEINLINE
void assign_complex(Complex &dest,
                    const FPT1 src_real,
                    const FPT2 src_imag)
{
    real(dest) = src_real;
    imag(dest) = src_imag;
}

/**
 * Overwrite \c dest_real with Re <tt>src</tt> and \c dest_imag with Re
 * <tt>src_imag</tt>.
 *
 * @param dest_real destination real part
 * @param dest_imag destination imag part
 * @param src source
 */
template<typename FPT, class Complex>
SUZERAIN_FORCEINLINE
void assign_components(FPT &dest_real,
                       FPT &dest_imag,
                       const Complex &src)
{
    dest_real = real(src);
    dest_imag = imag(src);
}

/**
 * Overwrite \c dest with <tt>alpha*src</tt>.
 *
 * @param dest destination
 * @param src source
 * @param alpha multiplicative real scaling factor
 */
template<class Complex1, class Complex2, typename FPT>
SUZERAIN_FORCEINLINE
void assign_complex_scaled(Complex1 &dest,
                           const Complex2 &src,
                           const FPT alpha)
{
    real(dest) = alpha*real(src);
    imag(dest) = alpha*imag(src);
}

/**
 * Overwrite \c dest with <tt>alpha*src*I^ipower</tt> where
 * \c I is the imaginary unit.
 *
 * @param dest destination
 * @param src source
 * @param alpha multiplicative real scaling factor
 * @param ipower exponent on the imaginary unit to include in the scaling
 */
template<class Complex1, class Complex2, typename FPT>
SUZERAIN_FORCEINLINE
void assign_complex_scaled_ipower(Complex1 &dest,
                                  const Complex2 &src,
                                  const FPT alpha,
                                  const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            real(dest) =  alpha*imag(src);
            imag(dest) = -alpha*real(src);
            break;
        case 2: // I^2 = -1 = I^-2
            real(dest) = -alpha*real(src);
            imag(dest) = -alpha*imag(src);
            break;
        case 1: // I^1 = I = I^-3
            real(dest) = -alpha*imag(src);
            imag(dest) =  alpha*real(src);
            break;
        case 0: // I^0 = 1
            real(dest) =  alpha*real(src);
            imag(dest) =  alpha*imag(src);
            break;
    }
}

} // namespace complex

namespace functional {

/**
 * A specialization of the assign functor to properly handle the case where \c
 * Target is a recognized complex type according to
 * ::suzerain::complex::traits::is_complex.  It uses
 * ::suzerain::complex::assign_complex to perform the assignment, and therefore
 * supports all types that \c assign_complex does.
 */
template<class Target, class Source>
struct assign<
    Target,
    Source,
    typename boost::enable_if<
        ::suzerain::complex::traits::is_complex<Target>
    >::type >
{
    /**
     * Create an instance which assigns \c s when applied.
     *
     * @param s source of assignment operation occurring via
     *          <tt>operator()</tt>.
     */
    assign(const Source &s) : s_(s) {};

    /**
     * Assign the value provided at construction to \c t.
     *
     * @param t to be assigned.
     */
    void operator()(Target& t) const {
        ::suzerain::complex::assign_complex(t, s_);
    }

private:
    const Source &s_; /**< Source for assignment operations */
};

} // namespace functional

} // namespace suzerain

#endif // __SUZERAIN_COMPLEX_HPP
