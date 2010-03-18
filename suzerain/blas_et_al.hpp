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
 * blas_et_al.hpp: C++ wrappers around blas_et_al.h's functionality
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_BLAS_ET_AL_HPP__
#define __SUZERAIN_BLAS_ET_AL_HPP__

#include <suzerain/blas_et_al.h>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>
#include <suzerain/allocator.hpp>
#include <suzerain/complex.hpp>

/** @file
 * C++ wrappers for external BLAS and LAPACK routines necessary for Suzerain.
 * Provides function overloading atop the pure C interfaces in blas_et_al.h.
 */
namespace suzerain {

/**
 * Provides function templates for BLAS level 1, 2, and 3 functionality.  Also
 * includes functionality from the <a
 * href="http://www.netlib.org/blas/blast-forum/">BLAS Technical Forum
 * Standard</a>.
 */
namespace blas {

/*! \name Memory allocation
 * @{
 */

/*! @copydoc suzerain_blas_malloc(size_t) */
inline void * malloc(size_t size)
{
    return suzerain_blas_malloc(size);
}

/*! @copydoc suzerain_blas_calloc(size_t,size_t) */
inline void * calloc(size_t nmemb, size_t size)
{
    return suzerain_blas_calloc(nmemb, size);
}

/*! @copydoc suzerain_blas_free(void*) */
inline void free(void *ptr)
{
    return suzerain_blas_free(ptr);
}

/*! @} */

/*! \name BLAS level 1 operations
 * @{
 */

/*! @copydoc suzerain_blas_sswap */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void swap(
        const Integer1 n,
        float *x,
        const Integer2 incx,
        float *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_sswap(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx),
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sswap */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void swap(
        const Integer1 n,
        double *x,
        const Integer2 incx,
        double *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_dswap(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx),
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sswap */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
> >::type swap(
        const Integer1 n,
        Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_cswap(boost::numeric_cast<int>(n),
                               reinterpret_cast<float (*)[2]>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<float (*)[2]>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sswap */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
> >::type swap(
        const Integer1 n,
        Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_zswap(boost::numeric_cast<int>(n),
                               reinterpret_cast<double (*)[2]>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<double (*)[2]>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_scopy */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void copy(
        const Integer1 n,
        const float *x,
        const Integer2 incx,
        float *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_scopy(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx),
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_scopy */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void copy(
        const Integer1 n,
        const double *x,
        const Integer2 incx,
        double *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_dcopy(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx),
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sdot */
template< typename Integer1, typename Integer2, typename Integer3 >
inline float dot(
        const Integer1 n,
        const float *x,
        const Integer2 incx,
        const float *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_sdot(boost::numeric_cast<int>(n),
                              x,
                              boost::numeric_cast<int>(incx),
                              y,
                              boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sdot */
template< typename Integer1, typename Integer2, typename Integer3 >
inline double dot(
        const Integer1 n,
        const double *x,
        const Integer2 incx,
        const double *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_ddot(boost::numeric_cast<int>(n),
                              x,
                              boost::numeric_cast<int>(incx),
                              y,
                              boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sasum */
template< typename Integer1, typename Integer2 >
inline float asum(
        const Integer1 n,
        const float *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_sasum(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_sasum */
template< typename Integer1, typename Integer2 >
inline double asum(
        const Integer1 n,
        const double *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_dasum(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_saxpy */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void axpy(
        const Integer1 n,
        const float alpha,
        const float *x,
        const Integer2 incx,
        float *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_saxpy(boost::numeric_cast<int>(n),
                               alpha,
                               x,
                               boost::numeric_cast<int>(incx),
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_saxpy */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void axpy(
        const Integer1 n,
        const double alpha,
        const double *x,
        const Integer2 incx,
        double *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_daxpy(boost::numeric_cast<int>(n),
                               alpha,
                               x,
                               boost::numeric_cast<int>(incx),
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_saxpy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2, typename Complex3 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>,
    suzerain::complex::traits::is_complex_float<Complex3>
> >::type axpy(
        const Integer1 n,
        const Complex1 alpha,
        const Complex2 *x,
        const Integer2 incx,
        Complex3 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_caxpy(boost::numeric_cast<int>(n),
                               reinterpret_cast<const float (*)[2]>(alpha),
                               reinterpret_cast<const float (*)[2]>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<float (*)[2]>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_saxpy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
> >::type axpy(
        const Integer1 n,
        const float alpha,
        const Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    const float alpha_complex[2] = { alpha, 0 };
    return suzerain_blas_caxpy(boost::numeric_cast<int>(n),
                               alpha_complex,
                               reinterpret_cast<const float (*)[2]>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<float (*)[2]>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_saxpy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2, typename Complex3 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>,
    suzerain::complex::traits::is_complex_double<Complex3>
> >::type axpy(
        const Integer1 n,
        const Complex1 alpha,
        const Complex2 *x,
        const Integer2 incx,
        Complex3 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_zaxpy(boost::numeric_cast<int>(n),
                               reinterpret_cast<const double (*)[2]>(alpha),
                               reinterpret_cast<const double (*)[2]>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<double (*)[2]>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_saxpy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
> >::type axpy(
        const Integer1 n,
        const double alpha,
        const Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    const double alpha_complex[2] = { alpha, 0 };
    return suzerain_blas_zaxpy(boost::numeric_cast<int>(n),
                               alpha_complex,
                               reinterpret_cast<const double (*)[2]>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<double (*)[2]>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_saxpby */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void axpby(
        const Integer1 n,
        const float alpha,
        const float *x,
        const Integer2 incx,
        const float beta,
        float *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_saxpby(boost::numeric_cast<int>(n),
                                alpha,
                                x,
                                boost::numeric_cast<int>(incx),
                                beta,
                                y,
                                boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_saxpby */
template< typename Integer1, typename Integer2, typename Integer3 >
inline void axpby(
        const Integer1 n,
        const double alpha,
        const double *x,
        const Integer2 incx,
        const double beta,
        double *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_daxpby(boost::numeric_cast<int>(n),
                                alpha,
                                x,
                                boost::numeric_cast<int>(incx),
                                beta,
                                y,
                                boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_swaxpby */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4 >
inline void waxpby(
        const Integer1 n,
        const float alpha,
        const float *x,
        const Integer2 incx,
        const float beta,
        const float *y,
        const Integer3 incy,
        float *w,
        const Integer4 incw)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    return suzerain_blas_swaxpby(boost::numeric_cast<int>(n),
                                 alpha,
                                 x,
                                 boost::numeric_cast<int>(incx),
                                 beta,
                                 y,
                                 boost::numeric_cast<int>(incy),
                                 w,
                                 boost::numeric_cast<int>(incw));
}

/*! @copydoc suzerain_blas_swaxpby */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4 >
inline void waxpby(
        const int n,
        const double alpha,
        const double *x,
        const int incx,
        const double beta,
        const double *y,
        const int incy,
        double *w,
        const int incw)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    return suzerain_blas_dwaxpby(boost::numeric_cast<int>(n),
                                 alpha,
                                 x,
                                 boost::numeric_cast<int>(incx),
                                 beta,
                                 y,
                                 boost::numeric_cast<int>(incy),
                                 w,
                                 boost::numeric_cast<int>(incw));
}

/*! @copydoc suzerain_blas_sscal */
template< typename Integer1, typename Integer2 >
inline void scal(
        const Integer1 n,
        const float alpha,
        float *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_sscal(boost::numeric_cast<int>(n),
                               alpha,
                               x,
                               boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_sscal */
template< typename Integer1, typename Integer2 >
inline void scal(
        const Integer1 n,
        const double alpha,
        double *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_dscal(boost::numeric_cast<int>(n),
                               alpha,
                               x,
                               boost::numeric_cast<int>(incx));
}

/*! @} */

/*! \name BLAS level 2 operations
 * @{
 */

/*! @copydoc suzerain_blas_sgbmv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6,
          typename Integer7 >
inline void gbmv(
        const char trans,
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const float alpha,
        const float *a,
        const Integer5 lda,
        const float *x,
        const Integer6 incx,
        const float beta,
        float *y,
        const Integer7 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer7>::value);
    return suzerain_blas_sgbmv(trans,
                               boost::numeric_cast<int>(m),
                               boost::numeric_cast<int>(n),
                               boost::numeric_cast<int>(kl),
                               boost::numeric_cast<int>(ku),
                               alpha,
                               a,
                               boost::numeric_cast<int>(lda),
                               x,
                               boost::numeric_cast<int>(incx),
                               beta,
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sgbmv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6,
          typename Integer7 >
inline void gbmv(
        const char trans,
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const double alpha,
        const double *a,
        const Integer5 lda,
        const double *x,
        const Integer6 incx,
        const double beta,
        double *y,
        const Integer7 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer7>::value);
    return suzerain_blas_dgbmv(trans,
                               boost::numeric_cast<int>(m),
                               boost::numeric_cast<int>(n),
                               boost::numeric_cast<int>(kl),
                               boost::numeric_cast<int>(ku),
                               alpha,
                               a,
                               boost::numeric_cast<int>(lda),
                               x,
                               boost::numeric_cast<int>(incx),
                               beta,
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_sgb_acc */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6 >
inline void gb_acc(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const float alpha,
        const float *a,
        const Integer5 lda,
        const float beta,
        float *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_blas_sgb_acc(boost::numeric_cast<int>(m),
                                 boost::numeric_cast<int>(n),
                                 boost::numeric_cast<int>(kl),
                                 boost::numeric_cast<int>(ku),
                                 alpha,
                                 a,
                                 boost::numeric_cast<int>(lda),
                                 beta,
                                 b,
                                 boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_blas_sgb_acc */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6 >
inline void gb_acc(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const double alpha,
        const double *a,
        const Integer5 lda,
        const double beta,
        double *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_blas_dgb_acc(boost::numeric_cast<int>(m),
                                 boost::numeric_cast<int>(n),
                                 boost::numeric_cast<int>(kl),
                                 boost::numeric_cast<int>(ku),
                                 alpha,
                                 a,
                                 boost::numeric_cast<int>(lda),
                                 beta,
                                 b,
                                 boost::numeric_cast<int>(ldb));
}

/*! @} */

/*! \name BLAS level 3 operations
 * @{
 */

/*! @} */

/**
 * An aligned, heap-based allocation policy using <tt>suzerain::blas::malloc</tt>
 * and <tt>free</tt>.
 *
 * @see Lai Shiaw San Kent's article <a
 * href="http://www.codeproject.com/KB/cpp/allocator.aspx">C++ Standard
 * Allocator, An Introduction and Implementation</a> for more information.
 */
template<typename T>
class allocator_blas_policy {
public:
    // Typedefs following std::allocator contract
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    /** Rebind to obtain allocator_blas_policy<U> */
    template<typename U>
    struct rebind {
        typedef allocator_blas_policy<U> other;
    };

    // Make explicit all constructors, including copy constructors
    explicit allocator_blas_policy() {}
    ~allocator_blas_policy() {}
    explicit allocator_blas_policy(allocator_blas_policy const&) {}
    template <typename U>
    explicit allocator_blas_policy(allocator_blas_policy<U> const&) {}

    /** Memory allocation uses ::suzerain::blas::malloc and std::free */
    //@{
    /**
     * Allocate space for \c cnt instances of type \c value_type.
     * @throws std::bad_alloc on allocation failure
     */
    pointer allocate(size_type cnt,
                     typename std::allocator<void>::const_pointer = 0)
    {
        if (SUZERAIN_UNLIKELY(cnt == 0)) return NULL;

        pointer p = reinterpret_cast<pointer>(
                ::suzerain::blas::malloc(cnt * sizeof(T)));

        if (SUZERAIN_UNLIKELY(!p)) throw std::bad_alloc();

        return p;
    }
    void deallocate(pointer p, size_type) { ::suzerain::blas::free(p); }
    //@}

    /** max_size method following std::allocator contract */
    size_type max_size() const
    {
        return std::numeric_limits<size_type>::max() / sizeof(T);
    }
};

/** Provide information about compatibility between different allocators. */
//@{
/**
 * Memory can be deallocated from other allocator_blas_policy
 * instantiations.
 **/
template<typename T, typename T2>
inline bool operator==(allocator_blas_policy<T> const&,
                       allocator_blas_policy<T2> const&)
{
    return true;
}

/** Memory cannot be deallocated by other allocator types */
template<typename T, typename OtherAllocator>
inline bool operator==(allocator_blas_policy<T> const&,
                       OtherAllocator const&)
{
    return false;
}
//@}

/**
 * Simulate a template typedef for a BLAS-aligned STL allocator.
 * Combines suzerain::allocator with allocator_blas_policy.
 */
template<class T>
struct allocator
{
    typedef typename ::suzerain::allocator<T, allocator_blas_policy<T> > type;
};

} // namespace blas

/**
 * Provides function templates for LAPACK functionality.
 */
namespace lapack {

/*! @copydoc suzerain_lapack_sgbtrf */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5 >
inline int gbtrf(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        float *ab,
        const Integer5 ldab,
        int *ipiv)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_lapack_sgbtrf(boost::numeric_cast<int>(m),
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  ab,
                                  boost::numeric_cast<int>(ldab),
                                  ipiv);
}

/*! @copydoc suzerain_lapack_sgbtrf */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5 >
inline int gbtrf(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        double *ab,
        const Integer5 ldab,
        int *ipiv)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_lapack_dgbtrf(boost::numeric_cast<int>(m),
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  ab,
                                  boost::numeric_cast<int>(ldab),
                                  ipiv);
}

/*! @copydoc suzerain_lapack_sgbtrs */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6 >
inline int gbtrs(
        const char trans,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        const float *ab,
        const Integer5 ldab,
        const int *ipiv,
        float *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_sgbtrs(trans,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  boost::numeric_cast<int>(nrhs),
                                  ab,
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  b,
                                  boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_lapack_sgbtrs */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6 >
inline int gbtrs(
        const char trans,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        const double *ab,
        const Integer5 ldab,
        const int *ipiv,
        double *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_dgbtrs(trans,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  boost::numeric_cast<int>(nrhs),
                                  ab,
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  b,
                                  boost::numeric_cast<int>(ldb));
}

} // namespace lapack

} // namespace suzerain

#endif /* __SUZERAIN_BLAS_ET_AL_HPP__ */
