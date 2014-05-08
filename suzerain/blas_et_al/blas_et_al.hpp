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

#ifndef SUZERAIN_BLAS_ET_AL_HPP
#define SUZERAIN_BLAS_ET_AL_HPP

/** @file
 * C++ wrappers around functionality within blas_et_al.h
 */

#include <boost/numeric/conversion/cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>

#include <suzerain/allocator.hpp>
#include <suzerain/blas_et_al/blas_et_al.h>
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

/*! @copydoc SUZERAIN_BLAS_ALIGNMENT */
static const size_t alignment = SUZERAIN_BLAS_ALIGNMENT;

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

/*! calloc \c nmemb elements of size <tt>sizeof(T)</tt>. */
template<typename T>
inline T * calloc_as(size_t nmemb)
{
    return reinterpret_cast<T *>(calloc(nmemb, sizeof(T)));
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

/*! @copydoc suzerain_blas_dswap */
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

/*! @copydoc suzerain_blas_cswap */
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
                               reinterpret_cast<std::complex<float> *>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<std::complex<float> *>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_zswap */
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
                               reinterpret_cast<complex_double *>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<complex_double *>(y),
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

/*! @copydoc suzerain_blas_dcopy */
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

/*! @copydoc suzerain_blas_ccopy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
> >::type copy(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_ccopy(boost::numeric_cast<int>(n),
                               reinterpret_cast<const std::complex<float> *>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<std::complex<float> *>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_zcopy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
> >::type copy(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    return suzerain_blas_zcopy(boost::numeric_cast<int>(n),
                               reinterpret_cast<const complex_double *>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<complex_double *>(y),
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

/*! @copydoc suzerain_blas_ddot */
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

/*!
 * \brief Compute \f$ x \cdot{} y \f$ using BLAS's dotc.
 *
 * \param n Number of elements in \c x and \c y.
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param y Second source vector.
 * \param incy Second source vector stride.
 *
 * \see A BLAS reference for more details.
 */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
>, std::complex<float> >::type dot(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx,
        const Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    std::complex<float> result;
    suzerain_blas_cdotc(boost::numeric_cast<int>(n),
                        reinterpret_cast<const std::complex<float> *>(x),
                        boost::numeric_cast<int>(incx),
                        reinterpret_cast<const std::complex<float> *>(y),
                        boost::numeric_cast<int>(incy),
                        &result);
    return result;
}

/*!
 * \brief Compute \f$ x \cdot{} y \f$ using BLAS's dotc.
 *
 * \param n Number of elements in \c x and \c y.
 * \param x First source vector.
 * \param incx First source vector stride.
 * \param y Second source vector.
 * \param incy Second source vector stride.
 *
 * \see A BLAS reference for more details.
 */
template< typename Integer1, typename Integer2, typename Integer3,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
>, std::complex<double> >::type dot(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx,
        const Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    std::complex<double> result;
    suzerain_blas_zdotc(boost::numeric_cast<int>(n),
                        reinterpret_cast<const complex_double *>(x),
                        boost::numeric_cast<int>(incx),
                        reinterpret_cast<const complex_double *>(y),
                        boost::numeric_cast<int>(incy),
                        &result);
    return result;
}

/*! @copydoc suzerain_blas_snrm2 */
template< typename Integer1, typename Integer2 >
inline float nrm2(
        const Integer1 n,
        const float *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_snrm2(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_dnrm2 */
template< typename Integer1, typename Integer2 >
inline double nrm2(
        const Integer1 n,
        const double *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_dnrm2(boost::numeric_cast<int>(n),
                               x,
                               boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_scnrm2 */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_float<Complex1>, float
>::type nrm2(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_scnrm2(boost::numeric_cast<int>(n),
                                reinterpret_cast<const std::complex<float> *>(x),
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_dznrm2 */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex1>, double
>::type nrm2(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_dznrm2(boost::numeric_cast<int>(n),
                                reinterpret_cast<const complex_double *>(x),
                                boost::numeric_cast<int>(incx));
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

/*! @copydoc suzerain_blas_dasum */
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

/*! @copydoc suzerain_blas_scasum */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_float<Complex1>, float
>::type asum(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_scasum(boost::numeric_cast<int>(n),
                                reinterpret_cast<const std::complex<float> *>(x),
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_dzasum */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex1>, double
>::type asum(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_dzasum(boost::numeric_cast<int>(n),
                                reinterpret_cast<const complex_double *>(x),
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_isamax */
template< typename Integer1, typename Integer2 >
inline int iamax(
        const Integer1 n,
        const float *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_isamax(boost::numeric_cast<int>(n),
                                x,
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_idamax */
template< typename Integer1, typename Integer2 >
inline int iamax(
        const Integer1 n,
        const double *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_idamax(boost::numeric_cast<int>(n),
                                x,
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_icamax */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_float<Complex1>, int
>::type iamax(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_icamax(boost::numeric_cast<int>(n),
                                reinterpret_cast<const std::complex<float> *>(x),
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_izamax */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex1>, int
>::type iamax(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_izamax(boost::numeric_cast<int>(n),
                                reinterpret_cast<const complex_double *>(x),
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_isamin */
template< typename Integer1, typename Integer2 >
inline int iamin(
        const Integer1 n,
        const float *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_isamin(boost::numeric_cast<int>(n),
                                x,
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_idamin */
template< typename Integer1, typename Integer2 >
inline int iamin(
        const Integer1 n,
        const double *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_idamin(boost::numeric_cast<int>(n),
                                x,
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_icamin */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_float<Complex1>, int
>::type iamin(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_icamin(boost::numeric_cast<int>(n),
                                reinterpret_cast<const std::complex<float> *>(x),
                                boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_izamin */
template< typename Integer1, typename Integer2, typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex1>, int
>::type iamin(
        const Integer1 n,
        const Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    return suzerain_blas_izamin(boost::numeric_cast<int>(n),
                                reinterpret_cast<const complex_double *>(x),
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

/*! @copydoc suzerain_blas_daxpy */
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

/*! @copydoc suzerain_blas_caxpy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename AlphaType,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
> >::type axpy(
        const Integer1 n,
        const AlphaType &alpha,
        const Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    std::complex<float> alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    return suzerain_blas_caxpy(boost::numeric_cast<int>(n),
                               alpha_complex,
                               reinterpret_cast<const std::complex<float> *>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<std::complex<float> *>(y),
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_zaxpy */
template< typename Integer1, typename Integer2, typename Integer3,
          typename AlphaType,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
> >::type axpy(
        const Integer1 n,
        const AlphaType &alpha,
        const Complex1 *x,
        const Integer2 incx,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    complex_double alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    return suzerain_blas_zaxpy(boost::numeric_cast<int>(n),
                               alpha_complex,
                               reinterpret_cast<const complex_double *>(x),
                               boost::numeric_cast<int>(incx),
                               reinterpret_cast<complex_double *>(y),
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

/*! @copydoc suzerain_blas_daxpby */
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

/*! @copydoc suzerain_blas_caxpby */
template< typename Integer1, typename Integer2, typename Integer3,
          typename AlphaType, typename BetaType,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
> >::type axpby(
        const Integer1 n,
        const AlphaType &alpha,
        const Complex1 *x,
        const Integer2 incx,
        const BetaType &beta,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    std::complex<float> alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    std::complex<float> beta_complex;
    suzerain::complex::assign_complex(beta_complex, beta);
    return suzerain_blas_caxpby(boost::numeric_cast<int>(n),
                                alpha,
                                reinterpret_cast<const std::complex<float> *>(x),
                                boost::numeric_cast<int>(incx),
                                beta,
                                reinterpret_cast<std::complex<float> *>(y),
                                boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_zaxpby */
template< typename Integer1, typename Integer2, typename Integer3,
          typename AlphaType, typename BetaType,
          typename Complex1, typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
> >::type axpby(
        const Integer1 n,
        const AlphaType &alpha,
        const Complex1 *x,
        const Integer2 incx,
        const BetaType &beta,
        Complex2 *y,
        const Integer3 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    complex_double alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    complex_double beta_complex;
    suzerain::complex::assign_complex(beta_complex, beta);
    return suzerain_blas_zaxpby(boost::numeric_cast<int>(n),
                                alpha_complex,
                                reinterpret_cast<const complex_double *>(x),
                                boost::numeric_cast<int>(incx),
                                beta_complex,
                                reinterpret_cast<complex_double *>(y),
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

/*! @copydoc suzerain_blas_dwaxpby */
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

/*! @copydoc suzerain_blas_dscal */
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

/*! @copydoc suzerain_blas_cscal */
template< typename Integer1, typename Integer2,
          typename AlphaType, typename Complex1  >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_float<Complex1>
>::type scal(
        const Integer1 n,
        const AlphaType &alpha,
        Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    std::complex<float> alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    return suzerain_blas_cscal(boost::numeric_cast<int>(n),
                               alpha_complex,
                               reinterpret_cast<std::complex<float> *>(x),
                               boost::numeric_cast<int>(incx));
}

/*! @copydoc suzerain_blas_zscal */
template< typename Integer1, typename Integer2,
          typename AlphaType, typename Complex1  >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex1>
>::type scal(
        const Integer1 n,
        const AlphaType &alpha,
        Complex1 *x,
        const Integer2 incx)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    complex_double alpha_complex;
    suzerain::complex::assign_complex(alpha_complex, alpha);
    return suzerain_blas_zscal(boost::numeric_cast<int>(n),
                               alpha_complex,
                               reinterpret_cast<complex_double *>(x),
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

/*! @copydoc suzerain_blas_dgbmv */
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

/*! @copydoc suzerain_blas_ssbmv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5 >
inline void sbmv(
        const char uplo,
        const Integer1 n,
        const Integer2 k,
        const float alpha,
        const float *a,
        const Integer3 lda,
        const float *x,
        const Integer4 incx,
        const float beta,
        float *y,
        const Integer5 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_blas_ssbmv(uplo,
                               boost::numeric_cast<int>(n),
                               boost::numeric_cast<int>(k),
                               alpha,
                               a,
                               boost::numeric_cast<int>(lda),
                               x,
                               boost::numeric_cast<int>(incx),
                               beta,
                               y,
                               boost::numeric_cast<int>(incy));
}

/*! @copydoc suzerain_blas_dsbmv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5 >
inline void sbmv(
        const char uplo,
        const Integer1 n,
        const Integer2 k,
        const double alpha,
        const double *a,
        const Integer3 lda,
        const double *x,
        const Integer4 incx,
        const double beta,
        double *y,
        const Integer5 incy)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_blas_dsbmv(uplo,
                               boost::numeric_cast<int>(n),
                               boost::numeric_cast<int>(k),
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

/*! @copydoc suzerain_blas_dgb_acc */
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

    /** Memory allocation uses suzerain::blas::malloc and std::free */
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
                suzerain::blas::malloc(cnt * sizeof(T)));

        if (SUZERAIN_UNLIKELY(!p)) throw std::bad_alloc();

        return p;
    }
    void deallocate(pointer p, size_type) { suzerain::blas::free(p); }
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
    typedef typename suzerain::allocator<T, allocator_blas_policy<T> > type;
};


/**
 * Provides functors for BLAS operations.  These are sometimes useful when
 * mapping operations with common or nearly common signatures across data
 * fields, especially for higher dimensional cases.
 */
namespace functor {

/*! \name BLAS level 1 functors
 * @{
 */

/** A functor performing suzerain::blas::swap when invoked. */
class swap {
public:
    template< typename TN,
              typename TX, typename TINCX,
              typename TY, typename TINCY >
    void operator()(TN n, TX x, TINCX incx, TY y, TINCY incy) const
    {
        return suzerain::blas::swap(n, x, incx, y, incy);
    }
};

/** A functor performing suzerain::blas::copy when invoked. */
class copy {
public:
    template< typename TN,
              typename TX, typename TINCX,
              typename TY, typename TINCY >
    void operator()(TN n, TX x, TINCX incx, TY y, TINCY incy) const
    {
        return suzerain::blas::copy(n, x, incx, y, incy);
    }
};

/**
 * A functor performing suzerain::blas::axpy when invoked.
 * The \c alpha parameter is set at construction time.
 **/
template< typename Element >
class axpy {
public:
    axpy(const Element& alpha) : alpha_(alpha) {}

    template< typename TN,
              typename TX, typename TINCX,
              typename TY, typename TINCY >
    void operator()(TN n, TX x, TINCX incx, TY y, TINCY incy) const
    {
        return suzerain::blas::axpy(n, alpha_, x, incx, y, incy);
    }

private:
    const Element& alpha_;
};

/*! @} */

} // namespace functor

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

/*! @copydoc suzerain_lapack_dgbtrf */
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

/*! @copydoc suzerain_lapack_cgbtrf */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_float<Complex1>
, int>::type gbtrf(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        Complex1 *ab,
        const Integer5 ldab,
        int *ipiv)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_lapack_cgbtrf(boost::numeric_cast<int>(m),
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  reinterpret_cast<std::complex<float> *>(ab),
                                  boost::numeric_cast<int>(ldab),
                                  ipiv);
}

/*! @copydoc suzerain_lapack_zgbtrf */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex1>
, int>::type gbtrf(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        Complex1 *ab,
        const Integer5 ldab,
        int *ipiv)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_lapack_zgbtrf(boost::numeric_cast<int>(m),
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  reinterpret_cast<complex_double *>(ab),
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

/*! @copydoc suzerain_lapack_dgbtrs */
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

/*! @copydoc suzerain_lapack_cgbtrs */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6,
          typename Complex1,
          typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
>, int >::type gbtrs(
        const char trans,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        const Complex1 *ab,
        const Integer5 ldab,
        const int *ipiv,
        Complex2 *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_cgbtrs(trans,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  boost::numeric_cast<int>(nrhs),
                                  reinterpret_cast<const std::complex<float> *>(ab),
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  reinterpret_cast<std::complex<float> *>(b),
                                  boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_lapack_zgbtrs */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6,
          typename Complex1,
          typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
>, int >::type gbtrs(
        const char trans,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        const Complex1 *ab,
        const Integer5 ldab,
        const int *ipiv,
        Complex2 *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_zgbtrs(trans,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  boost::numeric_cast<int>(nrhs),
                                  reinterpret_cast<const complex_double *>(ab),
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  reinterpret_cast<complex_double *>(b),
                                  boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_lapack_sgbsv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6 >
inline int gbsv(
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        float *ab,
        const Integer5 ldab,
        int *ipiv,
        float *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_sgbsv(boost::numeric_cast<int>(n),
                                 boost::numeric_cast<int>(kl),
                                 boost::numeric_cast<int>(ku),
                                 boost::numeric_cast<int>(nrhs),
                                 ab,
                                 boost::numeric_cast<int>(ldab),
                                 ipiv,
                                 b,
                                 boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_lapack_dgbsv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6 >
inline int gbsv(
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        double *ab,
        const Integer5 ldab,
        int *ipiv,
        double *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_dgbsv(boost::numeric_cast<int>(n),
                                 boost::numeric_cast<int>(kl),
                                 boost::numeric_cast<int>(ku),
                                 boost::numeric_cast<int>(nrhs),
                                 ab,
                                 boost::numeric_cast<int>(ldab),
                                 ipiv,
                                 b,
                                 boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_lapack_cgbsv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6,
          typename Complex1,
          typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
>, int >::type gbsv(
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        Complex1 *ab,
        const Integer5 ldab,
        int *ipiv,
        Complex2 *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_cgbsv(boost::numeric_cast<int>(n),
                                 boost::numeric_cast<int>(kl),
                                 boost::numeric_cast<int>(ku),
                                 boost::numeric_cast<int>(nrhs),
                                 reinterpret_cast<std::complex<float> *>(ab),
                                 boost::numeric_cast<int>(ldab),
                                 ipiv,
                                 reinterpret_cast<std::complex<float> *>(b),
                                 boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_lapack_zgbsv */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Integer6,
          typename Complex1,
          typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
>, int >::type gbsv(
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Integer4 nrhs,
        Complex1 *ab,
        const Integer5 ldab,
        int *ipiv,
        Complex2 *b,
        const Integer6 ldb)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer6>::value);
    return suzerain_lapack_zgbsv(boost::numeric_cast<int>(n),
                                 boost::numeric_cast<int>(kl),
                                 boost::numeric_cast<int>(ku),
                                 boost::numeric_cast<int>(nrhs),
                                 reinterpret_cast<complex_double *>(ab),
                                 boost::numeric_cast<int>(ldab),
                                 ipiv,
                                 reinterpret_cast<complex_double *>(b),
                                 boost::numeric_cast<int>(ldb));
}

/*! @copydoc suzerain_lapack_sgbcon */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4 >
inline int gbcon(
        const char norm,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const float *ab,
        const Integer4 ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        float *work,
        int *iwork)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    return suzerain_lapack_sgbcon(norm,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  ab,
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  anorm,
                                  rcond,
                                  work,
                                  iwork);
}

/*! @copydoc suzerain_lapack_dgbcon */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4 >
inline int gbcon(
        const char norm,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const double *ab,
        const Integer4 ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        double *work,
        int *iwork)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    return suzerain_lapack_dgbcon(norm,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  ab,
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  anorm,
                                  rcond,
                                  work,
                                  iwork);
}

/*! @copydoc suzerain_lapack_cgbcon */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Complex1,
          typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_float<Complex1>,
    suzerain::complex::traits::is_complex_float<Complex2>
> >::type gbcon(
        const char norm,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Complex1 *ab,
        const Integer4 ldab,
        const int *ipiv,
        const float anorm,
        float *rcond,
        Complex2 *work,
        float *rwork)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    return suzerain_lapack_cgbcon(norm,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  reinterpret_cast<const std::complex<float> *>(ab),
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  anorm,
                                  rcond,
                                  reinterpret_cast<std::complex<float> *>(work),
                                  rwork);
}

/*! @copydoc suzerain_lapack_zgbcon */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Complex1,
          typename Complex2 >
inline typename boost::enable_if<boost::mpl::and_<
    suzerain::complex::traits::is_complex_double<Complex1>,
    suzerain::complex::traits::is_complex_double<Complex2>
> >::type gbcon(
        const char norm,
        const Integer1 n,
        const Integer2 kl,
        const Integer3 ku,
        const Complex1 *ab,
        const Integer4 ldab,
        const int *ipiv,
        const double anorm,
        double *rcond,
        Complex2 *work,
        double *rwork)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    return suzerain_lapack_zgbcon(norm,
                                  boost::numeric_cast<int>(n),
                                  boost::numeric_cast<int>(kl),
                                  boost::numeric_cast<int>(ku),
                                  reinterpret_cast<const complex_double *>(ab),
                                  boost::numeric_cast<int>(ldab),
                                  ipiv,
                                  anorm,
                                  rcond,
                                  reinterpret_cast<complex_double *>(work),
                                  rwork);
}

} // namespace lapack

/**
 * Provides function templates for BLAS-like extensions
 */
namespace blasext {

/*! @copydoc suzerain_blasext_sgbnorm1 */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5 >
inline int gbnorm1(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const float *a,
        const Integer5 lda,
        float *norm1)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_blasext_sgbnorm1(boost::numeric_cast<int>(m),
                                     boost::numeric_cast<int>(n),
                                     boost::numeric_cast<int>(kl),
                                     boost::numeric_cast<int>(ku),
                                     a,
                                     boost::numeric_cast<int>(lda),
                                     norm1);
}

/*! @copydoc suzerain_blasext_dgbnorm1 */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5 >
inline int gbnorm1(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const double *a,
        const Integer5 lda,
        double *norm1)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_blasext_dgbnorm1(boost::numeric_cast<int>(m),
                                     boost::numeric_cast<int>(n),
                                     boost::numeric_cast<int>(kl),
                                     boost::numeric_cast<int>(ku),
                                     a,
                                     boost::numeric_cast<int>(lda),
                                     norm1);
}

/*! @copydoc suzerain_blasext_cgbnorm1 */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_float<Complex1>, int
>::type gbnorm1(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const Complex1 *a,
        const Integer5 lda,
        float *norm1)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_blasext_cgbnorm1(boost::numeric_cast<int>(m),
                                     boost::numeric_cast<int>(n),
                                     boost::numeric_cast<int>(kl),
                                     boost::numeric_cast<int>(ku),
                                     reinterpret_cast<const std::complex<float> *>(a),
                                     boost::numeric_cast<int>(lda),
                                     norm1);
}

/*! @copydoc suzerain_blasext_zgbnorm1 */
template< typename Integer1,
          typename Integer2,
          typename Integer3,
          typename Integer4,
          typename Integer5,
          typename Complex1 >
inline typename boost::enable_if<
    suzerain::complex::traits::is_complex_double<Complex1>, int
>::type gbnorm1(
        const Integer1 m,
        const Integer2 n,
        const Integer3 kl,
        const Integer4 ku,
        const Complex1 *a,
        const Integer5 lda,
        double *norm1)
{
    BOOST_STATIC_ASSERT(boost::is_integral<Integer1>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer2>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer3>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer4>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<Integer5>::value);
    return suzerain_blasext_zgbnorm1(boost::numeric_cast<int>(m),
                                     boost::numeric_cast<int>(n),
                                     boost::numeric_cast<int>(kl),
                                     boost::numeric_cast<int>(ku),
                                     reinterpret_cast<const complex_double *>(a),
                                     boost::numeric_cast<int>(lda),
                                     norm1);
}

} // namespace blasext

} // namespace suzerain

#endif /* SUZERAIN_BLAS_ET_AL_HPP */
