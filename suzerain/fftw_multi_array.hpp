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
 * fftw_multi_array.hpp: perform FFTs atop boost::multi_array using FFTW
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef PECOS_SUZERAIN_FFTW_MULTI_ARRAY_HPP
#define PECOS_SUZERAIN_FFTW_MULTI_ARRAY_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <boost/array.hpp>
#include <boost/integer_traits.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/decay.hpp>
#include <boost/type_traits/integral_constant.hpp>
#include <boost/type_traits/is_array.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/typeof/typeof.hpp>
#include <fftw3.h>

namespace pecos { namespace suzerain {

/**
  * TODO Namespace documentation
  */
namespace fftw_multi_array {

/**
 * Implementation details underneath pecos::suzerain::fftw_multi_array
 * @internal
 */
namespace detail {


template<std::size_t NumDims,
         typename IndexType,
         typename MaxIndexType,
         class InputIterator>
bool increment(IndexType &index,
               const MaxIndexType &max_index,
               InputIterator index_order)
{
    using boost::integer_traits;

    typedef BOOST_TYPEOF_TPL(index[0])              index_element_type;
    typedef BOOST_TYPEOF_TPL(max_index[0])          max_index_element_type;
    typedef BOOST_TYPEOF_TPL(index[0]/max_index[0]) element_division_type;

    // Assert compile time algorithm preconditions valid when in debug mode
    BOOST_STATIC_ASSERT(NumDims > 0);
    BOOST_STATIC_ASSERT(!boost::is_const<index_element_type>::value);
    BOOST_STATIC_ASSERT(boost::is_integral<element_division_type>::value);

    // Because we have 4*NumDims integral operations in the fully unrolled loop
    // and NumDims is often small, do not break out of loop when overflow == 0.
    // The overflow == 0 condition causes an effective NOP after occurring.
    bool overflow = 1;
    for (std::size_t i = 0; i < NumDims; ++i) {
        const typename InputIterator::value_type n = *index_order++;

        // Assert runtime algorithm preconditions valid when in debug mode
        assert(1 <= max_index[n]);
        assert(boost::is_unsigned<index_element_type>::value || 0 <= index[n]);
        assert(static_cast<max_index_element_type>(index[n]) < max_index[n]);
        assert(index[n]< integer_traits<index_element_type>::const_max - 1);

        index[n] += overflow;                 // Handle incoming overflow
        overflow  = index[n]/max_index[n];    // Check outgoing overflow
        index[n] *= !overflow;                // Set to zero on outgoing
    }
    return !overflow;
}

template<std::size_t NumDims,
         typename IndexType,
         typename MaxIndexType>
bool increment_impl(IndexType &index,
                    const MaxIndexType &max_index,
                    const boost::false_type &)
{
    typedef typename IndexType::value_type value_type;
    typedef boost::counting_iterator<value_type> counting_iterator;

    return increment<NumDims, IndexType, MaxIndexType, counting_iterator>(
            index, max_index, counting_iterator(0));
}

template<std::size_t NumDims,
         typename IndexType,
         typename MaxIndexType>
bool increment_impl(IndexType &index,
                    const MaxIndexType &max_index,
                    const boost::true_type &)
{
    typedef typename std::iterator_traits<
            typename boost::decay<IndexType>::type
        >::value_type value_type;
    typedef boost::counting_iterator<value_type> counting_iterator;

    return increment<NumDims, IndexType, MaxIndexType, counting_iterator>(
            index, max_index, counting_iterator(0));
}

template<std::size_t NumDims,
         typename IndexType,
         typename MaxIndexType>
bool increment(IndexType &index, const MaxIndexType &max_index)
{
    return increment_impl<NumDims, IndexType, MaxIndexType>(
            index, max_index, boost::is_array<IndexType>());
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
template<typename FPT, typename Integral>
FPT integer_power(FPT x, Integral n)
{
    using std::numeric_limits;

    // Avoid shooting ourselves by accidentally requesting a negative power for
    // an integer input.  Long lines to ensure messages appear in compiler
    // error output when used improperly.
    BOOST_STATIC_ASSERT(numeric_limits<Integral>::is_integer);
    BOOST_STATIC_ASSERT(!numeric_limits<Integral>::is_signed || !numeric_limits<FPT>::is_integer);

    FPT retval = 1.0;
    // Convert all requests into one involving a positive power
    if (n < 0) {
        x = ((FPT) 1)/x;
        n = -n;
    }
    // Repeated squaring method.  Returns 0.0^0 = 1.0; continuous in x
    do {
        if (n & 1) retval *= x;  /* for n odd */
        n >>= 1;
        x *= x;
    } while (n);

    return retval;
}

/**
 * Overwrite \c dest with <tt>src</tt>
 *
 * @param dest destination
 * @param src source
 */
template<typename FPT>
void assign_complex(fftw_complex &dest,
                    const std::complex<FPT> &src)
{
    dest[0] = src.real();
    dest[1] = src.imag();
}

/**
 * Overwrite \c dest with <tt>src</tt>
 *
 * @param dest destination
 * @param src source
 */
void assign_complex(fftw_complex &dest,
                    const fftw_complex &src)
{
    dest[0] = src[0];
    dest[1] = src[1];
}

/**
 * Overwrite \c dest with <tt>src_real</tt> + \f$\sqrt{-1}\f$ <tt>src_imag</tt>
 *
 * @param dest destination
 * @param src_real real part of the source
 * @param src_imag imag part of the source
 */
template<typename FPT>
void assign_complex(fftw_complex &dest,
                    const FPT &src_real,
                    const FPT &src_imag)
{
    dest[0] = src_real;
    dest[1] = src_imag;
}

/**
 * Overwrite \c dest with <tt>src_real</tt> + \f$\sqrt{-1}\f$ <tt>src_imag</tt>
 *
 * @param dest destination
 * @param src_real real part of the source
 * @param src_imag imag part of the source
 */
template<typename FPT1, typename FPT2>
void assign_complex(std::complex<FPT1> &dest,
                    const FPT2 &src_real,
                    const FPT2 &src_imag)
{
    dest.real() = src_real;
    dest.imag() = src_imag;
}

/**
 * Overwrite \c dest_real with Re <tt>src</tt> and \c dest_imag with Re
 * <tt>src_imag</tt>
 *
 * @param dest_real destination real part
 * @param dest_imag destination imag part
 * @param src source
 */
template<typename FPT>
void assign_components(FPT &dest_real,
                       FPT &dest_imag,
                       const std::complex<FPT> &src)
{
    dest_real = src.real();
    dest_imag = src.imag();
}

/**
 * Overwrite \c dest_real with Re <tt>src</tt> and \c dest_imag with Re
 * <tt>src_imag</tt>
 *
 * @param dest_real destination real part
 * @param dest_imag destination imag part
 * @param src source
 */
template<typename FPT>
void assign_components(FPT &dest_real,
                       FPT &dest_imag,
                       const fftw_complex &src)
{
    dest_real = src[0];
    dest_imag = src[1];
}

/**
 * Overwrite \c dest with <tt>alpha*src</tt>.
 *
 * @param dest destination
 * @param src source
 * @param alpha multiplicative real scaling factor
 */
template<typename FPT1, typename FPT2>
void assign_complex_scaled(fftw_complex &dest,
                           const std::complex<FPT1> &src,
                           const FPT2 alpha)
{
    dest[0] = alpha*src.real();
    dest[1] = alpha*src.imag();
}

/**
 * Overwrite \c dest with <tt>alpha*src</tt>.
 *
 * @param dest destination
 * @param src source
 * @param alpha multiplicative real scaling factor
 */
template<typename FPT>
void assign_complex_scaled(fftw_complex &dest,
                           const fftw_complex &src,
                           const FPT alpha)
{
    dest[0] = alpha*src[0];
    dest[1] = alpha*src[1];
}

/**
 * Overwrite \c dest with <tt>alpha*src_real + I*alpha*src_imag*I</tt> where
 * \c I is the imaginary unit.
 *
 * @param dest destination
 * @param src_real real part of the source
 * @param src_imag imag part of the source
 * @param alpha multiplicative real scaling factor
 */
template<typename FPT1, typename FPT2>
void assign_complex_scaled(fftw_complex &dest,
                           const FPT1 &src_real,
                           const FPT1 &src_imag,
                           const FPT2 alpha)

{
    dest[0] = alpha*src_real;
    dest[1] = alpha*src_imag;
}

/**
 * Overwrite \c dest with <tt>alpha*src_real + I*alpha*src_imag</tt> where \c I
 * is the imaginary unit.
 *
 * @param dest destination
 * @param src_real real part of the source
 * @param src_imag imag part of the source
 * @param alpha multiplicative real scaling factor
 */
template<typename FPT1, typename FPT2, typename FPT3>
void assign_complex_scaled(std::complex<FPT1> &dest,
                           const FPT2 &src_real,
                           const FPT2 &src_imag,
                           const FPT3 alpha)
{
    dest.real() = alpha*src_real;
    dest.imag() = alpha*src_imag;
}

/**
 * Overwrite \c dest_real with Re <tt>alpha*src</tt> and \c dest_imag
 * with Im <tt>alpha*src</tt>.
 *
 * @param dest_real destination real part
 * @param dest_imag destination imag part
 * @param src source
 * @param alpha multiplicative real scaling factor
 */
template<typename FPT1, typename FPT2>
void assign_components_scaled(FPT1 &dest_real,
                              FPT2 &dest_imag,
                              const std::complex<FPT1> &src,
                              const FPT2 alpha)
{
    dest_real = alpha*src.real();
    dest_imag = alpha*src.imag();
}

/**
 * Overwrite \c dest_real with Re <tt>alpha*src</tt> and \c dest_imag
 * with Im <tt>alpha*src_imag</tt>.
 *
 * @param dest_real destination real part
 * @param dest_imag destination imag part
 * @param src source
 */
template<typename FPT1, typename FPT2>
void assign_components_scaled(FPT1 &dest_real,
                              FPT1 &dest_imag,
                              const fftw_complex &src,
                              const FPT2 alpha)
{
    dest_real = alpha*src[0];
    dest_imag = alpha*src[1];
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
template<typename FPT1, typename FPT2>
void assign_complex_scaled_ipower(fftw_complex &dest,
                                  const std::complex<FPT1> &src,
                                  const FPT2 alpha,
                                  const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            dest[0] =  alpha*src.imag();
            dest[1] = -alpha*src.real();
            break;
        case 2: // I^2 = -1 = I^-2
            dest[0] = -alpha*src.real();
            dest[1] = -alpha*src.imag();
            break;
        case 1: // I^1 = I = I^-3
            dest[0] = -alpha*src.imag();
            dest[1] =  alpha*src.real();
            break;
        case 0: // I^0 = 1
            dest[0] =  alpha*src.real();
            dest[1] =  alpha*src.imag();
            break;
    }
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
template<typename FPT>
void assign_complex_scaled_ipower(fftw_complex &dest,
                                  const fftw_complex &src,
                                  const FPT alpha,
                                  const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            dest[0] =  alpha*src[1];
            dest[1] = -alpha*src[0];
            break;
        case 2: // I^2 = -1 = I^-2
            dest[0] = -alpha*src[0];
            dest[1] = -alpha*src[1];
            break;
        case 1: // I^1 = I = I^-3
            dest[0] = -alpha*src[1];
            dest[1] =  alpha*src[0];
            break;
        case 0: // I^0 = 1
            dest[0] =  alpha*src[0];
            dest[1] =  alpha*src[1];
            break;
    }
}

/**
 * Overwrite \c dest with <tt>(alpha*src_real +
 * I*alpha*src_imag)*I^ipower</tt> where \c I is the imaginary unit.
 *
 * @param dest destination
 * @param src_real real part of the source
 * @param src_imag imag part of the source
 * @param alpha multiplicative real scaling factor
 * @param ipower exponent on the imaginary unit to include in the scaling
 */
template<typename FPT1, typename FPT2>
void assign_complex_scaled_ipower(fftw_complex &dest,
                                  const FPT1 &src_real,
                                  const FPT1 &src_imag,
                                  const FPT2 alpha,
                                  const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            dest[0] =  alpha*src_imag;
            dest[1] = -alpha*src_real;
            break;
        case 2: // I^2 = -1 = I^-2
            dest[0] = -alpha*src_real;
            dest[1] = -alpha*src_imag;
            break;
        case 1: // I^1 = I = I^-3
            dest[0] = -alpha*src_imag;
            dest[1] =  alpha*src_real;
            break;
        case 0: // I^0 = 1
            dest[0] =  alpha*src_real;
            dest[1] =  alpha*src_imag;
            break;
    }
}

/**
 * Overwrite \c dest with <tt>(alpha*src_real +
 * I*alpha*src_imag)*I^ipower</tt> where \c I is the imaginary unit.
 *
 * @param dest destination
 * @param src_real real part of the source
 * @param src_imag imag part of the source
 * @param alpha multiplicative real scaling factor
 * @param ipower exponent on the imaginary unit to include in the scaling
 */
template<typename FPT1, typename FPT2, typename FPT3>
void assign_complex_scaled_ipower(std::complex<FPT1> &dest,
                                  const FPT2 &src_real,
                                  const FPT2 &src_imag,
                                  const FPT3 alpha,
                                  const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            dest.real() =  alpha*src_imag;
            dest.imag() = -alpha*src_real;
            break;
        case 2: // I^2 = -1 = I^-2
            dest.real() = -alpha*src_real;
            dest.imag() = -alpha*src_imag;
            break;
        case 1: // I^1 = I = I^-3
            dest.real() = -alpha*src_imag;
            dest.imag() =  alpha*src_real;
            break;
        case 0: // I^0 = 1
            dest.real() =  alpha*src_real;
            dest.imag() =  alpha*src_imag;
            break;
    }
}

/**
 * Overwrite \c dest_real with Re <tt>alpha*src*I^ipower</tt> and \c dest_imag
 * with Im <tt>alpha*src_imag*I^ipower</tt> where \c I is the imaginary unit.
 *
 * @param dest_real destination real part
 * @param dest_imag destination imag part
 * @param src source
 * @param alpha multiplicative real scaling factor
 * @param ipower exponent on the imaginary unit to include in the scaling
 */
template<typename FPT1, typename FPT2>
void assign_components_scaled_ipower(FPT1 &dest_real,
                                     FPT2 &dest_imag,
                                     const std::complex<FPT1> &src,
                                     const FPT2 alpha,
                                     const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            dest_real =  alpha*src.imag();
            dest_imag = -alpha*src.real();
            break;
        case 2: // I^2 = -1 = I^-2
            dest_real = -alpha*src.real();
            dest_imag = -alpha*src.imag();
            break;
        case 1: // I^1 = I = I^-3
            dest_real = -alpha*src.imag();
            dest_imag =  alpha*src.real();
            break;
        case 0: // I^0 = 1
            dest_real =  alpha*src.real();
            dest_imag =  alpha*src.imag();
            break;
    }
}

/**
 * Overwrite \c dest_real with Re <tt>alpha*src*I^ipower</tt> and \c dest_imag
 * with Re <tt>alpha*src_imag*I^ipower</tt> where \c I is the imaginary unit.
 *
 * @param dest_real destination real part
 * @param dest_imag destination imag part
 * @param src source
 * @param alpha multiplicative scaling factor
 */
template<typename FPT1, typename FPT2>
void assign_components_scaled_ipower(FPT1 &dest_real,
                                     FPT1 &dest_imag,
                                     const fftw_complex &src,
                                     const FPT2 alpha,
                                     const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            dest_imag = -alpha*src[0];
            dest_real =  alpha*src[1];
            break;
        case 2: // I^2 = -1 = I^-2
            dest_real = -alpha*src[0];
            dest_imag = -alpha*src[1];
            break;
        case 1: // I^1 = I = I^-3
            dest_imag =  alpha*src[0];
            dest_real = -alpha*src[1];
            break;
        case 0: // I^0 = 1
            dest_real =  alpha*src[0];
            dest_imag =  alpha*src[1];
            break;
    }
}

} // namespace detail

template<class ComplexMultiArray1, class ComplexMultiArray2>
void transform_c2c(const size_t transform_dim,
                   const ComplexMultiArray1 &in,
                   ComplexMultiArray2 &out,
                   const int fftw_sign,
                   const unsigned fftw_flags = 0)
{
    using boost::integer_traits;

    // Typedefs fixed separately each ComplexMultiArray template parameters
    typedef typename ComplexMultiArray1::element element1;
    typedef typename ComplexMultiArray2::element element2;
    BOOST_STATIC_ASSERT(sizeof(element1) == sizeof(element2));

    // Typedefs expected to be consistent across both template parameters
    BOOST_STATIC_ASSERT((boost::is_same<
                typename ComplexMultiArray1::index,
                typename ComplexMultiArray2::index
            >::value));
    typedef typename ComplexMultiArray1::index index;
    BOOST_STATIC_ASSERT((boost::is_same<
                typename ComplexMultiArray1::size_type,
                typename ComplexMultiArray2::size_type
            >::value));
    typedef typename ComplexMultiArray1::size_type size_type;
    BOOST_STATIC_ASSERT(    ComplexMultiArray1::dimensionality
                         == ComplexMultiArray2::dimensionality);
    const size_type dimensionality = ComplexMultiArray1::dimensionality;

    // Typedefs due to implementation choices
    typedef boost::array<index,dimensionality>       index_array;
    typedef int                                      shape_type; // Per FFTW
    typedef boost::array<shape_type,dimensionality>  shape_array;

    // Ensure we transform a dimension that exists in the data
    assert(0 <= transform_dim && transform_dim < dimensionality);
    // Ensure we are operating on a complex-valued array
    // C99 _Complex may require additional handling just below
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element1>::value)
           || (boost::is_same<element1, fftw_complex>::value));
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element2>::value)
           || (boost::is_same<element2, fftw_complex>::value));
    // Copy all shape information into integers well-suited for FFTW
    shape_array shape_in, shape_out;
    {
        const size_type * const p_shape_in  = in.shape();
        const size_type * const p_shape_out = out.shape();
        for (size_type n = 0; n < dimensionality; ++n) {
            // Ensure out shape at least as large as in shape for
            // all non-transformed directions
            assert(n == transform_dim || p_shape_in[n] <= p_shape_out[n]);
            // Ensure we won't accidentally truncate the shape values
            // Also checks that FFTW can handle the transform size
            assert(p_shape_in[n]  <= integer_traits<shape_type>::const_max);
            assert(p_shape_out[n] <= integer_traits<shape_type>::const_max);

            shape_in[n]  = p_shape_in[n];
            shape_out[n] = p_shape_out[n];
        }
    }
    // Ensure transformation direction parameter is sane
    assert(fftw_sign == FFTW_FORWARD || fftw_sign == FFTW_BACKWARD);
    // Transform size determined from the physical space dimension
    const shape_type transform_n = (fftw_sign == FFTW_FORWARD)
        ? shape_in[transform_dim] : shape_out[transform_dim];

    // We choose to always use an intermediate buffer for the transform:
    //  1) Avoids nuking in or out during FFTW planning
    //  2) Allows us to always enforce FFTW memory alignment recommendations
    //  3) Allows us to repeatedly apply the same, simple FFTW plan
    //  4) Minimizes the costs of potentially non-stride-1 access
    //  5) Always gives us in-place transform performance for the FFT
    //  6) Greatly simplifies transform, dealiasing, and differentiation code
    //  7) Simplifies moving to other FFT libraries in the future, e.g. ESSL
    boost::shared_array<fftw_complex> buffer(
        static_cast<fftw_complex *>(
            fftw_malloc(sizeof(fftw_complex)*transform_n)),
        std::ptr_fun(fftw_free));
    assert(buffer);

    // Construct the FFTW in-place plan for the transform buffer
    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan(
            fftw_plan_dft_1d(transform_n,
                             buffer.get(),
                             buffer.get(),
                             fftw_sign,
                             fftw_flags | FFTW_DESTROY_INPUT),
            std::ptr_fun(fftw_destroy_plan));
    assert(plan);

    // Dereference all constant parameters outside main processing loop
    // Pulls array information into our stack frame
    index_array index_bases_in, index_bases_out;
    {
        const index * const p_index_bases_in  = in.index_bases();
        const index * const p_index_bases_out = out.index_bases();
        for (size_type n = 0; n < dimensionality; ++n) {
            index_bases_in[n]  = p_index_bases_in[n];
            index_bases_out[n] = p_index_bases_out[n];
        }
    }
    const index      stride_in_transform_dim  = in.strides()[transform_dim];
    const index      stride_out_transform_dim = out.strides()[transform_dim];
    const shape_type shape_in_transform_dim   = shape_in[transform_dim];
    const shape_type shape_out_transform_dim  = shape_out[transform_dim];
    // Parameters that control dealiasing copy operations
    const shape_type last_n_pos_in         = shape_in_transform_dim/2;
    const shape_type last_n_pos_out        = shape_out_transform_dim/2;
    const shape_type last_n_pos_transform  = transform_n/2;
    const shape_type last_n_pos_copyin
        = std::min(last_n_pos_in, last_n_pos_transform);
    const shape_type last_n_pos_copyout
        = std::min(last_n_pos_out, last_n_pos_transform);
    const shape_type first_n_neg_in        = -(shape_in_transform_dim-1)/2;
    const shape_type first_n_neg_out       = -(shape_out_transform_dim-1)/2;
    const shape_type first_n_neg_transform = -(transform_n-1)/2;
    // Must scale data to account for possible differences in grid sizes
    // Must additionally normalize after backwards transform completes
    typedef BOOST_TYPEOF(buffer[0][0]) fftw_real;
    const fftw_real input_scale_factor
        = ((fftw_real) transform_n) / shape_in_transform_dim;
    const fftw_real output_scale_factor
        = ((fftw_real) shape_out_transform_dim) / transform_n
        * (fftw_sign == FFTW_FORWARD ? 1.0 : ((fftw_real) 1.0)/transform_n);

    // Prepare per-pencil outer loop index and loop bounds
    shape_array loop_shape(shape_in);   // Iterate over all dimensions...
    loop_shape[transform_dim] = 1;      // ...except the transformed one
    index_array loop_index = {{   }};   // Initialize to default value
    for (size_type n = 0; n < dimensionality; ++n) {
        assert(loop_index[n] == 0);     // Check initialization correct
    }
    index_array dereference_index;      // To be adjusted by index_bases
    fftw_complex * p_buffer;            // Used to walk the buffer
    shape_type n;                       // Used to track current wavenumber

    // TODO Walk fastest dimensions first in increment routine

    // Process each of the transform_dim pencils in turn
    do {
        // Obtain pointer to this input pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_in.begin(), dereference_index.begin(),
                       std::plus<index>());
        const element1 * p_pencil_in = &(in(dereference_index));

        // Copy input into transform buffer and pad any excess with zeros
        // TODO differentiate prior to FFTW_BACKWARD if requested
        p_buffer = buffer.get();
        for (n = 0; n <= last_n_pos_copyin; ++n) {
            detail::assign_complex_scaled(
                    *(p_buffer++), *p_pencil_in, input_scale_factor);
            p_pencil_in += stride_in_transform_dim;
        }
        if (shape_in_transform_dim > transform_n) {
            for (/* init from above */; n <= last_n_pos_in; ++n) {
                p_pencil_in += stride_in_transform_dim;
            }
            for (n = first_n_neg_in; n < first_n_neg_transform; ++n) {
                p_pencil_in += stride_in_transform_dim;
            }
        } else {
            for (/* init from above */; n <= last_n_pos_transform; ++n) {
                detail::assign_complex(*(p_buffer++), 0, 0);
            }
            for (n = first_n_neg_transform; n < first_n_neg_in; ++n) {
                detail::assign_complex(*(p_buffer++), 0, 0);
            }
        }
        for (/* init from above */; n <= -1; ++n) {
            detail::assign_complex_scaled(
                    *(p_buffer++), *p_pencil_in, input_scale_factor);
            p_pencil_in += stride_in_transform_dim;
        }

        fftw_execute(plan.get()); // Pull the strings!  Pull the strings!

        // Obtain pointer to this output pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_out.begin(), dereference_index.begin(),
                       std::plus<index>());
        element2 * p_pencil_out = &(out(dereference_index));

        // Copy transform buffer into output truncating auxiliary modes
        // TODO differentiate after FFTW_FORWARD if requested
        p_buffer = buffer.get();
        for (n = 0; n <= last_n_pos_copyout; ++n) {
            detail::assign_complex_scaled( *p_pencil_out,
                    (*p_buffer)[0], (*p_buffer)[1], output_scale_factor);
            ++p_buffer;
            p_pencil_out += stride_out_transform_dim;
        }
        if (shape_out_transform_dim > transform_n) {
            for (/* init from above */; n <= last_n_pos_out; ++n) {
                detail::assign_complex(*p_pencil_out, 0, 0);
                p_pencil_out += stride_out_transform_dim;
            }
            for (n = first_n_neg_out; n < first_n_neg_transform; ++n) {
                detail::assign_complex(*p_pencil_out, 0, 0);
                p_pencil_out += stride_out_transform_dim;
            }
        } else {
            for (/* init from above */; n <= last_n_pos_transform; ++n) {
                ++p_buffer;
            }
            for (n = first_n_neg_transform; n < first_n_neg_out; ++n) {
                ++p_buffer;
            }
        }
        for (/* init from above */; n <= -1; ++n) {
            detail::assign_complex_scaled(*p_pencil_out,
                        (*p_buffer)[0], (*p_buffer)[1], output_scale_factor);
            ++p_buffer;
            p_pencil_out += stride_out_transform_dim;
        }

    } while (detail::increment<dimensionality>(loop_index, loop_shape));

} /* transform_c2c */

} /* fftw_multi_array */ } /* suzerain */ } /* pecos */

#endif // PECOS_SUZERAIN_FFTW_MULTI_ARRAY_HPP
