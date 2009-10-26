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
#include <boost/concept/assert.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/numeric/conversion/converter.hpp>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_pointer.hpp>
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

/**
 * Increment the next appropriate index in \c indices according to the upper
 * bounds given in \c max_indices.  Both \c indices and \c max_indices must
 * have at least \c NumDims elements.
 *
 * @param indices contains current values on input and one incremented
 *                value if method returns true.
 * @param max_indices contains the upper bounds on each index.
 * @param index_order contains the order in which the indices should be
 *                incremented.
 *
 * @pre  <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 * @pre  <tt>index_order[i] < NumDims</tt> for <tt>0<=i<NumDims</tt>
 * @post <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 *
 * @return true if one of \c indices was incremented and false otherwise.
 */
template<std::size_t NumDims,
         typename Mutable_RandomAccessIterator,
         typename RandomAccessIterator,
         typename InputIterator>
bool increment(Mutable_RandomAccessIterator indices,
               RandomAccessIterator max_indices,
               InputIterator index_order)
{
    BOOST_STATIC_ASSERT(NumDims > 0);
    BOOST_CONCEPT_ASSERT((boost::Mutable_RandomAccessIterator<Mutable_RandomAccessIterator>));
    BOOST_CONCEPT_ASSERT((boost::RandomAccessIterator<RandomAccessIterator>));
    BOOST_CONCEPT_ASSERT((boost::InputIterator<InputIterator>));

    typedef typename std::iterator_traits<
        Mutable_RandomAccessIterator>::value_type index_type;
    typedef typename std::iterator_traits<
        RandomAccessIterator>::value_type max_index_type;
    typedef typename std::iterator_traits<
        InputIterator>::value_type index_order_type;

    BOOST_CONCEPT_ASSERT((boost::Integer<index_type>));
    BOOST_CONCEPT_ASSERT((boost::Integer<max_index_type>));
    BOOST_CONCEPT_ASSERT((boost::Integer<index_order_type>));

    // Because we have 4*NumDims integral operations in the fully unrolled loop
    // and NumDims is often small, do not break out of loop when overflow == 0.
    // The overflow == 0 condition causes an effective NOP after occurring.
    bool overflow = 1;
    for (std::size_t i = 0; i < NumDims; ++i, ++index_order) {

        index_type           &index     = indices[*index_order];
        const max_index_type &max_index = max_indices[*index_order];

        assert(1 <= max_index);
        assert(0 <= index);
        assert(index < max_index);

        index    += overflow;           // Handle incoming overflow
        overflow  = index/max_index;    // Check outgoing overflow
        index    *= !overflow;          // Set to zero on outgoing
    }
    return !overflow;
}

/**
 * Increment the next appropriate index in \c indices according to the upper
 * bounds given in \c max_indices.  Indices are incremented according to
 * their position within \c indices with the 0th index being the fastest.
 *
 * @param indices contains current values on input and one incremented
 *                value if method returns true.
 * @param max_indices contains the upper bounds on each index.
 *
 * @return true if one of \c indices was incremented and false otherwise.
 *
 * @see increment(Mutable_RandomAccessIterator, RandomAccessIterator, InputIterator) for more details.
 */
template<std::size_t NumDims,
         typename Mutable_RandomAccessIterator,
         typename RandomAccessIterator>
bool increment(Mutable_RandomAccessIterator indices,
               RandomAccessIterator max_indices)
{
    return increment<NumDims>(
            indices, max_indices, boost::make_counting_iterator(0));
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
template<typename FPT>
void assign_complex(std::complex<FPT> &dest,
                    const fftw_complex &src)
{
    dest.real() = src[0];
    dest.imag() = src[1];
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
template<typename FPT1, typename FPT2>
void assign_complex_scaled(std::complex<FPT1> &dest,
                           const fftw_complex &src,
                           const FPT2 alpha)
{
     dest.real() = alpha*src[0];
     dest.imag() = alpha*src[1];
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
 * Overwrite \c dest with <tt>alpha*src</tt>.
 *
 * @param dest destination
 * @param src source
 * @param alpha multiplicative real scaling factor
 */
template<typename FPT1, typename FPT2, typename FPT3>
void assign_complex_scaled(std::complex<FPT1> &dest,
                           const std::complex<FPT2> &src,
                           const FPT3 alpha)
{
    dest.real() = alpha*src.real();
    dest.imag() = alpha*src.imag();
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
 * @warning \c dest and \c src must not refer to the same data
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
 * Overwrite \c dest with <tt>alpha*src*I^ipower</tt> where
 * \c I is the imaginary unit.
 *
 * @param dest destination
 * @param src source
 * @param alpha multiplicative real scaling factor
 * @param ipower exponent on the imaginary unit to include in the scaling
 */
template<typename FPT1, typename FPT2>
void assign_complex_scaled_ipower(std::complex<FPT1> &dest,
                                  const fftw_complex &src,
                                  const FPT2 alpha,
                                  const int ipower)
{
    switch (ipower & 3) { // Modulo-four-like operation for 2s complement
        case 3: // I^3 = -I = I^-1
            dest.real() =  alpha*src[1];
            dest.imag() = -alpha*src[0];
            break;
        case 2: // I^2 = -1 = I^-2
            dest.real() = -alpha*src[0];
            dest.imag() = -alpha*src[1];
            break;
        case 1: // I^1 = I = I^-3
            dest.real() = -alpha*src[1];
            dest.imag() =  alpha*src[0];
            break;
        case 0: // I^0 = 1
            dest.real() =  alpha*src[0];
            dest.imag() =  alpha*src[1];
            break;
    }
}

template<typename T>
struct complex_traits {};

template<typename FPT>
struct complex_traits<std::complex<FPT> > {
    typedef typename std::complex<FPT>::value_type value_type;
};

template<>
struct complex_traits<fftwf_complex> {
    typedef float value_type;
};

template<>
struct complex_traits<fftw_complex> {
    typedef double value_type;
};

template<>
struct complex_traits<fftwl_complex> {
    typedef long double value_type;
};

struct copy_complex {
    template<class ComplexDestination,
             class ComplexSource,
             typename SignedInteger>
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& /* Don't care */) const
    {
        assign_complex(dest, src);
    }
};

template<class ComplexDestination>
struct copy_complex_scaled {
    const typename complex_traits<ComplexDestination>::value_type alpha_;

    copy_complex_scaled(
        const typename complex_traits<ComplexDestination>::value_type alpha)
        : alpha_(alpha) {}

    template<class ComplexSource, typename SignedInteger>
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& /* Don't care */) const
    {
        assign_complex_scaled(dest, src, alpha_);
    }
};

template<class ComplexDestination>
struct copy_complex_differentiate {
    typedef typename complex_traits<ComplexDestination>::value_type scalar;

    const int ipower_;
    const scalar twopioverlength_;

    copy_complex_differentiate(
        const int ipower,
        const scalar length = 2.0*M_PI)
        : ipower_(ipower), twopioverlength_(2.0*M_PI/length) {}

    template<class ComplexSource, typename SignedInteger>
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& n) const
    {
        assign_complex_scaled_ipower(
            dest, src, integer_power(twopioverlength_*n, ipower_), n);
    }
};

template<class ComplexDestination>
struct copy_complex_scaled_differentiate {
    typedef typename complex_traits<ComplexDestination>::value_type scalar;

    const scalar alpha_;
    const int ipower_;
    const scalar twopioverlength_;

    copy_complex_scaled_differentiate(
        const scalar alpha,
        const int ipower,
        const scalar length = 2.0*M_PI)
        : alpha_(alpha), ipower_(ipower), twopioverlength_(2.0*M_PI/length) {}

    template<class ComplexSource, typename SignedInteger>
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& n) const
    {
        assign_complex_scaled_ipower(
            dest, src, alpha_ * integer_power(twopioverlength_*n, ipower_), n);
    }
};

template<class InputIterator,
         class OutputIterator,
         typename SizeType,
         typename StrideType,
         class C2CBufferProcessor>
void c2c_buffer_process(OutputIterator out,
                        const SizeType size_out,
                        const StrideType stride_out,
                        InputIterator  in,
                        const SizeType size_in,
                        const StrideType stride_in,
                        const C2CBufferProcessor c2c_buffer_processor)
{
    // Highest positive wavenumber
    const SizeType last_n_pos_out  = size_out / 2;
    const SizeType last_n_pos_in   = size_in  / 2;
    const SizeType last_n_pos_copy = std::min(last_n_pos_out, last_n_pos_in);
    // Lowest negative wavenumber
    const SizeType first_n_neg_out = -(size_out - 1)/2;
    const SizeType first_n_neg_in  = -(size_in  - 1)/2;

    SizeType n = 0; // Always tracks current wavenumber during iteration
    for (n = 0; n <= last_n_pos_copy; ++n) {
        c2c_buffer_processor(*out, *in, n);
        std::advance(in,  stride_in);
        std::advance(out, stride_out);
    }
    if (size_in > size_out) {
        for (/* init from above */; n <= last_n_pos_in; ++n) {
            std::advance(in, stride_in);
        }
        for (n = first_n_neg_in; n < first_n_neg_out; ++n) {
            std::advance(in, stride_in);
        }
    } else {
        for (/* init from above */; n <= last_n_pos_out; ++n) {
            detail::assign_complex(*out, 0, 0);
            std::advance(out, stride_out);
        }
        for (n = first_n_neg_out; n < first_n_neg_in; ++n) {
            detail::assign_complex(*out, 0, 0);
            std::advance(out, stride_out);
        }
    }
    for (/* init from above */; n <= -1; ++n) {
        c2c_buffer_processor(*out, *in, n);
        std::advance(in,  stride_in);
        std::advance(out, stride_out);
    }
}


} // namespace detail

template<class ComplexMultiArray1, class ComplexMultiArray2>
void transform_c2c(const size_t transform_dim,
                   const ComplexMultiArray1 &in,
                   ComplexMultiArray2 &out,
                   const int fftw_sign,
                   const int derivative = 0,
                   const unsigned fftw_flags = 0)
{
    // Typedefs fixed separately each ComplexMultiArray template parameters
    typedef typename ComplexMultiArray1::element element1;
    typedef typename ComplexMultiArray2::element element2;
    // Ensure we are operating on complex-valued arrays
    // C99 _Complex may require additional handling
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element1>::value)
           || (boost::is_same<element1, fftw_complex>::value));
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element2>::value)
           || (boost::is_same<element2, fftw_complex>::value));
    // Ensure element types occupy the same storage
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
    // Copy all shape information into integers well-suited for FFTW
    shape_array shape_in, shape_out;
    {
        typedef boost::numeric::converter<shape_type,size_type> Size2ShapeType;
        std::transform(in.shape(), in.shape() + dimensionality,
                    shape_in.begin(), Size2ShapeType());
        std::transform(out.shape(), out.shape() + dimensionality,
                    shape_out.begin(), Size2ShapeType());
    }
    // Ensure out shape at least as large as in shape for
    // all non-transformed directions
    for (size_type n = 0; n < dimensionality; ++n) {
        assert(n == transform_dim || shape_in[n] <= shape_out[n]);
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
    std::copy(in.index_bases(), in.index_bases() + dimensionality,
              index_bases_in.begin());
    std::copy(out.index_bases(), out.index_bases() + dimensionality,
              index_bases_out.begin());
    const index      stride_in_transform_dim  = in.strides()[transform_dim];
    const index      stride_out_transform_dim = out.strides()[transform_dim];
    const shape_type shape_in_transform_dim   = shape_in[transform_dim];
    const shape_type shape_out_transform_dim  = shape_out[transform_dim];

    // Prepare per-pencil outer loop index and loop bounds
    shape_array loop_shape(shape_in);   // Iterate over all dimensions...
    loop_shape[transform_dim] = 1;      // ...except the transformed one
    index_array loop_index = {{   }};   // Initialize to default value
    for (size_type n = 0; n < dimensionality; ++n) {
        assert(loop_index[n] == 0);     // Check initialization correct
    }
    index_array dereference_index;      // To be adjusted by index_bases

    // TODO Walk fastest dimensions first in increment routine

    // Process each of the transform_dim pencils in turn
    do {
        // Obtain index for the current input pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_in.begin(), dereference_index.begin(),
                       std::plus<index>());

        // Copy input into transform buffer performing any needed scaling, etc.
        if (fftw_sign == FFTW_BACKWARD && derivative != 0) {
            detail::c2c_buffer_process(
                buffer.get(), transform_n, index(1),
                &in(dereference_index),
                shape_in_transform_dim,
                stride_in_transform_dim,
                detail::copy_complex_differentiate<fftw_complex>(derivative));
        } else {
            detail::c2c_buffer_process(
                buffer.get(), transform_n, index(1),
                &in(dereference_index),
                shape_in_transform_dim,
                stride_in_transform_dim,
                detail::copy_complex());
        }

        // Pull the strings!  Pull the strings!
        fftw_execute(plan.get());

        // Obtain index for the current output pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_out.begin(), dereference_index.begin(),
                       std::plus<index>());

        // Copy transform buffer into output performing any needed scaling, etc.
        if (fftw_sign == FFTW_FORWARD) {
            typedef typename detail::complex_traits<element2>::value_type scalar;
            const scalar normalization = scalar(1.0)/transform_n;
            if (derivative == 0) {
                detail::c2c_buffer_process(
                    &out(dereference_index),
                    shape_out_transform_dim,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index(1),
                    detail::copy_complex_scaled<element2>(normalization));
            } else {
                detail::c2c_buffer_process(
                    &out(dereference_index),
                    shape_out_transform_dim,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index(1),
                    detail::copy_complex_scaled_differentiate<element2>(
                        normalization, derivative));
            }
        } else {
            detail::c2c_buffer_process(
                &out(dereference_index),
                shape_out_transform_dim,
                stride_out_transform_dim,
                buffer.get(), transform_n, index(1),
                detail::copy_complex());
        }

    } while (detail::increment<dimensionality>(loop_index.begin(),
                                               loop_shape.begin()));

} /* transform_c2c */

} /* fftw_multi_array */ } /* suzerain */ } /* pecos */

#endif // PECOS_SUZERAIN_FFTW_MULTI_ARRAY_HPP
