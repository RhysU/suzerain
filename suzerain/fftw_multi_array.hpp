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

// TODO Broken details::assign_* if FFTW3 discovers the C99 _Complex type

namespace pecos { namespace suzerain {

/**
 * Provides routines for performing FFTs atop the Boost.MultiArray concept.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
 *      Boost.MultiArray</a> for more information on the MultiArray concept.
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
 * Compares indices based on element comparison from some underlying
 * RandomAccessIterator.
 */
template<class RandomAccessIterator>
struct indexed_element_comparator {

    BOOST_CONCEPT_ASSERT((boost::RandomAccessIterator<RandomAccessIterator>));

    /** Used to lookup elements during comparison */
    const RandomAccessIterator table_;


    /**
     * Create an instance that compares indices based on \c table.
     *
     * @param table to use for comparisons
     */
    indexed_element_comparator(const RandomAccessIterator table)
        : table_(table) {};


    /**
     * Performs the comparison of two indices.
     *
     * @param xi left index to compare
     * @param yi right index to compare
     *
     * @return true if <tt>table_[xi] < table_[yi]</tt> and false otherwise.
     */
    bool operator()(const std::size_t &xi, const std::size_t &yi) const
    {
        return table_[xi] < table_[yi];
    }
};

/**
 * Constructs an indexed_element_comparator from \c table.
 *
 * @param table to use for indexed comparison.
 *
 * @return an indexed_element_comparison using \c table for element lookup.
 */
template<class RandomAccessIterator>
indexed_element_comparator<RandomAccessIterator>
make_indexed_element_comparator(const RandomAccessIterator table)
{
    return indexed_element_comparator<RandomAccessIterator>(table);
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
 * Obtain the real part of a complex number.
 *
 * @param z complex number stored as a two-element array.
 *
 * @return <tt>Re(z)</tt>
 */
template<typename FPT>
inline
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
inline
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
inline
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
inline
const FPT& imag(const FPT (&z)[2]) {
    return z[1];
}

/**
 * Overwrite \c dest with \c src.
 *
 * @param dest destination
 * @param src source
 */
template<class Complex1, class Complex2>
inline
void assign_complex(Complex1 &dest, const Complex2 &src)
{
    real(dest) = real(src);
    imag(dest) = imag(src);
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
inline
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
inline
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
inline
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
inline
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

/**
 * Primary template declaration for FFT traits.
 * Must be specialized for all types of interest.
 **/
template<typename Complex>
struct transform_traits;

/** FFT traits specialized for \c fftwf_complex */
template<>
struct transform_traits<fftwf_complex> {

    /** Complex number type */
    typedef fftwf_complex complex_type;

    /** Real and imaginary components type */
    typedef float real_type;

    /** Corresponding precision FFTW complex type */
    typedef fftwf_complex fftw_complex_type;

    /** Corresponding precision FFTW plan type */
    typedef fftwf_plan fftw_plan_type;

    /** Corresponding FFTW C2C planning function */
    static fftw_plan_type (* const plan_c2c_1d)(
                int, fftw_complex_type*, fftw_complex_type*, int, unsigned
            ) = &fftwf_plan_dft_1d;

    /** Corresponding FFTW execute plan function */
    static void (* const execute_plan)(fftw_plan_type) = &fftwf_execute;

    /** Corresponding FFTW destroy plan function */
    static void (* const destroy_plan)(fftw_plan_type) = &fftwf_destroy_plan;
};

/** FFT traits specialized for \c fftw_complex */
template<>
struct transform_traits<fftw_complex> {

    /** Complex number type */
    typedef fftw_complex complex_type;

    /** Real and imaginary components type */
    typedef double real_type;

    /** Corresponding precision FFTW complex type */
    typedef fftw_complex fftw_complex_type;

    /** Corresponding precision FFTW plan type */
    typedef fftw_plan fftw_plan_type;

    /** Corresponding FFTW C2C planning function */
    static fftw_plan_type (* const plan_c2c_1d)(
                int, fftw_complex_type*, fftw_complex_type*, int, unsigned
            ) = &fftw_plan_dft_1d;

    /** Corresponding FFTW execute plan function */
    static void (* const execute_plan)(fftw_plan_type) = &fftw_execute;

    /** Corresponding FFTW destroy plan function */
    static void (* const destroy_plan)(fftw_plan_type) = &fftw_destroy_plan;
};

/** FFT traits specialized for \c fftwl_complex */
template<>
struct transform_traits<fftwl_complex> {

    /** Complex number type */
    typedef fftwl_complex complex_type;

    /** Real and imaginary components type */
    typedef long double real_type;

    /** Corresponding precision FFTW complex type */
    typedef fftwl_complex fftw_complex_type;

    /** Corresponding precision FFTW plan type */
    typedef fftwl_plan fftw_plan_type;

    /** Corresponding FFTW C2C planning function */
    static fftw_plan_type (* const plan_c2c_1d)(
                int, fftw_complex_type*, fftw_complex_type*, int, unsigned
            ) = &fftwl_plan_dft_1d;

    /** Corresponding FFTW execute plan function */
    static void (* const execute_plan)(fftw_plan_type) = &fftwl_execute;

    /** Corresponding FFTW destroy plan function */
    static void (* const destroy_plan)(fftw_plan_type) = &fftwl_destroy_plan;
};

/** FFT traits specialized for <tt>std::complex<float></tt> */
template<>
struct transform_traits<std::complex<float> >
    : transform_traits<fftwf_complex> {

    /** Complex number type */
    typedef typename std::complex<float> complex_type;

};

/** FFT traits specialized for <tt>std::complex<double></tt> */
template<>
struct transform_traits<std::complex<double> >
    : transform_traits<fftw_complex> {

    /** Complex number type */
    typedef typename std::complex<double> complex_type;

};

/** FFT traits specialized for <tt>std::complex<long double></tt> */
template<>
struct transform_traits<std::complex<long double> >
    : transform_traits<fftwl_complex> {

    /** Complex number type */
    typedef typename std::complex<long double> complex_type;
};

/** A copy-only functor for manipulating complex values */
struct complex_copy {
    /**
     * Copies \c src to \c dest
     *
     * @param dest destination
     * @param src source
     * @param dontcare ignored within this functor
     */
    template<class ComplexDestination,
             class ComplexSource,
             typename SignedInteger>
    inline
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& dontcare) const
    {
        assign_complex(dest, src);
    }
};

/** A copy-and-scale functor for manipulating complex values */
template<class ComplexDestination>
struct complex_copy_scale {
    /** Scaling factor applied by the functor */
    const typename transform_traits<ComplexDestination>::real_type alpha_;

    /**
     * Create a functor instance with the supplied scaling factor
     *
     * @param alpha scaling factor to use within operator()
     */
    complex_copy_scale(
        const typename transform_traits<ComplexDestination>::real_type alpha)
        : alpha_(alpha) {}

    /** Copies <tt>alpha*src</tt> to \c dest.
     *
     * @param dest destination
     * @param src source
     * @param dontcare ignored within this functor
     */
    template<class ComplexSource, typename SignedInteger>
    inline
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& dontcare) const
    {
        assign_complex_scaled(dest, src, alpha_);
    }
};

/**
 * A copy-and-differentiate functor for manipulating complex values
 * in the context of a wavenumber index.
 **/
template<class ComplexDestination>
struct complex_copy_differentiate {
    /** Scaling factor type */
    typedef typename transform_traits<ComplexDestination>::real_type scalar;

    /** Derivative order to take within operator() */
    const int derivative_;

    /** Domain length-based factor used to find wavenumber from index */
    const scalar twopioverlength_;

    /**
     * Create a functor instance applying the given derivative operator
     * for a domain of given length.
     *
     * @param derivative derivative to take
     * @param length length of the domain
     */
    complex_copy_differentiate(
        const int derivative,
        const scalar length)
        : derivative_(derivative), twopioverlength_(2.0*M_PI/length) {}

    /**
     * Copies <tt>(2*pi*n*I/length)^derivative * src</tt> to \c dest
     * where \c I is the imaginary unit.
     *
     * @param dest destination
     * @param src source
     * @param n wavenumber index to use
     */
    template<class ComplexSource, typename SignedInteger>
    inline
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& n) const
    {
        assign_complex_scaled_ipower(
                dest,
                src,
                integer_power(twopioverlength_*n, derivative_),
                derivative_);
    }
};

/**
 * A copy-scale-and-differentiate functor for manipulating complex values
 * in the context of a wavenumber index.
 **/
template<class ComplexDestination>
struct complex_copy_scale_differentiate {
    /** Scaling factor type */
    typedef typename transform_traits<ComplexDestination>::real_type scalar;

    /** Scaling factor applied by the functor */
    const scalar alpha_;

    /** Derivative order to take within operator() */
    const int derivative_;

    /** Domain length-based factor used to find wavenumber from index */
    const scalar twopioverlength_;

    /**
     * Create a functor instance applying the given derivative operator
     * for a domain of given length and then scaling the result.
     *
     * @param alpha scaling factor to apply
     * @param derivative derivative to take
     * @param length length of the domain
     */
    complex_copy_scale_differentiate(
        const scalar alpha,
        const int derivative,
        const scalar length)
        : alpha_(alpha),
          derivative_(derivative),
          twopioverlength_(2.0*M_PI/length) {}

    /**
     * Copies <tt>(2*pi*n*I/length)^derivative * alpha * src</tt> to \c dest
     * where \c I is the imaginary unit.
     *
     * @param dest destination
     * @param src source
     * @param n wavenumber index to use
     */
    template<class ComplexSource, typename SignedInteger>
    inline
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& n) const
    {
        assign_complex_scaled_ipower(
                dest,
                src,
                alpha_ * integer_power(twopioverlength_*n, derivative_),
                derivative_);
    }
};

/**
 * Transforms \c in to \c out taking into account FFTW's complex-to-complex
 * storage order and dealiasing considerations.  If \c size_in and \c size_out
 * are identical, then this routine performs a simple strided transformation.
 *
 * @param out start of the output storage
 * @param size_out number of logical elements in the output storage
 * @param stride_out stride between logical elements in the output storage
 * @param in start of the input storage
 * @param size_in number of logical elements in the input storage
 * @param stride_in stride between logical elements in the input storage
 * @param c2c_buffer_processor a functor to apply at each location.
 *
 * @see complex_copy, complex_copy_scale, complex_copy_differentiate, and
 *      complex_copy_scale_differentiate for valid \c c2c_buffer_processor
 *      values and a better idea of the \c c2c_buffer_processor concept
 * @see the <a href="http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html">FFTW documentation</a> for the expected wavenumber order.
 */
template<class OutputIterator,
         class InputIterator,
         typename SizeType,
         typename StrideType,
         class C2CBufferProcessor>
void c2c_buffer_process(OutputIterator out,
                        const SizeType size_out,
                        const StrideType stride_out,
                        InputIterator  in,
                        const SizeType size_in,
                        const StrideType stride_in,
                        const C2CBufferProcessor &c2c_buffer_processor)
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

/**
 * Perform a FFT on each 1D "pencil" of \c in storing the result in \c out.  If
 * desired, the data can be differentiated.  Processing is done in a way that
 * attempts to minimize memory access under the assumption that either \c in or
 * \c out will require non-stride-one access.  If \c in and \c out are two
 * distinct data sets, then \c in is not destroyed in the computation.  Data is
 * normalized in wave space in a transform-size-independent way.
 *
 * @param transform_dim zero-indexed dimension indicating the pencils to
 *                      be transformed.
 * @param in an instance modeling the MultiArray concept containing the
 *           complex input data to transform.
 * @param out an instance modeling the MultiArray concept to contain the
 *            complex output data from the transform.
 * @param fftw_sign either \c FFTW_FORWARD or \c FFTW_BACKWARD to transform
 *                  from physical to wave space or from wave space to
 *                  physical space, respectively.
 * @param domain_length Used when differentiating the data during the
 *                   transformation.
 * @param derivative If nonzero, the data is differentiated during the transform
 *                   process.
 * @param fftw_flags FFTW planner flags to use when computing the transform.
 *                   For example, \c FFTW_MEASURE or \c FFTW_PATIENT.
 *
 * @see The FFTW documentation for <a href="http://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html">
 *      the expected wavenumber order</a> in each 1D wave space pencil.
 * @see The FFTW documentation for information on the
 *      <a href="http://www.fftw.org/fftw3_doc/Planner-Flags.html">
 *      available planner flags</a>.
 *
 * @note We choose to always use an intermediate buffer for the transform:
 *    @li Avoids nuking in or out during FFTW planning
 *    @li Allows us to always enforce FFTW memory alignment recommendations
 *    @li Allows us to repeatedly apply the same, simple FFTW plan
 *    @li Minimizes the costs of potentially non-stride-1 access
 *    @li Always gives us in-place transform performance for the FFT
 *    @li Greatly simplifies transform, dealiasing, and differentiation code
 *    @li Simplifies moving to other FFT libraries in the future, e.g. ESSL
 */
template<class TransformTraits,
         class ComplexMultiArray1,
         class ComplexMultiArray2,
         typename DomainLengthType>
void transform_c2c(
    const size_t transform_dim,
    const ComplexMultiArray1 &in,
    ComplexMultiArray2 &out,
    const int fftw_sign,
    const DomainLengthType domain_length,
    const int derivative,
    const unsigned fftw_flags)
{
    // Typedefs fixed separately by ComplexMultiArray template parameters
    typedef typename ComplexMultiArray1::element element1;
    typedef typename ComplexMultiArray2::element element2;
    typedef typename ComplexMultiArray1::index index1;
    typedef typename ComplexMultiArray2::index index2;
    typedef typename ComplexMultiArray1::size_type size_type1;
    typedef typename ComplexMultiArray2::size_type size_type2;

    // Ensure we are operating on complex-valued arrays
    // If available, C99 _Complex will be typedefed to fftw_complex, etc.
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element1>::value)
           || (boost::is_same<element1, fftw_complex>::value));
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element2>::value)
           || (boost::is_same<element2, fftw_complex>::value));
    // Ensure element types occupy the same storage
    BOOST_STATIC_ASSERT(sizeof(element1) == sizeof(element2));
    // Ensure dimension of both in and out is consistent
    BOOST_STATIC_ASSERT(    ComplexMultiArray1::dimensionality
                         == ComplexMultiArray2::dimensionality);
    const size_type1 dimensionality = ComplexMultiArray1::dimensionality;

    // Typedefs due to implementation choices
    typedef boost::array<index1,dimensionality>      index_array1;
    typedef boost::array<index2,dimensionality>      index_array2;
    typedef int                                      shape_type; // Per FFTW
    typedef boost::array<shape_type,dimensionality>  shape_array;

    // Ensure we transform a dimension that exists in the data
    assert(0 <= transform_dim && transform_dim < dimensionality);
    // Copy all shape information into integers well-suited for FFTW
    shape_array shape_in, shape_out;
    {
        using boost::numeric::converter;
        typedef converter<shape_type,size_type1> Size1ToShapeType;
        std::transform(in.shape(), in.shape() + dimensionality,
                       shape_in.begin(), Size1ToShapeType());
        typedef converter<shape_type,size_type2> Size2ToShapeType;
        std::transform(out.shape(), out.shape() + dimensionality,
                       shape_out.begin(), Size2ToShapeType());
    }
    // Ensure out shape at least as large as in shape for
    // all non-transformed directions
    for (size_type1 n = 0; n < dimensionality; ++n) {
        assert(n == transform_dim || shape_in[n] <= shape_out[n]);
    }
    // Ensure transformation direction parameter is sane
    assert(fftw_sign == FFTW_FORWARD || fftw_sign == FFTW_BACKWARD);
    // Transform size determined from the physical space dimension
    const shape_type transform_n = (fftw_sign == FFTW_FORWARD)
        ? shape_in[transform_dim] : shape_out[transform_dim];

    // Prepare the in-place transform buffer and construct the FFTW plan
    typedef typename TransformTraits::fftw_complex_type fftw_complex_type;
    boost::shared_array<fftw_complex_type> buffer(
        static_cast<fftw_complex_type *>(
            fftw_malloc(sizeof(fftw_complex_type)*transform_n)),
        std::ptr_fun(fftw_free));
    assert(buffer);
    typedef typename TransformTraits::fftw_plan_type fftw_plan_type;
    boost::shared_ptr<
            typename boost::remove_pointer<fftw_plan_type>::type
        > plan(TransformTraits::plan_c2c_1d(transform_n,
                                            buffer.get(),
                                            buffer.get(),
                                            fftw_sign,
                                            fftw_flags | FFTW_DESTROY_INPUT),
               std::ptr_fun(TransformTraits::destroy_plan));
    assert(plan);

    // Dereference all constant parameters outside main processing loop
    // Pulls array information into our stack frame
    index_array1 index_bases_in;
    std::copy(in.index_bases(), in.index_bases() + dimensionality,
              index_bases_in.begin());
    index_array2 index_bases_out;
    std::copy(out.index_bases(), out.index_bases() + dimensionality,
              index_bases_out.begin());
    const index1     stride_in_transform_dim  = in.strides()[transform_dim];
    const index2     stride_out_transform_dim = out.strides()[transform_dim];
    const shape_type shape_in_transform_dim   = shape_in[transform_dim];
    const shape_type shape_out_transform_dim  = shape_out[transform_dim];

    // Prepare per-pencil outer loop index and loop bounds
    shape_array loop_shape(shape_in);   // Iterate over all dimensions...
    loop_shape[transform_dim] = 1;      // ...except the transformed one
    index_array1 loop_index = {{   }};  // Initialize to default value
    for (size_type1 n = 0; n < dimensionality; ++n) {
        assert(loop_index[n] == 0);     // Check initialization correct
    }
    index_array1 dereference_index1;     // To be adjusted by index_bases
    index_array2 dereference_index2;     // To be adjusted by index_bases

    // Walk fastest dimensions first when incrementing across pencils
    index_array1 increment_order;
    for (index1 n = 0; n < dimensionality; ++n) { increment_order[n] = n; }
    std::stable_sort(increment_order.begin(), increment_order.end(),
                     detail::make_indexed_element_comparator(in.strides()));

    // Process each of the transform_dim pencils in turn
    do {
        // Obtain index for the current input pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_in.begin(), dereference_index1.begin(),
                       std::plus<index1>());

        // Copy input into transform buffer performing any needed scaling, etc.
        if (fftw_sign == FFTW_BACKWARD && derivative != 0) {
            detail::c2c_buffer_process(
                buffer.get(), transform_n, index1(1),
                &in(dereference_index1),
                shape_in_transform_dim,
                stride_in_transform_dim,
                detail::complex_copy_differentiate<fftw_complex>(
                    derivative, domain_length));
        } else {
            detail::c2c_buffer_process(
                buffer.get(), transform_n, index1(1),
                &in(dereference_index1),
                shape_in_transform_dim,
                stride_in_transform_dim,
                detail::complex_copy());
        }

        // Pull the strings!  Pull the strings!
        TransformTraits::execute_plan(plan.get());

        // Obtain index for the current output pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_out.begin(), dereference_index2.begin(),
                       std::plus<index2>());

        // Copy transform buffer into output performing any needed scaling, etc.
        if (fftw_sign == FFTW_FORWARD) {
            typedef typename
                detail::transform_traits<element2>::real_type scalar;
            const scalar normalization = scalar(1.0)/transform_n;
            if (derivative == 0) {
                detail::c2c_buffer_process(
                    &out(dereference_index2),
                    shape_out_transform_dim,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index2(1),
                    detail::complex_copy_scale<element2>(normalization));
            } else {
                detail::c2c_buffer_process(
                    &out(dereference_index2),
                    shape_out_transform_dim,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index2(1),
                    detail::complex_copy_scale_differentiate<element2>(
                        normalization, derivative, domain_length));
            }
        } else {
            detail::c2c_buffer_process(
                &out(dereference_index2),
                shape_out_transform_dim,
                stride_out_transform_dim,
                buffer.get(), transform_n, index2(1),
                detail::complex_copy());
        }

    } while (detail::increment<dimensionality>(loop_index.begin(),
                                               loop_shape.begin(),
                                               increment_order.begin()));
} /* transform_c2c */

} // namespace detail

/**
 * Perform a forward FFT on each 1D "pencil" of \c in storing the result in \c
 * out.  That is, data is transformed from physical space to wave space.  If
 * desired, the data can be differentiated.
 *
 * @param transform_dim zero-indexed dimension indicating the pencils to
 *                      be transformed.
 * @param in an instance modeling the MultiArray concept containing the
 *           complex physical space input data to transform.
 * @param out an instance modeling the MultiArray concept to contain the
 *            complex wave space output data from the transform.
 * @param domain_length Used when differentiating the data during the
 *                   transformation.
 * @param derivative If nonzero, the data is differentiated during the transform
 *                   process.
 * @param fftw_flags FFTW planner flags to use when computing the transform.
 *                   For example, \c FFTW_MEASURE or \c FFTW_PATIENT.
 *
 * @see detail::c2c_transform for more details on the transform process.
 */
template<class ComplexMultiArray1,
         class ComplexMultiArray2>
void forward_c2c(
    const size_t transform_dim,
    const ComplexMultiArray1 &in,
    ComplexMultiArray2 &out,
    const typename detail::transform_traits<
            typename ComplexMultiArray2::element         // wave space scalar
        >::real_type domain_length = 2.0*M_PI,
    const int derivative = 0,
    const unsigned fftw_flags = 0)
{
    return detail::transform_c2c<
            // Transform traits based on physical space types
            detail::transform_traits<typename ComplexMultiArray1::element>,
            ComplexMultiArray1,
            ComplexMultiArray2,
            typename detail::transform_traits<
                    typename ComplexMultiArray2::element // wave space scalar
                >::real_type
        >(
            transform_dim,
            in,
            out,
            FFTW_FORWARD,
            domain_length,
            derivative,
            fftw_flags
         );
}

/**
 * Perform a backward FFT on each 1D "pencil" of \c in storing the result in \c
 * out.  That is, data is transformed from wave space to physical space.  If
 * desired, the data can be differentiated.
 *
 * @param transform_dim zero-indexed dimension indicating the pencils to
 *                      be transformed.
 * @param in an instance modeling the MultiArray concept containing the
 *           complex wave space input data to transform.
 * @param out an instance modeling the MultiArray concept to contain the
 *            complex physical space output data from transform.
 * @param domain_length Used when differentiating the data during the
 *                   transformation.
 * @param derivative If nonzero, the data is differentiated during the transform
 *                   process.
 * @param fftw_flags FFTW planner flags to use when computing the transform.
 *                   For example, \c FFTW_MEASURE or \c FFTW_PATIENT.
 *
 * @see detail::c2c_transform for more details on the transform process.
 */
template<class ComplexMultiArray1,
         class ComplexMultiArray2>
void backward_c2c(
    const size_t transform_dim,
    const ComplexMultiArray1 &in,
    ComplexMultiArray2 &out,
    const typename detail::transform_traits<
            typename ComplexMultiArray1::element         // wave space scalar
        >::real_type domain_length = 2.0*M_PI,
    const int derivative = 0,
    const unsigned fftw_flags = 0)
{
    return detail::transform_c2c<
            // Transform traits based on physical space types
            detail::transform_traits<typename ComplexMultiArray2::element>,
            ComplexMultiArray1,
            ComplexMultiArray2,
            typename detail::transform_traits<
                    typename ComplexMultiArray1::element // wave space scalar
                >::real_type
        >(
            transform_dim,
            in,
            out,
            FFTW_BACKWARD,
            domain_length,
            derivative,
            fftw_flags
         );
}

} /* fftw_multi_array */ } /* suzerain */ } /* pecos */

#endif // PECOS_SUZERAIN_FFTW_MULTI_ARRAY_HPP
