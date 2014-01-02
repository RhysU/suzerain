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

#ifndef SUZERAIN_PENCILFFT_HPP
#define SUZERAIN_PENCILFFT_HPP

/** @file
 * Perform FFTs atop boost::multi_array using FFTW
 */

#include <fftw3.h>

#include <suzerain/common.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/math.hpp>
#include <suzerain/utility.hpp>

// TODO Much of transform_c2c, forward_r2c, and backward_c2r is boilerplate
// TODO Bomb when encountering in-place transform with mismatched strides

namespace suzerain {

/**
 * Provides routines for performing FFTs along a "pencil", one direction of a
 * potentially multidimensional array.   In particular, provides such routines
 * atop the Boost.MultiArray concept.
 *
 * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
 *      Boost.MultiArray</a> for more information on the MultiArray concept.
 */
namespace pencilfft {

/**
 * Implementation details underneath suzerain::pencilfft
 * @internal
 */
namespace internal {

/**
 * Increment the next appropriate index in \c indices according to
 * the upper bounds given in \c max_indices.  All parameters must have
 * at least \c NumDims elements.
 *
 * @param indices contains current values on input and at least one
 *                incremented value if method returns true.
 * @param max_indices contains the upper bounds on each index.
 * @param index_order contains the order in which the indices should be
 *                incremented.
 *
 * @pre  <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 * @pre  <tt>index_order[i] < NumDims</tt> for <tt>0<=i<NumDims</tt>
 * @post <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 *
 * @return true if at least one of \c indices was incremented
 *         and false otherwise.  The contents of \c indices are undefined
 *         when false is returned.
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
    // and NumDims is often small, do not break out of loop when overflow is
    // false.  The overflow == false condition causes an effective NOP.
    bool overflow = true;
    #pragma unroll
    for (std::size_t i = 0; i < NumDims; ++i, ++index_order) {

        index_type           &index     = indices[*index_order];
        const max_index_type &max_index = max_indices[*index_order];

        assert(1 <= max_index);
        assert(suzerain::is_nonnegative(index));
        assert(boost::numeric_cast<max_index_type>(index) < max_index);

        index    += overflow;           // Handle incoming overflow
        overflow  = index/max_index;    // Check outgoing overflow
        index    *= !overflow;          // Set to zero on outgoing
    }
    return !overflow;
}

/**
 * Decrement the next appropriate index in \c indices according to
 * the upper bounds given in \c max_indices.  All parameters must
 * have at least \c NumDims elements.
 *
 * @param indices contains current values on input and at least
 *                one decremented value if method returns true.
 * @param max_indices contains the upper bounds on each index.
 * @param index_order contains the order in which the indices should be
 *                decremented.
 *
 * @pre  <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 * @pre  <tt>index_order[i] < NumDims</tt> for <tt>0<=i<NumDims</tt>
 * @post <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 *
 * @return true if at least one of \c indices was decremented
 *         and false otherwise.  The contents of \c indices are undefined
 *         when false is returned.
 */
template<std::size_t NumDims,
         typename Mutable_RandomAccessIterator,
         typename RandomAccessIterator,
         typename InputIterator>
bool decrement(Mutable_RandomAccessIterator indices,
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
    // and NumDims is often small, do not break out of loop when underflow is
    // false.  The underflow == false condition causes an effective NOP.
    bool underflow = true;
    #pragma unroll
    for (std::size_t i = 0; i < NumDims; ++i, ++index_order) {

        index_type           &index     = indices[*index_order];
        const max_index_type &max_index = max_indices[*index_order];

        assert(1 <= max_index);
        assert(suzerain::is_nonnegative(index));
        assert(boost::numeric_cast<max_index_type>(index) < max_index);

        index     -= underflow;                  // Handle incoming underflow
        underflow  = (index == index_type(-1));  // Check outgoing underflow
        index     += underflow * max_index;      // Set max_index-1 on outgoing
    }
    return !underflow;
}

/**
 * Increment or decrement the next appropriate index in \c indices according
 * to the directions in \c indices_increasing and the upper bounds given in
 * \c max_indices.  All parameters must have at least \c NumDims elements.
 *
 * @param indices contains current values on input and at least one
 *                modified value if method returns true.
 * @param indices_increasing contains the direction in which each index
 *                should be modified.  Elements must be coercible to type bool.
 *                Values in \c indices are incremented if the corresponding
 *                element in \c indices_increasing is true.
 * @param max_indices contains the upper bounds on each index.
 * @param index_order contains the order in which the indices should be
 *                incremented or decremented.
 *
 * @pre  <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 * @pre  <tt>index_order[i] < NumDims</tt> for <tt>0<=i<NumDims</tt>
 * @post <tt>0 <= indices[i] < max_indices[i]</tt> for <tt>0<=i<NumDims</tt>
 *
 * @return true if at least one of \c indices was modified
 *         and false otherwise.  The contents of \c indices are undefined
 *         when false is returned.
 */
template<std::size_t NumDims,
         typename Mutable_RandomAccessIterator,
         typename RandomAccessIterator1,
         typename RandomAccessIterator2,
         typename InputIterator>
bool crement(Mutable_RandomAccessIterator indices,
             RandomAccessIterator1 indices_increasing,
             RandomAccessIterator2 max_indices,
             InputIterator index_order)
{
    BOOST_STATIC_ASSERT(NumDims > 0);
    BOOST_CONCEPT_ASSERT((boost::Mutable_RandomAccessIterator<Mutable_RandomAccessIterator>));
    BOOST_CONCEPT_ASSERT((boost::RandomAccessIterator<RandomAccessIterator1>));
    BOOST_CONCEPT_ASSERT((boost::RandomAccessIterator<RandomAccessIterator2>));
    BOOST_CONCEPT_ASSERT((boost::InputIterator<InputIterator>));

    typedef typename std::iterator_traits<
        Mutable_RandomAccessIterator>::value_type index_type;
    typedef typename std::iterator_traits<
        RandomAccessIterator2>::value_type max_index_type;
    typedef typename std::iterator_traits<
        InputIterator>::value_type index_order_type;

    BOOST_CONCEPT_ASSERT((boost::Integer<index_type>));
    BOOST_CONCEPT_ASSERT((boost::Integer<max_index_type>));
    BOOST_CONCEPT_ASSERT((boost::Integer<index_order_type>));

    // Because we have O(NumDims) integral operations in the fully unrolled
    // loop and NumDims is O(3), do not break out of loop when carry_bit is
    // false.  The carry_bit == false condition causes an effective NOP.
    bool carry_bit = true;
    #pragma unroll
    for (std::size_t i = 0; i < NumDims; ++i, ++index_order) {

        index_type           &index        = indices[*index_order];
        const bool           is_increasing = indices_increasing[*index_order];
        const max_index_type &max_index    = max_indices[*index_order];

        assert(1 <= max_index);
        assert(suzerain::is_nonnegative(index));
        assert(boost::numeric_cast<max_index_type>(index) < max_index);

        index += is_increasing*carry_bit - !is_increasing*carry_bit;
        const bool overflow  = is_increasing*(index/max_index);
        const bool underflow = !is_increasing*(index == index_type(-1));
        index    *= !overflow;             // Set to zero if new overflow
        index    += underflow * max_index; // Set max_index-1 if new underflow
        carry_bit = overflow + underflow;
    }
    return !carry_bit;
}

/**
 * Increment the next appropriate index in \c indices according to the upper
 * bounds given in \c max_indices.  Indices are incremented according to
 * their position within \c indices with the 0th index being the fastest.
 *
 * @param indices contains current values on input and at least
 *                one incremented value if method returns true.
 * @param max_indices contains the upper bounds on each index.
 *
 * @return true if at least one of \c indices was incremented
 *         and false otherwise.  The contents of \c indices are undefined
 *         when false is returned.
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
 * Decrement the next appropriate index in \c indices according to the upper
 * bounds given in \c max_indices.  Indices are incremented according to
 * their position within \c indices with the 0th index being the fastest.
 *
 * @param indices contains current values on input and at least
 *                one decremented value if method returns true.
 * @param max_indices contains the upper bounds on each index.
 *
 * @return true if at least one of \c indices was decremented
 *         and false otherwise.  The contents of \c indices are undefined
 *         when false is returned.
 *
 * @see decrement(Mutable_RandomAccessIterator, RandomAccessIterator, InputIterator) for more details.
 */
template<std::size_t NumDims,
         typename Mutable_RandomAccessIterator,
         typename RandomAccessIterator>
bool decrement(Mutable_RandomAccessIterator indices,
               RandomAccessIterator max_indices)
{
    return decrement<NumDims>(
            indices, max_indices, boost::make_counting_iterator(0));
}

/**
 * Increment or decrement the next appropriate index in \c indices according
 * to the directions in \c indices_increasing and the upper bounds given in
 * \c max_indices.  Indices are incremented according to
 * their position within \c indices with the 0th index being the fastest.
 *
 * @param indices contains current values on input and at least
 *                one modified value if method returns true.
 * @param indices_increasing contains the direction in which each index
 *                should be modified.  Elements must be coercible to type bool.
 *                Values in \c indices are incremented if the corresponding
 *                element in \c indices_increasing is true.
 * @param max_indices contains the upper bounds on each index.
 *
 * @return true if at least one of \c indices was modified and false
 *         otherwise.  The contents of \c indices are undefined when
 *         false is returned.
 *
 * @see crement(Mutable_RandomAccessIterator, RandomAccessIterator1, RandomAccessIterator2, InputIterator) for more details.
 */
template<std::size_t NumDims,
         typename Mutable_RandomAccessIterator,
         typename RandomAccessIterator1,
         typename RandomAccessIterator2>
bool crement(Mutable_RandomAccessIterator indices,
             RandomAccessIterator1 indices_increasing,
             RandomAccessIterator2 max_indices)
{
    return crement<NumDims>(indices,
                            indices_increasing,
                            max_indices,
                            boost::make_counting_iterator(0));
}

/**
 * Compares indices based on element magnitude comparison from some underlying
 * RandomAccessIterator.  Element magnitude is calculated using
 * <tt>std::abs</tt>.
 */
template<class RandomAccessIterator>
struct indexed_element_magnitude_comparator {

    BOOST_CONCEPT_ASSERT((boost::RandomAccessIterator<RandomAccessIterator>));

    /** Used to lookup elements during comparison */
    const RandomAccessIterator table_;

    /**
     * Create an instance that compares indices based on \c table.
     *
     * @param table to use for comparisons
     */
    indexed_element_magnitude_comparator(const RandomAccessIterator table)
        : table_(table) {};

    /**
     * Performs the comparison of two indices.
     *
     * @param xi left index to compare
     * @param yi right index to compare
     *
     * @return true if <tt>std::abs(table_[xi]) < std::abs(table_[yi])</tt>
     *         and false otherwise.
     */
    bool operator()(const std::size_t &xi, const std::size_t &yi) const
    {
        return std::abs(table_[xi]) < std::abs(table_[yi]);
    }
};

/**
 * Constructs an indexed_element_magnitude_comparator from \c table.
 *
 * @param table to use for indexed comparison.
 *
 * @return an indexed_element_magnitude_comparison using \c table for element
 *         lookup.
 */
template<class RandomAccessIterator>
inline
indexed_element_magnitude_comparator<RandomAccessIterator>
make_indexed_element_magnitude_comparator(const RandomAccessIterator table)
{
    return indexed_element_magnitude_comparator<RandomAccessIterator>(table);
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
    typedef suzerain::complex::traits::real<fftwf_complex>::type real_type;

    /** Corresponding precision FFTW complex type */
    typedef fftwf_complex fftw_complex_type;

    /** Corresponding precision FFTW plan type */
    typedef fftwf_plan fftw_plan_type;

    /** FFTW complex-to-complex planning */
    static fftw_plan_type plan_c2c_1d(int n,
                                      fftw_complex_type *in,
                                      fftw_complex_type *out,
                                      int sign,
                                      unsigned flags) {
        return fftwf_plan_dft_1d(n, in, out, sign, flags);
    }

    /** FFTW real-to-complex planning */
    static fftw_plan_type plan_r2c_1d(int n,
                                      real_type *in,
                                      fftw_complex_type *out,
                                      unsigned flags) {
        return fftwf_plan_dft_r2c_1d(n, in, out, flags);
    }

    /** FFTW complex-to-real planning */
    static fftw_plan_type plan_c2r_1d(int n,
                                      fftw_complex_type *in,
                                      real_type *out,
                                      unsigned flags) {
        return fftwf_plan_dft_c2r_1d(n, in, out, flags);
    }

    /** FFTW execute plan */
    static void execute_plan(fftw_plan_type p) {
        return fftwf_execute(p);
    }

    /** FFTW destroy plan */
    static void destroy_plan(fftw_plan_type p) {
        return fftwf_destroy_plan(p);
    }
};

/** FFT traits specialized for \c fftw_complex */
template<>
struct transform_traits<fftw_complex> {

    /** Complex number type */
    typedef fftw_complex complex_type;

    /** Real and imaginary components type */
    typedef suzerain::complex::traits::real<fftw_complex>::type real_type;

    /** Corresponding precision FFTW complex type */
    typedef fftw_complex fftw_complex_type;

    /** Corresponding precision FFTW plan type */
    typedef fftw_plan fftw_plan_type;

    /** FFTW complex-to-complex planning */
    static fftw_plan_type plan_c2c_1d(int n,
                                      fftw_complex_type *in,
                                      fftw_complex_type *out,
                                      int sign,
                                      unsigned flags) {
        return fftw_plan_dft_1d(n, in, out, sign, flags);
    }

    /** FFTW real-to-complex planning */
    static fftw_plan_type plan_r2c_1d(int n,
                                      real_type *in,
                                      fftw_complex_type *out,
                                      unsigned flags) {
        return fftw_plan_dft_r2c_1d(n, in, out, flags);
    }

    /** FFTW complex-to-real planning */
    static fftw_plan_type plan_c2r_1d(int n,
                                      fftw_complex_type *in,
                                      real_type *out,
                                      unsigned flags) {
        return fftw_plan_dft_c2r_1d(n, in, out, flags);
    }

    /** FFTW execute plan */
    static void execute_plan(fftw_plan_type p) {
        return fftw_execute(p);
    }

    /** FFTW destroy plan */
    static void destroy_plan(fftw_plan_type p) {
        return fftw_destroy_plan(p);
    }
};

/** FFT traits specialized for \c fftwl_complex */
template<>
struct transform_traits<fftwl_complex> {

    /** Complex number type */
    typedef fftwl_complex complex_type;

    /** Real and imaginary components type */
    typedef suzerain::complex::traits::real<fftwl_complex>::type real_type;

    /** Corresponding precision FFTW complex type */
    typedef fftwl_complex fftw_complex_type;

    /** Corresponding precision FFTW plan type */
    typedef fftwl_plan fftw_plan_type;

    /** FFTW complex-to-complex planning */
    static fftw_plan_type plan_c2c_1d(int n,
                                      fftw_complex_type *in,
                                      fftw_complex_type *out,
                                      int sign,
                                      unsigned flags) {
        return fftwl_plan_dft_1d(n, in, out, sign, flags);
    }

    /** FFTW real-to-complex planning */
    static fftw_plan_type plan_r2c_1d(int n,
                                      real_type *in,
                                      fftw_complex_type *out,
                                      unsigned flags) {
        return fftwl_plan_dft_r2c_1d(n, in, out, flags);
    }

    /** FFTW complex-to-real planning */
    static fftw_plan_type plan_c2r_1d(int n,
                                      fftw_complex_type *in,
                                      real_type *out,
                                      unsigned flags) {
        return fftwl_plan_dft_c2r_1d(n, in, out, flags);
    }

    /** FFTW execute plan */
    static void execute_plan(fftw_plan_type p) {
        return fftwl_execute(p);
    }

    /** FFTW destroy plan */
    static void destroy_plan(fftw_plan_type p) {
        return fftwl_destroy_plan(p);
    }
};

/** FFT traits specialized for <tt>std::complex<float></tt> */
template<>
struct transform_traits<std::complex<float> >
    : transform_traits<fftwf_complex> {

    /** Complex number type */
    typedef std::complex<float> complex_type;

};

/** FFT traits specialized for <tt>std::complex<double></tt> */
template<>
struct transform_traits<std::complex<double> >
    : transform_traits<fftw_complex> {

    /** Complex number type */
    typedef std::complex<double> complex_type;

};

/** FFT traits specialized for <tt>std::complex<long double></tt> */
template<>
struct transform_traits<std::complex<long double> >
    : transform_traits<fftwl_complex> {

    /** Complex number type */
    typedef std::complex<long double> complex_type;
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
    SUZERAIN_FORCEINLINE
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& dontcare) const
    {
        SUZERAIN_UNUSED(dontcare);
        suzerain::complex::assign_complex(dest, src);
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
    SUZERAIN_FORCEINLINE
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& dontcare) const
    {
        SUZERAIN_UNUSED(dontcare);
        suzerain::complex::assign_complex_scaled(dest, src, alpha_);
    }
};

/**
 * A copy-and-differentiate functor for manipulating complex values
 * in the context of a wavenumber index.
 **/
template<class ComplexDestination>
struct complex_copy_differentiate {
    /** Scaling factor type */
    typedef typename transform_traits<ComplexDestination>::real_type real_type;

    /** Derivative order to take within operator() */
    const int derivative_;

    /** Domain length-based factor used to find wavenumber from index */
    const real_type twopioverlength_;

    /**
     * Create a functor instance applying the given derivative operator
     * for a domain of given length.
     *
     * @param derivative derivative to take
     * @param length length of the domain
     */
    complex_copy_differentiate(
        const int derivative,
        const real_type length)
        : derivative_(derivative),
          twopioverlength_(2*boost::math::constants::pi<real_type>()/length) {}

    /**
     * Copies <tt>(2*pi*n*I/length)^derivative * src</tt> to \c dest
     * where \c I is the imaginary unit.
     *
     * @param dest destination
     * @param src source
     * @param n wavenumber index to use
     */
    template<class ComplexSource, typename SignedInteger>
    SUZERAIN_FORCEINLINE
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& n) const
    {
        suzerain::complex::assign_complex_scaled_ipower(
                dest,
                src,
                suzerain::math::integer_power(twopioverlength_*n, derivative_),
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
    typedef typename transform_traits<ComplexDestination>::real_type real_type;

    /** Scaling factor applied by the functor */
    const real_type alpha_;

    /** Derivative order to take within operator() */
    const int derivative_;

    /** Domain length-based factor used to find wavenumber from index */
    const real_type twopioverlength_;

    /**
     * Create a functor instance applying the given derivative operator
     * for a domain of given length and then scaling the result.
     *
     * @param alpha scaling factor to apply
     * @param derivative derivative to take
     * @param length length of the domain
     */
    complex_copy_scale_differentiate(
        const real_type alpha,
        const int derivative,
        const real_type length)
        : alpha_(alpha),
          derivative_(derivative),
          twopioverlength_(2*boost::math::constants::pi<real_type>()/length) {}

    /**
     * Copies <tt>(2*pi*n*I/length)^derivative * alpha * src</tt> to \c dest
     * where \c I is the imaginary unit.
     *
     * @param dest destination
     * @param src source
     * @param n wavenumber index to use
     */
    template<class ComplexSource, typename SignedInteger>
    SUZERAIN_FORCEINLINE
    void operator()(ComplexDestination &dest,
                    const ComplexSource &src,
                    const SignedInteger& n) const
    {
        suzerain::complex::assign_complex_scaled_ipower(
                dest,
                src,
                alpha_ * suzerain::math::integer_power(twopioverlength_*n,
                                                       derivative_),
                derivative_);
    }
};

/**
 * Transforms \c in to \c out taking into account FFTW's complex-to-complex
 * storage order and dealiasing considerations.  If \c size_in and \c size_out
 * are identical, then this routine performs a simple strided transformation.
 * Below, the number of "logical elements" means the number of physical space
 * complex-valued data points or wave space complex-valued wavenumbers.
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
void c2c_fullbuffer_process(OutputIterator out,
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

    // Code relies on loop conditions being checked prior to any loop entry
    SizeType n; // Always tracks current wavenumber during iteration
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
            suzerain::complex::assign_complex(*out, 0, 0);
            std::advance(out, stride_out);
        }
        for (n = first_n_neg_out; n < first_n_neg_in; ++n) {
            suzerain::complex::assign_complex(*out, 0, 0);
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
 * Transforms \c in to \c out taking into account FFTW's real-to-complex and
 * complex-to-real storage order and dealiasing considerations.  If \c size_in
 * and \c size_out are identical, then this routine performs a simple strided
 * transformation.  Below, the number of "logical elements" means the
 * corresponding number of physical space real-valued data points.
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
void c2c_halfbuffer_process(OutputIterator out,
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

    // Code relies on loop conditions being checked prior to any loop entry
    SizeType n; // Always tracks current wavenumber during iteration
    for (n = 0; n <= last_n_pos_copy; ++n) {
        c2c_buffer_processor(*out, *in, n);
        std::advance(in,  stride_in);
        std::advance(out, stride_out);
    }
    for (/* init from above */; n <= last_n_pos_out; ++n) {
        suzerain::complex::assign_complex(*out, 0, 0);
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
    using suzerain::is_nonnegative;

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
           || (boost::is_same<element1, fftwf_complex>::value)
           || (boost::is_same<element1, fftw_complex>::value)
           || (boost::is_same<element1, fftwl_complex>::value));
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element2>::value)
           || (boost::is_same<element2, fftwf_complex>::value)
           || (boost::is_same<element2, fftw_complex>::value)
           || (boost::is_same<element2, fftwl_complex>::value));
    // Ensure dimension of both in and out is consistent
    BOOST_STATIC_ASSERT(    ComplexMultiArray1::dimensionality
                         == ComplexMultiArray2::dimensionality);
    const size_type1 dimensionality = ComplexMultiArray1::dimensionality;

    // Typedefs due to implementation choices
    typedef array<index1,dimensionality>      index_array1;
    typedef array<index2,dimensionality>      index_array2;
    typedef int                               shape_type; // Per FFTW
    typedef array<shape_type,dimensionality>  shape_array;

    // Ensure we transform a dimension that exists in the data
    assert(is_nonnegative(transform_dim) && transform_dim < dimensionality);
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
    shared_array<fftw_complex_type> buffer(
        static_cast<fftw_complex_type *>(
            fftw_malloc(sizeof(fftw_complex_type)*transform_n)),
        std::ptr_fun(fftw_free));
    assert(buffer);
    typedef typename TransformTraits::fftw_plan_type fftw_plan_type;
    shared_ptr<
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
    for (size_type1 n = 0; n < dimensionality; ++n) { increment_order[n] = n; }
    std::sort(
        increment_order.begin(), increment_order.end(),
        internal::make_indexed_element_magnitude_comparator(in.strides()));

    // Process each of the transform_dim pencils in turn
    do {
        // Obtain index for the current input pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_in.begin(), dereference_index1.begin(),
                       std::plus<index1>());

        // Copy input into transform buffer performing any needed scaling, etc.
        if (fftw_sign == FFTW_BACKWARD && derivative != 0) {
            internal::c2c_fullbuffer_process(
                buffer.get(), transform_n, index1(1),
                &in(dereference_index1),
                shape_in_transform_dim,
                stride_in_transform_dim,
                internal::complex_copy_differentiate<fftw_complex_type>(
                    derivative, domain_length));
        } else {
            internal::c2c_fullbuffer_process(
                buffer.get(), transform_n, index1(1),
                &in(dereference_index1),
                shape_in_transform_dim,
                stride_in_transform_dim,
                internal::complex_copy());
        }

        // Pull the string!  Pull the string!
        TransformTraits::execute_plan(plan.get());

        // Obtain index for the current output pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_out.begin(), dereference_index2.begin(),
                       std::plus<index2>());

        // Copy transform buffer into output performing any needed scaling, etc.
        if (fftw_sign == FFTW_FORWARD) {
            typedef typename
                internal::transform_traits<element2>::real_type real_type;
            const real_type normalization = real_type(1)/transform_n;
            if (derivative == 0) {
                internal::c2c_fullbuffer_process(
                    &out(dereference_index2),
                    shape_out_transform_dim,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index2(1),
                    internal::complex_copy_scale<element2>(normalization));
            } else {
                internal::c2c_fullbuffer_process(
                    &out(dereference_index2),
                    shape_out_transform_dim,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index2(1),
                    internal::complex_copy_scale_differentiate<element2>(
                        normalization, derivative, domain_length));
            }
        } else {
            internal::c2c_fullbuffer_process(
                &out(dereference_index2),
                shape_out_transform_dim,
                stride_out_transform_dim,
                buffer.get(), transform_n, index2(1),
                internal::complex_copy());
        }

    } while (internal::increment<dimensionality>(loop_index.begin(),
                                                 loop_shape.begin(),
                                                 increment_order.begin()));
} /* transform_c2c */

} // namespace internal

/**
 * Perform a forward complex-to-complex FFT on each 1D "pencil" of \c in
 * storing the result in \c out.  Data is transformed from physical
 * space to wave space.  If desired, the data can be differentiated.
 * For the transform to proceed without any dealiasing, it must be true that
 * <tt>in.shape()[transform_dim] == out.shape()[transform_dim]</tt>.
 *
 * @param transform_dim zero-indexed dimension indicating the pencils to
 *                      be transformed.
 * @param in an instance modeling the MultiArray concept containing the
 *           complex-valued physical space input data to transform.
 * @param out an instance modeling the MultiArray concept to contain the
 *            complex-valued wave space output data from the transform.
 * @param domain_length Used when differentiating the data during the
 *                   transformation.
 * @param derivative If nonzero, the data is differentiated during the
 *                   transform process.
 * @param fftw_flags FFTW planner flags to use when computing the transform.
 *                   For example, \c FFTW_MEASURE or \c FFTW_PATIENT.
 *
 * @see internal::c2c_transform for more details on the transform process.
 */
template<class ComplexMultiArray1,
         class ComplexMultiArray2>
void forward_c2c(
    const size_t transform_dim,
    const ComplexMultiArray1 &in,
    ComplexMultiArray2 &out,
    const typename internal::transform_traits<
            typename ComplexMultiArray2::element         // wave space scalar
        >::real_type domain_length = 2.0*M_PI,
    const int derivative = 0,
    const unsigned fftw_flags = 0)
{
    return internal::transform_c2c<
            // Transform traits based on physical space types
            internal::transform_traits<typename ComplexMultiArray1::element>,
            ComplexMultiArray1,
            ComplexMultiArray2,
            typename internal::transform_traits<
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
 * Perform a backward complex-to-complex FFT on each 1D "pencil" of \c in
 * storing the result in \c out.  Data is transformed from wave space
 * to physical space.  If desired, the data can be differentiated.
 * For the transform to proceed without any dealiasing, it must be true that
 * <tt>in.shape()[transform_dim] == out.shape()[transform_dim]</tt>.
 *
 * @param transform_dim zero-indexed dimension indicating the pencils to
 *                      be transformed.
 * @param in an instance modeling the MultiArray concept containing the
 *           complex-valued wave space input data to transform.
 * @param out an instance modeling the MultiArray concept to contain the
 *            complex-valued physical space output data from transform.
 * @param domain_length Used when differentiating the data during the
 *                   transformation.
 * @param derivative If nonzero, the data is differentiated during the
 *                   transform process.
 * @param fftw_flags FFTW planner flags to use when computing the transform.
 *                   For example, \c FFTW_MEASURE or \c FFTW_PATIENT.
 *
 * @see internal::c2c_transform for more details on the transform process.
 */
template<class ComplexMultiArray1,
         class ComplexMultiArray2>
void backward_c2c(
    const size_t transform_dim,
    const ComplexMultiArray1 &in,
    ComplexMultiArray2 &out,
    const typename internal::transform_traits<
            typename ComplexMultiArray1::element         // wave space scalar
        >::real_type domain_length = 2.0*M_PI,
    const int derivative = 0,
    const unsigned fftw_flags = 0)
{
    return internal::transform_c2c<
            // Transform traits based on physical space types
            internal::transform_traits<typename ComplexMultiArray2::element>,
            ComplexMultiArray1,
            ComplexMultiArray2,
            typename internal::transform_traits<
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

/**
 * Perform a forward real-to-complex FFT on each 1D "pencil" of \c in storing
 * the result in \c out.  Data is transformed from physical space to
 * wave space.  If desired, the data can be differentiated.
 * For the transform to proceed without any dealiasing, it must be true that
 * <tt>in.shape()[transform_dim] == out.shape()[transform_dim]/2 + 1</tt>.
 *
 * @param transform_dim zero-indexed dimension indicating the pencils to
 *                      be transformed.
 * @param in an instance modeling the MultiArray concept containing the
 *           real-valued physical space input data to transform.
 * @param out an instance modeling the MultiArray concept to contain the
 *            complex-valued wave space output data from the transform.
 *            Only the nonnegative wavenumber indices are stored because the
 *            forward transform of real-valued data has conjugate symmetry.
 * @param domain_length Used when differentiating the data during the
 *                   transformation.
 * @param derivative If nonzero, the data is differentiated during the
 *                   transform process.
 * @param fftw_flags FFTW planner flags to use when computing the transform.
 *                   For example, \c FFTW_MEASURE or \c FFTW_PATIENT.
 *
 * @see internal::c2c_transform for more details on the transform process.
 */
template<class RealMultiArray,
         class ComplexMultiArray>
void forward_r2c(
    const size_t transform_dim,
    const RealMultiArray &in,
    ComplexMultiArray &out,
    const typename internal::transform_traits<
            typename ComplexMultiArray::element          // wave space scalar
        >::real_type domain_length = 2.0*M_PI,
    const int derivative      = 0,
    const unsigned fftw_flags = 0)
{
    using suzerain::is_nonnegative;

    // TransformTraits fixed by the wave space type
    typedef internal::transform_traits<
            typename ComplexMultiArray::element> transform_traits;

    // Typedefs fixed separately by MultiArray template parameters
    typedef typename RealMultiArray::element      element1;
    typedef typename ComplexMultiArray::element   element2;
    typedef typename RealMultiArray::index        index1;
    typedef typename ComplexMultiArray::index     index2;
    typedef typename RealMultiArray::size_type    size_type1;
    typedef typename ComplexMultiArray::size_type size_type2;

    // Ensure we are operating on real-valued input and complex-valued output
    // If available, C99 _Complex will be typedefed to fftw_complex, etc.
    BOOST_STATIC_ASSERT(
              (! boost::is_complex<element1>::value)
           && (! boost::is_same<element1, fftwf_complex>::value)
           && (! boost::is_same<element1, fftw_complex>::value)
           && (! boost::is_same<element1, fftwl_complex>::value));
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element2>::value)
           || (boost::is_same<element2, fftwf_complex>::value)
           || (boost::is_same<element2, fftw_complex>::value)
           || (boost::is_same<element2, fftwl_complex>::value));
    // Ensure dimension of both in and out is consistent
    BOOST_STATIC_ASSERT(    RealMultiArray::dimensionality
                         == ComplexMultiArray::dimensionality);
    const size_type1 dimensionality = RealMultiArray::dimensionality;

    // Typedefs due to implementation choices
    typedef array<index1,dimensionality>      index_array1;
    typedef array<index2,dimensionality>      index_array2;
    typedef int                               shape_type; // Per FFTW
    typedef array<shape_type,dimensionality>  shape_array;

    // Ensure we transform a dimension that exists in the data
    assert(is_nonnegative(transform_dim) && transform_dim < dimensionality);
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
    // Transform size determined from the physical space dimension
    const shape_type transform_n = shape_in[transform_dim];

    // Prepare the in-place transform buffer and construct the FFTW plan
    typedef typename transform_traits::fftw_complex_type fftw_complex_type;
    shared_array<fftw_complex_type> buffer(
        static_cast<fftw_complex_type *>(
            fftw_malloc(sizeof(fftw_complex_type)*(transform_n/2+1))),
        std::ptr_fun(fftw_free));
    assert(buffer);
    typedef typename transform_traits::fftw_plan_type fftw_plan_type;
    shared_ptr<
            typename boost::remove_pointer<fftw_plan_type>::type
        > plan(transform_traits::plan_r2c_1d(
                    transform_n,
                    reinterpret_cast<typename transform_traits::real_type *>(
                            buffer.get()),
                    buffer.get(),
                    fftw_flags | FFTW_DESTROY_INPUT),
               std::ptr_fun(transform_traits::destroy_plan));
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
    const shape_type shape_out_transform_dim  = shape_out[transform_dim];

    // Prepare per-pencil outer loop index and loop bounds
    shape_array loop_shape(shape_in);   // Iterate over all dimensions...
    loop_shape[transform_dim] = 1;      // ...except the transformed one
    index_array1 loop_index;            // Start at the highest index
    std::transform(loop_shape.begin(), loop_shape.end(), loop_index.begin(),
                   std::bind2nd(std::minus<shape_type>(), shape_type(1)));
    index_array1 dereference_index1;     // To be adjusted by index_bases
    index_array2 dereference_index2;     // To be adjusted by index_bases

    // Walk fastest dimensions first when decrementing across pencils
    index_array1 decrement_order;
    for (size_type1 n = 0; n < dimensionality; ++n) { decrement_order[n] = n; }
    std::sort(
        decrement_order.begin(), decrement_order.end(),
        internal::make_indexed_element_magnitude_comparator(in.strides()));

    // Process each of the transform_dim pencils in turn
    do {
        // Obtain index for the current input pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_in.begin(), dereference_index1.begin(),
                       std::plus<index1>());

        // Copy strided real input into stride one transform buffer
        {
            const element1 * p_first = &in(dereference_index1);
            const element1 * const p_last
                = p_first + transform_n * stride_in_transform_dim;
            typename transform_traits::real_type * p_result
                = reinterpret_cast<typename transform_traits::real_type *>(
                        buffer.get());
            while (p_first != p_last) {
                *p_result++ = *p_first;
                p_first += stride_in_transform_dim;
            }
        }

        // Pull the string!  Pull the string!
        transform_traits::execute_plan(plan.get());

        // Obtain index for the current output pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_out.begin(), dereference_index2.begin(),
                       std::plus<index2>());

        // Copy complex buffer into output performing any needed scaling, etc.
        {
            typedef typename
                internal::transform_traits<element2>::real_type real_type;
            const real_type normalization = real_type(1)/transform_n;
            if (derivative == 0) {
                internal::c2c_halfbuffer_process(
                    &out(dereference_index2),
                    (shape_out_transform_dim - 1)*2 /* logical */,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index2(1),
                    internal::complex_copy_scale<element2>(normalization));
            } else {
                internal::c2c_halfbuffer_process(
                    &out(dereference_index2),
                    (shape_out_transform_dim - 1)*2 /* logical */,
                    stride_out_transform_dim,
                    buffer.get(), transform_n, index2(1),
                    internal::complex_copy_scale_differentiate<element2>(
                        normalization, derivative, domain_length));
            }
        }

    } while (internal::decrement<dimensionality>(loop_index.begin(),
                                                 loop_shape.begin(),
                                                 decrement_order.begin()));
} /* forward_r2c */

/**
 * Perform a backward complex-to-real FFT on each 1D "pencil" of \c in storing
 * the result in \c out.  Data is transformed from wave space to
 * physical space.  If desired, the data can be differentiated.
 * For the transform to proceed without any dealiasing, it must be true that
 * <tt>in.shape()[transform_dim]/2 + 1 == out.shape()[transform_dim]</tt>.
 *
 * @param transform_dim zero-indexed dimension indicating the pencils to
 *                      be transformed.
 * @param in an instance modeling the MultiArray concept containing the
 *           complex-valued wave space input data to transform.
 *           The input consists of only nonnegative wavenumber indices.
 * @param out an instance modeling the MultiArray concept to contain the
 *            real-valued physical space output data from transform.
 * @param domain_length Used when differentiating the data during the
 *                   transformation.
 * @param derivative If nonzero, the data is differentiated during the
 *                   transform process.
 * @param fftw_flags FFTW planner flags to use when computing the transform.
 *                   For example, \c FFTW_MEASURE or \c FFTW_PATIENT.
 *
 * @see internal::c2c_transform for more details on the transform process.
 */
template<class ComplexMultiArray,
         class RealMultiArray>
void backward_c2r(
    const size_t transform_dim,
    const ComplexMultiArray &in,
    RealMultiArray &out,
    const typename internal::transform_traits<
            typename ComplexMultiArray::element          // wave space scalar
        >::real_type domain_length = 2.0*M_PI,
    const int derivative      = 0,
    const unsigned fftw_flags = 0)
{
    using suzerain::is_nonnegative;

    // TransformTraits fixed by the wave space type
    typedef internal::transform_traits<
            typename ComplexMultiArray::element> transform_traits;

    // Typedefs fixed separately by MultiArray template parameters
    typedef typename ComplexMultiArray::element   element1;
    typedef typename RealMultiArray::element      element2;
    typedef typename ComplexMultiArray::index     index1;
    typedef typename RealMultiArray::index        index2;
    typedef typename ComplexMultiArray::size_type size_type1;
    typedef typename RealMultiArray::size_type    size_type2;

    // Ensure we are operating on real-valued input and complex-valued output
    // If available, C99 _Complex will be typedefed to fftw_complex, etc.
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element1>::value)
           || (boost::is_same<element1, fftwf_complex>::value)
           || (boost::is_same<element1, fftw_complex>::value)
           || (boost::is_same<element1, fftwl_complex>::value));
    BOOST_STATIC_ASSERT(
              (! boost::is_complex<element2>::value)
           && (! boost::is_same<element2, fftwf_complex>::value)
           && (! boost::is_same<element2, fftw_complex>::value)
           && (! boost::is_same<element2, fftwl_complex>::value));
    // Ensure dimension of both in and out is consistent
    BOOST_STATIC_ASSERT(    ComplexMultiArray::dimensionality
                         == RealMultiArray::dimensionality);
    const size_type1 dimensionality = ComplexMultiArray::dimensionality;

    // Typedefs due to implementation choices
    typedef array<index1,dimensionality>      index_array1;
    typedef array<index2,dimensionality>      index_array2;
    typedef int                               shape_type; // Per FFTW
    typedef array<shape_type,dimensionality>  shape_array;

    // Ensure we transform a dimension that exists in the data
    assert(is_nonnegative(transform_dim) && transform_dim < dimensionality);
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
    // Transform size determined from the physical space dimension
    const shape_type transform_n = shape_out[transform_dim];

    // Prepare the in-place transform buffer and construct the FFTW plan
    typedef typename transform_traits::fftw_complex_type fftw_complex_type;
    shared_array<fftw_complex_type> buffer(
        static_cast<fftw_complex_type *>(
            fftw_malloc(sizeof(fftw_complex_type)*(transform_n/2+1))),
        std::ptr_fun(fftw_free));
    assert(buffer);
    typedef typename transform_traits::fftw_plan_type fftw_plan_type;
    shared_ptr<
            typename boost::remove_pointer<fftw_plan_type>::type
        > plan(transform_traits::plan_c2r_1d(
                    transform_n,
                    buffer.get(),
                    reinterpret_cast<typename transform_traits::real_type *>(
                            buffer.get()),
                    fftw_flags | FFTW_DESTROY_INPUT),
               std::ptr_fun(transform_traits::destroy_plan));
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
    for (size_type1 n = 0; n < dimensionality; ++n) { increment_order[n] = n; }
    std::sort(
        increment_order.begin(), increment_order.end(),
        internal::make_indexed_element_magnitude_comparator(in.strides()));

    // Process each of the transform_dim pencils in turn
    do {
        // Obtain index for the current input pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_in.begin(), dereference_index1.begin(),
                       std::plus<index1>());

        // Copy complex input into buffer performing any needed differentiation
        {
            if (derivative == 0) {
                internal::c2c_halfbuffer_process(
                    buffer.get(), transform_n, index1(1),
                    &in(dereference_index1),
                    (shape_in_transform_dim - 1)*2 /* logical */,
                    stride_in_transform_dim,
                    internal::complex_copy());
            } else {
                internal::c2c_halfbuffer_process(
                    buffer.get(), transform_n, index1(1),
                    &in(dereference_index1),
                    (shape_in_transform_dim - 1)*2 /* logical */,
                    stride_in_transform_dim,
                    internal::complex_copy_differentiate<fftw_complex_type>(
                        derivative, domain_length));
            }
        }

        // Pull the string!  Pull the string!
        transform_traits::execute_plan(plan.get());

        // Obtain index for the current output pencil's starting position
        std::transform(loop_index.begin(), loop_index.end(),
                       index_bases_out.begin(), dereference_index2.begin(),
                       std::plus<index2>());

        // Copy stride one transform buffer to strided real output
        {
            typedef typename transform_traits::real_type real_type;
            const real_type * p_first
                = reinterpret_cast<real_type *>(buffer.get());
            const real_type * const p_last = p_first + transform_n;
            element2 * p_result = &out(dereference_index2);
            while (p_first != p_last) {
                *p_result = *p_first++;
                p_result += stride_out_transform_dim;
            }
        }

    } while (internal::increment<dimensionality>(loop_index.begin(),
                                                 loop_shape.begin(),
                                                 increment_order.begin()));
} /* backward_c2r */

} /* pencilfft */ } /* suzerain */

#endif // SUZERAIN_PENCILFFT_HPP
