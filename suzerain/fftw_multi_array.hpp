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
#include <fftw3.h>
#include <boost/array.hpp>
#include <boost/integer_traits.hpp>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/typeof/typeof.hpp>

namespace pecos { namespace suzerain { namespace fftw_multi_array {

template<std::size_t NumDims, typename IndexType, typename MaxIndexType>
bool increment(IndexType &index, const MaxIndexType &max_index)
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
    for (std::size_t n = 0; n < NumDims; ++n) {
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

/**
 * Implementation details underneath pecos::suzerain::fftw_multi_array
 * @internal
 */
namespace detail {

/**
 * Overwrite \c dest with \c src
 *
 * @param dest destination
 * @param src source
 */
template<typename FPT>
void assign_complex(fftw_complex &dest, const std::complex<FPT> &src)
{
    dest[0] = src.real();
    dest[1] = src.imag();
}

/**
 * Overwrite \c dest with \c src
 *
 * @param dest destination
 * @param src source
 */
void assign_complex(fftw_complex &dest, const fftw_complex &src)
{
    dest[0] = src[0];
    dest[1] = src[1];
}

/**
 * Overwrite \c dest with \c src_real + \f$\sqrt{-1}\f$ \c src_imag
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
 * Overwrite \c dest with \c src_real + \f$\sqrt{-1}\f$ \c src_imag
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
 * Overwrite \c dest_real with Re \c src and \c dest_imag with Re \c src_imag
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
 * Overwrite \c dest_real with Re \c src and \c dest_imag with Re \c src_imag
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

} // namespace detail

template<class ComplexMultiArray1, class ComplexMultiArray2>
void c2c_transform(const size_t transform_dim,
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
    const shape_type last_kn_pos_in         = shape_in_transform_dim/2;
    const shape_type last_kn_pos_out        = shape_out_transform_dim/2;
    const shape_type last_kn_pos_transform  = transform_n/2;
    const shape_type last_kn_pos_copyin
        = std::min(last_kn_pos_in, last_kn_pos_transform);
    const shape_type last_kn_pos_copyout
        = std::min(last_kn_pos_out, last_kn_pos_transform);
    const shape_type first_kn_neg_in        = -(shape_in_transform_dim-1)/2;
    const shape_type first_kn_neg_out       = -(shape_out_transform_dim-1)/2;
    const shape_type first_kn_neg_transform = -(transform_n-1)/2;
    // Normalization only required during backwards transform
    const double normalization_factor
        = (fftw_sign == FFTW_BACKWARD) ? 1.0/transform_n : 1.0;

    // Prepare per-pencil outer loop index and loop bounds
    shape_array loop_shape(shape_in);   // Iterate over all dimensions...
    loop_shape[transform_dim] = 1;      // ...except the transformed one
    index_array loop_index = {{   }};   // Initialize to default value
    for (size_type n = 0; n < dimensionality; ++n) {
        assert(loop_index[n] == 0);     // Check initialization correct
    }
    index_array dereference_index;      // To be adjusted by index_bases
    fftw_complex * p_buffer;            // Used to walk the buffer
    shape_type kn;                      // Used to track current wavenumber

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
        kn = 0;
        for (/* from above */; kn <= last_kn_pos_copyin; ++kn) {
            detail::assign_complex(*(p_buffer++), *p_pencil_in);
            p_pencil_in += stride_in_transform_dim;
        }
        if (shape_in_transform_dim > transform_n) {
            for (/* from above */; kn <= last_kn_pos_in; ++kn) {
                p_pencil_in += stride_in_transform_dim;
            }
            for (kn = first_kn_neg_in; kn < first_kn_neg_transform; ++kn) {
                p_pencil_in += stride_in_transform_dim;
            }
        } else {
            for (/* from above */; kn <= last_kn_pos_transform; ++kn) {
                detail::assign_complex(*(p_buffer++), 0, 0);
            }
            for (kn = first_kn_neg_transform; kn < first_kn_neg_in; ++kn) {
                detail::assign_complex(*(p_buffer++), 0, 0);
            }
        }
        for (/* from above */; kn <= -1; ++kn) {
            detail::assign_complex(*(p_buffer++), *p_pencil_in);
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
        kn = 0;
        for (/* from above */; kn <= last_kn_pos_copyout; ++kn) {
            detail::assign_complex(*p_pencil_out,
                        (*p_buffer)[0] * normalization_factor,
                        (*p_buffer)[1] * normalization_factor);
            ++p_buffer;
            p_pencil_out += stride_out_transform_dim;
        }
        if (shape_out_transform_dim > transform_n) {
            for (/* from above */; kn <= last_kn_pos_out; ++kn) {
                detail::assign_complex(*p_pencil_out, 0, 0);
                p_pencil_out += stride_out_transform_dim;
            }
            for (kn = first_kn_neg_out; kn < first_kn_neg_transform; ++kn) {
                detail::assign_complex(*p_pencil_out, 0, 0);
                p_pencil_out += stride_out_transform_dim;
            }
        } else {
            for (/* from above */; kn <= last_kn_pos_transform; ++kn) {
                ++p_buffer;
            }
            for (kn = first_kn_neg_transform; kn < first_kn_neg_out; ++kn) {
                ++p_buffer;
            }
        }
        for (/* from above */; kn <= -1; ++kn) {
            detail::assign_complex(*p_pencil_out,
                        (*p_buffer)[0] * normalization_factor,
                        (*p_buffer)[1] * normalization_factor);
            ++p_buffer;
            p_pencil_out += stride_out_transform_dim;
        }

    } while (increment<dimensionality>(loop_index, loop_shape));

} /* c2c_transform */

} /* fftw_multi_array */ } /* suzerain */ } /* pecos */

#endif // PECOS_SUZERAIN_FFTW_MULTI_ARRAY_HPP
