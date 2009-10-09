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

#include <cassert>
#include <cstddef>
#include <vector>
#include <fftw3.h>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_unsigned.hpp>
#include <boost/typeof/typeof.hpp>

/* DEBUG */
#include <iostream>
#include <iterator>

namespace pecos { namespace suzerain { namespace fftw_multi_array {

template<std::size_t NumDims, typename IndexType, typename MaxIndexType>
bool increment(IndexType &index, const MaxIndexType &max_index)
{
    typedef BOOST_TYPEOF_TPL(index[0])          index_element_type;
    typedef BOOST_TYPEOF_TPL(max_index[0])      max_index_element_type;
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
        assert(index[n] < std::numeric_limits<index_element_type>::max() - 1);

        index[n] += overflow;                 // Handle incoming overflow
        overflow  = index[n]/max_index[n];    // Check outgoing overflow
        index[n] *= !overflow;                // Set to zero on outgoing
    }
    return !overflow;
}

template<class ComplexMultiArray>
void c2c_transform(const size_t transform_dim,
                   ComplexMultiArray &in,
                   ComplexMultiArray &out,
                   const int fftw_sign,
                   const unsigned fftw_flags,
                   const double dealias_by = 1.0)
{
    typedef typename ComplexMultiArray::element                   element;
    typedef typename ComplexMultiArray::index                     index;
    typedef typename ComplexMultiArray::size_type                 size_type;
    typedef typename ComplexMultiArray::index_range               index_range;
    typedef boost::array<index,ComplexMultiArray::dimensionality> index_array;
    const size_type dimensionality = ComplexMultiArray::dimensionality;

    /* Ensure we transform a dimension that exists in the data */
    assert(transform_dim < dimensionality);
    /* Ensure we are operating on a complex-valued array */
    BOOST_STATIC_ASSERT(
              (boost::is_complex<element>::value)
           || (boost::is_same<element, fftw_complex>::value));
    /* Ensure transformation direction and dealiasing choice are sane */
    assert(   (fftw_sign == FFTW_FORWARD  && dealias_by == 1.0)
           || (fftw_sign == FFTW_BACKWARD && dealias_by >= 1.0));
    /* Ensure the in and out arrays have the same shape */
    const size_type * const shape = in.shape();
    for (size_type n = 0; n < dimensionality; ++n) {
        assert(shape[n] == out.shape()[n]);
    }

    /* Determine FFTW advanced plan parameters */
    const int rank = 1;
    const int n[1] = { shape[transform_dim] };
    assert( shape[transform_dim] <= static_cast<size_type>(n[0])); /* Sanity */
//    const int howmany;

    index_array out_strides, in_strides;
    index_array out_indices, in_indices;
    index_array in_last_indices;
    for (size_type n = 0; n < dimensionality; ++n) {
        out_strides[n] = out.strides()[n];
        in_strides[n]  = in.strides()[n];
        out_indices[n] = out.index_bases()[n];
        in_indices[n] = in.index_bases()[n];
        in_last_indices[n] = in_indices[n] + shape[n];
    }

//    std::copy(in_indices.begin(), in_indices.end(),
//            std::ostream_iterator<index>(std::cout, ",")); /* DEBUG */

//    std::copy(in_indices.begin(), in_indices.end(),
//            std::ostream_iterator<index>(std::cout, ",")); /* DEBUG */




//    int howmany,
//    fftw_complex *in,
//    const int *inembed,
//    int istride,
//    int idist,
//    fftw_complex *out,
//    const int *onembed,
//    int ostride,
//    int odist,
//    int sign,
//    unsigned flags

    /* Construct the FFTW transform plan and plan parameters */
//    boost::shared_ptr<boost::remove_pointer<fftw_plan>::type> plan(
//        fftw_plan_dft_c2r_1d(
//                        NR,
//                        reinterpret_cast<fftw_complex*>(c.data()),
//                        r.data(),
//                        FFTW_ESTIMATE),
//                    std::ptr_fun(fftw_destroy_plan))

//         fftw_plan_many_dft(int rank,
//                            const int *n,
//                            int howmany,
//                            fftw_complex *in,
//                            const int *inembed,
//                            int istride,
//                            int idist,
//                            fftw_complex *out,
//                            const int *onembed,
//                            int ostride,
//                            int odist,
//                            int sign,
//                            unsigned flags),
//        std::ptr_fun(fftw_destroy_plan));

//    const fftw_plan plan = fftw_plan_dft(1, )
} /* c2c_transform */

} /* fftw_multi_array */ } /* suzerain */ } /* pecos */

#endif // PECOS_SUZERAIN_FFTW_MULTI_ARRAY_HPP
