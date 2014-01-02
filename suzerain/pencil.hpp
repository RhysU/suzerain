//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
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

#ifndef SUZERAIN_PENCIL_HPP
#define SUZERAIN_PENCIL_HPP

/** @file
 * Addressing support for pencil_grid's physical and wave space data.
 */

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/contiguous_memory.hpp>
#include <suzerain/multi_array.hpp>
#include <suzerain/pencil_grid.hpp>
#include <suzerain/storage.hpp>

namespace suzerain {

/**
 * Encapsulates a pencil-storage-based scalar field to be in-place transformed
 * between physical and from wave space by \ref pencil_grid.  Handles details
 * like maintaining size information, indexing, appropriate storage order.
 * Here \e physical indicates real-valued storage of pointwise values and \e
 * wave indicates complex-valued storage of Fourier coefficients.
 *
 * The \f$ x \f$, \f$ y \f$, and \f$ z \f$ directions are designed to be the
 * streamwise, wall-normal, and spanwise directions respectively.  The storage
 * is arranged so that the wall-normal direction is stride one in wave space
 * while the streamwise direction is stride one in physical space.
 */
template<
    typename FPT       = real_t,
    typename Allocator = typename blas::allocator<FPT>::type
    >
class pencil
    : public  boost::noncopyable,
      private contiguous_memory<FPT,Allocator>
{
public:

    // Index- and size-related types
    typedef boost::multi_array_types::index           index;
    typedef boost::multi_array_types::size_type       size_type;
    typedef boost::multi_array_types::difference_type difference_type;

    // Types specifying storage ordering
    typedef typename storage::general<
        boost::mpl::vector_c<std::size_t,0,2,1> > physical_storage_order_type;
    typedef typename storage::general<
        boost::mpl::vector_c<std::size_t,1,0,2> > wave_storage_order_type;

    // Types for contained data
    typedef FPT real_type;
    typedef typename std::complex<FPT> complex_type;
    typedef typename multi_array::ref<real_type,3>    physical_type;
    typedef typename multi_array::ref<complex_type,3> wave_type;

    /**
     * Construct a scalar pencil with the given characteristics.
     * The pencil is filled with zeros at construction time.
     *
     * @warning Note that unlike P3DFFT's \c get_dims return values, \c pstart
     * and \c wstart are to be zero-indexed.
     *
     * @param pstart zero-indexed starting location in physical space within
     *               global grid.
     * @param psize  size of the pencil in physical space.
     * @param wstart zero-indexed starting location in wave space within the
     *               global grid.
     * @param wsize  size of the pencil in wave space.
     *
     * @see <a href="http://www.sgi.com/tech/stl/RandomAccessContainer.html">
     *      Random Access Container</a> for details on the contract expected
     *      of the arguments.
     */
    template<typename RandomAccessContainer1,
             typename RandomAccessContainer2,
             typename RandomAccessContainer3,
             typename RandomAccessContainer4>
    pencil(const RandomAccessContainer1& pstart,
           const RandomAccessContainer2& psize,
           const RandomAccessContainer3& wstart,
           const RandomAccessContainer4& wsize);

    /**
     * Construct a local pencil matching the characteristics of a global
     * pencil_grid.  The pencil is filled with zeros at construction time.
     *
     * @param pg The pencil_grid to match.
     */
    pencil(const pencil_grid& pg);

    /**
     * A real-valued, local physical-space view of the storage stored row-major
     * XZY with X being the fastest direction.
     */
    physical_type physical;

    /**
     * A real-valued, global physical-space view of the storage stored
     * row-major XZY with X being the fastest direction.  The
     * <tt>index_bases()</tt> are set such that iterating between
     * <tt>index_bases(</tt>) and <tt>index_bases() + shape()</tt> will provide
     * global indexes.
     */
    physical_type global_physical;

    /**
     * A complex-valued, local wave-space view of the storage stored row-major
     * YXZ with Y being the fastest direction.
     */
    wave_type wave;

    /**
     * A wave-valued, global wave-space view of the storage stored row-major
     * YXZ with Y being the fastest direction.  The <tt>index_bases()</tt> are
     * set such that iterating between <tt>index_bases(</tt>) and
     * <tt>index_bases() + shape()</tt> will provide global indexes.
     */
    wave_type global_wave;

    /**
     * @name Iterator-based access to all physical and wave-space data.
     * Iteration access is linear across the underlying storage.
     * @{ */

    /**
     * Retrieve an iterator to the beginning of this pencil's storage.
     */
    typename Allocator::pointer begin() {
        return this->memory_begin();
    }

    /**
     * Retrieve a constant iterator to the beginning of this pencil's storage.
     */
    typename Allocator::const_pointer begin() const {
        return this->memory_begin();
    }

    /**
     * Retrieve an iterator to the end of this pencil's storage.
     */
    typename Allocator::pointer end() {
        return this->memory_end();
    }

    /**
     * Retrieve a constant iterator to the end of this pencil's storage.
     */
    typename Allocator::const_pointer end() const {
        return this->memory_end();
    }
    /**  @} */

private:

    template<typename ForwardIterator>
    static
    array<typename std::iterator_traits<ForwardIterator>::value_type,3>
    make_collection(ForwardIterator iter)
    {
        typedef typename std::iterator_traits<ForwardIterator>::value_type T;
        array<T,3> retval;
        retval[0] = *iter;
        retval[1] = *++iter;
        retval[2] = *++iter;
        return retval;
    }
};

template<typename FPT, typename Allocator>
template<typename RandomAccessContainer1,
         typename RandomAccessContainer2,
         typename RandomAccessContainer3,
         typename RandomAccessContainer4 >
pencil<FPT,Allocator>::pencil(const RandomAccessContainer1& pstart,
                              const RandomAccessContainer2& psize,
                              const RandomAccessContainer3& wstart,
                              const RandomAccessContainer4& wsize)
    : contiguous_memory<FPT,Allocator>(std::max<typename Allocator::size_type>(
              psize[0]*psize[1]*psize[2],
                (   sizeof(typename wave_type::element)
                  / sizeof(typename physical_type::element))
              * wsize[0]*wsize[1]*wsize[2])),
      physical(
              this->memory_begin(),
              psize,
              physical_storage_order_type()),
      global_physical(
              this->memory_begin(),
              typename physical_type::extent_gen()
                    [typename physical_type::extent_range(
                        pstart[0], pstart[0] + psize[0])]
                    [typename physical_type::extent_range(
                        pstart[1], pstart[1] + psize[1])]
                    [typename physical_type::extent_range(
                        pstart[2], pstart[2] + psize[2])],
              physical_storage_order_type()),
      wave(
              reinterpret_cast<typename wave_type::element *>(
                  this->memory_begin()),
              wsize,
              wave_storage_order_type()),
      global_wave(
              reinterpret_cast<typename wave_type::element *>(
                  this->memory_begin()),
              typename wave_type::extent_gen()
                    [typename wave_type::extent_range(
                        wstart[0], wstart[0] + wsize[0])]
                    [typename wave_type::extent_range(
                        wstart[1], wstart[1] + wsize[1])]
                    [typename wave_type::extent_range(
                        wstart[2], wstart[2] + wsize[2])],
              wave_storage_order_type())
{
    std::fill(this->memory_begin(), this->memory_end(), 0);  // Fill with zeros
}

template<typename FPT, typename Allocator>
pencil<FPT,Allocator>::pencil(const pencil_grid& pg)
    : contiguous_memory<FPT,Allocator>(std::max<typename Allocator::size_type>(
                pg.local_physical_extent.prod(),
                  (   sizeof(typename wave_type::element)
                    / sizeof(typename physical_type::element))
                * pg.local_wave_extent.prod())),
      physical(
              this->memory_begin(),
              make_collection(pg.local_physical_extent.data()),
              physical_storage_order_type()),
      global_physical(
              this->memory_begin(),
              typename physical_type::extent_gen()
                    [typename physical_type::extent_range(
                            pg.local_physical_start[0],
                            pg.local_physical_end  [0])]
                    [typename physical_type::extent_range(
                            pg.local_physical_start[1],
                            pg.local_physical_end  [1])]
                    [typename physical_type::extent_range(
                            pg.local_physical_start[2],
                            pg.local_physical_end  [2])],
              physical_storage_order_type()),
      wave(
              reinterpret_cast<typename wave_type::element *>(
                  this->memory_begin()),
              make_collection(pg.local_wave_extent.data()),
              wave_storage_order_type()),
      global_wave(
              reinterpret_cast<typename wave_type::element *>(
                  this->memory_begin()),
              typename wave_type::extent_gen()
                    [typename wave_type::extent_range(
                            pg.local_wave_start[0],
                            pg.local_wave_end  [0])]
                    [typename wave_type::extent_range(
                            pg.local_wave_start[1],
                            pg.local_wave_end  [1])]
                    [typename wave_type::extent_range(
                            pg.local_wave_start[2],
                            pg.local_wave_end  [2])],
              wave_storage_order_type())
{
    std::fill(this->memory_begin(), this->memory_end(), 0);  // Fill with zeros
}

} // namespace suzerain

#endif // SUZERAIN_PENCIL_HPP
