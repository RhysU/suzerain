/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
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
 * pencil.hpp: Class to support P3DFFT's physical and wave space data layout
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */
#ifndef __SUZERAIN_PENCIL_H
#define __SUZERAIN_PENCIL_H

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/pencil_grid.hpp>

namespace suzerain
{

/** Encapsulates a pencil-storage-based scalar field to be in-place transformed
 * between physical and from wave space by P3DFFT.  Handles details like
 * maintaining size information, indexing, appropriate storage order.  Here \e
 * physical indicates real-valued storage of pointwise values and \e wave
 * indicates complex-valued storage of Fourier coefficients.
 *
 * The \f$ x \f$, \f$ y \f$, and \f$ z \f$ directions are designed to be the
 * streamwise, wall-normal, and spanwise directions respectively.  The storage
 * is arranged so that the wall-normal direction is stride one in wave space
 * while the streamwise direction is stride one in physical space.
 */
template<
    typename FPT       = double,
    typename Allocator = typename suzerain::blas::allocator<FPT>::type
    >
class pencil
    : public boost::noncopyable,
      private Allocator
{
public:

    /**
     * @name Index and index offset types
     * Type choices mimic those found in <tt>boost::multi_array_types</tt>.
     * @see <a href="http://www.boost.org/doc/libs/release/libs/multi_array">
     *      Boost.MultiArray</a> for more information on
     *      <tt>boost::multi_array_types</tt>.
     * @{ */
    /**
     * This is a signed integral type used for indexing into pencils.
     * It is also used to represent strides and index bases.
     */
    typedef boost::multi_array_types::index index;

    /**
     * This is an unsigned integral type.
     * It is primarily used to specify pencil shapes.
     */
    typedef boost::multi_array_types::size_type size_type;

    /**
     * This is a signed integral type used to represent the distance
     * between two iterators.  It is the same type as
     * <tt>std::iterator_traits<iterator>::difference_type</tt>.
     **/
    typedef boost::multi_array_types::difference_type difference_type;
    /**  @} */


    /**
     * @name Real-valued types
     * Used for physical space values and the real and imaginary portions
     * of wave space coefficients.
     * @{ */
    typedef FPT                real_type;
    typedef real_type&         real_reference;
    typedef const real_type&   const_real_reference;
    typedef real_type*         real_pointer;
    typedef const real_type*   const_real_pointer;
    typedef real_pointer       real_iterator;
    typedef const_real_pointer const_real_iterator;
    /**  @} */

    /**
     * @name Complex-valued types
     * Used for wave space values.
     * @{ */
    typedef typename std::complex<FPT> complex_type;
    typedef complex_type&              complex_reference;
    typedef const complex_type&        const_complex_reference;
    typedef complex_type*              complex_pointer;
    typedef const complex_type*        const_complex_pointer;
    typedef complex_pointer            complex_iterator;
    typedef const_complex_pointer      const_complex_iterator;
    /**  @} */

private: // Declared above public members to enforce initialization order

    /** Total amount of real_type data stored within the pencil, including both
     *  physical and wave space storage requirements.
     */
    const size_type data_nelem_;

    /** Sufficient in size to simultaneously house physical and wave space
     * data.  Contained physical and wave instances store information within
     * this array.
     */
    real_pointer const data_;

public:

    class wave_space; // Forward declaration

    /**
     * Provides access to the real-valued representation of the field in
     * physical space, assuming P3DFFT's \c p3dfft_btran_c2r has been applied
     * to the pencil.
     *
     * The underlying storage is column major in (X,Z,Y) index order.  For
     * stride reasons, three loops iterating across physical_space should
     * resemble
     * \code
     *  // p an instance of pencil<>::physical_space
     *  for (pencil<>::index j = 0; j < p.size_y; ++j)
     *      for (pencil<>::index k = 0; k < p.size_z; ++k)
     *          for (pencil<>::index i = 0; i < p.size_x; ++i)
     *              // Access p(i,j,k) here
     * \endcode
     */
    //TODO Add comments about half storage of last dimension
    class physical_space : boost::noncopyable
    {
    public:
        /**
         *  @name Starting offsets within the global pencil_grid
         *  Inclusive index.
         * @{ */
        const index start_x; /**< Starting streamwise offset */
        const index start_y; /**< Starting wall-normal offset */
        const index start_z; /**< Starting spanwise offset */
        /**  @} */

        /**
         *  @name Ending offsets within the global pencil_grid
         *  Exclusive index.
         * @{ */
        const index end_x; /**< Ending streamwise offset */
        const index end_y; /**< Ending wall-normal offset */
        const index end_z; /**< Ending spanwise offset */
        /**  @} */

        /**
         * @name Size of the physical_space data within the pencil.
         * @{ */
        const size_type size_x; /**< Size in streamwise direction */
        const size_type size_y; /**< Size in wall-normal direction */
        const size_type size_z; /**< Size in spanwise direction */
        const size_type size;   /**< Products of sizes in three dimensions */
        /**  @} */

        /**
         * @name Index-based access to physical space data
         * @{ */
        /** Mutable access to physical space data at given offset */
        real_reference operator()(
            const index x, const index y, const index z);

        /** Immutable access to physical space data at given offset */
        const_real_reference operator()(
            const index x, const index y, const index z) const;
        /**  @} */

        /**
         * @name Iterator-based access to physical space data
         * Iteration access is linear across the underlying storage.
         * @{ */
        real_iterator       begin();
        const_real_iterator begin() const;
        real_iterator       end();
        const_real_iterator end() const;
        /**  @} */

        /**
         * @name Offset information methods for physical space data
         * @{ */
        /** Compute the linear offset to the (\c x, \c y, \c z) element within
         * the physical_space data.  This is the "original" orientation per the
         * P3DFFT manual page 4.  The layout is column major (Fortran) storage
         * in (X,Z,Y) order.
         *
         * @param x desired entry offset in streamwise direction
         * @param y desired entry offset in wall-normal direction
         * @param z desired entry offset in spanwise direction
         *
         * @return the linear, 1D offset where (\c x, \c y, \c z) is stored.
         */
        index offset(
            const index x,
            const index y,
            const index z) const;

        /** Compute the (\c x, \c y, \c z) physical space indices associated
         * with offset \c i.  Useful for turning an offset into a tuple that
         * can be logged or displayed.  Returned indices are for the local
         * pencil data, not the global grid.
         *
         * @param[in] i offset found per \c offset method
         * @param[out] x index in the streamwise direction
         * @param[out] y index in the wall-normal direction
         * @param[out] z index in the spanwise direction
         */
        void inverse_offset(
            const index  i,
            index &x,
            index &y,
            index &z) const;

        /** Compute the (\c x, \c y, \c z) global physical space indices
         * associated with offset \c i.  Useful for turning an offset into a
         * tuple that can be logged or displayed.  Returned indices are for the
         * global grid.
         *
         * @param[in]  i offset found per \c offset method
         * @param[out] x index in the streamwise direction
         * @param[out] y index in the wall-normal direction
         * @param[out] z index in the spanwise direction
         */
        void inverse_global_offset(
            const index  i,
            index &x,
            index &y,
            index &z) const;
         /* @} */

    private:
        /** Allows pencil to construct instances of physical_space */
        friend class pencil;

        /** Intended for use by pencil only.  The containing pencil
         * instance is required to perform all memory allocation and
         * deallocation.
         *
         * @param start starting location within the global pencil grid
         * @param size  sizes of the data stored within this instance
         * @param data  location where coefficients are to be found,
         *              must be sufficiently large to hold all elements.
         */
        template<typename RandomAccessContainer>
        physical_space(
            const RandomAccessContainer &start,
            const RandomAccessContainer &size,
            real_pointer data);

        /** Raw real_type data where coefficients are stored. */
        real_pointer const data_;

        /** Precomputed size_x * size_z */
        const index size_xz_;
    };

    /**
     * Provides access to the complex-valued representation of the field in
     * wave space, assuming P3DFFT's \c p3dfft_ftran_r2c has been applied
     * to the pencil.
     *
     * The underlying storage is column major in (Y,X,Z) index order.  For
     * stride reasons, three loops iterating across physical_space should
     * resemble
     * \code
     *  // p an instance of pencil::wave_space
     *  for (pencil<>::index k = 0; k < p.size_z; ++k)
     *      for (pencil<>::index i = 0; i < p.size_x; ++i)
     *          for (pencil<>::index j = 0; j < p.size_y; ++j)
     *              // Access p(i,j,k) here
     * \endcode
     */
    class wave_space : boost::noncopyable
    {
    public:
        /**
         *  @name Starting offsets within the global pencil_grid
         *  Inclusive index.
         * @{ */
        const index start_x; /**< Starting streamwise offset */
        const index start_y; /**< Starting wall-normal offset */
        const index start_z; /**< Starting spanwise offset */
        /**  @} */

        /**
         *  @name Ending offsets within the global pencil_grid
         *  Exclusive index.
         * @{ */
        const index end_x; /**< Ending streamwise offset */
        const index end_y; /**< Ending wall-normal offset */
        const index end_z; /**< Ending spanwise offset */
        /**  @} */

        /**
         * @name Size of the wave_space data within the pencil.
         * @{ */
        const size_type size_x; /**< Size in streamwise direction */
        const size_type size_y; /**< Size in wall-normal direction */
        const size_type size_z; /**< Size in spanwise direction */
        const size_type size;   /**< Products of sizes in three dimensions */
        /**  @} */

        /**
         * @name Index-based access to wave space data
         * @{ */
        /** Mutable access to wave space data at given offset */
        complex_reference operator()(
            const index x, const index y, const index z);

        /** Immutable access to wave space data at given offset */
        const_complex_reference operator()(
            const index x, const index y, const index z) const;

        /** Mutable access to real coefficients at given offset */
        real_reference real(
            const index x, const index y, const index z);

        /** Immutable access to space real coefficients at given offset */
        const_real_reference real(
            const index x, const index y, const index z) const;

        /** Mutable access to imaginary coefficients at given offset */
        real_reference imag(
            const index x, const index y, const index z);

        /** Immutable access to imaginary coefficients at given offset */
        const_real_reference imag(
            const index x, const index y, const index z) const;
        /**  @} */

        /**
         * @name Iterator-based access to wave space data
         * Iteration access is linear across the underlying storage.
         * @{ */
        complex_iterator       begin();
        const_complex_iterator begin() const;
        complex_iterator       end();
        const_complex_iterator end() const;
        /**  @} */

        /**
         * @name Offset information methods for wave space data
         * @{ */
        /** Compute the linear offset to the (\c x, \c y, \c z) element within
         * the wave_space data.  This is the "transposed" orientation per the
         * P3DFFT manual page 4.  The layout is column major (Fortran) storage
         * in (Y,X,Z) order.
         *
         * @param x desired entry offset in streamwise direction
         * @param y desired entry offset in wall-normal direction
         * @param z desired entry offset in spanwise direction
         *
         * @return the linear, 1D offset where (\c x, \c y, \c z) is stored.
         */
        index offset(
            const index x,
            const index y,
            const index z) const;

        /** Compute the (\c x, \c y, \c z) wave space indices associated with
         * offset \c i.  Useful for turning an offset into a tuple that can be
         * logged or displayed.  Returned indices are for the local pencil
         * data, not the global grid.
         *
         * @param[in]  i offset found per \c offset method
         * @param[out] x index in the streamwise direction
         * @param[out] y index in the wall-normal direction
         * @param[out] z index in the spanwise direction
         */
        void inverse_offset(
            const index  i,
            index &x,
            index &y,
            index &z) const;

        /** Compute the (\c x, \c y, \c z) global wave space indices associated
         * with offset \c i.  Useful for turning an offset into a tuple that
         * can be logged or displayed.  Returned indices are for the global
         * grid.
         *
         * @param[in]  i offset found per \c offset method
         * @param[out] x index in the streamwise direction
         * @param[out] y index in the wall-normal direction
         * @param[out] z index in the spanwise direction
         */
        void inverse_global_offset(
            const index  i,
            index &x,
            index &y,
            index &z) const;
        /**  @} */

    private:
        /** Allows pencil to construct instances of physical_space */
        friend class pencil;

        /** Intended for use by pencil only.  The containing pencil
         * instance is required to perform all memory allocation and
         * deallocation.
         *
         * @param start starting location within the global pencil grid
         * @param size  sizes of the data stored within this instance
         * @param data  location where coefficients are to be found,
         *              must be sufficiently large to hold all elements.
         */
        template<typename RandomAccessContainer>
        wave_space(
            const RandomAccessContainer &start,
            const RandomAccessContainer &size,
            real_pointer data);

        /** Raw complex_type data where coefficients are stored. */
        complex_pointer const data_complex_;

        /** Raw real_type data where coefficients are stored. */
        real_pointer    const data_real_;

        /** Precomputed size_x * size_y */
        const index size_xy_;
    };

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
     * @throw std::invalid_argument if any index is negative.
     */
    template<typename RandomAccessContainer>
    pencil(const RandomAccessContainer &pstart,
           const RandomAccessContainer &psize,
           const RandomAccessContainer &wstart,
           const RandomAccessContainer &wsize)
    throw(std::invalid_argument);

    /**
     * Construct a local pencil matching the characteristics of a global
     * pencil_grid.  The pencil is filled with zeros at construction time.
     *
     * @param pg The pencil_grid to match.
     */
    pencil(const pencil_grid &pg);

    /** Virtual destructor */
    virtual ~pencil();

    /**
     * @name Iterator-based access to all physical and wave space data.
     * Iteration access is linear across the underlying storage.
     * @{ */
    real_iterator       begin();
    const_real_iterator begin() const;
    real_iterator       end();
    const_real_iterator end() const;
    /**  @} */


    /** Access the raw underlying storage directly, intended for use with
     * P3DFFT's \c p3dfft_ftran_r2c and \c p3dfft_btran_c2r methods.
     */
    real_pointer data();

    /** Use to access all physical_space data for this instance */
    physical_space physical;

    /** Use to access all wave_space data for this instance */
    wave_space wave;
};

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_iterator
pencil<FPT,Allocator>::begin()
{
    return data_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_real_iterator
pencil<FPT,Allocator>::begin() const
{
    return data_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_iterator
pencil<FPT,Allocator>::end()
{
    return data_ + data_nelem_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_real_iterator
pencil<FPT,Allocator>::end() const
{
    return data_ + data_nelem_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_pointer
pencil<FPT,Allocator>::data()
{
    return data_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::index
pencil<FPT,Allocator>::physical_space::offset(
    const index x,
    const index y,
    const index z) const
{
    return x + z*size_x + y*size_xz_;
}

template<typename FPT, typename Allocator>
inline
void
pencil<FPT,Allocator>::physical_space::inverse_offset(
    const index  i,
    index &x,
    index &y,
    index &z) const
{
    y = i / (size_xz_);
    z = i / size_x - y*size_z;
    x = i - y*size_xz_ - z*size_x;
}

template<typename FPT, typename Allocator>
inline
void
pencil<FPT,Allocator>::physical_space::inverse_global_offset(
    const index  i,
    index &x,
    index &y,
    index &z) const
{
    inverse_offset(i, x, y, z);
    x += start_x;
    y += start_y;
    z += start_z;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::index
pencil<FPT,Allocator>::wave_space::offset(
    const index x,
    const index y,
    const index z) const
{
    // TODO Assert STRIDE1 specified during P3DFFT compilation
    return y + x*size_y + z*size_xy_;
}

template<typename FPT, typename Allocator>
inline
void
pencil<FPT,Allocator>::wave_space::inverse_offset(
    const index  i,
    index &x,
    index &y,
    index &z) const
{
    z = i / (size_xy_);
    x = i / size_y - z*size_x;
    y = i - z*size_xy_ - x*size_y;
}

template<typename FPT, typename Allocator>
inline
void
pencil<FPT,Allocator>::wave_space::inverse_global_offset(
    const index  i,
    index &x,
    index &y,
    index &z) const
{
    inverse_offset(i, x, y, z);
    x += start_x;
    y += start_y;
    z += start_z;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_reference
pencil<FPT,Allocator>::physical_space::operator()(
    const index x, const index y, const index z)
{
    return data_[offset(x, y, z)];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_real_reference
pencil<FPT,Allocator>::physical_space::operator()(
    const index x, const index y, const index z) const
{
    return data_[offset(x, y, z)];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::complex_reference
pencil<FPT,Allocator>::wave_space::operator()(
    const index x, const index y, const index z)
{
    return data_complex_[offset(x, y, z)];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_complex_reference
pencil<FPT,Allocator>::wave_space::operator()(
    const index x, const index y, const index z) const
{
    return data_complex_[offset(x, y, z)];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_reference
pencil<FPT,Allocator>::wave_space::real(
    const index x, const index y, const index z)
{
    return data_real_[2*offset(x, y, z)];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_real_reference
pencil<FPT,Allocator>::wave_space::real(
    const index x, const index y, const index z) const
{
    return data_real_[2*offset(x, y, z)];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_reference
pencil<FPT,Allocator>::wave_space::imag(
    const index x, const index y, const index z)
{
    return data_real_[2*offset(x, y, z) + 1];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_real_reference
pencil<FPT,Allocator>::wave_space::imag(
    const index x, const index y, const index z) const
{
    return data_real_[2*offset(x, y, z) + 1];
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_iterator
pencil<FPT,Allocator>::physical_space::begin()
{
    return data_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_real_iterator
pencil<FPT,Allocator>::physical_space::begin() const
{
    return data_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::real_iterator
pencil<FPT,Allocator>::physical_space::end()
{
    return data_ + size;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_real_iterator
pencil<FPT,Allocator>::physical_space::end() const
{
    return data_ + size;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::complex_iterator
pencil<FPT,Allocator>::wave_space::begin()
{
    return data_complex_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_complex_iterator
pencil<FPT,Allocator>::wave_space::begin() const
{
    return data_complex_;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::complex_iterator
pencil<FPT,Allocator>::wave_space::end()
{
    return data_complex_ + size;
}

template<typename FPT, typename Allocator>
inline
typename pencil<FPT,Allocator>::const_complex_iterator
pencil<FPT,Allocator>::wave_space::end() const
{
    return data_complex_ + size;
}

template<typename FPT, typename Allocator>
template< typename RandomAccessContainer >
pencil<FPT,Allocator>::pencil(
    const RandomAccessContainer &pstart, const RandomAccessContainer &psize,
    const RandomAccessContainer &wstart, const RandomAccessContainer &wsize)
throw(std::invalid_argument)
        : data_nelem_(std::max(
              psize[0]*psize[1]*psize[2],
              2*wsize[0]*wsize[1]*wsize[2])),
        data_(Allocator::allocate(data_nelem_)),
        physical(pstart, psize, data_),
        wave(wstart, wsize, data_)
{
    std::fill_n(data_, data_nelem_, real_type(0)); // Fill with zeros
}

template<typename FPT, typename Allocator>
pencil<FPT,Allocator>::pencil(const pencil_grid &pg)
        : data_nelem_(std::max(
                 pg.local_physical_extent[0]
                    *pg.local_physical_extent[1]
                    *pg.local_physical_extent[2],
                 2*pg.local_wave_extent[0]
                    *pg.local_wave_extent[1]
                    *pg.local_wave_extent[2])),
        data_(Allocator::allocate(data_nelem_)),
        physical(pg.local_physical_start,
                 pg.local_physical_extent, data_),
        wave(pg.local_wave_start,
             pg.local_wave_extent, data_)
{
    std::fill_n(data_, data_nelem_, real_type(0)); // Fill with zeros
}

template<typename FPT, typename Allocator>
pencil<FPT,Allocator>::~pencil()
{
    Allocator::deallocate(data_, data_nelem_);
}

template<typename FPT, typename Allocator>
template<typename RandomAccessContainer>
pencil<FPT,Allocator>::physical_space::physical_space(
    const RandomAccessContainer &start,
    const RandomAccessContainer &size,
    real_pointer data)
    :
    start_x(start[0]), start_y(start[1]), start_z(start[2]),
    end_x(start[0]+size[0]), end_y(start[1]+size[1]), end_z(start[2]+size[2]),
    size_x(size[0]), size_y(size[1]), size_z(size[2]),
    size(size_x*size_y*size_z),
    data_(data),
    size_xz_(size_x*size_z)
{
    // NOP
}

template<typename FPT, typename Allocator>
template<typename RandomAccessContainer>
pencil<FPT,Allocator>::wave_space::wave_space(
    const RandomAccessContainer &start,
    const RandomAccessContainer &size,
    real_pointer data)
    :
    start_x(start[0]), start_y(start[1]), start_z(start[2]),
    end_x(start[0]+size[0]), end_y(start[1]+size[1]), end_z(start[2]+size[2]),
    size_x(size[0]), size_y(size[1]), size_z(size[2]),
    size(size_x*size_y*size_z),
    data_complex_(reinterpret_cast<complex_pointer>(data)),
    data_real_(data),
    size_xy_(size_x*size_y)
{
    // NOP
}

} // namespace suzerain

#endif // __SUZERAIN_PENCIL_H
