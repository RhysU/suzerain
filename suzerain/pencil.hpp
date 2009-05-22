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
#ifndef PECOS_SUZERAIN_PENCIL_H
#define PECOS_SUZERAIN_PENCIL_H

#include <algorithm>
#include <boost/noncopyable.hpp>
#include <boost/scoped_array.hpp>
#include <boost/static_assert.hpp>
#include <complex>
#include <cstddef>

#include <suzerain/pencil_grid.hpp>

namespace pecos
{

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
template < typename T = double, typename G = pencil_grid<> >
class pencil : boost::noncopyable
{
public:
    /**
     * @name Real-valued types
     * Used for physical space values and the real and imaginary portions
     * of wave space coefficients.
     * @{ */
    typedef T                  real_type;
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
    typedef std::complex<T>       complex_type;
    typedef complex_type&         complex_reference;
    typedef const complex_type&   const_complex_reference;
    typedef complex_type*         complex_pointer;
    typedef const complex_type*   const_complex_pointer;
    typedef complex_pointer       complex_iterator;
    typedef const_complex_pointer const_complex_iterator;
    /**  @} */

    /** Dimension type used to describe the underlying pencil_grid */
    typedef typename G::dim_type dim_type;

    /** Offset type used to access coefficients within the %pencil */
    typedef std::size_t size_type;

    /** Relative offset type used to describe differences between offsets */
    typedef std::ptrdiff_t difference_type;

private:
    // Ensure interface design assumptions valid when instantiated
    BOOST_STATIC_ASSERT(2*sizeof(real_type) == sizeof(complex_type));

    // Declared above public members to enforce correct initialization order
    /** Total amount of real_type data stored within the pencil, including both
     *  physical and wave space storage requirements.
     */
    const size_type                      data_nelem_;

    /** Sufficient in size to simultaneously house physical and wave space
     * data.  Contained physical and wave instances store information within
     * this array.
     */
    const boost::scoped_array<real_type> data_;

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
     *  // p an instance of pencil<T,G>::physical
     *  for (pencil<T,G>::size_type j = 0; j < p.size_y; ++j)
     *      for (pencil<T,G>::size_type k = 0; k < p.size_z; ++k)
     *          for (pencil<T,G>::size_type i = 0; i < p.size_x; ++i)
     *              // Access p(i,j,k) here
     * \endcode
     */
    class physical_space : boost::noncopyable
    {
    public:
        /**
         *  @name Starting offsets within the global pencil_grid
         *  Inclusive index.
         * @{ */
        const dim_type start_x; /**< Starting streamwise offset */
        const dim_type start_y; /**< Starting wall-normal offset */
        const dim_type start_z; /**< Starting spanwise offset */
        /**  @} */

        /**
         *  @name Ending offsets within the global pencil_grid
         *  Exclusive index.
         * @{ */
        const dim_type end_x; /**< Ending streamwise offset */
        const dim_type end_y; /**< Ending wall-normal offset */
        const dim_type end_z; /**< Ending spanwise offset */
        /**  @} */

        /**
         * @name Size of the physical_space data within the pencil.
         * @{ */
        const dim_type size_x; /**< Size in streamwise direction */
        const dim_type size_y; /**< Size in wall-normal direction */
        const dim_type size_z; /**< Size in spanwise direction */
        const dim_type size;   /**< Combined size in three dimensions */
        /**  @} */

        /**
         * @name Index-based access to physical space data
         * @{ */
        /** Mutable access to physical space data at given offset */
        real_reference operator()(
            const size_type x, const size_type y, const size_type z);

        /** Immutable access to physical space data at given offset */
        const_real_reference operator()(
            const size_type x, const size_type y, const size_type z) const;
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

    private:
        friend class pencil<T,G>;

        physical_space(
            const dim_type start[3], const dim_type size[3], real_pointer data)
        throw(domain_error);

        size_type offset(
            const size_type x,
            const size_type y,
            const size_type z) const;

        real_pointer const data_;
    };

    class wave_space : boost::noncopyable
    {
    public:
        const dim_type start_x;
        const dim_type start_y;
        const dim_type start_z;
        const dim_type end_x;
        const dim_type end_y;
        const dim_type end_z;
        const dim_type size_x;
        const dim_type size_y;
        const dim_type size_z;
        const dim_type size;

        complex_reference operator()(
            const size_type x, const size_type y, const size_type z);
        const_complex_reference operator()(
            const size_type x, const size_type y, const size_type z) const;

        real_reference real(
            const size_type x, const size_type y, const size_type z);
        const_real_reference real(
            const size_type x, const size_type y, const size_type z) const;
        real_reference imag(
            const size_type x, const size_type y, const size_type z);
        const_real_reference imag(
            const size_type x, const size_type y, const size_type z) const;

        complex_iterator       begin();
        const_complex_iterator begin() const;
        complex_iterator       end();
        const_complex_iterator end() const;

    private:
        friend class pencil<T,G>;

        wave_space(
            const dim_type start[3],
            const dim_type size[3],
            real_pointer data)
        throw(domain_error);

        size_type offset(
            const size_type x,
            const size_type y,
            const size_type z) const;

        complex_pointer const data_complex_;
        real_pointer    const data_real_;
    };

    pencil(const dim_type pstart[3], const dim_type psize[3],
           const dim_type wstart[3], const dim_type wsize[3])
    throw(domain_error);

    real_iterator       begin();
    const_real_iterator begin() const;
    real_iterator       end();
    const_real_iterator end() const;

    physical_space physical;
    wave_space     wave;
};

template<typename T, typename G>
pencil<T, G>::pencil(
    const dim_type pstart[3], const dim_type psize[3],
    const dim_type wstart[3], const dim_type wsize[3])
throw(domain_error)
        : data_nelem_(std::max(
              psize[0]*psize[1]*psize[2],
            2*wsize[0]*wsize[1]*wsize[2])),
        data_(new real_type[data_nelem_]),
        physical(pstart, psize, data_.get()),
        wave(wstart, wsize, data_.get())
{
    // NOP
}

template<typename T, typename G>
inline
pencil<T, G>::real_iterator
pencil<T, G>::begin()
{
    return data_.get();
}

template<typename T, typename G>
inline
pencil<T, G>::const_real_iterator
pencil<T, G>::begin() const
{
    return data_.get();
}

template<typename T, typename G>
inline
pencil<T, G>::real_iterator
pencil<T, G>::end()
{
    return data_.get() + data_nelem_;
}

template<typename T, typename G>
inline
pencil<T, G>::const_real_iterator
pencil<T, G>::end() const
{
    return data_.get() + data_nelem_;
}

template<typename T, typename G>
pencil<T, G>::physical_space::physical_space(
    const dim_type start[3], const dim_type size[3], real_pointer data)
throw(domain_error)
    :
    start_x(start[0]), start_y(start[1]), start_z(start[2]),
    end_x(start[0]+size[0]), end_y(start[1]+size[1]), end_z(start[2]+size[2]),
    size_x(size[0]), size_y(size[1]), size_z(size[2]),
    size(size_x*size_y*size_z),
    data_(data)
{
    if (start_x < 0) throw domain_error();

    if (start_y < 0) throw domain_error();

    if (start_z < 0) throw domain_error();

    if (size_x  < 0) throw domain_error();

    if (size_y  < 0) throw domain_error();

    if (size_z  < 0) throw domain_error();
}

template<typename T, typename G>
pencil<T, G>::wave_space::wave_space(
    const dim_type start[3],
    const dim_type size[3],
    real_pointer data)
throw(domain_error)
    :
    start_x(start[0]), start_y(start[1]), start_z(start[2]),
    end_x(start[0]+size[0]), end_y(start[1]+size[1]), end_z(start[2]+size[2]),
    size_x(size[0]), size_y(size[1]), size_z(size[2]),
    size(size_x*size_y*size_z),
    data_complex_(reinterpret_cast<complex_pointer>(data)),
    data_real_(data)
{
    if (start_x < 0) throw domain_error();

    if (start_y < 0) throw domain_error();

    if (start_z < 0) throw domain_error();

    if (size_x  < 0) throw domain_error();

    if (size_y  < 0) throw domain_error();

    if (size_z  < 0) throw domain_error();
}


template<typename T, typename G>
inline
pencil<T, G>::size_type
pencil<T, G>::physical_space::offset(
    const size_type x,
    const size_type y,
    const size_type z) const
{
    // "original" orientation per the P3DFFT manual page 4
    // Intended for X streamwise, Z spanwise, and Y wall-normal
    // Column major (Fortran) storage in (X,Z,Y) order
    return x + z*size_x + y*size_x*size_z;
}

template<typename T, typename G>
inline
pencil<T, G>::size_type
pencil<T, G>::wave_space::offset(
    const size_type x,
    const size_type y,
    const size_type z) const
{
    // "transposed" orientation per the P3DFFT manual page 4
    // Intended for X streamwise, Z spanwise, and Y wall-normal
    // Column major (Fortran) storage in (Y,Z,X) order
    // TODO Assert STRIDE1 not specified during P3DFFT compilation
    return y + z*size_y + x*size_y*size_z;
}

template<typename T, typename G>
inline
pencil<T, G>::real_reference
pencil<T, G>::physical_space::operator()(
    const size_type x, const size_type y, const size_type z)
{
    return data_[offset(x, y, z)];
}

template<typename T, typename G>
inline
pencil<T, G>::const_real_reference
pencil<T, G>::physical_space::operator()(
    const size_type x, const size_type y, const size_type z) const
{
    return data_[offset(x, y, z)];
}

template<typename T, typename G>
inline
pencil<T, G>::complex_reference
pencil<T, G>::wave_space::operator()(
    const size_type x, const size_type y, const size_type z)
{
    return data_complex_[offset(x, y, z)];
}

template<typename T, typename G>
inline
pencil<T, G>::const_complex_reference
pencil<T, G>::wave_space::operator()(
    const size_type x, const size_type y, const size_type z) const
{
    return data_complex_[offset(x, y, z)];
}

template<typename T, typename G>
inline
pencil<T, G>::real_reference
pencil<T, G>::wave_space::real(
    const size_type x, const size_type y, const size_type z)
{
    return data_real_[2*offset(x, y, z)];
}

template<typename T, typename G>
inline
pencil<T, G>::const_real_reference
pencil<T, G>::wave_space::real(
    const size_type x, const size_type y, const size_type z) const
{
    return data_real_[2*offset(x, y, z)];
}

template<typename T, typename G>
inline
pencil<T, G>::real_reference
pencil<T, G>::wave_space::imag(
    const size_type x, const size_type y, const size_type z)
{
    return data_real_[2*offset(x, y, z) + 1];
}

template<typename T, typename G>
inline
pencil<T, G>::const_real_reference
pencil<T, G>::wave_space::imag(
    const size_type x, const size_type y, const size_type z) const
{
    return data_real_[2*offset(x, y, z) + 1];
}

template<typename T, typename G>
inline
pencil<T, G>::real_iterator
pencil<T, G>::physical_space::begin()
{
    return data_;
}

template<typename T, typename G>
inline
pencil<T, G>::const_real_iterator
pencil<T, G>::physical_space::begin() const
{
    return data_;
}

template<typename T, typename G>
inline
pencil<T, G>::real_iterator
pencil<T, G>::physical_space::end()
{
    return data_ + size;
}

template<typename T, typename G>
inline
pencil<T, G>::const_real_iterator
pencil<T, G>::physical_space::end() const
{
    return data_ + size;
}

template<typename T, typename G>
inline
pencil<T, G>::complex_iterator
pencil<T, G>::wave_space::begin()
{
    return data_complex_;
}

template<typename T, typename G>
inline
pencil<T, G>::const_complex_iterator
pencil<T, G>::wave_space::begin() const
{
    return data_complex_;
}

template<typename T, typename G>
inline
pencil<T, G>::complex_iterator
pencil<T, G>::wave_space::end()
{
    return data_complex_ + size;
}

template<typename T, typename G>
inline
pencil<T, G>::const_complex_iterator
pencil<T, G>::wave_space::end() const
{
    return data_complex_ + size;
}

} // namespace suzerain

} // namespace pecos

#endif // PECOS_SUZERAIN_PENCIL_H
