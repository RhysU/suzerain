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
#ifndef PECOS_SUZERAIN_PENCIL
#define PECOS_SUZERAIN_PENCIL

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/scoped_array.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>
#include <complex>

#include "pencil_grid.hpp"

namespace pecos
{

namespace suzerain
{

template < typename T = double, typename G = pencil_grid<> >
class pencil : boost::noncopyable
{
private:
    typedef T               real_type;
    typedef std::complex<T> complex_type;

    typedef typename
        boost::numeric::ublas::shallow_array_adaptor<real_type>
        physical_space_adaptor_type;
    typedef typename
        boost::numeric::ublas::vector<T, physical_space_adaptor_type>
        physical_space_vector_type;

    typedef typename
        boost::numeric::ublas::shallow_array_adaptor<complex_type>
        wave_space_adaptor_type;
    typedef typename
        boost::numeric::ublas::vector<complex_type, wave_space_adaptor_type>
        wave_space_vector_type;

    // Ensure interface design assumptions valid when instantiated
    BOOST_STATIC_ASSERT(2*sizeof(real_type) == sizeof(complex_type));
    BOOST_STATIC_ASSERT((boost::is_same <
        typename physical_space_vector_type::size_type,
        typename wave_space_vector_type::size_type >::value));
    BOOST_STATIC_ASSERT((boost::is_same <
        typename physical_space_vector_type::difference_type,
        typename wave_space_vector_type::difference_type >::value));
    BOOST_STATIC_ASSERT((boost::is_same <
        real_type,
        typename physical_space_vector_type::value_type >::value));
    BOOST_STATIC_ASSERT((boost::is_same <
        complex_type,
        typename wave_space_vector_type::value_type >::value));

public:
    typedef typename G::dim_type dim_type;

    // Prior static assertions make these identical with wave_space varieties
    typedef typename
        physical_space_vector_type::size_type       size_type;
    typedef typename
        physical_space_vector_type::difference_type difference_type;

private:
    // Declared above public members to enforce correct initialization order
    const size_type                      array_nelem_;
    const boost::scoped_array<real_type> array_;

public:
    class wave_space; // Forward declaration

    class physical_space : boost::noncopyable
    {
    public:
        typedef typename
            physical_space_vector_type::const_iterator  const_iterator;
        typedef typename
            physical_space_vector_type::const_pointer   const_pointer;
        typedef typename
            physical_space_vector_type::const_reference const_reference;
        typedef typename
            physical_space_vector_type::iterator        iterator;
        typedef typename
            physical_space_vector_type::pointer         pointer;
        typedef typename
            physical_space_vector_type::reference       reference;
        typedef typename
            physical_space_vector_type::value_type      value_type;

        const dim_type start_x;
        const dim_type start_y;
        const dim_type start_z;
        const dim_type size_x;
        const dim_type size_y;
        const dim_type size_z;
        const dim_type size;

        reference operator()(
            const size_type x, const size_type y, const size_type z);
        const_reference operator()(
            const size_type x, const size_type y, const size_type z) const;

        iterator       begin();
        const_iterator begin() const;
        iterator       end();
        const_iterator end() const;

    private:
        friend class pencil<T,G>;

        physical_space(
            const dim_type start[3], const dim_type size[3], pointer data)
        throw(domain_error);

        size_type offset(
            const size_type x,
            const size_type y,
            const size_type z) const;

        physical_space_adaptor_type adaptor_;
        physical_space_vector_type  vector_;
    };

    class wave_space : boost::noncopyable
    {
    public:
        typedef typename
            wave_space_vector_type::const_iterator  const_iterator;
        typedef typename
            wave_space_vector_type::const_pointer   const_pointer;
        typedef typename
            wave_space_vector_type::const_reference const_reference;
        typedef typename
            wave_space_vector_type::iterator        iterator;
        typedef typename
            wave_space_vector_type::pointer         pointer;
        typedef typename
            wave_space_vector_type::reference       reference;
        typedef typename
            wave_space_vector_type::value_type      value_type;

        const dim_type start_x;
        const dim_type start_y;
        const dim_type start_z;
        const dim_type size_x;
        const dim_type size_y;
        const dim_type size_z;
        const dim_type size;

        reference operator()(
            const size_type x, const size_type y, const size_type z);
        const_reference operator()(
            const size_type x, const size_type y, const size_type z) const;

        physical_space::reference real(
            const size_type x, const size_type y, const size_type z);
        physical_space::const_reference real(
            const size_type x, const size_type y, const size_type z) const;
        physical_space::reference imag(
            const size_type x, const size_type y, const size_type z);
        physical_space::const_reference imag(
            const size_type x, const size_type y, const size_type z) const;

        iterator       begin();
        const_iterator begin() const;
        iterator       end();
        const_iterator end() const;

    private:
        friend class pencil<T,G>;

        wave_space(
            const dim_type start[3],
            const dim_type size[3],
            physical_space::pointer data)
        throw(domain_error);

        size_type offset(
            const size_type x,
            const size_type y,
            const size_type z) const;

        wave_space_adaptor_type adaptor_;
        wave_space_vector_type  vector_;
        physical_space_adaptor_type adaptor_components_;
        physical_space_vector_type  vector_components_;
    };

    pencil(const dim_type pstart[3], const dim_type psize[3],
           const dim_type wstart[3], const dim_type wsize[3])
    throw(domain_error);

    typedef typename physical_space::const_iterator const_iterator;
    typedef typename physical_space::iterator       iterator;
    iterator       begin();
    const_iterator begin() const;
    iterator       end();
    const_iterator end() const;

    physical_space physical;
    wave_space     wave;
};

template<typename T, typename G>
pencil<T, G>::pencil(
    const dim_type pstart[3], const dim_type psize[3],
    const dim_type wstart[3], const dim_type wsize[3])
throw(domain_error)
        : array_nelem_(std::max(
              psize[0]*psize[1]*psize[2],
            2*wsize[0]*wsize[1]*wsize[2])),
        array_(new real_type[array_nelem_]),
        physical(pstart, psize, array_.get()),
        wave(wstart, wsize, array_.get())
{
    // NOP
}

template<typename T, typename G>
pencil<T, G>::iterator pencil<T, G>::begin()
{
    return array_.get();
}

template<typename T, typename G>
pencil<T, G>::const_iterator pencil<T, G>::begin() const
{
    return array_.get();
}

template<typename T, typename G>
pencil<T, G>::iterator pencil<T, G>::end()
{
    return array_.get() + array_nelem_;
}

template<typename T, typename G>
pencil<T, G>::const_iterator pencil<T, G>::end() const
{
    return array_.get() + array_nelem_;
}

template<typename T, typename G>
pencil<T, G>::physical_space::physical_space(
    const dim_type start[3], const dim_type size[3], pointer data)
throw(domain_error)
        :
        start_x(start[0]), start_y(start[1]), start_z(start[2]),
        size_x(size[0]), size_y(size[1]), size_z(size[2]),
        size(size_x*size_y*size_z),
        adaptor_(this->size, data),
        vector_(this->size, adaptor_)
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
    physical_space::pointer data)
throw(domain_error)
        :
        start_x(start[0]), start_y(start[1]), start_z(start[2]),
        size_x(size[0]), size_y(size[1]), size_z(size[2]),
        size(size_x*size_y*size_z),
        adaptor_(this->size,reinterpret_cast<wave_space::pointer>(data)),
        vector_(this->size, adaptor_),
        adaptor_components_(2*this->size, data),
        vector_components_(2*this->size, adaptor_components_)
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
pencil<T, G>::physical_space::reference
pencil<T, G>::physical_space::operator()(
    const size_type x, const size_type y, const size_type z)
{
    return vector_(physical_space::offset(x, y, z));
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::const_reference
pencil<T, G>::physical_space::operator()(
    const size_type x, const size_type y, const size_type z) const
{
    return vector_(physical_space::offset(x, y, z));
}

template<typename T, typename G>
inline
pencil<T, G>::wave_space::reference
pencil<T, G>::wave_space::operator()(
    const size_type x, const size_type y, const size_type z)
{
    return vector_(wave_space::offset(x, y, z));
}

template<typename T, typename G>
inline
pencil<T, G>::wave_space::const_reference
pencil<T, G>::wave_space::operator()(
    const size_type x, const size_type y, const size_type z) const
{
    return vector_(wave_space::offset(x, y, z));
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::reference
pencil<T, G>::wave_space::real(
    const size_type x, const size_type y, const size_type z)
{
    return vector_components_(2*wave_space::offset(x, y, z));
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::const_reference
pencil<T, G>::wave_space::real(
    const size_type x, const size_type y, const size_type z) const
{
    return vector_components_(2*wave_space::offset(x, y, z));
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::reference
pencil<T, G>::wave_space::imag(
    const size_type x, const size_type y, const size_type z)
{
    return vector_components_(2*wave_space::offset(x, y, z) + 1);
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::const_reference
pencil<T, G>::wave_space::imag(
    const size_type x, const size_type y, const size_type z) const
{
    return vector_components_(2*wave_space::offset(x, y, z) + 1);
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::iterator
pencil<T, G>::physical_space::begin()
{
    return vector_.begin();
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::const_iterator
pencil<T, G>::physical_space::begin() const
{
    return vector_.begin();
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::iterator
pencil<T, G>::physical_space::end()
{
    return vector_.end();
}

template<typename T, typename G>
inline
pencil<T, G>::physical_space::const_iterator
pencil<T, G>::physical_space::end() const
{
    return vector_.end();
}

template<typename T, typename G>
inline
pencil<T, G>::wave_space::iterator
pencil<T, G>::wave_space::begin()
{
    return vector_.begin();
}

template<typename T, typename G>
inline
pencil<T, G>::wave_space::const_iterator
pencil<T, G>::wave_space::begin() const
{
    return vector_.begin();
}

template<typename T, typename G>
inline
pencil<T, G>::wave_space::iterator
pencil<T, G>::wave_space::end()
{
    return vector_.end();
}

template<typename T, typename G>
inline
pencil<T, G>::wave_space::const_iterator
pencil<T, G>::wave_space::end() const
{
    return vector_.end();
}

}

}

#endif // PECOS_SUZERAIN_PENCIL
