/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2010 The PECOS Development Team
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
 * mpi_datatype.hpp: lookup functionality for MPI_Datatype constants
 *
 * $Id$
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __SUZERAIN_MPI_DATATYPE_HPP
#define __SUZERAIN_MPI_DATATYPE_HPP

#include <suzerain/common.hpp>
#include <suzerain/mpi.h>

namespace suzerain {

namespace mpi {

/** Provides template-friendly lookup of MPI_Datatype constants */
template< typename T, class Enable = void > struct datatype {
    BOOST_MPL_ASSERT_MSG(
        sizeof(T) == 0, SANITY_FAILURE_BASE_TEMPLATE_INSTANTIATED, (T));
};

/** Remove any pointer modifiers during MPI_Datatype lookups */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_pointer<T> >::type>
    : public datatype<typename boost::remove_pointer<T>::type> {};

/** Remove any reference modifiers during MPI_Datatype lookups */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_reference<T> >::type>
    : public datatype<typename boost::remove_reference<T>::type> {};

/** Remove any extent modifiers during MPI_Datatype lookups */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_array<T> >::type>
    : public datatype<typename boost::remove_extent<T>::type> {};

/** Provides MPI_Datatype lookup for <tt>char</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    char
> >::type> {
    operator MPI_Datatype () const { return MPI_CHAR; }
};

/** Provides MPI_Datatype lookup for <tt>signed char</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed char
> >::type> {
    operator MPI_Datatype () const { return MPI_SIGNED_CHAR; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned char</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned char
> >::type> {
    operator MPI_Datatype () const { return MPI_UNSIGNED_CHAR; }
};

/** Provides MPI_Datatype lookup for <tt>signed short int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed short int
> >::type> {
    operator MPI_Datatype () const { return MPI_SHORT; }
};

/** Provides MPI_Datatype lookup for <tt>signed int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed int
> >::type> {
    operator MPI_Datatype () const { return MPI_INT; }
};

/** Provides MPI_Datatype lookup for <tt>signed long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed long int
> >::type> {
    operator MPI_Datatype () const { return MPI_LONG; }
};

/** Provides MPI_Datatype lookup for <tt>signed long long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed long long int
> >::type> {
    operator MPI_Datatype () const { return MPI_LONG_LONG; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned short int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned short int
> >::type> {
    operator MPI_Datatype () const { return MPI_UNSIGNED_SHORT; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned int
> >::type> {
    operator MPI_Datatype () const { return MPI_UNSIGNED; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned long int
> >::type> {
    operator MPI_Datatype () const { return MPI_UNSIGNED_LONG; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned long long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned long long int
> >::type> {
    operator MPI_Datatype () const { return MPI_UNSIGNED_LONG_LONG; }
};

/** Provides MPI_Datatype lookup for <tt>float</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    float
> >::type> {
    operator MPI_Datatype () const { return MPI_FLOAT; }
};

/** Provides MPI_Datatype lookup for <tt>double</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    double
> >::type> {
    operator MPI_Datatype () const { return MPI_DOUBLE; }
};

/** Provides MPI_Datatype lookup for <tt>long double</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    long double
> >::type> {
    operator MPI_Datatype () const { return MPI_LONG_DOUBLE; }
};

/** Provides MPI_Datatype lookup for <tt>long double</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    wchar_t
> >::type> {
    operator MPI_Datatype () const { return MPI_WCHAR; }
};

/**
 * Perform MPI_Datatype lookup based on an argument.
 *
 * @param t Used to determine the appropriate MPI_Datatype.
 *
 * @return the MPI_Datatype corresponding to the type of \c t.
 */
template<typename T>
MPI_Datatype datatype_of(const T& t)
{
    SUZERAIN_UNUSED(t);
    return datatype<T>();
}

} // namespace mpi

} // namespace suzerain

#endif // __SUZERAIN_MPI_DATATYPE_HPP
