//--------------------------------------------------------------------------
//
// Copyright (C) 2010, 2011, 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------
// mpi_datatype.hpp: lookup functionality for MPI_Datatype constants
// $Id$

#ifndef SUZERAIN_MPI_DATATYPE_HPP
#define SUZERAIN_MPI_DATATYPE_HPP

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
    static const MPI_Datatype value = MPI_CHAR;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>signed char</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed char
> >::type> {
    static const MPI_Datatype value = MPI_SIGNED_CHAR;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned char</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned char
> >::type> {
    static const MPI_Datatype value = MPI_UNSIGNED_CHAR;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>signed short int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed short int
> >::type> {
    static const MPI_Datatype value = MPI_SHORT;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>signed int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed int
> >::type> {
    static const MPI_Datatype value = MPI_INT;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>signed long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed long int
> >::type> {
    static const MPI_Datatype value = MPI_LONG;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>signed long long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    signed long long int
> >::type> {
    static const MPI_Datatype value = MPI_LONG_LONG;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned short int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned short int
> >::type> {
    static const MPI_Datatype value = MPI_UNSIGNED_SHORT;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned int
> >::type> {
    static const MPI_Datatype value = MPI_UNSIGNED;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned long int
> >::type> {
    static const MPI_Datatype value = MPI_UNSIGNED_LONG;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>unsigned long long int</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    unsigned long long int
> >::type> {
    static const MPI_Datatype value = MPI_UNSIGNED_LONG_LONG;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>float</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    float
> >::type> {
    static const MPI_Datatype value = MPI_FLOAT;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>double</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    double
> >::type> {
    static const MPI_Datatype value = MPI_DOUBLE;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>long double</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    long double
> >::type> {
    static const MPI_Datatype value = MPI_LONG_DOUBLE;
    operator MPI_Datatype () const { return value; }
};

/** Provides MPI_Datatype lookup for <tt>long double</tt>s */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_same<
    typename boost::remove_cv<T>::type,
    wchar_t
> >::type> {
    static const MPI_Datatype value = MPI_WCHAR;
    operator MPI_Datatype () const { return value; }
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

#endif // SUZERAIN_MPI_DATATYPE_HPP
