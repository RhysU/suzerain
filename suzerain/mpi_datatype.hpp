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

#ifndef SUZERAIN_MPI_DATATYPE_HPP
#define SUZERAIN_MPI_DATATYPE_HPP

/** @file
 * Automatic lookup of type-specific \c MPI_Datatype constants.
 */

#include <suzerain/common.hpp>
#include <suzerain/mpi.h>

namespace suzerain {

namespace mpi {

// Strategy is to remove array, reference, and pointer modifiers to
// obtain a possibly cv-qualified fundamental or complex type followed
// by dispatching on MPI-compatible fundamental or complex types.

/** Provides template-friendly lookup of MPI_Datatype constants */
template< typename T, class Enable = void > struct datatype {
    BOOST_MPL_ASSERT_MSG(
        sizeof(T) == 0, SANITY_FAILURE_BASE_TEMPLATE_INSTANTIATED, (T));
};

/** Removes any extent modifiers during MPI_Datatype lookups */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_array<T> >::type>
    : public datatype<typename boost::remove_extent<T>::type> {};

/** Removes any reference modifiers during MPI_Datatype lookups */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_reference<T> >::type>
    : public datatype<typename boost::remove_reference<T>::type> {};

/** Removes any pointer modifiers during MPI_Datatype lookups */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_pointer<T> >::type>
    : public datatype<typename boost::remove_pointer<T>::type> {};

/** Removes any const or volatile modifiers on a fundamental MPI_Datatype */
template<typename T>
struct datatype<T, typename boost::enable_if<boost::is_fundamental<T> >::type>
    : public datatype<typename boost::remove_cv<T>::type> {};

/** Removes any const modifier on a complex MPI_Datatype */
template<typename T>
struct datatype<const T, typename boost::enable_if<boost::mpl::and_<
    boost::is_complex<T>,
    boost::mpl::not_<boost::is_array<T> >,
    boost::mpl::not_<boost::is_reference<T> >,
    boost::mpl::not_<boost::is_pointer<T> >
> >::type>
    : public datatype<T> {};

/** Removes any volatile modifier on a complex MPI_Datatype */
template<typename T>
struct datatype<volatile T, typename boost::enable_if<boost::mpl::and_<
    boost::is_complex<T>,
    boost::mpl::not_<boost::is_array<T> >,
    boost::mpl::not_<boost::is_reference<T> >,
    boost::mpl::not_<boost::is_pointer<T> >
> >::type>
    : public datatype<T> {};

/** Specialized MPI_Datatype case for <tt>char</tt>s */
template<>
struct datatype<char> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<char>::value = MPI_CHAR;

/** Specialized MPI_Datatype case for <tt>signed char</tt>s */
template<>
struct datatype<signed char> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<signed char>::value = MPI_SIGNED_CHAR;

/** Specialized MPI_Datatype case for <tt>unsigned char</tt>s */
template<>
struct datatype<unsigned char> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<unsigned char>::value = MPI_UNSIGNED_CHAR;

/** Specialized MPI_Datatype case for <tt>signed short int</tt>s */
template<>
struct datatype<signed short int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<signed short int>::value = MPI_SHORT;

/** Specialized MPI_Datatype case for <tt>signed int</tt>s */
template<>
struct datatype<signed int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<signed int>::value = MPI_INT;

/** Specialized MPI_Datatype case for <tt>signed long int</tt>s */
template<>
struct datatype<signed long int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<signed long int>::value = MPI_LONG;

/** Specialized MPI_Datatype case for <tt>signed long long int</tt>s */
template<>
struct datatype<signed long long int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<signed long long int>::value = MPI_LONG_LONG;

/** Specialized MPI_Datatype case for <tt>unsigned short int</tt>s */
template<>
struct datatype<unsigned short int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<unsigned short int>::value = MPI_UNSIGNED_SHORT;

/** Specialized MPI_Datatype case for <tt>unsigned int</tt>s */
template<>
struct datatype<unsigned int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<unsigned int>::value = MPI_UNSIGNED;

/** Specialized MPI_Datatype case for <tt>unsigned long int</tt>s */
template<>
struct datatype<unsigned long int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<unsigned long int>::value = MPI_UNSIGNED_LONG;

/** Specialized MPI_Datatype case for <tt>unsigned long long int</tt>s */
template<>
struct datatype<unsigned long long int> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<unsigned long long int>::value = MPI_UNSIGNED_LONG_LONG;

/** Specialized MPI_Datatype case for <tt>float</tt>s */
template<>
struct datatype<float> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<float>::value = MPI_FLOAT;

/** Specialized MPI_Datatype case for <tt>double</tt>s */
template<>
struct datatype<double> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<double>::value = MPI_DOUBLE;

/** Specialized MPI_Datatype case for <tt>long double</tt>s */
template<>
struct datatype<long double> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<long double>::value = MPI_LONG_DOUBLE;

/** Specialized MPI_Datatype case for <tt>std::complex<float></tt>s */
template<>
struct datatype<std::complex<float> > {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<std::complex<float> >::value = MPI_C_FLOAT_COMPLEX;

/** Specialized MPI_Datatype case for <tt>std::complex<double></tt>s */
template<>
struct datatype<std::complex<double> > {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<std::complex<double> >::value = MPI_C_DOUBLE_COMPLEX;

/** Specialized MPI_Datatype case for <tt>std::complex<long double></tt>s */
template<>
struct datatype<std::complex<long double> > {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<std::complex<long double> >::value = MPI_C_LONG_DOUBLE_COMPLEX;

/** Specialized MPI_Datatype case for <tt>wchar_t</tt>s */
template<>
struct datatype<wchar_t> {
    static const MPI_Datatype value;
    operator MPI_Datatype () const { return value; }
};

const MPI_Datatype datatype<wchar_t>::value = MPI_WCHAR;

/**
 * Perform MPI_Datatype lookup based on an argument.
 *
 * @return the MPI_Datatype corresponding to the type of \c t.
 */
template<typename T>
MPI_Datatype datatype_of(const T&)
{
    return datatype<T>();
}

} // namespace mpi

} // namespace suzerain

#endif // SUZERAIN_MPI_DATATYPE_HPP
