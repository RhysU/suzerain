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

// Boilerplate common to each template specialization
// (static const member initialized in matching .cpp)
#define SPECIALIZE(T,C)                                  \
    template<> struct datatype< T > {                    \
        static const MPI_Datatype value;                 \
        operator MPI_Datatype () const { return value; } \
    }

// Template specializations translating from C++ to MPI_Datatypes
SPECIALIZE(char,                      MPI_CHAR);
SPECIALIZE(signed char,               MPI_SIGNED_CHAR);
SPECIALIZE(unsigned char,             MPI_UNSIGNED_CHAR);
SPECIALIZE(wchar_t,                   MPI_WCHAR);
SPECIALIZE(signed short int,          MPI_SHORT);
SPECIALIZE(signed int,                MPI_INT);
SPECIALIZE(signed long int,           MPI_LONG);
SPECIALIZE(signed long long int,      MPI_LONG_LONG);
SPECIALIZE(unsigned short int,        MPI_UNSIGNED_SHORT);
SPECIALIZE(unsigned int,              MPI_UNSIGNED);
SPECIALIZE(unsigned long int,         MPI_UNSIGNED_LONG);
SPECIALIZE(unsigned long long int,    MPI_UNSIGNED_LONG_LONG);
SPECIALIZE(float,                     MPI_FLOAT);
SPECIALIZE(double,                    MPI_DOUBLE);
SPECIALIZE(long double,               MPI_LONG_DOUBLE);
SPECIALIZE(std::complex<float>,       MPI_C_FLOAT_COMPLEX);
SPECIALIZE(std::complex<double>,      MPI_C_DOUBLE_COMPLEX);
SPECIALIZE(std::complex<long double>, MPI_C_LONG_DOUBLE_COMPLEX);

#undef SPECIALIZE

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
