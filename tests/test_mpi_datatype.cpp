//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#include <suzerain/mpi_datatype.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

#pragma warning(disable:1572)

using suzerain::mpi::datatype;

BOOST_AUTO_TEST_CASE( datatype_char )
{
    const MPI_Datatype e = MPI_CHAR;
    BOOST_CHECK(e == datatype<char>());
    BOOST_CHECK(e == datatype<char*>());
    BOOST_CHECK(e == datatype<char&>());
    BOOST_CHECK(e == datatype<char[3]>());
    BOOST_CHECK(e == datatype<char const>());
    BOOST_CHECK(e == datatype<char const*>());
    BOOST_CHECK(e == datatype<char const&>());
    BOOST_CHECK(e == datatype<char const[4]>());
    BOOST_CHECK(e == datatype<char volatile>());
    BOOST_CHECK(e == datatype<char volatile*>());
    BOOST_CHECK(e == datatype<char volatile&>());
    BOOST_CHECK(e == datatype<char volatile[]>());
    BOOST_CHECK(e == datatype<char volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_short )
{
    const MPI_Datatype e = MPI_SHORT;

    BOOST_CHECK(e == datatype<short>());
    BOOST_CHECK(e == datatype<short*>());
    BOOST_CHECK(e == datatype<short&>());
    BOOST_CHECK(e == datatype<short[3]>());
    BOOST_CHECK(e == datatype<short const>());
    BOOST_CHECK(e == datatype<short const*>());
    BOOST_CHECK(e == datatype<short const&>());
    BOOST_CHECK(e == datatype<short const[4]>());
    BOOST_CHECK(e == datatype<short volatile>());
    BOOST_CHECK(e == datatype<short volatile*>());
    BOOST_CHECK(e == datatype<short volatile&>());
    BOOST_CHECK(e == datatype<short volatile[]>());
    BOOST_CHECK(e == datatype<short volatile const>());

    BOOST_CHECK(e == datatype<signed short>());
    BOOST_CHECK(e == datatype<signed short*>());
    BOOST_CHECK(e == datatype<signed short&>());
    BOOST_CHECK(e == datatype<signed short[3]>());
    BOOST_CHECK(e == datatype<signed short const>());
    BOOST_CHECK(e == datatype<signed short const*>());
    BOOST_CHECK(e == datatype<signed short const&>());
    BOOST_CHECK(e == datatype<signed short const[4]>());
    BOOST_CHECK(e == datatype<signed short volatile>());
    BOOST_CHECK(e == datatype<signed short volatile*>());
    BOOST_CHECK(e == datatype<signed short volatile&>());
    BOOST_CHECK(e == datatype<signed short volatile[]>());
    BOOST_CHECK(e == datatype<signed short volatile const>());

    BOOST_CHECK(e == datatype<short int>());
    BOOST_CHECK(e == datatype<short int*>());
    BOOST_CHECK(e == datatype<short int&>());
    BOOST_CHECK(e == datatype<short int[3]>());
    BOOST_CHECK(e == datatype<short int const>());
    BOOST_CHECK(e == datatype<short int const*>());
    BOOST_CHECK(e == datatype<short int const&>());
    BOOST_CHECK(e == datatype<short int const[4]>());
    BOOST_CHECK(e == datatype<short int volatile>());
    BOOST_CHECK(e == datatype<short int volatile*>());
    BOOST_CHECK(e == datatype<short int volatile&>());
    BOOST_CHECK(e == datatype<short int volatile[]>());
    BOOST_CHECK(e == datatype<short int volatile const>());

    BOOST_CHECK(e == datatype<signed short int>());
    BOOST_CHECK(e == datatype<signed short int*>());
    BOOST_CHECK(e == datatype<signed short int&>());
    BOOST_CHECK(e == datatype<signed short int[3]>());
    BOOST_CHECK(e == datatype<signed short int const>());
    BOOST_CHECK(e == datatype<signed short int const*>());
    BOOST_CHECK(e == datatype<signed short int const&>());
    BOOST_CHECK(e == datatype<signed short int const[4]>());
    BOOST_CHECK(e == datatype<signed short int volatile>());
    BOOST_CHECK(e == datatype<signed short int volatile*>());
    BOOST_CHECK(e == datatype<signed short int volatile&>());
    BOOST_CHECK(e == datatype<signed short int volatile[]>());
    BOOST_CHECK(e == datatype<signed short int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_int )
{
    const MPI_Datatype e = MPI_INT;

    BOOST_CHECK(e == datatype<int>());
    BOOST_CHECK(e == datatype<int*>());
    BOOST_CHECK(e == datatype<int&>());
    BOOST_CHECK(e == datatype<int[3]>());
    BOOST_CHECK(e == datatype<int const>());
    BOOST_CHECK(e == datatype<int const*>());
    BOOST_CHECK(e == datatype<int const&>());
    BOOST_CHECK(e == datatype<int const[4]>());
    BOOST_CHECK(e == datatype<int volatile>());
    BOOST_CHECK(e == datatype<int volatile*>());
    BOOST_CHECK(e == datatype<int volatile&>());
    BOOST_CHECK(e == datatype<int volatile[]>());
    BOOST_CHECK(e == datatype<int volatile const>());

    BOOST_CHECK(e == datatype<signed int>());
    BOOST_CHECK(e == datatype<signed int*>());
    BOOST_CHECK(e == datatype<signed int&>());
    BOOST_CHECK(e == datatype<signed int[3]>());
    BOOST_CHECK(e == datatype<signed int const>());
    BOOST_CHECK(e == datatype<signed int const*>());
    BOOST_CHECK(e == datatype<signed int const&>());
    BOOST_CHECK(e == datatype<signed int const[4]>());
    BOOST_CHECK(e == datatype<signed int volatile>());
    BOOST_CHECK(e == datatype<signed int volatile*>());
    BOOST_CHECK(e == datatype<signed int volatile&>());
    BOOST_CHECK(e == datatype<signed int volatile[]>());
    BOOST_CHECK(e == datatype<signed int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_long )
{
    const MPI_Datatype e = MPI_LONG;

    BOOST_CHECK(e == datatype<long>());
    BOOST_CHECK(e == datatype<long*>());
    BOOST_CHECK(e == datatype<long&>());
    BOOST_CHECK(e == datatype<long[3]>());
    BOOST_CHECK(e == datatype<long const>());
    BOOST_CHECK(e == datatype<long const*>());
    BOOST_CHECK(e == datatype<long const&>());
    BOOST_CHECK(e == datatype<long const[4]>());
    BOOST_CHECK(e == datatype<long volatile>());
    BOOST_CHECK(e == datatype<long volatile*>());
    BOOST_CHECK(e == datatype<long volatile&>());
    BOOST_CHECK(e == datatype<long volatile[]>());
    BOOST_CHECK(e == datatype<long volatile const>());

    BOOST_CHECK(e == datatype<signed long>());
    BOOST_CHECK(e == datatype<signed long*>());
    BOOST_CHECK(e == datatype<signed long&>());
    BOOST_CHECK(e == datatype<signed long[3]>());
    BOOST_CHECK(e == datatype<signed long const>());
    BOOST_CHECK(e == datatype<signed long const*>());
    BOOST_CHECK(e == datatype<signed long const&>());
    BOOST_CHECK(e == datatype<signed long const[4]>());
    BOOST_CHECK(e == datatype<signed long volatile>());
    BOOST_CHECK(e == datatype<signed long volatile*>());
    BOOST_CHECK(e == datatype<signed long volatile&>());
    BOOST_CHECK(e == datatype<signed long volatile[]>());
    BOOST_CHECK(e == datatype<signed long volatile const>());

    BOOST_CHECK(e == datatype<long int>());
    BOOST_CHECK(e == datatype<long int*>());
    BOOST_CHECK(e == datatype<long int&>());
    BOOST_CHECK(e == datatype<long int[3]>());
    BOOST_CHECK(e == datatype<long int const>());
    BOOST_CHECK(e == datatype<long int const*>());
    BOOST_CHECK(e == datatype<long int const&>());
    BOOST_CHECK(e == datatype<long int const[4]>());
    BOOST_CHECK(e == datatype<long int volatile>());
    BOOST_CHECK(e == datatype<long int volatile*>());
    BOOST_CHECK(e == datatype<long int volatile&>());
    BOOST_CHECK(e == datatype<long int volatile[]>());
    BOOST_CHECK(e == datatype<long int volatile const>());

    BOOST_CHECK(e == datatype<signed long int>());
    BOOST_CHECK(e == datatype<signed long int*>());
    BOOST_CHECK(e == datatype<signed long int&>());
    BOOST_CHECK(e == datatype<signed long int[3]>());
    BOOST_CHECK(e == datatype<signed long int const>());
    BOOST_CHECK(e == datatype<signed long int const*>());
    BOOST_CHECK(e == datatype<signed long int const&>());
    BOOST_CHECK(e == datatype<signed long int const[4]>());
    BOOST_CHECK(e == datatype<signed long int volatile>());
    BOOST_CHECK(e == datatype<signed long int volatile*>());
    BOOST_CHECK(e == datatype<signed long int volatile&>());
    BOOST_CHECK(e == datatype<signed long int volatile[]>());
    BOOST_CHECK(e == datatype<signed long int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_long_long )
{
    const MPI_Datatype e = MPI_LONG_LONG;
    BOOST_CHECK_EQUAL(MPI_LONG_LONG_INT, MPI_LONG_LONG);

    BOOST_CHECK(e == datatype<long long>());
    BOOST_CHECK(e == datatype<long long*>());
    BOOST_CHECK(e == datatype<long long&>());
    BOOST_CHECK(e == datatype<long long[3]>());
    BOOST_CHECK(e == datatype<long long const>());
    BOOST_CHECK(e == datatype<long long const*>());
    BOOST_CHECK(e == datatype<long long const&>());
    BOOST_CHECK(e == datatype<long long const[4]>());
    BOOST_CHECK(e == datatype<long long volatile>());
    BOOST_CHECK(e == datatype<long long volatile*>());
    BOOST_CHECK(e == datatype<long long volatile&>());
    BOOST_CHECK(e == datatype<long long volatile[]>());
    BOOST_CHECK(e == datatype<long long volatile const>());

    BOOST_CHECK(e == datatype<signed long long>());
    BOOST_CHECK(e == datatype<signed long long*>());
    BOOST_CHECK(e == datatype<signed long long&>());
    BOOST_CHECK(e == datatype<signed long long[3]>());
    BOOST_CHECK(e == datatype<signed long long const>());
    BOOST_CHECK(e == datatype<signed long long const*>());
    BOOST_CHECK(e == datatype<signed long long const&>());
    BOOST_CHECK(e == datatype<signed long long const[4]>());
    BOOST_CHECK(e == datatype<signed long long volatile>());
    BOOST_CHECK(e == datatype<signed long long volatile*>());
    BOOST_CHECK(e == datatype<signed long long volatile&>());
    BOOST_CHECK(e == datatype<signed long long volatile[]>());
    BOOST_CHECK(e == datatype<signed long long volatile const>());

    BOOST_CHECK(e == datatype<long long int>());
    BOOST_CHECK(e == datatype<long long int*>());
    BOOST_CHECK(e == datatype<long long int&>());
    BOOST_CHECK(e == datatype<long long int[3]>());
    BOOST_CHECK(e == datatype<long long int const>());
    BOOST_CHECK(e == datatype<long long int const*>());
    BOOST_CHECK(e == datatype<long long int const&>());
    BOOST_CHECK(e == datatype<long long int const[4]>());
    BOOST_CHECK(e == datatype<long long int volatile>());
    BOOST_CHECK(e == datatype<long long int volatile*>());
    BOOST_CHECK(e == datatype<long long int volatile&>());
    BOOST_CHECK(e == datatype<long long int volatile[]>());
    BOOST_CHECK(e == datatype<long long int volatile const>());

    BOOST_CHECK(e == datatype<signed long long int>());
    BOOST_CHECK(e == datatype<signed long long int*>());
    BOOST_CHECK(e == datatype<signed long long int&>());
    BOOST_CHECK(e == datatype<signed long long int[3]>());
    BOOST_CHECK(e == datatype<signed long long int const>());
    BOOST_CHECK(e == datatype<signed long long int const*>());
    BOOST_CHECK(e == datatype<signed long long int const&>());
    BOOST_CHECK(e == datatype<signed long long int const[4]>());
    BOOST_CHECK(e == datatype<signed long long int volatile>());
    BOOST_CHECK(e == datatype<signed long long int volatile*>());
    BOOST_CHECK(e == datatype<signed long long int volatile&>());
    BOOST_CHECK(e == datatype<signed long long int volatile[]>());
    BOOST_CHECK(e == datatype<signed long long int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_signed_char )
{
    const MPI_Datatype e = MPI_SIGNED_CHAR;
    BOOST_CHECK(e == datatype<signed char>());
    BOOST_CHECK(e == datatype<signed char*>());
    BOOST_CHECK(e == datatype<signed char&>());
    BOOST_CHECK(e == datatype<signed char[3]>());
    BOOST_CHECK(e == datatype<signed char const>());
    BOOST_CHECK(e == datatype<signed char const*>());
    BOOST_CHECK(e == datatype<signed char const&>());
    BOOST_CHECK(e == datatype<signed char const[4]>());
    BOOST_CHECK(e == datatype<signed char volatile>());
    BOOST_CHECK(e == datatype<signed char volatile*>());
    BOOST_CHECK(e == datatype<signed char volatile&>());
    BOOST_CHECK(e == datatype<signed char volatile[]>());
    BOOST_CHECK(e == datatype<signed char volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_unsigned_char )
{
    const MPI_Datatype e = MPI_UNSIGNED_CHAR;
    BOOST_CHECK(e == datatype<unsigned char>());
    BOOST_CHECK(e == datatype<unsigned char*>());
    BOOST_CHECK(e == datatype<unsigned char&>());
    BOOST_CHECK(e == datatype<unsigned char[3]>());
    BOOST_CHECK(e == datatype<unsigned char const>());
    BOOST_CHECK(e == datatype<unsigned char const*>());
    BOOST_CHECK(e == datatype<unsigned char const&>());
    BOOST_CHECK(e == datatype<unsigned char const[4]>());
    BOOST_CHECK(e == datatype<unsigned char volatile>());
    BOOST_CHECK(e == datatype<unsigned char volatile*>());
    BOOST_CHECK(e == datatype<unsigned char volatile&>());
    BOOST_CHECK(e == datatype<unsigned char volatile[]>());
    BOOST_CHECK(e == datatype<unsigned char volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_unsigned_short )
{
    const MPI_Datatype e = MPI_UNSIGNED_SHORT;

    BOOST_CHECK(e == datatype<unsigned short>());
    BOOST_CHECK(e == datatype<unsigned short*>());
    BOOST_CHECK(e == datatype<unsigned short&>());
    BOOST_CHECK(e == datatype<unsigned short[3]>());
    BOOST_CHECK(e == datatype<unsigned short const>());
    BOOST_CHECK(e == datatype<unsigned short const*>());
    BOOST_CHECK(e == datatype<unsigned short const&>());
    BOOST_CHECK(e == datatype<unsigned short const[4]>());
    BOOST_CHECK(e == datatype<unsigned short volatile>());
    BOOST_CHECK(e == datatype<unsigned short volatile*>());
    BOOST_CHECK(e == datatype<unsigned short volatile&>());
    BOOST_CHECK(e == datatype<unsigned short volatile[]>());
    BOOST_CHECK(e == datatype<unsigned short volatile const>());

    BOOST_CHECK(e == datatype<unsigned short int>());
    BOOST_CHECK(e == datatype<unsigned short int*>());
    BOOST_CHECK(e == datatype<unsigned short int&>());
    BOOST_CHECK(e == datatype<unsigned short int[3]>());
    BOOST_CHECK(e == datatype<unsigned short int const>());
    BOOST_CHECK(e == datatype<unsigned short int const*>());
    BOOST_CHECK(e == datatype<unsigned short int const&>());
    BOOST_CHECK(e == datatype<unsigned short int const[4]>());
    BOOST_CHECK(e == datatype<unsigned short int volatile>());
    BOOST_CHECK(e == datatype<unsigned short int volatile*>());
    BOOST_CHECK(e == datatype<unsigned short int volatile&>());
    BOOST_CHECK(e == datatype<unsigned short int volatile[]>());
    BOOST_CHECK(e == datatype<unsigned short int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_unsigned )
{
    const MPI_Datatype e = MPI_UNSIGNED;

    BOOST_CHECK(e == datatype<unsigned>());
    BOOST_CHECK(e == datatype<unsigned*>());
    BOOST_CHECK(e == datatype<unsigned&>());
    BOOST_CHECK(e == datatype<unsigned[3]>());
    BOOST_CHECK(e == datatype<unsigned const>());
    BOOST_CHECK(e == datatype<unsigned const*>());
    BOOST_CHECK(e == datatype<unsigned const&>());
    BOOST_CHECK(e == datatype<unsigned const[4]>());
    BOOST_CHECK(e == datatype<unsigned volatile>());
    BOOST_CHECK(e == datatype<unsigned volatile*>());
    BOOST_CHECK(e == datatype<unsigned volatile&>());
    BOOST_CHECK(e == datatype<unsigned volatile[]>());
    BOOST_CHECK(e == datatype<unsigned volatile const>());

    BOOST_CHECK(e == datatype<unsigned int>());
    BOOST_CHECK(e == datatype<unsigned int*>());
    BOOST_CHECK(e == datatype<unsigned int&>());
    BOOST_CHECK(e == datatype<unsigned int[3]>());
    BOOST_CHECK(e == datatype<unsigned int const>());
    BOOST_CHECK(e == datatype<unsigned int const*>());
    BOOST_CHECK(e == datatype<unsigned int const&>());
    BOOST_CHECK(e == datatype<unsigned int const[4]>());
    BOOST_CHECK(e == datatype<unsigned int volatile>());
    BOOST_CHECK(e == datatype<unsigned int volatile*>());
    BOOST_CHECK(e == datatype<unsigned int volatile&>());
    BOOST_CHECK(e == datatype<unsigned int volatile[]>());
    BOOST_CHECK(e == datatype<unsigned int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_unsigned_long )
{
    const MPI_Datatype e = MPI_UNSIGNED_LONG;

    BOOST_CHECK(e == datatype<unsigned long>());
    BOOST_CHECK(e == datatype<unsigned long*>());
    BOOST_CHECK(e == datatype<unsigned long&>());
    BOOST_CHECK(e == datatype<unsigned long[3]>());
    BOOST_CHECK(e == datatype<unsigned long const>());
    BOOST_CHECK(e == datatype<unsigned long const*>());
    BOOST_CHECK(e == datatype<unsigned long const&>());
    BOOST_CHECK(e == datatype<unsigned long const[4]>());
    BOOST_CHECK(e == datatype<unsigned long volatile>());
    BOOST_CHECK(e == datatype<unsigned long volatile*>());
    BOOST_CHECK(e == datatype<unsigned long volatile&>());
    BOOST_CHECK(e == datatype<unsigned long volatile[]>());
    BOOST_CHECK(e == datatype<unsigned long volatile const>());

    BOOST_CHECK(e == datatype<unsigned long int>());
    BOOST_CHECK(e == datatype<unsigned long int*>());
    BOOST_CHECK(e == datatype<unsigned long int&>());
    BOOST_CHECK(e == datatype<unsigned long int[3]>());
    BOOST_CHECK(e == datatype<unsigned long int const>());
    BOOST_CHECK(e == datatype<unsigned long int const*>());
    BOOST_CHECK(e == datatype<unsigned long int const&>());
    BOOST_CHECK(e == datatype<unsigned long int const[4]>());
    BOOST_CHECK(e == datatype<unsigned long int volatile>());
    BOOST_CHECK(e == datatype<unsigned long int volatile*>());
    BOOST_CHECK(e == datatype<unsigned long int volatile&>());
    BOOST_CHECK(e == datatype<unsigned long int volatile[]>());
    BOOST_CHECK(e == datatype<unsigned long int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_unsigned_long_long )
{
    const MPI_Datatype e = MPI_UNSIGNED_LONG_LONG;

    BOOST_CHECK(e == datatype<unsigned long long>());
    BOOST_CHECK(e == datatype<unsigned long long*>());
    BOOST_CHECK(e == datatype<unsigned long long&>());
    BOOST_CHECK(e == datatype<unsigned long long[3]>());
    BOOST_CHECK(e == datatype<unsigned long long const>());
    BOOST_CHECK(e == datatype<unsigned long long const*>());
    BOOST_CHECK(e == datatype<unsigned long long const&>());
    BOOST_CHECK(e == datatype<unsigned long long const[4]>());
    BOOST_CHECK(e == datatype<unsigned long long volatile>());
    BOOST_CHECK(e == datatype<unsigned long long volatile*>());
    BOOST_CHECK(e == datatype<unsigned long long volatile&>());
    BOOST_CHECK(e == datatype<unsigned long long volatile[]>());
    BOOST_CHECK(e == datatype<unsigned long long volatile const>());

    BOOST_CHECK(e == datatype<unsigned long long int>());
    BOOST_CHECK(e == datatype<unsigned long long int*>());
    BOOST_CHECK(e == datatype<unsigned long long int&>());
    BOOST_CHECK(e == datatype<unsigned long long int[3]>());
    BOOST_CHECK(e == datatype<unsigned long long int const>());
    BOOST_CHECK(e == datatype<unsigned long long int const*>());
    BOOST_CHECK(e == datatype<unsigned long long int const&>());
    BOOST_CHECK(e == datatype<unsigned long long int const[4]>());
    BOOST_CHECK(e == datatype<unsigned long long int volatile>());
    BOOST_CHECK(e == datatype<unsigned long long int volatile*>());
    BOOST_CHECK(e == datatype<unsigned long long int volatile&>());
    BOOST_CHECK(e == datatype<unsigned long long int volatile[]>());
    BOOST_CHECK(e == datatype<unsigned long long int volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_float )
{
    const MPI_Datatype e = MPI_FLOAT;
    BOOST_CHECK(e == datatype<float>());
    BOOST_CHECK(e == datatype<float*>());
    BOOST_CHECK(e == datatype<float&>());
    BOOST_CHECK(e == datatype<float[3]>());
    BOOST_CHECK(e == datatype<float const>());
    BOOST_CHECK(e == datatype<float const*>());
    BOOST_CHECK(e == datatype<float const&>());
    BOOST_CHECK(e == datatype<float const[4]>());
    BOOST_CHECK(e == datatype<float volatile>());
    BOOST_CHECK(e == datatype<float volatile*>());
    BOOST_CHECK(e == datatype<float volatile&>());
    BOOST_CHECK(e == datatype<float volatile[]>());
    BOOST_CHECK(e == datatype<float volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_double )
{
    const MPI_Datatype e = MPI_DOUBLE;
    BOOST_CHECK(e == datatype<double>());
    BOOST_CHECK(e == datatype<double*>());
    BOOST_CHECK(e == datatype<double&>());
    BOOST_CHECK(e == datatype<double[3]>());
    BOOST_CHECK(e == datatype<double const>());
    BOOST_CHECK(e == datatype<double const*>());
    BOOST_CHECK(e == datatype<double const&>());
    BOOST_CHECK(e == datatype<double const[4]>());
    BOOST_CHECK(e == datatype<double volatile>());
    BOOST_CHECK(e == datatype<double volatile*>());
    BOOST_CHECK(e == datatype<double volatile&>());
    BOOST_CHECK(e == datatype<double volatile[]>());
    BOOST_CHECK(e == datatype<double volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_long_double )
{
    const MPI_Datatype e = MPI_LONG_DOUBLE;
    BOOST_CHECK(e == datatype<long double>());
    BOOST_CHECK(e == datatype<long double*>());
    BOOST_CHECK(e == datatype<long double&>());
    BOOST_CHECK(e == datatype<long double[3]>());
    BOOST_CHECK(e == datatype<long double const>());
    BOOST_CHECK(e == datatype<long double const*>());
    BOOST_CHECK(e == datatype<long double const&>());
    BOOST_CHECK(e == datatype<long double const[4]>());
    BOOST_CHECK(e == datatype<long double volatile>());
    BOOST_CHECK(e == datatype<long double volatile*>());
    BOOST_CHECK(e == datatype<long double volatile&>());
    BOOST_CHECK(e == datatype<long double volatile[]>());
    BOOST_CHECK(e == datatype<long double volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_complex_float )
{
    const MPI_Datatype e = MPI_C_FLOAT_COMPLEX;
    typedef std::complex<float> complex_type;
    BOOST_CHECK(e == datatype< complex_type >());
    BOOST_CHECK(e == datatype< complex_type *>());
    BOOST_CHECK(e == datatype< complex_type &>());
    BOOST_CHECK(e == datatype< complex_type [3]>());
    BOOST_CHECK(e == datatype< complex_type  const>());
    BOOST_CHECK(e == datatype< complex_type  const*>());
    BOOST_CHECK(e == datatype< complex_type  const&>());
    BOOST_CHECK(e == datatype< complex_type  const[4]>());
    BOOST_CHECK(e == datatype< complex_type  volatile>());
    BOOST_CHECK(e == datatype< complex_type  volatile*>());
    BOOST_CHECK(e == datatype< complex_type  volatile&>());
    BOOST_CHECK(e == datatype< complex_type  volatile[]>());
    BOOST_CHECK(e == datatype< complex_type  volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_complex_double )
{
    const MPI_Datatype e = MPI_C_DOUBLE_COMPLEX;
    typedef std::complex<double> complex_type;
    BOOST_CHECK(e == datatype< complex_type >());
    BOOST_CHECK(e == datatype< complex_type *>());
    BOOST_CHECK(e == datatype< complex_type &>());
    BOOST_CHECK(e == datatype< complex_type [3]>());
    BOOST_CHECK(e == datatype< complex_type  const>());
    BOOST_CHECK(e == datatype< complex_type  const*>());
    BOOST_CHECK(e == datatype< complex_type  const&>());
    BOOST_CHECK(e == datatype< complex_type  const[4]>());
    BOOST_CHECK(e == datatype< complex_type  volatile>());
    BOOST_CHECK(e == datatype< complex_type  volatile*>());
    BOOST_CHECK(e == datatype< complex_type  volatile&>());
    BOOST_CHECK(e == datatype< complex_type  volatile[]>());
    BOOST_CHECK(e == datatype< complex_type  volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_complex_long_double )
{
    const MPI_Datatype e = MPI_C_LONG_DOUBLE_COMPLEX;
    typedef std::complex<long double> complex_type;
    BOOST_CHECK(e == datatype< complex_type >());
    BOOST_CHECK(e == datatype< complex_type *>());
    BOOST_CHECK(e == datatype< complex_type &>());
    BOOST_CHECK(e == datatype< complex_type [3]>());
    BOOST_CHECK(e == datatype< complex_type  const>());
    BOOST_CHECK(e == datatype< complex_type  const*>());
    BOOST_CHECK(e == datatype< complex_type  const&>());
    BOOST_CHECK(e == datatype< complex_type  const[4]>());
    BOOST_CHECK(e == datatype< complex_type  volatile>());
    BOOST_CHECK(e == datatype< complex_type  volatile*>());
    BOOST_CHECK(e == datatype< complex_type  volatile&>());
    BOOST_CHECK(e == datatype< complex_type  volatile[]>());
    BOOST_CHECK(e == datatype< complex_type  volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_wchar_t )
{
    const MPI_Datatype e = MPI_WCHAR;
    BOOST_CHECK(e == datatype<wchar_t>());
    BOOST_CHECK(e == datatype<wchar_t*>());
    BOOST_CHECK(e == datatype<wchar_t&>());
    BOOST_CHECK(e == datatype<wchar_t[3]>());
    BOOST_CHECK(e == datatype<wchar_t const>());
    BOOST_CHECK(e == datatype<wchar_t const*>());
    BOOST_CHECK(e == datatype<wchar_t const&>());
    BOOST_CHECK(e == datatype<wchar_t const[4]>());
    BOOST_CHECK(e == datatype<wchar_t volatile>());
    BOOST_CHECK(e == datatype<wchar_t volatile*>());
    BOOST_CHECK(e == datatype<wchar_t volatile&>());
    BOOST_CHECK(e == datatype<wchar_t volatile[]>());
    BOOST_CHECK(e == datatype<wchar_t volatile const>());
}

BOOST_AUTO_TEST_CASE( datatype_of_test )
{
    using suzerain::mpi::datatype_of;

    BOOST_CHECK_EQUAL(MPI_CHAR,   datatype_of('c') );
    BOOST_CHECK_EQUAL(MPI_INT,    datatype_of(5)   );
    BOOST_CHECK_EQUAL(MPI_LONG,   datatype_of(5l)  );
    BOOST_CHECK_EQUAL(MPI_FLOAT,  datatype_of(1.2f));
    BOOST_CHECK_EQUAL(MPI_DOUBLE, datatype_of(1.2) );

    int x = 5;
    int &y = x;
    const int &z = x;
    BOOST_CHECK_EQUAL(MPI_INT, datatype_of(y));
    BOOST_CHECK_EQUAL(MPI_INT, datatype_of(z));
}
