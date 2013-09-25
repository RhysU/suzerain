//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
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

#include <suzerain/math.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

#include "test_tools.hpp"

typedef boost::mpl::list< double, float > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( minnan, T, test_types )
{
    using boost::math::isnan;
    using std::numeric_limits;
    using suzerain::math::minnan;

    const T a   = 1;
    const T b   = 2;
    const T nan = std::numeric_limits<T>::quiet_NaN();

    BOOST_CHECK_EQUAL(a, minnan(a,b));
    BOOST_CHECK_EQUAL(a, minnan(b,a));
    BOOST_CHECK((isnan)(minnan(nan,a)));
    BOOST_CHECK((isnan)(minnan(nan,b)));
    BOOST_CHECK((isnan)(minnan(a,nan)));
    BOOST_CHECK((isnan)(minnan(b,nan)));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( maxnan, T, test_types )
{
    using boost::math::isnan;
    using std::numeric_limits;
    using suzerain::math::maxnan;

    const T a   = 2;
    const T b   = 1;
    const T nan = std::numeric_limits<T>::quiet_NaN();

    BOOST_CHECK_EQUAL(a, maxnan(a,b));
    BOOST_CHECK_EQUAL(a, maxnan(b,a));
    BOOST_CHECK((isnan)(maxnan(nan,a)));
    BOOST_CHECK((isnan)(maxnan(nan,b)));
    BOOST_CHECK((isnan)(maxnan(a,nan)));
    BOOST_CHECK((isnan)(maxnan(b,nan)));
}

BOOST_AUTO_TEST_CASE( integer_power )
{
    using suzerain::math::integer_power;

    BOOST_CHECK_EQUAL( 1u, integer_power(2u, 0u));
    BOOST_CHECK_EQUAL( 2u, integer_power(2u, 1u));
    BOOST_CHECK_EQUAL( 4u, integer_power(2u, 2u));
    BOOST_CHECK_EQUAL( 8u, integer_power(2u, 3u));
    BOOST_CHECK_EQUAL(16u, integer_power(2u, 4u));

    BOOST_CHECK_EQUAL( 1ul, integer_power(2ul, 0ul));
    BOOST_CHECK_EQUAL( 2ul, integer_power(2ul, 1ul));
    BOOST_CHECK_EQUAL( 4ul, integer_power(2ul, 2ul));
    BOOST_CHECK_EQUAL( 8ul, integer_power(2ul, 3ul));
    BOOST_CHECK_EQUAL(16ul, integer_power(2ul, 4ul));

    BOOST_CHECK_EQUAL( 1.0, integer_power(2.0, 0));
    BOOST_CHECK_EQUAL( 2.0, integer_power(2.0, 1));
    BOOST_CHECK_EQUAL( 4.0, integer_power(2.0, 2));
    BOOST_CHECK_EQUAL( 8.0, integer_power(2.0, 3));
    BOOST_CHECK_EQUAL(16.0, integer_power(2.0, 4));

    BOOST_CHECK_EQUAL( 1.0f, integer_power(2.0f, 0));
    BOOST_CHECK_EQUAL( 2.0f, integer_power(2.0f, 1));
    BOOST_CHECK_EQUAL( 4.0f, integer_power(2.0f, 2));
    BOOST_CHECK_EQUAL( 8.0f, integer_power(2.0f, 3));
    BOOST_CHECK_EQUAL(16.0f, integer_power(2.0f, 4));

    BOOST_CHECK_EQUAL( 1.0,    integer_power(2.0, -0));
    BOOST_CHECK_EQUAL( 0.5,    integer_power(2.0, -1));
    BOOST_CHECK_EQUAL( 0.25,   integer_power(2.0, -2));
    BOOST_CHECK_EQUAL( 0.125,  integer_power(2.0, -3));
    BOOST_CHECK_EQUAL( 0.0625, integer_power(2.0, -4));

    BOOST_CHECK_EQUAL( 1.0f,    integer_power(2.0f, -0));
    BOOST_CHECK_EQUAL( 0.5f,    integer_power(2.0f, -1));
    BOOST_CHECK_EQUAL( 0.25f,   integer_power(2.0f, -2));
    BOOST_CHECK_EQUAL( 0.125f,  integer_power(2.0f, -3));
    BOOST_CHECK_EQUAL( 0.0625f, integer_power(2.0f, -4));

    BOOST_CHECK_EQUAL( 1.0, integer_power(0.0, 0));
    BOOST_CHECK_EQUAL( 0.0, integer_power(0.0, 1));

    BOOST_CHECK_EQUAL( 1.0f, integer_power(0.0f, 0));
    BOOST_CHECK_EQUAL( 0.0f, integer_power(0.0f, 1));
}

BOOST_AUTO_TEST_CASE( fixed_integer_power )
{
    using suzerain::math::fixed_integer_power;

    BOOST_CHECK_EQUAL( 1u, fixed_integer_power<0>(2u));
    BOOST_CHECK_EQUAL( 2u, fixed_integer_power<1>(2u));
    BOOST_CHECK_EQUAL( 4u, fixed_integer_power<2>(2u));
    BOOST_CHECK_EQUAL( 8u, fixed_integer_power<3>(2u));
    BOOST_CHECK_EQUAL(16u, fixed_integer_power<4>(2u));

    BOOST_CHECK_EQUAL( 1ul, fixed_integer_power<0>(2ul));
    BOOST_CHECK_EQUAL( 2ul, fixed_integer_power<1>(2ul));
    BOOST_CHECK_EQUAL( 4ul, fixed_integer_power<2>(2ul));
    BOOST_CHECK_EQUAL( 8ul, fixed_integer_power<3>(2ul));
    BOOST_CHECK_EQUAL(16ul, fixed_integer_power<4>(2ul));

    BOOST_CHECK_EQUAL( 1.0, fixed_integer_power<0>(2.0));
    BOOST_CHECK_EQUAL( 2.0, fixed_integer_power<1>(2.0));
    BOOST_CHECK_EQUAL( 4.0, fixed_integer_power<2>(2.0));
    BOOST_CHECK_EQUAL( 8.0, fixed_integer_power<3>(2.0));
    BOOST_CHECK_EQUAL(16.0, fixed_integer_power<4>(2.0));

    BOOST_CHECK_EQUAL( 1.0f, fixed_integer_power<0>(2.0f));
    BOOST_CHECK_EQUAL( 2.0f, fixed_integer_power<1>(2.0f));
    BOOST_CHECK_EQUAL( 4.0f, fixed_integer_power<2>(2.0f));
    BOOST_CHECK_EQUAL( 8.0f, fixed_integer_power<3>(2.0f));
    BOOST_CHECK_EQUAL(16.0f, fixed_integer_power<4>(2.0f));

    BOOST_CHECK_EQUAL( 1.0, fixed_integer_power<0>(0.0));
    BOOST_CHECK_EQUAL( 0.0, fixed_integer_power<1>(0.0));

    BOOST_CHECK_EQUAL( 1.0f, fixed_integer_power<0>(0.0f));
    BOOST_CHECK_EQUAL( 0.0f, fixed_integer_power<1>(0.0f));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( linspace, T, test_types )
{
    using suzerain::math::linspace;

    BOOST_CHECK_THROW(linspace(T(0),T(10),0,(T*)NULL), std::invalid_argument);

    {
        T result;

        BOOST_CHECK_EQUAL(linspace(T(567),T(567),1,&result), &result+1);
        BOOST_CHECK_EQUAL(result, T(567));

        BOOST_CHECK_EQUAL(linspace(T(-123),T(-123),1,&result), &result+1);
        BOOST_CHECK_EQUAL(result, T(-123));

        BOOST_CHECK_THROW(
                linspace(T(0),T(10),1,(T *)NULL), std::invalid_argument);
    }

    {
        T result[2];

        // Strictly positive, increasing
        BOOST_CHECK_EQUAL(linspace(T(234),T(567),2,result),result+2);
        BOOST_CHECK_EQUAL(result[0],T(234));
        BOOST_CHECK_EQUAL(result[1],T(567));

        // Spans zero, increasing
        BOOST_CHECK_EQUAL(linspace(T(-234),T(567),2,result),result+2);
        BOOST_CHECK_EQUAL(result[0],T(-234));
        BOOST_CHECK_EQUAL(result[1],T( 567));

        // Strictly negative, increasing
        BOOST_CHECK_EQUAL(linspace(T(-567),T(-234),2,result),result+2);
        BOOST_CHECK_EQUAL(result[0],T(-567));
        BOOST_CHECK_EQUAL(result[1],T(-234));

        // Strictly negative, decreasing
        BOOST_CHECK_EQUAL(linspace(T(-234),T(-567),2,result),result+2);
        BOOST_CHECK_EQUAL(result[0],T(-234));
        BOOST_CHECK_EQUAL(result[1],T(-567));
    }

    {
        const T close_enough = std::numeric_limits<T>::epsilon();
        T result[3];

        // Strictly positive, increasing
        BOOST_CHECK_EQUAL(linspace(T(234),T(567),3,result),result+3);
        BOOST_CHECK_EQUAL(result[0],T(234));
        BOOST_CHECK_CLOSE(result[1],T(234+567)/2,close_enough);
        BOOST_CHECK_EQUAL(result[2],T(567));

        // Spans zero, increasing
        BOOST_CHECK_EQUAL(linspace(T(-234),T(567),3,result),result+3);
        BOOST_CHECK_EQUAL(result[0],T(-234));
        BOOST_CHECK_CLOSE(result[1],T(-234+567)/2,close_enough);
        BOOST_CHECK_EQUAL(result[2],T( 567));

        // Strictly negative, increasing
        BOOST_CHECK_EQUAL(linspace(T(-567),T(-234),3,result),result+3);
        BOOST_CHECK_EQUAL(result[0],T(-567));
        BOOST_CHECK_CLOSE(result[1],T(-567-234)/2,close_enough);
        BOOST_CHECK_EQUAL(result[2],T(-234));

        // Strictly negative, decreasing
        BOOST_CHECK_EQUAL(linspace(T(-234),T(-567),3,result),result+3);
        BOOST_CHECK_EQUAL(result[0],T(-234));
        BOOST_CHECK_CLOSE(result[1],T(-234-567)/2,close_enough);
        BOOST_CHECK_EQUAL(result[2],T(-567));
    }

    {
        const T close_enough = std::numeric_limits<T>::epsilon();
        T result[5];
        BOOST_CHECK_EQUAL(linspace(T(1),T(5),5,result),result+5);
        BOOST_CHECK_EQUAL(result[0],T(1));
        BOOST_CHECK_CLOSE(result[1],T(2),close_enough);
        BOOST_CHECK_CLOSE(result[2],T(3),close_enough);
        BOOST_CHECK_CLOSE(result[3],T(4),close_enough);
        BOOST_CHECK_EQUAL(result[4],T(5));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( logspace, T, test_types )
{
    using suzerain::math::logspace;

    BOOST_CHECK_THROW(logspace(T(0),T(10),0,(T*)NULL), std::invalid_argument);

    {
        T result;

        BOOST_CHECK_EQUAL(logspace(T(5),T(5),1,&result,T(5)), &result+1);
        BOOST_CHECK_EQUAL(result, std::pow(T(5),T(5)));

        BOOST_CHECK_EQUAL(logspace(T(-123),T(-123),1,&result,T(7)), &result+1);
        BOOST_CHECK_EQUAL(result, std::pow(T(7),T(-123)));

        BOOST_CHECK_THROW(
                logspace(T(0),T(10),1,(T *)NULL), std::invalid_argument);
    }

    {
        T result[2];

        // Strictly positive, increasing
        BOOST_CHECK_EQUAL(logspace(T(2),T(5),2,result),result+2);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(10),T(2)));
        BOOST_CHECK_EQUAL(result[1],std::pow(T(10),T(5)));

        // Spans zero, increasing
        BOOST_CHECK_EQUAL(logspace(T(-2),T(5),2,result,T(8)),result+2);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(8),T(-2)));
        BOOST_CHECK_EQUAL(result[1],std::pow(T(8),T( 5)));

        // Strictly negative, increasing
        BOOST_CHECK_EQUAL(logspace(T(-5),T(-2),2,result,T(3)),result+2);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(3),T(-5)));
        BOOST_CHECK_EQUAL(result[1],std::pow(T(3),T(-2)));

        // Strictly negative, decreasing
        BOOST_CHECK_EQUAL(logspace(T(-2),T(-5),2,result),result+2);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(10),T(-2)));
        BOOST_CHECK_EQUAL(result[1],std::pow(T(10),T(-5)));
    }

    {
        const T close_enough = std::numeric_limits<T>::epsilon()*125;
        T result[3];

        // Strictly positive, increasing
        BOOST_CHECK_EQUAL(logspace(T(2),T(5),3,result),result+3);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(10),T(2)));
        BOOST_CHECK_CLOSE(result[1],std::pow(T(10),T(2+5)/2),close_enough);
        BOOST_CHECK_EQUAL(result[2],std::pow(T(10),T(5)));

        // Spans zero, increasing
        BOOST_CHECK_EQUAL(logspace(T(-2),T(5),3,result,T(6)),result+3);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(6),T(-2)));
        BOOST_CHECK_CLOSE(result[1],std::pow(T(6),T(-2+5)/2),close_enough);
        BOOST_CHECK_EQUAL(result[2],std::pow(T(6),T( 5)));

        // Strictly negative, increasing
        BOOST_CHECK_EQUAL(logspace(T(-5),T(-2),3,result),result+3);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(10),T(-5)));
        BOOST_CHECK_CLOSE(result[1],std::pow(T(10),T(-5-2)/2),close_enough);
        BOOST_CHECK_EQUAL(result[2],std::pow(T(10),T(-2)));

        // Strictly negative, decreasing
        BOOST_CHECK_EQUAL(logspace(T(-2),T(-5),3,result,T(4)),result+3);
        BOOST_CHECK_EQUAL(result[0],std::pow(T(4),T(-2)));
        BOOST_CHECK_CLOSE(result[1],std::pow(T(4),T(-2-5)/2),close_enough);
        BOOST_CHECK_EQUAL(result[2],std::pow(T(4),T(-5)));
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( stretchspace, T, test_types )
{
    using suzerain::math::stretchspace;

    BOOST_CHECK_THROW(
            stretchspace(T(0),T(10),0,T(1),(T*)NULL), std::invalid_argument);
    BOOST_CHECK_THROW(
            stretchspace(T(0),T(10),1,T(1),(T*)NULL), std::invalid_argument);
    BOOST_CHECK_THROW(
            stretchspace(T(0),T(10),2,T(1),(T*)NULL), std::invalid_argument);

    {
        const std::size_t N = 3;
        T result[N];
        const T close_enough = std::numeric_limits<T>::epsilon()*1000;

        // Strictly positive, increasing, alpha = 1
        BOOST_CHECK_EQUAL(stretchspace(T(1),T(5),N,T(1),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(1));
        BOOST_CHECK_EQUAL(result[1],T(3));
        BOOST_CHECK_EQUAL(result[2],T(5));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(1),close_enough);

        // Spans zero, increasing, alpha = 1
        BOOST_CHECK_EQUAL(stretchspace(T(-1),T(3),N,T(1),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-1));
        BOOST_CHECK_EQUAL(result[1],T( 1));
        BOOST_CHECK_EQUAL(result[2],T( 3));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(1),close_enough);

        // Strictly negative, increasing, alpha = 1
        BOOST_CHECK_EQUAL(stretchspace(T(-5),T(-1),N,T(1),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-5));
        BOOST_CHECK_EQUAL(result[1],T(-3));
        BOOST_CHECK_EQUAL(result[2],T(-1));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(1),close_enough);

        // Strictly negative, decreasing, alpha = 1
        BOOST_CHECK_EQUAL(stretchspace(T(-1),T(-5),N,T(1),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-1));
        BOOST_CHECK_EQUAL(result[1],T(-3));
        BOOST_CHECK_EQUAL(result[2],T(-5));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(1),close_enough);

        // Strictly positive, increasing, alpha = 2
        BOOST_CHECK_EQUAL(stretchspace(T(1),T(5),N,T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(1));
        BOOST_CHECK_CLOSE(result[1],T(7)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[2],T(5));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(2),close_enough);

        // Spans zero, increasing, alpha = 2
        BOOST_CHECK_EQUAL(stretchspace(T(-1),T(3),N,T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-1));
        BOOST_CHECK_CLOSE(result[1],T( 1)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[2],T( 3));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(2),close_enough);

        // Strictly negative, increasing, alpha = 2
        BOOST_CHECK_EQUAL(stretchspace(T(-5),T(-1),N,T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-5));
        BOOST_CHECK_CLOSE(result[1],T(-11)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[2],T(-1));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(2),close_enough);

        // Strictly negative, decreasing, alpha = 2
        BOOST_CHECK_EQUAL(stretchspace(T(-1),T(-5),N,T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-1));
        BOOST_CHECK_CLOSE(result[1],T(-7)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[2],T(-5));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),T(2),close_enough);

        // Strictly positive, increasing, alpha = 1/2
        BOOST_CHECK_EQUAL(stretchspace(T(1),T(5),N,1/T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(1));
        BOOST_CHECK_CLOSE(result[1],T(11)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[2],T(5));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),1/T(2),close_enough);

        // Spans zero, increasing, alpha = 1/2
        BOOST_CHECK_EQUAL(stretchspace(T(-1),T(3),N,1/T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-1));
        BOOST_CHECK_CLOSE(result[1],T( 5)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[2],T( 3));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),1/T(2),close_enough);

        // Strictly negative, increasing, alpha = 1/2
        BOOST_CHECK_EQUAL(stretchspace(T(-5),T(-1),N,1/T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-5));
        BOOST_CHECK_CLOSE(result[1],T(-7)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[2],T(-1));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),1/T(2),close_enough);

        // Strictly negative, decreasing, alpha = 1/2
        BOOST_CHECK_EQUAL(stretchspace(T(-1),T(-5),N,1/T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(-1));
        BOOST_CHECK_EQUAL(result[1],T(-11)/T(3));
        BOOST_CHECK_EQUAL(result[2],T(-5));
        BOOST_CHECK_CLOSE(
            (result[2]-result[1])/(result[1]-result[0]),1/T(2),close_enough);
    }

    {
        const std::size_t N = 5;
        T result[N];
        const T close_enough = std::numeric_limits<T>::epsilon()*1000;

        BOOST_CHECK_EQUAL(stretchspace(T(1),T(5),N,T(2),result),result+N);
        BOOST_CHECK_EQUAL(result[0],T(1));
        BOOST_CHECK_CLOSE(result[1],T(5)/T(3),close_enough);
        BOOST_CHECK_CLOSE(result[2],T(23)/T(9),close_enough);
        BOOST_CHECK_CLOSE(result[3],T(11)/T(3),close_enough);
        BOOST_CHECK_EQUAL(result[4],T(5));
        BOOST_CHECK_CLOSE(
            (result[N-1]-result[N-2])/(result[1]-result[0]),T(2),close_enough);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( bump_classic, T, test_types )
{
    using namespace suzerain::math::bump;
    const T eps   = std::sqrt(std::numeric_limits<T>::epsilon());

    BOOST_CHECK_EQUAL(classic<T>(-1.1), 0);
    BOOST_CHECK_EQUAL(classic<T>(-1.0), 0);
    BOOST_CHECK_GT   (classic<T>(-0.9), 0);
    BOOST_CHECK_GT   (classic<T>(-eps), 0);
    BOOST_CHECK_GT   (classic<T>( 0.0), classic<T>(-eps));
    BOOST_CHECK_GT   (classic<T>( 0.0), classic<T>( eps));
    BOOST_CHECK_GT   (classic<T>( eps), 0);
    BOOST_CHECK_GT   (classic<T>(+0.9), 0);
    BOOST_CHECK_EQUAL(classic<T>(+1.0), 0);
    BOOST_CHECK_EQUAL(classic<T>(+1.1), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( bump_scaled, T, test_types )
{
    using namespace suzerain::math::bump;
    const T eps   = std::sqrt(std::numeric_limits<T>::epsilon());

    // One argument
    BOOST_CHECK_EQUAL(scaled<T>(-1.1), 0);
    BOOST_CHECK_EQUAL(scaled<T>(-1.0), 0);
    BOOST_CHECK_GT   (scaled<T>(-0.9), 0);
    BOOST_CHECK_GT   (scaled<T>(-eps), 0);
    BOOST_CHECK_GT   (scaled<T>( 0.0), scaled<T>(-eps));
    BOOST_CHECK_EQUAL(scaled<T>( 0.0), 1);
    BOOST_CHECK_GT   (scaled<T>( 0.0), scaled<T>( eps));
    BOOST_CHECK_GT   (scaled<T>( eps), 0);
    BOOST_CHECK_GT   (scaled<T>(+0.9), 0);
    BOOST_CHECK_EQUAL(scaled<T>(+1.0), 0);
    BOOST_CHECK_EQUAL(scaled<T>(+1.1), 0);

    // Two arguments
    BOOST_CHECK_EQUAL(scaled<T>(-1.1, 2), 0);
    BOOST_CHECK_EQUAL(scaled<T>(-1.0, 2), 0);
    BOOST_CHECK_GT   (scaled<T>(-0.9, 2), 0);
    BOOST_CHECK_GT   (scaled<T>(-eps, 2), 0);
    BOOST_CHECK_GT   (scaled<T>( 0.0, 2), scaled<T>(-eps, 2));
    BOOST_CHECK_EQUAL(scaled<T>( 0.0, 2), 2);
    BOOST_CHECK_GT   (scaled<T>( 0.0, 2), scaled<T>( eps, 2));
    BOOST_CHECK_GT   (scaled<T>( eps, 2), 0);
    BOOST_CHECK_GT   (scaled<T>(+0.9, 2), 0);
    BOOST_CHECK_EQUAL(scaled<T>(+1.0, 2), 0);
    BOOST_CHECK_EQUAL(scaled<T>(+1.1, 2), 0);

    // Three arguments
    BOOST_CHECK_EQUAL(scaled<T>(-1.1, 2, 1), 1);
    BOOST_CHECK_EQUAL(scaled<T>(-1.0, 2, 1), 1);
    BOOST_CHECK_GT   (scaled<T>(-0.9, 2, 1), 1);
    BOOST_CHECK_GT   (scaled<T>(-eps, 2, 1), 1);
    BOOST_CHECK_GT   (scaled<T>( 0.0, 2, 1), scaled<T>(-eps, 2, 1));
    BOOST_CHECK_EQUAL(scaled<T>( 0.0, 2, 1), 2);
    BOOST_CHECK_GT   (scaled<T>( 0.0, 2, 1), scaled<T>( eps, 2, 1));
    BOOST_CHECK_GT   (scaled<T>( eps, 2, 1), 1);
    BOOST_CHECK_GT   (scaled<T>(+0.9, 2, 1), 1);
    BOOST_CHECK_EQUAL(scaled<T>(+1.0, 2, 1), 1);
    BOOST_CHECK_EQUAL(scaled<T>(+1.1, 2, 1), 1);

    // Four arguments
    BOOST_CHECK_EQUAL(scaled<T>(-1.1, 2, 1, 2), 1);
    BOOST_CHECK_EQUAL(scaled<T>(-1.0, 2, 1, 2), 1);
    BOOST_CHECK_GT   (scaled<T>(-0.9, 2, 1, 2), 1);
    BOOST_CHECK_GT   (scaled<T>(-eps, 2, 1, 2), 1);
    BOOST_CHECK_GT   (scaled<T>( 0.0, 2, 1, 2), scaled<T>(-eps, 2, 1, 2));
    BOOST_CHECK_EQUAL(scaled<T>( 0.0, 2, 1, 2), 2);
    BOOST_CHECK_GT   (scaled<T>( 0.0, 2, 1, 2), scaled<T>( eps, 2, 1, 2));
    BOOST_CHECK_GT   (scaled<T>( eps, 2, 1, 2), 1);
    BOOST_CHECK_GT   (scaled<T>(+0.9, 2, 1, 2), 1);
    BOOST_CHECK_EQUAL(scaled<T>(+1.0, 2, 1, 2), 1);
    BOOST_CHECK_EQUAL(scaled<T>(+1.1, 2, 1, 2), 1);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( bump_shifted, T, test_types )
{
    using namespace suzerain::math::bump;
    const T l = 1, r = 3, m = (l + r) / 2, eps = (r - l) / 10;

    // Three arguments
    BOOST_CHECK_EQUAL(shifted<T>(l-eps, l, r), 0);
    BOOST_CHECK_EQUAL(shifted<T>(l    , l, r), 0);
    BOOST_CHECK_GT   (shifted<T>(l+eps, l, r), 0);
    BOOST_CHECK_GT   (shifted<T>(m-eps, l, r), 0);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r), shifted<T>(m-eps, l, r));
    BOOST_CHECK_EQUAL(shifted<T>(m    , l, r), 1);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r), shifted<T>(m+eps, l, r));
    BOOST_CHECK_GT   (shifted<T>(m+eps, l, r), 0);
    BOOST_CHECK_GT   (shifted<T>(r-eps, l, r), 0);
    BOOST_CHECK_EQUAL(shifted<T>(r    , l, r), 0);
    BOOST_CHECK_EQUAL(shifted<T>(r+eps, l, r), 0);

    // Four arguments
    BOOST_CHECK_EQUAL(shifted<T>(l-eps, l, r, 2), 0);
    BOOST_CHECK_EQUAL(shifted<T>(l    , l, r, 2), 0);
    BOOST_CHECK_GT   (shifted<T>(l+eps, l, r, 2), 0);
    BOOST_CHECK_GT   (shifted<T>(m-eps, l, r, 2), 0);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r, 2), shifted<T>(m-eps, l, r, 2));
    BOOST_CHECK_EQUAL(shifted<T>(m    , l, r, 2), 2);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r, 2), shifted<T>(m+eps, l, r, 2));
    BOOST_CHECK_GT   (shifted<T>(m+eps, l, r, 2), 0);
    BOOST_CHECK_GT   (shifted<T>(r-eps, l, r, 2), 0);
    BOOST_CHECK_EQUAL(shifted<T>(r    , l, r, 2), 0);
    BOOST_CHECK_EQUAL(shifted<T>(r+eps, l, r, 2), 0);

    // Five arguments
    BOOST_CHECK_EQUAL(shifted<T>(l-eps, l, r, 2, 1), 1);
    BOOST_CHECK_EQUAL(shifted<T>(l    , l, r, 2, 1), 1);
    BOOST_CHECK_GT   (shifted<T>(l+eps, l, r, 2, 1), 1);
    BOOST_CHECK_GT   (shifted<T>(m-eps, l, r, 2, 1), 1);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r, 2, 1), shifted<T>(m-eps, l, r, 2, 1));
    BOOST_CHECK_EQUAL(shifted<T>(m    , l, r, 2, 1), 2);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r, 2, 1), shifted<T>(m+eps, l, r, 2, 1));
    BOOST_CHECK_GT   (shifted<T>(m+eps, l, r, 2, 1), 1);
    BOOST_CHECK_GT   (shifted<T>(r-eps, l, r, 2, 1), 1);
    BOOST_CHECK_EQUAL(shifted<T>(r    , l, r, 2, 1), 1);
    BOOST_CHECK_EQUAL(shifted<T>(r+eps, l, r, 2, 1), 1);

    // Six arguments
    BOOST_CHECK_EQUAL(shifted<T>(l-eps, l, r, 2, 1, 2), 1);
    BOOST_CHECK_EQUAL(shifted<T>(l    , l, r, 2, 1, 2), 1);
    BOOST_CHECK_GT   (shifted<T>(l+eps, l, r, 2, 1, 2), 1);
    BOOST_CHECK_GT   (shifted<T>(m-eps, l, r, 2, 1, 2), 1);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r, 2, 1, 2), shifted<T>(m-eps, l, r, 2, 1, 2));
    BOOST_CHECK_EQUAL(shifted<T>(m    , l, r, 2, 1, 2), 2);
    BOOST_CHECK_GT   (shifted<T>(m    , l, r, 2, 1, 2), shifted<T>(m+eps, l, r, 2, 1, 2));
    BOOST_CHECK_GT   (shifted<T>(m+eps, l, r, 2, 1, 2), 1);
    BOOST_CHECK_GT   (shifted<T>(r-eps, l, r, 2, 1, 2), 1);
    BOOST_CHECK_EQUAL(shifted<T>(r    , l, r, 2, 1, 2), 1);
    BOOST_CHECK_EQUAL(shifted<T>(r+eps, l, r, 2, 1, 2), 1);
}
