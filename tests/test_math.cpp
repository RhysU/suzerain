#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <suzerain/math.hpp>
#include "test_tools.hpp"

#pragma warning(disable:383)

typedef boost::mpl::list< double, float > test_types;

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
        const T close_enough = std::numeric_limits<T>::epsilon()*100;
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
