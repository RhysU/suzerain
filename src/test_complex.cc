#define BOOST_TEST_MODULE $Id$
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include <complex>
#include <log4cxx/logger.h>

log4cxx::LoggerPtr logger = log4cxx::Logger::getRootLogger();

BOOST_AUTO_TEST_CASE( check_complex_double_is_two_doubles )
{
  BOOST_REQUIRE_EQUAL(
    sizeof(std::complex<double>), 2*sizeof(double) );
}

//BOOST_AUTO_TEST_CASE( check_complex_double_is_two_doubles )
//{
//  BOOST_REQUIRE_EQUAL(
//    sizeof(std::complex<double>), 2*sizeof(double) );
//}
