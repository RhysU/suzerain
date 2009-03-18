#define BOOST_TEST_MODULE $Id$

#include "config.h"

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/foreach.hpp>
#include <complex>
#include <log4cxx/logger.h>
#include <typeinfo>

log4cxx::LoggerPtr logger = log4cxx::Logger::getRootLogger();

namespace ublas = boost::numeric::ublas;

BOOST_AUTO_TEST_CASE( check_complex_double_is_two_doubles )
{
  BOOST_REQUIRE_EQUAL( sizeof(std::complex<double>), 2*sizeof(double) );
}

BOOST_AUTO_TEST_CASE( shared_c_array )
{
  const std::size_t N = 6;
  double carray_double[N] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
  std::complex<double> *carray_complex
    = reinterpret_cast<std::complex<double> *>(carray_double);

  for (std::size_t i = 0; i < N/2; ++i)
    {
      BOOST_CHECK_EQUAL(carray_double[2*i],   carray_complex[i].real());
      BOOST_CHECK_EQUAL(carray_double[2*i+1], carray_complex[i].imag());
    }

  for (std::size_t i = 0; i < N/2; ++i)
    {
      carray_complex[i] += carray_complex[1];
    }

  for (std::size_t i = 0; i < N/2; ++i)
    {
      BOOST_CHECK_EQUAL(carray_double[2*i],   carray_complex[i].real());
      BOOST_CHECK_EQUAL(carray_double[2*i+1], carray_complex[i].imag());
    }
}

BOOST_AUTO_TEST_CASE( shallow_array_adaptor )
{
  const std::size_t N = 3;
  double raw[N];

  ublas::shallow_array_adaptor<double> adaptor(N, raw);
  ublas::vector<double, ublas::shallow_array_adaptor<double> > vec(N, adaptor);

  BOOST_CHECK_EQUAL( vec.size(), N );

  std::fill(&raw[0], &raw[N], 1.0);
  BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());

  std::fill(vec.begin(), vec.end(), 2.0);
  BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());
}

BOOST_AUTO_TEST_CASE( shallow_array_adaptor_complex )
{
  const std::size_t N = 3;
  std::complex<double> raw[N];

  ublas::shallow_array_adaptor<std::complex<double> > adaptor(N, raw);
  ublas::vector<std::complex<double>,
    ublas::shallow_array_adaptor<std::complex<double> > > vec(N, adaptor);

  BOOST_CHECK_EQUAL( vec.size(), N );

  std::fill(&raw[0], &raw[N], std::complex<double>(1.0,-1.0));
  BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());

  std::fill(vec.begin(), vec.end(), std::complex<double>(2.0,-2.0));
  BOOST_CHECK_EQUAL_COLLECTIONS(&raw[0], &raw[N], vec.begin(), vec.end());
}
