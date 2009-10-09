#define BOOST_TEST_MODULE $Id$
#include <complex>
#include <functional>
#include <numeric>
#include <vector>
#include <boost/array.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/multi_array.hpp>
#include <suzerain/fftw_multi_array.hpp>

#include "test_tools.hpp"

/* DEBUG */
#include <iostream>
#include <iterator>

using namespace pecos::suzerain;

BOOST_AUTO_TEST_CASE( increment_1d_degenerate )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 1 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_1d_normal )
{
    const int n        = 1;
    int       index[n] = { 0 };
    const int shape[n] = { 3 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_first )
{
    const int           n     = 2;
    boost::array<int,n> index = {{ 0, 0 }};
    boost::array<int,n> shape = {{ 3, 1 }};

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_second )
{
    const int                   n     = 2;
    boost::array<int, n>        index = {{ 0, 0 }};
    boost::array<std::size_t,n> shape = {{ 1, 3 }};

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_degenerate_all )
{
    const int  n        = 2;
    int        index[n] = { 0, 0 };
    const long shape[n] = { 1, 1 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_2d_normal )
{
    const int    n        = 2;
    unsigned int index[n] = { 0, 0 };
    unsigned int shape[n] = { 2, 2 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_all )
{
    const int         n        = 3;
    signed int        index[n] = { 0, 0, 0 };
    const signed long shape[n] = { 1, 1, 1 };

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_degenerate_middle )
{
    const int                   n        = 3;
    long                        index[n] =  { 0, 0, 0 };
    boost::array<std::size_t,n> shape    = {{ 3, 1, 3 }};

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 2);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 2);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}

BOOST_AUTO_TEST_CASE( increment_3d_normal )
{
    const int n = 3;
    std::vector<short> index(n, 0);
    std::vector<int> shape(n, 2);

    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 0);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 0);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 0);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), true);
    BOOST_CHECK_EQUAL(index[0], 1);
    BOOST_CHECK_EQUAL(index[1], 1);
    BOOST_CHECK_EQUAL(index[2], 1);
    BOOST_CHECK_EQUAL(fftw_multi_array::increment<n>(index,shape), false);
}


BOOST_AUTO_TEST_CASE( multi_array_sample )
{
  // Create a 3D array that is 3 x 4 x 2
  typedef boost::multi_array<double, 3> array_type;
  typedef array_type::index index;
  array_type A(boost::extents[3][4][2]);

  // Assign values to the elements
  int values = 0;
  for(index i = 0; i != 3; ++i)
    for(index j = 0; j != 4; ++j)
      for(index k = 0; k != 2; ++k)
        A[i][j][k] = values++;

  // Verify values
  int verify = 0;
  for(index i = 0; i != 3; ++i)
    for(index j = 0; j != 4; ++j)
      for(index k = 0; k != 2; ++k)
        BOOST_CHECK_EQUAL(A[i][j][k], verify++);
}


BOOST_AUTO_TEST_CASE( in_progress )
{
    typedef boost::multi_array<std::complex<double>,2> array_type;
    array_type in(boost::extents[3][2]);
    array_type out(boost::extents[3][2]);

    fftw_multi_array::c2c_transform(0, in, out, FFTW_FORWARD, FFTW_PRESERVE_INPUT);
}
