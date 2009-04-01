#define BOOST_TEST_MODULE $Id$

#include "config.h"

#include <boost/test/included/unit_test.hpp>
#include <boost/type_traits.hpp>

#include "pencil.hpp"

BOOST_AUTO_TEST_CASE( declare_pointer )
{
  using namespace pecos::suzerain;

  pencil<> *p = NULL;
}

BOOST_AUTO_TEST_CASE( constructor )
{
  using namespace pecos::suzerain;

  const pencil<>::dim_type pstart[] = { 0, 0,  0};
  const pencil<>::dim_type psize[]  = {16, 7,  4};
  const pencil<>::dim_type wstart[] = { 0, 0,  0};
  const pencil<>::dim_type wsize[]  = { 7, 4, 16};

  pencil<> p(pstart, psize, wstart, wsize);

  BOOST_CHECK_EQUAL( p.pstart_x, pstart[0] );
  BOOST_CHECK_EQUAL( p.pstart_y, pstart[1] );
  BOOST_CHECK_EQUAL( p.pstart_z, pstart[2] );
  BOOST_CHECK_EQUAL( p.psize_x,  psize[0]  );
  BOOST_CHECK_EQUAL( p.psize_y,  psize[1]  );
  BOOST_CHECK_EQUAL( p.psize_z,  psize[2]  );

  BOOST_CHECK_EQUAL( p.wstart_x, wstart[0] );
  BOOST_CHECK_EQUAL( p.wstart_y, wstart[1] );
  BOOST_CHECK_EQUAL( p.wstart_z, wstart[2] );
  BOOST_CHECK_EQUAL( p.wsize_x,  wsize[0]  );
  BOOST_CHECK_EQUAL( p.wsize_y,  wsize[1]  );
  BOOST_CHECK_EQUAL( p.wsize_z,  wsize[2]  );
}

BOOST_AUTO_TEST_CASE( real_access )
{
  using namespace pecos::suzerain;

  const pencil<>::dim_type pstart[] = { 0, 0,  0};
  const pencil<>::dim_type psize[]  = { 2, 2,  0};
  const pencil<>::dim_type wstart[] = { 0, 0,  0};
  const pencil<>::dim_type wsize[]  = { 2, 1,  0};

  pencil<> p(pstart, psize, wstart, wsize);

  // Z, Y, X loop order
  for (pencil<>::size_type k = 0; k < p.psize_z; ++k) {
    for (pencil<>::size_type j = 0; j < p.psize_y; ++j) {
      for (pencil<>::size_type i = 0; i < p.psize_x; ++i) {
        p.p(i,j,k) = (i+1)*(j+1)*(k+1);
      }
    }
  }

  // X, Y, Z loop order
  for (pencil<>::size_type i = 0; i < p.psize_x; ++i) {
    for (pencil<>::size_type j = 0; j < p.psize_y; ++j) {
      for (pencil<>::size_type k = 0; k < p.psize_z; ++k) {
        BOOST_CHECK_EQUAL(p.p(i,j,k), (i+1)*(j+1)*(k+1));
      }
    }
  }

}
