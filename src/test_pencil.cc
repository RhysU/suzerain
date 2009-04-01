#define BOOST_TEST_MODULE $Id$

#include "config.h"

#include <boost/test/included/unit_test.hpp>
#include <boost/type_traits.hpp>

#include "pencil.hpp"

BOOST_AUTO_TEST_CASE( declare_pointer )
{
  pecos::suzerain::pencil<> *p = NULL;
}

BOOST_AUTO_TEST_CASE( constructor )
{
  const int pstart[] = { 0, 0,  0};
  const int psize[]  = {16, 7,  4};
  const int wstart[] = { 0, 0,  0};
  const int wsize[]  = { 7, 4, 16};

  pecos::suzerain::pencil<> p(pstart, psize, wstart, wsize);

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
