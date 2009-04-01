#define BOOST_TEST_MODULE $Id$

#include "config.h"

#include <boost/test/included/unit_test.hpp>

#include "pencil_grid.h"
#include "utility.h"

BOOST_AUTO_TEST_CASE( declare_pointer )
{
  pecos::suzerain::pencil_grid<> *pg = NULL;
}

BOOST_AUTO_TEST_CASE( constructor_int )
{
  const int proc_dims[] = { 4, 2 };
  const int nx = 16, ny = 16, nz = 16;
  const int pstart[] = { 0, 0,  0};
  const int pend[]   = {15, 6,  3};
  const int psize[]  = {16, 7,  4};
  const int wstart[] = { 0, 0,  0};
  const int wend[]   = { 6, 3, 15};
  const int wsize[]  = { 7, 4, 16};

  pecos::suzerain::pencil_grid<int> pg(
      proc_dims, nx, ny, nz, pstart, pend, psize, wstart, wend, wsize);

  BOOST_CHECK_EQUAL( pg.pg1, proc_dims[0] );
  BOOST_CHECK_EQUAL( pg.pg2, proc_dims[1] );

  BOOST_CHECK_EQUAL( pg.nx, nx );
  BOOST_CHECK_EQUAL( pg.ny, ny );
  BOOST_CHECK_EQUAL( pg.nz, nz );

  BOOST_CHECK_EQUAL( pg.pstart_x, pstart[0] );
  BOOST_CHECK_EQUAL( pg.pstart_y, pstart[1] );
  BOOST_CHECK_EQUAL( pg.pstart_z, pstart[2] );
  BOOST_CHECK_EQUAL( pg.pend_x,   pend[0]   );
  BOOST_CHECK_EQUAL( pg.pend_y,   pend[1]   );
  BOOST_CHECK_EQUAL( pg.pend_z,   pend[2]   );
  BOOST_CHECK_EQUAL( pg.psize_x,  psize[0]  );
  BOOST_CHECK_EQUAL( pg.psize_y,  psize[1]  );
  BOOST_CHECK_EQUAL( pg.psize_z,  psize[2]  );

  BOOST_CHECK_EQUAL( pg.wstart_x, wstart[0] );
  BOOST_CHECK_EQUAL( pg.wstart_y, wstart[1] );
  BOOST_CHECK_EQUAL( pg.wstart_z, wstart[2] );
  BOOST_CHECK_EQUAL( pg.wend_x,   wend[0]   );
  BOOST_CHECK_EQUAL( pg.wend_y,   wend[1]   );
  BOOST_CHECK_EQUAL( pg.wend_z,   wend[2]   );
  BOOST_CHECK_EQUAL( pg.wsize_x,  wsize[0]  );
  BOOST_CHECK_EQUAL( pg.wsize_y,  wsize[1]  );
  BOOST_CHECK_EQUAL( pg.wsize_z,  wsize[2]  );

}

BOOST_AUTO_TEST_CASE( constructor_long )
{
  const long proc_dims[] = { 4, 2 };
  const long nx = 16, ny = 16, nz = 16;
  const long pstart[] = { 0, 0,  0};
  const long pend[]   = {15, 3,  3};
  const long psize[]  = {16, 4,  4};
  const long wstart[] = { 0, 0,  0};
  const long wend[]   = { 3, 3, 15};
  const long wsize[]  = { 4, 4, 16};

  pecos::suzerain::pencil_grid<long> pg(
      proc_dims, nx, ny, nz, pstart, pend, psize, wstart, wend, wsize);

  BOOST_CHECK_EQUAL( pg.pg1, proc_dims[0] );
  BOOST_CHECK_EQUAL( pg.pg2, proc_dims[1] );

  BOOST_CHECK_EQUAL( pg.nx, nx );
  BOOST_CHECK_EQUAL( pg.ny, ny );
  BOOST_CHECK_EQUAL( pg.nz, nz );

  BOOST_CHECK_EQUAL( pg.pstart_x, pstart[0] );
  BOOST_CHECK_EQUAL( pg.pstart_y, pstart[1] );
  BOOST_CHECK_EQUAL( pg.pstart_z, pstart[2] );
  BOOST_CHECK_EQUAL( pg.pend_x,   pend[0]   );
  BOOST_CHECK_EQUAL( pg.pend_y,   pend[1]   );
  BOOST_CHECK_EQUAL( pg.pend_z,   pend[2]   );
  BOOST_CHECK_EQUAL( pg.psize_x,  psize[0]  );
  BOOST_CHECK_EQUAL( pg.psize_y,  psize[1]  );
  BOOST_CHECK_EQUAL( pg.psize_z,  psize[2]  );

  BOOST_CHECK_EQUAL( pg.wstart_x, wstart[0] );
  BOOST_CHECK_EQUAL( pg.wstart_y, wstart[1] );
  BOOST_CHECK_EQUAL( pg.wstart_z, wstart[2] );
  BOOST_CHECK_EQUAL( pg.wend_x,   wend[0]   );
  BOOST_CHECK_EQUAL( pg.wend_y,   wend[1]   );
  BOOST_CHECK_EQUAL( pg.wend_z,   wend[2]   );
  BOOST_CHECK_EQUAL( pg.wsize_x,  wsize[0]  );
  BOOST_CHECK_EQUAL( pg.wsize_y,  wsize[1]  );
  BOOST_CHECK_EQUAL( pg.wsize_z,  wsize[2]  );

}
