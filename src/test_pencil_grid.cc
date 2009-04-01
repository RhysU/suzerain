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

  pecos::suzerain::pencil_grid<int> pg(proc_dims, nx, ny, nz);

  BOOST_CHECK_EQUAL( pg.pg1, proc_dims[0] );
  BOOST_CHECK_EQUAL( pg.pg2, proc_dims[1] );

  BOOST_CHECK_EQUAL( pg.nx, nx );
  BOOST_CHECK_EQUAL( pg.ny, ny );
  BOOST_CHECK_EQUAL( pg.nz, nz );
}

BOOST_AUTO_TEST_CASE( constructor_long )
{
  const long proc_dims[] = { 4, 2 };
  const long nx = 16, ny = 16, nz = 16;

  pecos::suzerain::pencil_grid<long> pg(proc_dims, nx, ny, nz);

  BOOST_CHECK_EQUAL( pg.pg1, proc_dims[0] );
  BOOST_CHECK_EQUAL( pg.pg2, proc_dims[1] );

  BOOST_CHECK_EQUAL( pg.nx, nx );
  BOOST_CHECK_EQUAL( pg.ny, ny );
  BOOST_CHECK_EQUAL( pg.nz, nz );

}
