#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/pencil_grid.hpp>
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( declare_pointer )
{
    pecos::suzerain::pencil_grid<> *pg = NULL;
}

BOOST_AUTO_TEST_CASE( constructor_default )
{

    using namespace pecos::suzerain;

    const pencil_grid<>::dim_type proc_dims[] = { 4, 2 };
    const pencil_grid<>::dim_type nx = 16, ny = 16, nz = 16;

    pencil_grid<> pg(proc_dims, nx, ny, nz);

    BOOST_CHECK_EQUAL( pg.pg1, proc_dims[0] );
    BOOST_CHECK_EQUAL( pg.pg2, proc_dims[1] );

    BOOST_CHECK_EQUAL( pg.nx, nx );
    BOOST_CHECK_EQUAL( pg.ny, ny );
    BOOST_CHECK_EQUAL( pg.nz, nz );
}

BOOST_AUTO_TEST_CASE( constructor_long )
{

    using namespace pecos::suzerain;

    const pencil_grid<long>::dim_type proc_dims[] = { 4, 2 };
    const pencil_grid<long>::dim_type nx = 16, ny = 16, nz = 16;

    pencil_grid<long> pg(proc_dims, nx, ny, nz);

    BOOST_CHECK_EQUAL( pg.pg1, proc_dims[0] );
    BOOST_CHECK_EQUAL( pg.pg2, proc_dims[1] );

    BOOST_CHECK_EQUAL( pg.nx, nx );
    BOOST_CHECK_EQUAL( pg.ny, ny );
    BOOST_CHECK_EQUAL( pg.nz, nz );

}
