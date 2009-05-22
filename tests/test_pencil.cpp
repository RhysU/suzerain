#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>

#include <suzerain/pencil.hpp>

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

    BOOST_CHECK_EQUAL( p.physical.start_x, pstart[0] );
    BOOST_CHECK_EQUAL( p.physical.start_y, pstart[1] );
    BOOST_CHECK_EQUAL( p.physical.start_z, pstart[2] );
    BOOST_CHECK_EQUAL( p.physical.end_x, pstart[0] + psize[0]);
    BOOST_CHECK_EQUAL( p.physical.end_y, pstart[1] + psize[1]);
    BOOST_CHECK_EQUAL( p.physical.end_z, pstart[2] + psize[2]);
    BOOST_CHECK_EQUAL( p.physical.size_x,  psize[0]  );
    BOOST_CHECK_EQUAL( p.physical.size_y,  psize[1]  );
    BOOST_CHECK_EQUAL( p.physical.size_z,  psize[2]  );

    BOOST_CHECK_EQUAL( p.wave.start_x, wstart[0] );
    BOOST_CHECK_EQUAL( p.wave.start_y, wstart[1] );
    BOOST_CHECK_EQUAL( p.wave.start_z, wstart[2] );
    BOOST_CHECK_EQUAL( p.wave.end_x, wstart[0] + wsize[0]);
    BOOST_CHECK_EQUAL( p.wave.end_y, wstart[1] + wsize[1]);
    BOOST_CHECK_EQUAL( p.wave.end_z, wstart[2] + wsize[2]);
    BOOST_CHECK_EQUAL( p.wave.size_x,  wsize[0]  );
    BOOST_CHECK_EQUAL( p.wave.size_y,  wsize[1]  );
    BOOST_CHECK_EQUAL( p.wave.size_z,  wsize[2]  );
}

BOOST_AUTO_TEST_CASE( storage_order )
{
    using namespace pecos::suzerain;

    const pencil<>::dim_type pstart[] = {  0,  0,  0};
    const pencil<>::dim_type psize[]  = { 11, 13, 17};
    const pencil<>::dim_type wstart[] = {  0,  0,  0};
    const pencil<>::dim_type wsize[]  = {  3,  5,  7};

    pencil<> p(pstart, psize, wstart, wsize);

    // Ensure memory locations are initialized for test below,
    // otherwise valgrind complains bitterly about uninitialized
    // locations being used for jumps
    p.physical(0,0,0) = 1;
    p.physical(1,0,0) = 1;
    p.physical(0,1,0) = 1;
    p.physical(0,0,1) = 1;
    p.wave(0,0,0) = 1;
    p.wave(1,0,0) = 1;
    p.wave(0,1,0) = 1;
    p.wave(0,0,1) = 1;

    // Check physical memory layout is (X,Z,Y) in column major storage
    BOOST_CHECK_EQUAL(
        &p.physical(0,0,0) + 1,                 &p.physical(1,0,0)); // x
    BOOST_CHECK_EQUAL(
        &p.physical(0,0,0) + psize[0]*psize[2], &p.physical(0,1,0)); // y
    BOOST_CHECK_EQUAL(
        &p.physical(0,0,0) + psize[0],          &p.physical(0,0,1)); // z

    // Check wavespace memory layout is (Y,Z,X) in column major storage
    BOOST_CHECK_EQUAL(
        &p.wave(0,0,0) + wsize[1]*wsize[2], &p.wave(1,0,0)); // x
    BOOST_CHECK_EQUAL(
        &p.wave(0,0,0) + 1,                 &p.wave(0,1,0)); // y
    BOOST_CHECK_EQUAL(
        &p.wave(0,0,0) + wsize[1],          &p.wave(0,0,1)); // z
}

BOOST_AUTO_TEST_CASE( real_access )
{

    using namespace pecos::suzerain;

    const pencil<>::dim_type pstart[] = { 0, 0,  0};
    const pencil<>::dim_type psize[]  = { 2, 2,  2};
    const pencil<>::dim_type wstart[] = { 0, 0,  0};
    const pencil<>::dim_type wsize[]  = { 2, 1,  1};

    pencil<> p(pstart, psize, wstart, wsize);

    // X, Z, Y loop order
    for (pencil<>::size_type i = 0; i < p.physical.size_x; ++i) {
        for (pencil<>::size_type k = 0; k < p.physical.size_z; ++k) {
            for (pencil<>::size_type j = 0; j < p.physical.size_y; ++j) {
                p.physical(i, j, k) = (i + 1) * (j + 1) * (k + 1);
            }
        }
    }

    // Y, Z, X loop order
    for (pencil<>::size_type j = 0; j < p.physical.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.physical.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.physical.size_x; ++i) {
                BOOST_CHECK_EQUAL(p.physical(i, j, k), (i + 1)*(j + 1)*(k + 1));
            }
        }
    }

    // Clear contents using iterator
    std::fill(p.physical.begin(), p.physical.end(), 0);
    // Check contents are clear
    for (pencil<>::size_type j = 0; j < p.physical.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.physical.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.physical.size_x; ++i) {
                BOOST_CHECK_EQUAL(p.physical(i, j, k), 0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( complex_access )
{

    using namespace pecos::suzerain;

    const pencil<>::dim_type pstart[] = { 1,  1,  1};
    const pencil<>::dim_type psize[]  = { 2,  3,  5};
    const pencil<>::dim_type wstart[] = { 2,  2,  2};
    const pencil<>::dim_type wsize[]  = { 7, 11, 13};

    pencil<> p(pstart, psize, wstart, wsize);

    // X, Z, Y loop order, assign complex value
    for (pencil<>::size_type i = 0; i < p.wave.size_x; ++i) {
        for (pencil<>::size_type k = 0; k < p.wave.size_z; ++k) {
            for (pencil<>::size_type j = 0; j < p.wave.size_y; ++j) {
                p.wave(i, j, k) = pencil<>::complex_type(
                    (i + 1) * (j + 1) * (k + 1),
                    (i - 1) * (j - 1) * (k - 1));
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (pencil<>::size_type j = 0; j < p.wave.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.wave.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.wave.size_x; ++i) {
                BOOST_CHECK_EQUAL(p.wave(i, j, k), pencil<>::complex_type(
                    (i + 1)*(j + 1)*(k + 1),
                    (i - 1)*(j - 1)*(k - 1)));
                BOOST_CHECK_EQUAL(p.wave.real(i, j, k), (i + 1)*(j + 1)*(k + 1));
                BOOST_CHECK_EQUAL(p.wave.imag(i, j, k), (i - 1)*(j - 1)*(k - 1));
            }
        }
    }

    // X, Z, Y loop order, assign real and imag values
    for (pencil<>::size_type i = 0; i < p.wave.size_x; ++i) {
        for (pencil<>::size_type k = 0; k < p.wave.size_z; ++k) {
            for (pencil<>::size_type j = 0; j < p.wave.size_y; ++j) {
                p.wave.real(i, j, k) = (i - 123) * (j - 123) * (k - 123);
                p.wave.imag(i, j, k) = (i + 123) * (j + 123) * (k + 123);
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (pencil<>::size_type j = 0; j < p.wave.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.wave.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.wave.size_x; ++i) {
                BOOST_CHECK_EQUAL(p.wave(i, j, k), pencil<>::complex_type(
                    (i - 123)*(j - 123)*(k - 123),
                    (i + 123)*(j + 123)*(k + 123)));
                BOOST_CHECK_EQUAL(
                    p.wave.real(i, j, k), (i - 123)*(j - 123)*(k - 123));
                BOOST_CHECK_EQUAL(
                    p.wave.imag(i, j, k), (i + 123)*(j + 123)*(k + 123));
            }
        }
    }

    // Clear contents using iterator
    std::fill(p.wave.begin(), p.wave.end(), pencil<>::complex_type(0));
    // Check contents are clear
    for (pencil<>::size_type j = 0; j < p.wave.size_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.wave.size_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.wave.size_x; ++i) {
                BOOST_CHECK_EQUAL(
                    p.wave(i, j, k), pencil<>::complex_type(0));
            }
        }
    }
}
