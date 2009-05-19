#define BOOST_TEST_MODULE $Id$

#include "config.h"

#include <algorithm>
#include <boost/test/included/unit_test.hpp>
#include <boost/typeof/typeof.hpp>
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

BOOST_AUTO_TEST_CASE( storage_order )
{
    using namespace pecos::suzerain;

    const pencil<>::dim_type pstart[] = {  0,  0,  0};
    const pencil<>::dim_type psize[]  = { 11, 13, 17};
    const pencil<>::dim_type wstart[] = {  0,  0,  0};
    const pencil<>::dim_type wsize[]  = {  3,  5,  7};

    pencil<> p(pstart, psize, wstart, wsize);

    // Check physical memory layout is (X,Z,Y) in column major storage
    BOOST_CHECK_EQUAL(&p.p(0,0,0) + 1,                 &p.p(1,0,0)); // x
    BOOST_CHECK_EQUAL(&p.p(0,0,0) + psize[0]*psize[2], &p.p(0,1,0)); // y
    BOOST_CHECK_EQUAL(&p.p(0,0,0) + psize[0],          &p.p(0,0,1)); // z

    // Check wavespace memory layout is (Y,Z,X) in column major storage
    BOOST_CHECK_EQUAL(&p.w(0,0,0) + wsize[1]*wsize[2], &p.w(1,0,0)); // x
    BOOST_CHECK_EQUAL(&p.w(0,0,0) + 1,                 &p.w(0,1,0)); // y
    BOOST_CHECK_EQUAL(&p.w(0,0,0) + wsize[1],          &p.w(0,0,1)); // z
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
    for (pencil<>::size_type i = 0; i < p.psize_x; ++i) {
        for (pencil<>::size_type k = 0; k < p.psize_z; ++k) {
            for (pencil<>::size_type j = 0; j < p.psize_y; ++j) {
                p.p(i, j, k) = (i + 1) * (j + 1) * (k + 1);
            }
        }
    }

    // Y, Z, X loop order
    for (pencil<>::size_type j = 0; j < p.psize_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.psize_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.psize_x; ++i) {
                BOOST_CHECK_EQUAL(p.p(i, j, k), (i + 1)*(j + 1)*(k + 1));
            }
        }
    }

    // Clear contents using iterator
    std::fill(p.pbegin(), p.pend(), 0);
    // Check contents are clear
    for (pencil<>::size_type j = 0; j < p.psize_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.psize_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.psize_x; ++i) {
                BOOST_CHECK_EQUAL(p.p(i, j, k), 0);
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
    for (pencil<>::size_type i = 0; i < p.wsize_x; ++i) {
        for (pencil<>::size_type k = 0; k < p.wsize_z; ++k) {
            for (pencil<>::size_type j = 0; j < p.wsize_y; ++j) {
                p.w(i, j, k) = pencil<>::wspace_value_type(
                    (i + 1) * (j + 1) * (k + 1),
                    (i - 1) * (j - 1) * (k - 1));
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (pencil<>::size_type j = 0; j < p.wsize_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.wsize_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.wsize_x; ++i) {
                BOOST_CHECK_EQUAL(p.w(i, j, k), pencil<>::wspace_value_type(
                    (i + 1)*(j + 1)*(k + 1),
                    (i - 1)*(j - 1)*(k - 1)));
                BOOST_CHECK_EQUAL(p.w_real(i, j, k), (i + 1)*(j + 1)*(k + 1));
                BOOST_CHECK_EQUAL(p.w_imag(i, j, k), (i - 1)*(j - 1)*(k - 1));
            }
        }
    }

    // X, Z, Y loop order, assign real and image values
    for (pencil<>::size_type i = 0; i < p.wsize_x; ++i) {
        for (pencil<>::size_type k = 0; k < p.wsize_z; ++k) {
            for (pencil<>::size_type j = 0; j < p.wsize_y; ++j) {
                p.w_real(i, j, k) = (i - 123) * (j - 123) * (k - 123);
                p.w_imag(i, j, k) = (i + 123) * (j + 123) * (k + 123);
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (pencil<>::size_type j = 0; j < p.wsize_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.wsize_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.wsize_x; ++i) {
                BOOST_CHECK_EQUAL(p.w(i, j, k), pencil<>::wspace_value_type(
                    (i - 123)*(j - 123)*(k - 123),
                    (i + 123)*(j + 123)*(k + 123)));
                BOOST_CHECK_EQUAL(
                    p.w_real(i, j, k), (i - 123)*(j - 123)*(k - 123));
                BOOST_CHECK_EQUAL(
                    p.w_imag(i, j, k), (i + 123)*(j + 123)*(k + 123));
            }
        }
    }

    // Clear contents using iterator
    std::fill(p.wbegin(), p.wend(), pencil<>::wspace_value_type(0));
    // Check contents are clear
    for (pencil<>::size_type j = 0; j < p.wsize_y; ++j) {
        for (pencil<>::size_type k = 0; k < p.wsize_z; ++k) {
            for (pencil<>::size_type i = 0; i < p.wsize_x; ++i) {
                BOOST_CHECK_EQUAL(
                    p.w(i, j, k), pencil<>::wspace_value_type(0));
            }
        }
    }
}
