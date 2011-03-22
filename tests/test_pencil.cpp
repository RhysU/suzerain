#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/pencil.hpp>
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( declare_pointer )
{

    using namespace suzerain;

    pencil<> *p;
    (void)p;
}

BOOST_AUTO_TEST_CASE( constructor )
{

    using namespace suzerain;

    const boost::array<pencil<>::index,3> pstart = {{ 0, 0,  0}};
    const boost::array<pencil<>::index,3> psize  = {{16, 7,  4}};
    const boost::array<pencil<>::index,3> wstart = {{ 0, 0,  0}};
    const boost::array<pencil<>::index,3> wsize  = {{ 7, 4, 16}};

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
    using namespace suzerain;

    const boost::array<pencil<>::index,3> pstart = {{  0,  0,  0}};
    const boost::array<pencil<>::index,3> psize  = {{ 11, 13, 17}};
    const boost::array<pencil<>::index,3> wstart = {{  0,  0,  0}};
    const boost::array<pencil<>::index,3> wsize  = {{  3,  5,  7}};

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

    // Check wavespace memory layout is (Y,X,Z) in column major storage
    BOOST_CHECK_EQUAL(
        &p.wave(0,0,0) + wsize[1],          &p.wave(1,0,0)); // x
    BOOST_CHECK_EQUAL(
        &p.wave(0,0,0) + 1,                 &p.wave(0,1,0)); // y
    BOOST_CHECK_EQUAL(
        &p.wave(0,0,0) + wsize[1]*wsize[0], &p.wave(0,0,1)); // z
}

BOOST_AUTO_TEST_CASE( offsets_and_inverse_offsets )
{
    using namespace suzerain;
    typedef pencil<>::index index;

    const boost::array<pencil<>::index,3> pstart = {{  5,  6,  7}};
    const boost::array<pencil<>::index,3> psize  = {{  2,  3,  5}};
    const boost::array<pencil<>::index,3> wstart = {{  1,  2,  3}};
    const boost::array<pencil<>::index,3> wsize  = {{  3,  5,  7}};

    pencil<> p(pstart, psize, wstart, wsize);

    // Check that we can invert physical space offsets
    for (index i = 0; i < (index) p.physical.size_x; ++i) {
        for (index j = 0; j < (index) p.physical.size_y; ++j) {
            for (index k = 0; k < (index) p.physical.size_z; ++k) {
                // Local offsets
                pencil<>::index x, y, z;
                p.physical.inverse_offset(
                    p.physical.offset(i,j,k), x, y, z);
                BOOST_CHECK_EQUAL(i, x);
                BOOST_CHECK_EQUAL(j, y);
                BOOST_CHECK_EQUAL(k, z);

                // Global offsets
                p.physical.inverse_global_offset(
                    p.physical.offset(i,j,k), x, y, z);
                BOOST_CHECK_EQUAL(i + pstart[0], x);
                BOOST_CHECK_EQUAL(j + pstart[1], y);
                BOOST_CHECK_EQUAL(k + pstart[2], z);
            }
        }
    }

    // Check that we can invert wave space offsets
    for (index i = 0; i < (index) p.wave.size_x; ++i) {
        for (index j = 0; j < (index) p.wave.size_y; ++j) {
            for (index k = 0; k < (index) p.wave.size_z; ++k) {
                // Local offsets
                pencil<>::index x, y, z;
                p.wave.inverse_offset(
                    p.wave.offset(i,j,k), x, y, z);
                BOOST_CHECK_EQUAL(i, x);
                BOOST_CHECK_EQUAL(j, y);
                BOOST_CHECK_EQUAL(k, z);

                // Global offsets
                p.wave.inverse_global_offset(
                    p.wave.offset(i,j,k), x, y, z);
                BOOST_CHECK_EQUAL(i + wstart[0], x);
                BOOST_CHECK_EQUAL(j + wstart[1], y);
                BOOST_CHECK_EQUAL(k + wstart[2], z);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( real_access )
{
    using namespace suzerain;
    typedef pencil<>::index index;

    const boost::array<pencil<>::index,3> pstart = {{ 0, 0,  0}};
    const boost::array<pencil<>::index,3> psize  = {{ 2, 2,  2}};
    const boost::array<pencil<>::index,3> wstart = {{ 0, 0,  0}};
    const boost::array<pencil<>::index,3> wsize  = {{ 2, 1,  1}};

    pencil<> p(pstart, psize, wstart, wsize);

    // X, Z, Y loop order
    for (index i = 0; i < (index) p.physical.size_x; ++i) {
        for (index k = 0; k < (index) p.physical.size_z; ++k) {
            for (index j = 0; j < (index) p.physical.size_y; ++j) {
#pragma warning(push,disable:810 2259)
                p.physical(i, j, k) = (i + 1) * (j + 1) * (k + 1);
#pragma warning(pop)
            }
        }
    }

    // Y, Z, X loop order
    for (index j = 0; j < (index) p.physical.size_y; ++j) {
        for (index k = 0; k < (index) p.physical.size_z; ++k) {
            for (index i = 0; i < (index) p.physical.size_x; ++i) {
                BOOST_CHECK_EQUAL(p.physical(i, j, k), (i + 1)*(j + 1)*(k + 1));
            }
        }
    }

    // Clear contents using iterator
    std::fill(p.physical.begin(), p.physical.end(), 0);
    // Check contents are clear
    for (index j = 0; j < (index) p.physical.size_y; ++j) {
        for (index k = 0; k < (index) p.physical.size_z; ++k) {
            for (index i = 0; i < (index) p.physical.size_x; ++i) {
                BOOST_CHECK_EQUAL(p.physical(i, j, k), 0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( complex_access )
{
    using namespace suzerain;
    typedef pencil<>::index index;

    const boost::array<pencil<>::index,3> pstart = {{ 1,  1,  1}};
    const boost::array<pencil<>::index,3> psize  = {{ 2,  3,  5}};
    const boost::array<pencil<>::index,3> wstart = {{ 2,  2,  2}};
    const boost::array<pencil<>::index,3> wsize  = {{ 7, 11, 13}};

    pencil<> p(pstart, psize, wstart, wsize);

    // X, Z, Y loop order, assign complex value
    for (index i = 0; i < (index) p.wave.size_x; ++i) {
        for (index k = 0; k < (index) p.wave.size_z; ++k) {
            for (index j = 0; j < (index) p.wave.size_y; ++j) {
                p.wave(i, j, k) = pencil<>::complex_type(
                    (i + 1) * (j + 1) * (k + 1),
                    (i - 1) * (j - 1) * (k - 1));
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (index j = 0; j < (index) p.wave.size_y; ++j) {
        for (index k = 0; k < (index) p.wave.size_z; ++k) {
            for (index i = 0; i < (index) p.wave.size_x; ++i) {
                BOOST_CHECK_EQUAL(p.wave(i, j, k), pencil<>::complex_type(
                    (i + 1)*(j + 1)*(k + 1),
                    (i - 1)*(j - 1)*(k - 1)));
                BOOST_CHECK_EQUAL(p.wave.real(i, j, k), (i + 1)*(j + 1)*(k + 1));
                BOOST_CHECK_EQUAL(p.wave.imag(i, j, k), (i - 1)*(j - 1)*(k - 1));
            }
        }
    }

    // X, Z, Y loop order, assign real and imag values
    for (index i = 0; i < (index) p.wave.size_x; ++i) {
        for (index k = 0; k < (index) p.wave.size_z; ++k) {
            for (index j = 0; j < (index) p.wave.size_y; ++j) {
#pragma warning(push,disable:810 2259)
                p.wave.real(i, j, k) = (i - 123) * (j - 123) * (k - 123);
                p.wave.imag(i, j, k) = (i + 123) * (j + 123) * (k + 123);
#pragma warning(pop)
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (index j = 0; j < (index) p.wave.size_y; ++j) {
        for (index k = 0; k < (index) p.wave.size_z; ++k) {
            for (index i = 0; i < (index) p.wave.size_x; ++i) {
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
    for (index j = 0; j < (index) p.wave.size_y; ++j) {
        for (index k = 0; k < (index) p.wave.size_z; ++k) {
            for (index i = 0; i < (index) p.wave.size_x; ++i) {
                BOOST_CHECK_EQUAL(
                    p.wave(i, j, k), pencil<>::complex_type(0));
            }
        }
    }
}
