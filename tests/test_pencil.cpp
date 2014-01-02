//--------------------------------------------------------------------------
//
// Copyright (C) 2008-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#include <suzerain/pencil.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

using suzerain::array;
using suzerain::pencil;

BOOST_AUTO_TEST_CASE( declare_pointer )
{
    pencil<> *p = NULL;
    (void)p;
    BOOST_CHECK(!p);
}

BOOST_AUTO_TEST_CASE( constructor )
{
    const array<pencil<>::index,3>     pstart = {{  1, 2,  3}};
    const array<pencil<>::size_type,3> psize  = {{ 16, 7,  4}};
    const array<pencil<>::index,3>     wstart = {{  3, 2,  1}};
    const array<pencil<>::size_type,3> wsize  = {{  7, 4, 16}};

    pencil<> p(pstart, psize, wstart, wsize);

    BOOST_CHECK_EQUAL( p.global_physical.index_bases()[0], pstart[0] );
    BOOST_CHECK_EQUAL( p.global_physical.index_bases()[1], pstart[1] );
    BOOST_CHECK_EQUAL( p.global_physical.index_bases()[2], pstart[2] );
    BOOST_CHECK_EQUAL( p.global_physical.shape()[0],       psize[0]  );
    BOOST_CHECK_EQUAL( p.global_physical.shape()[1],       psize[1]  );
    BOOST_CHECK_EQUAL( p.global_physical.shape()[2],       psize[2]  );

    BOOST_CHECK_EQUAL( p.physical.index_bases()[0], 0 );
    BOOST_CHECK_EQUAL( p.physical.index_bases()[1], 0 );
    BOOST_CHECK_EQUAL( p.physical.index_bases()[2], 0 );
    BOOST_CHECK_EQUAL( p.physical.shape()[0],       psize[0]  );
    BOOST_CHECK_EQUAL( p.physical.shape()[1],       psize[1]  );
    BOOST_CHECK_EQUAL( p.physical.shape()[2],       psize[2]  );

    BOOST_CHECK_EQUAL( p.global_wave.index_bases()[0], wstart[0] );
    BOOST_CHECK_EQUAL( p.global_wave.index_bases()[1], wstart[1] );
    BOOST_CHECK_EQUAL( p.global_wave.index_bases()[2], wstart[2] );
    BOOST_CHECK_EQUAL( p.global_wave.shape()[0],       wsize[0]  );
    BOOST_CHECK_EQUAL( p.global_wave.shape()[1],       wsize[1]  );
    BOOST_CHECK_EQUAL( p.global_wave.shape()[2],       wsize[2]  );

    BOOST_CHECK_EQUAL( p.wave.index_bases()[0], 0 );
    BOOST_CHECK_EQUAL( p.wave.index_bases()[1], 0 );
    BOOST_CHECK_EQUAL( p.wave.index_bases()[2], 0 );
    BOOST_CHECK_EQUAL( p.wave.shape()[0],       wsize[0]  );
    BOOST_CHECK_EQUAL( p.wave.shape()[1],       wsize[1]  );
    BOOST_CHECK_EQUAL( p.wave.shape()[2],       wsize[2]  );
}

BOOST_AUTO_TEST_CASE( storage_order )
{
    const array<pencil<>::index,3>     pstart = {{  0,  0,  0}};
    const array<pencil<>::size_type,3> psize  = {{ 11, 13, 17}};
    const array<pencil<>::index,3>     wstart = {{  0,  0,  0}};
    const array<pencil<>::size_type,3> wsize  = {{  3,  5,  7}};

    pencil<> p(pstart, psize, wstart, wsize);

    // Ensure memory locations are initialized for test below,
    // otherwise valgrind complains bitterly about uninitialized
    // locations being used for jumps
    p.physical[0][0][0] = 1;
    p.physical[1][0][0] = 1;
    p.physical[0][1][0] = 1;
    p.physical[0][0][1] = 1;
    p.wave[0][0][0] = 1;
    p.wave[1][0][0] = 1;
    p.wave[0][1][0] = 1;
    p.wave[0][0][1] = 1;

    // Check physical memory layout is (X,Z,Y) in column major storage
    BOOST_CHECK_EQUAL(
        &p.physical[0][0][0] + 1,                 &p.physical[1][0][0]); // x
    BOOST_CHECK_EQUAL(
        &p.physical[0][0][0] + psize[0]*psize[2], &p.physical[0][1][0]); // y
    BOOST_CHECK_EQUAL(
        &p.physical[0][0][0] + psize[0],          &p.physical[0][0][1]); // z

    // Check wavespace memory layout is (Y,X,Z) in column major storage
    BOOST_CHECK_EQUAL(
        &p.wave[0][0][0] + wsize[1],          &p.wave[1][0][0]); // x
    BOOST_CHECK_EQUAL(
        &p.wave[0][0][0] + 1,                 &p.wave[0][1][0]); // y
    BOOST_CHECK_EQUAL(
        &p.wave[0][0][0] + wsize[1]*wsize[0], &p.wave[0][0][1]); // z
}

BOOST_AUTO_TEST_CASE( real_access )
{
    typedef pencil<>::index     index;
    typedef pencil<>::size_type size_type;
    const array<index,3>     pstart = {{ 0, 0, 0}};
    const array<size_type,3> psize  = {{ 2, 2, 2}};
    const array<index,3>     wstart = {{ 0, 0, 0}};
    const array<size_type,3> wsize  = {{ 2, 1, 1}};

    pencil<> p(pstart, psize, wstart, wsize);

    // X, Z, Y loop order
    for (index i = 0; i < (index) p.physical.shape()[0]; ++i) {
        for (index k = 0; k < (index) p.physical.shape()[2]; ++k) {
            for (index j = 0; j < (index) p.physical.shape()[1]; ++j) {
#pragma warning(push,disable:810 2259)
                p.physical[i][j][k] = (i + 1) * (j + 1) * (k + 1);
#pragma warning(pop)
            }
        }
    }

    // Y, Z, X loop order
    for (index j = 0; j < (index) p.physical.shape()[1]; ++j) {
        for (index k = 0; k < (index) p.physical.shape()[2]; ++k) {
            for (index i = 0; i < (index) p.physical.shape()[0]; ++i) {
                BOOST_CHECK_EQUAL(p.physical[i][j][k],
                                  (i + 1)*(j + 1)*(k + 1));
            }
        }
    }

    // Clear contents using suzerain::multi_array::fill using Koenig lookup
    fill(p.physical, 0);

    // Check contents are clear
    for (index j = 0; j < (index) p.physical.shape()[1]; ++j) {
        for (index k = 0; k < (index) p.physical.shape()[2]; ++k) {
            for (index i = 0; i < (index) p.physical.shape()[0]; ++i) {
                BOOST_CHECK_EQUAL(p.physical[i][j][k], 0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( complex_access )
{
    typedef pencil<>::index     index;
    typedef pencil<>::size_type size_type;
    const array<index,3>     pstart = {{ 1, 1,  1}};
    const array<size_type,3> psize  = {{ 2, 3,  5}};
    const array<index,3>     wstart = {{ 2, 2,  2}};
    const array<size_type,3> wsize  = {{ 7, 11, 13}};

    pencil<> p(pstart, psize, wstart, wsize);

    // X, Z, Y loop order, assign complex value
    for (index i = 0; i < (index) p.wave.shape()[0]; ++i) {
        for (index k = 0; k < (index) p.wave.shape()[2]; ++k) {
            for (index j = 0; j < (index) p.wave.shape()[1]; ++j) {
                p.wave[i][j][k] = pencil<>::complex_type(
                    (i + 1) * (j + 1) * (k + 1),
                    (i - 1) * (j - 1) * (k - 1));
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (index j = 0; j < (index) p.wave.shape()[1]; ++j) {
        for (index k = 0; k < (index) p.wave.shape()[2]; ++k) {
            for (index i = 0; i < (index) p.wave.shape()[0]; ++i) {
                BOOST_CHECK_EQUAL(p.wave[i][j][k],
                                  pencil<>::complex_type(
                                        (i + 1)*(j + 1)*(k + 1),
                                        (i - 1)*(j - 1)*(k - 1)));
            }
        }
    }

    // X, Z, Y loop order, assign real and imag values
    for (index i = 0; i < (index) p.wave.shape()[0]; ++i) {
        for (index k = 0; k < (index) p.wave.shape()[2]; ++k) {
            for (index j = 0; j < (index) p.wave.shape()[1]; ++j) {
#pragma warning(push,disable:810 2259)
                p.wave[i][j][k].real() = (i - 123) * (j - 123) * (k - 123);
                p.wave[i][j][k].imag() = (i + 123) * (j + 123) * (k + 123);
#pragma warning(pop)
            }
        }
    }

    // Y, Z, X loop order, check values are correct
    for (index j = 0; j < (index) p.wave.shape()[1]; ++j) {
        for (index k = 0; k < (index) p.wave.shape()[2]; ++k) {
            for (index i = 0; i < (index) p.wave.shape()[0]; ++i) {
                BOOST_CHECK_EQUAL(p.wave[i][j][k],
                                  pencil<>::complex_type(
                                        (i - 123)*(j - 123)*(k - 123),
                                        (i + 123)*(j + 123)*(k + 123)));
            }
        }
    }

    // Clear contents using iterator
    fill(p.wave, 0);

    // Check contents are clear
    for (index j = 0; j < (index) p.wave.shape()[1]; ++j) {
        for (index k = 0; k < (index) p.wave.shape()[2]; ++k) {
            for (index i = 0; i < (index) p.wave.shape()[0]; ++i) {
                BOOST_CHECK_EQUAL(
                    p.wave[i][j][k], pencil<>::complex_type(0));
            }
        }
    }
}
