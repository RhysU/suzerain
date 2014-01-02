//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#include <suzerain/inorder.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>

#include "test_tools.hpp"

BOOST_AUTO_TEST_CASE( wavenumber_misc )
{
    BOOST_CHECK_EQUAL( 0, suzerain::inorder::wavenumber_min   (1));
    BOOST_CHECK_EQUAL( 0, suzerain::inorder::wavenumber_absmin(1));
    BOOST_CHECK_EQUAL( 0, suzerain::inorder::wavenumber_max   (1));

    BOOST_CHECK_EQUAL( 0, suzerain::inorder::wavenumber_min   (2));
    BOOST_CHECK_EQUAL( 0, suzerain::inorder::wavenumber_absmin(2));
    BOOST_CHECK_EQUAL( 1, suzerain::inorder::wavenumber_max   (2));

    BOOST_CHECK_EQUAL(-1, suzerain::inorder::wavenumber_min   (3));
    BOOST_CHECK_EQUAL( 1, suzerain::inorder::wavenumber_absmin(3));
    BOOST_CHECK_EQUAL( 1, suzerain::inorder::wavenumber_max   (3));

    BOOST_CHECK_EQUAL(-3, suzerain::inorder::wavenumber_min   (7));
    BOOST_CHECK_EQUAL( 3, suzerain::inorder::wavenumber_absmin(7));
    BOOST_CHECK_EQUAL( 3, suzerain::inorder::wavenumber_max   (7));

    BOOST_CHECK_EQUAL(-3, suzerain::inorder::wavenumber_min   (8));
    BOOST_CHECK_EQUAL( 3, suzerain::inorder::wavenumber_absmin(8));
    BOOST_CHECK_EQUAL( 4, suzerain::inorder::wavenumber_max   (8));
}

BOOST_AUTO_TEST_CASE( wavenumber )
{
    using suzerain::inorder::wavenumber;

    const int expected[][10] = {
        { 0                                    },
        { 0, 1                                 },
        { 0, 1,                             -1 },
        { 0, 1,  2,                         -1 },
        { 0, 1,  2,                     -2, -1 },
        { 0, 1,  2,  3,                 -2, -1 },
        { 0, 1,  2,  3,             -3, -2, -1 },
        { 0, 1,  2,  3,  4,         -3, -2, -1 },
        { 0, 1,  2,  3,  4,     -4, -3, -2, -1 },
        { 0, 1,  2,  3,  4,  5, -4, -3, -2, -1 }
    };
    for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
        int result[sizeof(expected[0])/sizeof(expected[0][0])];
        for (int j = 0; j < i+1; ++j) {
            result[j] = wavenumber(i + 1, j);
        }
        BOOST_CHECK_EQUAL_COLLECTIONS(
                expected[i], expected[i] + i + 1, result, result + i + 1);
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_valid )
{
    using suzerain::inorder::wavenumber_valid;

    BOOST_CHECK(!wavenumber_valid(1, -1));
    BOOST_CHECK( wavenumber_valid(1,  0));
    BOOST_CHECK(!wavenumber_valid(1,  1));

    BOOST_CHECK(!wavenumber_valid(2, -1));
    BOOST_CHECK( wavenumber_valid(2,  0));
    BOOST_CHECK( wavenumber_valid(2,  1));
    BOOST_CHECK(!wavenumber_valid(2,  2));

    BOOST_CHECK(!wavenumber_valid(3, -2));
    BOOST_CHECK( wavenumber_valid(3, -1));
    BOOST_CHECK( wavenumber_valid(3,  0));
    BOOST_CHECK( wavenumber_valid(3,  1));
    BOOST_CHECK(!wavenumber_valid(3,  2));

    BOOST_CHECK(!wavenumber_valid(4, -2));
    BOOST_CHECK( wavenumber_valid(4, -1));
    BOOST_CHECK( wavenumber_valid(4,  0));
    BOOST_CHECK( wavenumber_valid(4,  1));
    BOOST_CHECK( wavenumber_valid(4,  2));
    BOOST_CHECK(!wavenumber_valid(4,  3));

    BOOST_CHECK(!wavenumber_valid(5, -3));
    BOOST_CHECK( wavenumber_valid(5, -2));
    BOOST_CHECK( wavenumber_valid(5, -1));
    BOOST_CHECK( wavenumber_valid(5,  0));
    BOOST_CHECK( wavenumber_valid(5,  1));
    BOOST_CHECK( wavenumber_valid(5,  2));
    BOOST_CHECK(!wavenumber_valid(5,  3));

    BOOST_CHECK(!wavenumber_valid(6, -3));
    BOOST_CHECK( wavenumber_valid(6, -2));
    BOOST_CHECK( wavenumber_valid(6, -1));
    BOOST_CHECK( wavenumber_valid(6,  0));
    BOOST_CHECK( wavenumber_valid(6,  1));
    BOOST_CHECK( wavenumber_valid(6,  2));
    BOOST_CHECK( wavenumber_valid(6,  3));
    BOOST_CHECK(!wavenumber_valid(6,  4));
}

BOOST_AUTO_TEST_CASE( inorder_index )
{
    using suzerain::inorder::index;

    const int freq[][10] = {
        { 0                                    },
        { 0, 1                                 },
        { 0, 1,                             -1 },
        { 0, 1,  2,                         -1 },
        { 0, 1,  2,                     -2, -1 },
        { 0, 1,  2,  3,                 -2, -1 },
        { 0, 1,  2,  3,             -3, -2, -1 },
        { 0, 1,  2,  3,  4,         -3, -2, -1 },
        { 0, 1,  2,  3,  4,     -4, -3, -2, -1 },
        { 0, 1,  2,  3,  4,  5, -4, -3, -2, -1 }
    };
    const int expected[][10] = {
        { 0                                    },
        { 0, 1                                 },
        { 0, 1,                              2 },
        { 0, 1,  2,                          3 },
        { 0, 1,  2,                      3,  4 },
        { 0, 1,  2,  3,                  4,  5 },
        { 0, 1,  2,  3,              4,  5,  6 },
        { 0, 1,  2,  3,  4,          5,  6,  7 },
        { 0, 1,  2,  3,  4,      5,  6,  7,  8 },
        { 0, 1,  2,  3,  4,  5,  6,  7,  8,  9 }
    };
    for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
        int result[sizeof(expected[0])/sizeof(expected[0][0])];
        for (int j = 0; j < i+1; ++j) {
            result[j] = index(i + 1, freq[i][j]);
        }
        BOOST_CHECK_EQUAL_COLLECTIONS(
                expected[i], expected[i] + i + 1, result, result + i + 1);
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_abs )
{
    using suzerain::inorder::wavenumber_abs;

    const int expected[][10] = {
        { 0                                    },
        { 0, 1                                 },
        { 0, 1,                              1 },
        { 0, 1,  2,                          1 },
        { 0, 1,  2,                      2,  1 },
        { 0, 1,  2,  3,                  2,  1 },
        { 0, 1,  2,  3,              3,  2,  1 },
        { 0, 1,  2,  3,  4,          3,  2,  1 },
        { 0, 1,  2,  3,  4,      4,  3,  2,  1 },
        { 0, 1,  2,  3,  4,  5,  4,  3,  2,  1 }
    };
    for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
        int result[sizeof(expected[0])/sizeof(expected[0][0])];
        for (int j = 0; j < i+1; ++j) {
            result[j] = wavenumber_abs(i + 1, j);
        }
        BOOST_CHECK_EQUAL_COLLECTIONS(
                expected[i], expected[i] + i + 1, result, result + i + 1);
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_nyquist )
{
    using suzerain::inorder::wavenumber;
    using suzerain::inorder::wavenumber_nyquist;
    using suzerain::inorder::wavenumber_nyquist_index;

    const int expected[][10] = {
        { 0                                    },
        { 0, 1                                 },
        { 0, 0,                              0 },
        { 0, 0,  1,                          0 },
        { 0, 0,  0,                      0,  0 },
        { 0, 0,  0,  1,                  0,  0 },
        { 0, 0,  0,  0,              0,  0,  0 },
        { 0, 0,  0,  0,  1,          0,  0,  0 },
        { 0, 0,  0,  0,  0,      0,  0,  0,  0 },
        { 0, 0,  0,  0,  0,  1,  0,  0,  0,  0 }
    };
    for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
        int result[sizeof(expected[0])/sizeof(expected[0][0])];
        for (int j = 0; j < i+1; ++j) {
            const int w = wavenumber(i + 1, j);
            result[j]   = wavenumber_nyquist(i + 1, w);

            // Consistency between wavenumber_nyquist{,_index}
            BOOST_CHECK_EQUAL(result[j], wavenumber_nyquist_index(i + 1, j));
        }

        // Correctness of the results against expected
        BOOST_CHECK_EQUAL_COLLECTIONS(
                expected[i], expected[i] + i + 1, result, result + i + 1);
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_imagzero )
{
    using suzerain::inorder::wavenumber;
    using suzerain::inorder::wavenumber_imagzero;
    using suzerain::inorder::wavenumber_imagzero_index;

    const int expected[][10] = {
        { 1                                    },
        { 1, 1                                 },
        { 1, 0,                              0 },
        { 1, 0,  1,                          0 },
        { 1, 0,  0,                      0,  0 },
        { 1, 0,  0,  1,                  0,  0 },
        { 1, 0,  0,  0,              0,  0,  0 },
        { 1, 0,  0,  0,  1,          0,  0,  0 },
        { 1, 0,  0,  0,  0,      0,  0,  0,  0 },
        { 1, 0,  0,  0,  0,  1,  0,  0,  0,  0 }
    };
    for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
        int result[sizeof(expected[0])/sizeof(expected[0][0])];
        for (int j = 0; j < i+1; ++j) {
            const int w = wavenumber(i + 1, j);
            result[j]   = wavenumber_imagzero(i + 1, w);

            // Consistency between wavenumber_imagzero{,_index}
            BOOST_CHECK_EQUAL(result[j], wavenumber_imagzero_index(i + 1, j));
        }

        // Correctness of the results against expected
        BOOST_CHECK_EQUAL_COLLECTIONS(
                expected[i], expected[i] + i + 1, result, result + i + 1);
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_diff_nodealiasing )
{
    using suzerain::inorder::wavenumber_diff;

    const int expected[][8] = { {0},
                                {0,               0},
                                {0,1,            -1},
                                {0,1,    0,      -1},
                                {0,1,2,       -2,-1},
                                {0,1,2,  0,   -2,-1},
                                {0,1,2,3,  -3,-2,-1},
                                {0,1,2,3,0,-3,-2,-1} };

    for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
        const int N = i + 1, dN = N;
        for (int j = 0; j < dN; ++j) {
            BOOST_CHECK_EQUAL(expected[i][j], wavenumber_diff(N, dN, j));
        }
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_diff_dealiasing )
{
    using suzerain::inorder::wavenumber_diff;

    {
        const int expected[][12] = { {0,0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,1,0,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,0,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,2,0,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,0,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,3,0,0,0,0,0,-3,-2,-1},
                                     {0,1,2,3,0,0,0,0,0,-3,-2,-1} };

        for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
                BOOST_CHECK_EQUAL(expected[i][j], wavenumber_diff(N, dN, j));
            }
        }
    }

    {
        const int expected[][11] = { {0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,1,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,2,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,3,0,0,0,0,-3,-2,-1},
                                     {0,1,2,3,0,0,0,0,-3,-2,-1} };

        for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
                BOOST_CHECK_EQUAL(expected[i][j], wavenumber_diff(N, dN, j));
            }
        }
    }

    {
        const int expected[][9] = { {0,0,0,0,0, 0, 0, 0, 0},
                                    {0,0,0,0,0, 0, 0, 0, 0},
                                    {0,1,0,0,0, 0, 0, 0,-1},
                                    {0,1,0,0,0, 0, 0, 0,-1},
                                    {0,1,2,0,0, 0, 0,-2,-1},
                                    {0,1,2,0,0, 0, 0,-2,-1},
                                    {0,1,2,3,0, 0,-3,-2,-1},
                                    {0,1,2,3,0, 0,-3,-2,-1},
                                    {0,1,2,3,4,-4,-3,-2,-1} };

        for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
                BOOST_CHECK_EQUAL(expected[i][j], wavenumber_diff(N, dN, j));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( nondealiased_nodealiasing )
{
    using suzerain::inorder::wavenumber_translatable;

    const int expected[][8] = { {1},
                                {1,1},
                                {1,1,1},
                                {1,1,1,1},
                                {1,1,1,1,1},
                                {1,1,1,1,1,1},
                                {1,1,1,1,1,1,1},
                                {1,1,1,1,1,1,1,1} };

    for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
        const int N = i + 1, dN = N;
        for (int j = 0; j < dN; ++j) {
            BOOST_CHECK_EQUAL(expected[i][j], wavenumber_translatable(N, dN, j));
        }
    }
}

BOOST_AUTO_TEST_CASE( nondealiased_dealiasing )
{
    using suzerain::inorder::wavenumber_translatable;

    {
        const int expected[][12] = { {1,0,0,0,0,0,0,0,0,0,0,0},
                                     {1,1,0,0,0,0,0,0,0,0,0,0},
                                     {1,1,0,0,0,0,0,0,0,0,0,1},
                                     {1,1,1,0,0,0,0,0,0,0,0,1},
                                     {1,1,1,0,0,0,0,0,0,0,1,1},
                                     {1,1,1,1,0,0,0,0,0,0,1,1},
                                     {1,1,1,1,0,0,0,0,0,1,1,1},
                                     {1,1,1,1,1,0,0,0,0,1,1,1} };

        for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
//              BOOST_TEST_MESSAGE("Checking N = " << N  <<
//                                 ", dN = "       << dN <<
//                                 ", j = "        << j);
                BOOST_CHECK_EQUAL(expected[i][j], wavenumber_translatable(N, dN, j));
            }
        }
    }

    {
        const int expected[][11] = { {1,0,0,0,0,0,0,0,0,0,0},
                                     {1,1,0,0,0,0,0,0,0,0,0},
                                     {1,1,0,0,0,0,0,0,0,0,1},
                                     {1,1,1,0,0,0,0,0,0,0,1},
                                     {1,1,1,0,0,0,0,0,0,1,1},
                                     {1,1,1,1,0,0,0,0,0,1,1},
                                     {1,1,1,1,0,0,0,0,1,1,1},
                                     {1,1,1,1,1,0,0,0,1,1,1} };

        for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
//              BOOST_TEST_MESSAGE("Checking N = " << N  <<
//                                 ", dN = "       << dN <<
//                                 ", j = "        << j);
                BOOST_CHECK_EQUAL(expected[i][j], wavenumber_translatable(N, dN, j));
            }
        }
    }

    {
        const int expected[][9] = { {1,0,0,0,0,0,0,0,0},
                                    {1,1,0,0,0,0,0,0,0},
                                    {1,1,0,0,0,0,0,0,1},
                                    {1,1,1,0,0,0,0,0,1},
                                    {1,1,1,0,0,0,0,1,1},
                                    {1,1,1,1,0,0,0,1,1},
                                    {1,1,1,1,0,0,1,1,1},
                                    {1,1,1,1,1,0,1,1,1},
                                    {1,1,1,1,1,1,1,1,1} };

        for (int i = 0; i < (int) (sizeof(expected)/sizeof(expected[0])); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
//              BOOST_TEST_MESSAGE("Checking N = " << N  <<
//                                 ", dN = "       << dN <<
//                                 ", j = "        << j);
                BOOST_CHECK_EQUAL(expected[i][j], wavenumber_translatable(N, dN, j));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_translate_dN_equal_to_N )
{
    using suzerain::inorder::wavenumber_translate;

    // Obtain vectors as solutions to ease test case writing
    int dat[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
    int & kb1 = dat[0], & ke1 = dat[1], & kb2 = dat[2], & ke2 = dat[3];
    int &dkb1 = dat[4], &dke1 = dat[5], &dkb2 = dat[6], &dke2 = dat[7];

    // Sanity check that a range of N, dkb, and dke return NOPs
    for (int N = 4; N < 9; ++N) {
        for (int dkb = 0; dkb <= N; ++dkb) {
            for (int dke = dkb; dke <= N; ++dke) {
                wavenumber_translate(N,  N, dkb, dke,
                        kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
                BOOST_CHECK_EQUAL(dkb, kb1);
                BOOST_CHECK_EQUAL(dke, ke1);
                BOOST_CHECK_EQUAL(dke, kb2);
                BOOST_CHECK_EQUAL(dke, ke2);
                BOOST_CHECK_EQUAL(kb1, dkb1);
                BOOST_CHECK_EQUAL(ke1, dke1);
                BOOST_CHECK_EQUAL(kb2, dkb2);
                BOOST_CHECK_EQUAL(ke2, dke2);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_translate_dN_greater_than_N )
{
    using suzerain::inorder::wavenumber_translate;

    // Obtain vectors as solutions to ease test case writing
    int dat[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
    int & kb1 = dat[0], & ke1 = dat[1], & kb2 = dat[2], & ke2 = dat[3];
    int &dkb1 = dat[4], &dke1 = dat[5], &dkb2 = dat[6], &dke2 = dat[7];

    int N, dN;

    // Positive wavenumbers with no dealiasing impact
    dN = 13;
    N = 8;
    {
        const int dkb = 0, dke = 4;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 1, dke = 4;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 0, dke = 3;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    N = 7;
    {
        const int dkb = 0, dke = 4;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 1, dke = 4;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 0, dke = 3;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    dN = 8;
    N = 6;
    {
        const int dkb = 0, dke = 3;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    N = 5;
    {
        const int dkb = 0, dke = 3;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      dkb, dke, dke, dke,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }

    // Negative wavenumbers with no dealiasing impact
    dN = 13;
    N  = 8;
    {
        const int dkb = 10, dke = 13;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        5,   8,   8,   8,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 10, dke = 12;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        5,   7,   7,   7,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 11, dke = 13;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        6,   8,   8,   8,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    N = 7;
    {
        const int dkb = 10, dke = 13;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        4,   7,   7,   7,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 10, dke = 12;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        4,   6,   6,   6,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 11, dke = 13;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        5,   7,   7,   7,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    dN = 8;
    N = 6;
    {
        const int dkb = 6, dke = 8;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        4,   6,   6,   6,  dkb, dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    N = 5;
    {
        const int dkb = 6, dke = 8;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        3,   5,   5,   5,  dkb,  dke,  dke,  dke };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }

    //  Positive and negative wavenumbers
    dN = 13;
    N = 8;
    {
        const int dkb = 0, dke = 13;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      0,     5,   5,   8,    0,    5,   10,   13 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 1, dke = 12;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      1,     5,   5,   7,    1,    5,   10,   12 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    N = 7;
    {
        const int dkb = 0, dke = 13;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      0,     4,   4,   7,    0,    4,   10,   13 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 1, dke = 12;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      1,     4,   4,   6,    1,    4,   10,   12 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    dN = 8;
    N = 6;
    {
        const int dkb = 0, dke = 8;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      0,     4,   4,   6,    0,    4,    6,    8 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 1, dke = 7;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {      1,     4,   4,   5,    1,    4,    6,    7 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }

    // Some dealiased wavenumbers
    dN = 13;
    N = 8;
    {
        const int dkb = 2, dke = 8;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        2,   5,   5,   5,    2,    5,    8,    8 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 8, dke = 11;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        5,   6,   6,   6,   10,   11,   11,   11 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    N = 7;
    {
        const int dkb = 2, dke = 8;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        2,   4,   4,   4,    2,    4,    8,    8 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 8, dke = 11;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        4,   5,   5,   5,   10,   11,   11,   11 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    dN = 8;
    N = 6;
    {
        const int dkb = 2, dke = 5;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        2,   4,   4,   4,    2,    4,    5,    5 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    N = 5;
    {
        const int dkb = 2, dke = 5;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        2,   3,   3,   3,    2,    3,    5,    5 };
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
}

BOOST_AUTO_TEST_CASE( wavenumber_translate_N_greater_than_dN )
{
    using suzerain::inorder::wavenumber_translate;

    // Obtain vectors as solutions to ease test case writing
    int dat[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
    int & kb1 = dat[0], & ke1 = dat[1], & kb2 = dat[2], & ke2 = dat[3];
    int &dkb1 = dat[4], &dke1 = dat[5], &dkb2 = dat[6], &dke2 = dat[7];

    int N, dN;

    N = 13;
    dN = 8;
    {
        const int dkb = 0, dke = 8;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        0,   5,  10,  13,    0,    5,    5,    8};
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 1, dke = 7;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        1,   5,  10,  12,    1,    5,    5,    7};
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 2, dke = 6;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        2,   5,  10,  11,    2,    5,    5,    6};
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 3, dke = 5;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        3,   5,  10,  10,    3,    5,    5,    5};
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }

    N = 8;
    dN = 5;
    {
        const int dkb = 0, dke = 5;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        0,   3,   6,   8,    0,    3,    3,    5};
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
    {
        const int dkb = 1, dke = 4;
        wavenumber_translate(
                N, dN, dkb, dke, kb1, ke1, kb2, ke2, dkb1, dke1, dkb2, dke2);
        const int ex[8] = {        1,   3,   6,   7,    1,    3,    3,    4};
        BOOST_CHECK_EQUAL_COLLECTIONS(ex, ex + 8, dat, dat + 8);
    }
}
