//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#include <suzerain/gbmatrix.h>

#include <gsl/gsl_test.h>

#include <suzerain/common.h>

// Test Matrix A
//     [11,  12,  13,  14,   0,   0,   0,   0;
//      21,  22,  23,  24,  25,   0,   0,   0;
//      31,  32,  33,  34,  35,  36,   0,   0;
//       0,  42,  43,  44,  45,  46,  47,   0;
//       0,   0,  53,  54,  55,  56,  57,  58;
//       0,   0,   0,  64,  65,  66,  67,  68;
//       0,   0,   0,   0,  75,  76,  77,  78;
//       0,   0,   0,   0,   0,  86,  87,  88;
//       0,   0,   0,   0,   0,   0,  97,  98]
static const int am  = 9;
static const int an  = 8;
static const int al = 2;
static const int au = 3;
static const int la = 6;
static const int a[] = { 555, 555, 555,  11,  21,  31,
                         555, 555,  12,  22,  32,  42,
                         555,  13,  23,  33,  43,  53,
                          14,  24,  34,  44,  54,  64,
                          25,  35,  45,  55,  65,  75,
                          36,  46,  56,  66,  76,  86,
                          47,  57,  67,  77,  87,  97,
                          58,  68,  78,  88,  98, 555 };

// Test Matrix B
//      [11,  12,  13,  14,   0,   0,   0,   0,   0,
//       21,  22,  23,  24,  25,   0,   0,   0,   0,
//       31,  32,  33,  34,  35,  36,   0,   0,   0,
//        0,  42,  43,  44,  45,  46,  47,   0,   0,
//        0,   0,  53,  54,  55,  56,  57,  58,   0,
//        0,   0,   0,  64,  65,  66,  67,  68,  69,
//        0,   0,   0,   0,  75,  76,  77,  78,  79 ];
static const int bm  = 7;
static const int bn  = 9;
static const int bl = 2;
static const int bu = 3;
static const int lb = 6;
static const int b[] = { 555, 555, 555,  11,  21,  31,
                         555, 555,  12,  22,  32,  42,
                         555,  13,  23,  33,  43,  53,
                          14,  24,  34,  44,  54,  64,
                          25,  35,  45,  55,  65,  75,
                          36,  46,  56,  66,  76, 555,
                          47,  57,  67,  77, 555, 555,
                          58,  68,  78, 555, 555, 555,
                          69,  79, 555, 555, 555, 555 };

// Shorthand for routines to be tested
#define offset  suzerain_gbmatrix_offset
#define in_band suzerain_gbmatrix_in_band
#define col     suzerain_gbmatrix_col
#define row     suzerain_gbmatrix_row

static void test_suzerain_gbmatrix_offset()
{
    // Matrix A

    gsl_test_int(a[offset(la, al, au, 0, 0)], 11, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 0, 1)], 12, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 0, 2)], 13, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 0, 3)], 14, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 1, 0)], 21, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 1, 1)], 22, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 1, 2)], 23, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 1, 3)], 24, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 1, 4)], 25, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 2, 0)], 31, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 2, 1)], 32, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 2, 2)], 33, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 2, 3)], 34, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 2, 4)], 35, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 2, 5)], 36, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 3, 1)], 42, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 3, 2)], 43, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 3, 3)], 44, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 3, 4)], 45, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 3, 5)], 46, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 3, 6)], 47, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 4, 2)], 53, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 4, 3)], 54, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 4, 4)], 55, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 4, 5)], 56, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 4, 6)], 57, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 4, 7)], 58, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 5, 3)], 64, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 5, 4)], 65, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 5, 5)], 66, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 5, 6)], 67, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 5, 7)], 68, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 6, 4)], 75, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 6, 5)], 76, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 6, 6)], 77, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 6, 7)], 78, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 7, 5)], 86, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 7, 6)], 87, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 7, 7)], 88, "line %d", __LINE__);

    gsl_test_int(a[offset(la, al, au, 8, 6)], 97, "line %d", __LINE__);
    gsl_test_int(a[offset(la, al, au, 8, 7)], 98, "line %d", __LINE__);


    // Matrix B

    gsl_test_int(b[offset(lb, bl, bu, 0, 0)], 11, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 0, 1)], 12, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 0, 2)], 13, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 0, 3)], 14, "line %d", __LINE__);

    gsl_test_int(b[offset(lb, bl, bu, 1, 0)], 21, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 1, 1)], 22, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 1, 2)], 23, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 1, 3)], 24, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 1, 4)], 25, "line %d", __LINE__);

    gsl_test_int(b[offset(lb, bl, bu, 2, 0)], 31, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 2, 1)], 32, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 2, 2)], 33, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 2, 3)], 34, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 2, 4)], 35, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 2, 5)], 36, "line %d", __LINE__);

    gsl_test_int(b[offset(lb, bl, bu, 3, 1)], 42, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 3, 2)], 43, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 3, 3)], 44, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 3, 4)], 45, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 3, 5)], 46, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 3, 6)], 47, "line %d", __LINE__);

    gsl_test_int(b[offset(lb, bl, bu, 4, 2)], 53, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 4, 3)], 54, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 4, 4)], 55, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 4, 5)], 56, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 4, 6)], 57, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 4, 7)], 58, "line %d", __LINE__);

    gsl_test_int(b[offset(lb, bl, bu, 5, 3)], 64, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 5, 4)], 65, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 5, 5)], 66, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 5, 6)], 67, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 5, 7)], 68, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 5, 8)], 69, "line %d", __LINE__);

    gsl_test_int(b[offset(lb, bl, bu, 6, 4)], 75, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 6, 5)], 76, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 6, 6)], 77, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 6, 7)], 78, "line %d", __LINE__);
    gsl_test_int(b[offset(lb, bl, bu, 6, 8)], 79, "line %d", __LINE__);
}

static void test_suzerain_gbmatrix_in_band()
{
    // Matrix A

    gsl_test_int(in_band(al, au, 0, 0), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 0, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 0, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 0, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 0, 4), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 0, 5), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 0, 6), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 0, 7), 0, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 1, 0), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 1, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 1, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 1, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 1, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 1, 5), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 1, 6), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 1, 7), 0, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 2, 0), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 2, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 2, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 2, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 2, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 2, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 2, 6), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 2, 7), 0, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 3, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 3, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 3, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 3, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 3, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 3, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 3, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 3, 7), 0, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 4, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 4, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 4, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 4, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 4, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 4, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 4, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 4, 7), 1, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 5, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 5, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 5, 2), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 5, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 5, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 5, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 5, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 5, 7), 1, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 6, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 6, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 6, 2), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 6, 3), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 6, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 6, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 6, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 6, 7), 1, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 7, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 7, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 7, 2), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 7, 3), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 7, 4), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 7, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 7, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 7, 7), 1, "line %d", __LINE__);

    gsl_test_int(in_band(al, au, 8, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 8, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 8, 2), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 8, 3), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 8, 4), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 8, 5), 0, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 8, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(al, au, 8, 7), 1, "line %d", __LINE__);


    // Matrix B

    gsl_test_int(in_band(bl, bu, 0, 0), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 4), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 5), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 6), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 7), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 0, 8), 0, "line %d", __LINE__);

    gsl_test_int(in_band(bl, bu, 1, 0), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 5), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 6), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 7), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 1, 8), 0, "line %d", __LINE__);

    gsl_test_int(in_band(bl, bu, 2, 0), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 6), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 7), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 2, 8), 0, "line %d", __LINE__);

    gsl_test_int(in_band(bl, bu, 3, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 1), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 7), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 3, 8), 0, "line %d", __LINE__);

    gsl_test_int(in_band(bl, bu, 4, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 2), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 7), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 4, 8), 0, "line %d", __LINE__);

    gsl_test_int(in_band(bl, bu, 5, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 2), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 3), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 7), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 5, 8), 1, "line %d", __LINE__);

    gsl_test_int(in_band(bl, bu, 6, 0), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 1), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 2), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 3), 0, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 4), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 5), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 6), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 7), 1, "line %d", __LINE__);
    gsl_test_int(in_band(bl, bu, 6, 8), 1, "line %d", __LINE__);
}

static void test_suzerain_gbmatrix_col()
{
    int * c;
    int il, iu;

    // Matrix A

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 0, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 3, "line %d", __LINE__);
    gsl_test_int(c[0], 11, "line %d", __LINE__);
    gsl_test_int(c[1], 21, "line %d", __LINE__);
    gsl_test_int(c[2], 31, "line %d", __LINE__);

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 1, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 4, "line %d", __LINE__);
    gsl_test_int(c[0], 12, "line %d", __LINE__);
    gsl_test_int(c[1], 22, "line %d", __LINE__);
    gsl_test_int(c[2], 32, "line %d", __LINE__);
    gsl_test_int(c[3], 42, "line %d", __LINE__);

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 2, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 5, "line %d", __LINE__);
    gsl_test_int(c[0], 13, "line %d", __LINE__);
    gsl_test_int(c[1], 23, "line %d", __LINE__);
    gsl_test_int(c[2], 33, "line %d", __LINE__);
    gsl_test_int(c[3], 43, "line %d", __LINE__);
    gsl_test_int(c[4], 53, "line %d", __LINE__);

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 3, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 6, "line %d", __LINE__);
    gsl_test_int(c[0], 14, "line %d", __LINE__);
    gsl_test_int(c[1], 24, "line %d", __LINE__);
    gsl_test_int(c[2], 34, "line %d", __LINE__);
    gsl_test_int(c[3], 44, "line %d", __LINE__);
    gsl_test_int(c[4], 54, "line %d", __LINE__);
    gsl_test_int(c[5], 64, "line %d", __LINE__);

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 4, &il, &iu);
    gsl_test_int(il, 1, "line %d", __LINE__);
    gsl_test_int(iu, 7, "line %d", __LINE__);
    gsl_test_int(c[1], 25, "line %d", __LINE__);
    gsl_test_int(c[2], 35, "line %d", __LINE__);
    gsl_test_int(c[3], 45, "line %d", __LINE__);
    gsl_test_int(c[4], 55, "line %d", __LINE__);
    gsl_test_int(c[5], 65, "line %d", __LINE__);
    gsl_test_int(c[6], 75, "line %d", __LINE__);

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 5, &il, &iu);
    gsl_test_int(il, 2, "line %d", __LINE__);
    gsl_test_int(iu, 8, "line %d", __LINE__);
    gsl_test_int(c[2], 36, "line %d", __LINE__);
    gsl_test_int(c[3], 46, "line %d", __LINE__);
    gsl_test_int(c[4], 56, "line %d", __LINE__);
    gsl_test_int(c[5], 66, "line %d", __LINE__);
    gsl_test_int(c[6], 76, "line %d", __LINE__);
    gsl_test_int(c[7], 86, "line %d", __LINE__);

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 6, &il, &iu);
    gsl_test_int(il, 3, "line %d", __LINE__);
    gsl_test_int(iu, 9, "line %d", __LINE__);
    gsl_test_int(c[3], 47, "line %d", __LINE__);
    gsl_test_int(c[4], 57, "line %d", __LINE__);
    gsl_test_int(c[5], 67, "line %d", __LINE__);
    gsl_test_int(c[6], 77, "line %d", __LINE__);
    gsl_test_int(c[7], 87, "line %d", __LINE__);
    gsl_test_int(c[8], 97, "line %d", __LINE__);

    c = col(am, an, al, au, (void *)a, la, sizeof(int), 7, &il, &iu);
    gsl_test_int(il, 4, "line %d", __LINE__);
    gsl_test_int(iu, 9, "line %d", __LINE__);
    gsl_test_int(c[4], 58, "line %d", __LINE__);
    gsl_test_int(c[5], 68, "line %d", __LINE__);
    gsl_test_int(c[6], 78, "line %d", __LINE__);
    gsl_test_int(c[7], 88, "line %d", __LINE__);
    gsl_test_int(c[8], 98, "line %d", __LINE__);

    // Matrix B

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 0, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 3, "line %d", __LINE__);
    gsl_test_int(c[0], 11, "line %d", __LINE__);
    gsl_test_int(c[1], 21, "line %d", __LINE__);
    gsl_test_int(c[2], 31, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 1, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 4, "line %d", __LINE__);
    gsl_test_int(c[0], 12, "line %d", __LINE__);
    gsl_test_int(c[1], 22, "line %d", __LINE__);
    gsl_test_int(c[2], 32, "line %d", __LINE__);
    gsl_test_int(c[3], 42, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 2, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 5, "line %d", __LINE__);
    gsl_test_int(c[0], 13, "line %d", __LINE__);
    gsl_test_int(c[1], 23, "line %d", __LINE__);
    gsl_test_int(c[2], 33, "line %d", __LINE__);
    gsl_test_int(c[3], 43, "line %d", __LINE__);
    gsl_test_int(c[4], 53, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 3, &il, &iu);
    gsl_test_int(il, 0, "line %d", __LINE__);
    gsl_test_int(iu, 6, "line %d", __LINE__);
    gsl_test_int(c[0], 14, "line %d", __LINE__);
    gsl_test_int(c[1], 24, "line %d", __LINE__);
    gsl_test_int(c[2], 34, "line %d", __LINE__);
    gsl_test_int(c[3], 44, "line %d", __LINE__);
    gsl_test_int(c[4], 54, "line %d", __LINE__);
    gsl_test_int(c[5], 64, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 4, &il, &iu);
    gsl_test_int(il, 1, "line %d", __LINE__);
    gsl_test_int(iu, 7, "line %d", __LINE__);
    gsl_test_int(c[1], 25, "line %d", __LINE__);
    gsl_test_int(c[2], 35, "line %d", __LINE__);
    gsl_test_int(c[3], 45, "line %d", __LINE__);
    gsl_test_int(c[4], 55, "line %d", __LINE__);
    gsl_test_int(c[5], 65, "line %d", __LINE__);
    gsl_test_int(c[6], 75, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 5, &il, &iu);
    gsl_test_int(il, 2, "line %d", __LINE__);
    gsl_test_int(iu, 7, "line %d", __LINE__);
    gsl_test_int(c[2], 36, "line %d", __LINE__);
    gsl_test_int(c[3], 46, "line %d", __LINE__);
    gsl_test_int(c[4], 56, "line %d", __LINE__);
    gsl_test_int(c[5], 66, "line %d", __LINE__);
    gsl_test_int(c[6], 76, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 6, &il, &iu);
    gsl_test_int(il, 3, "line %d", __LINE__);
    gsl_test_int(iu, 7, "line %d", __LINE__);
    gsl_test_int(c[3], 47, "line %d", __LINE__);
    gsl_test_int(c[4], 57, "line %d", __LINE__);
    gsl_test_int(c[5], 67, "line %d", __LINE__);
    gsl_test_int(c[6], 77, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 7, &il, &iu);
    gsl_test_int(il, 4, "line %d", __LINE__);
    gsl_test_int(iu, 7, "line %d", __LINE__);
    gsl_test_int(c[4], 58, "line %d", __LINE__);
    gsl_test_int(c[5], 68, "line %d", __LINE__);
    gsl_test_int(c[6], 78, "line %d", __LINE__);

    c = col(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 8, &il, &iu);
    gsl_test_int(il, 5, "line %d", __LINE__);
    gsl_test_int(iu, 7, "line %d", __LINE__);
    gsl_test_int(c[5], 69, "line %d", __LINE__);
    gsl_test_int(c[6], 79, "line %d", __LINE__);
}

static void test_suzerain_gbmatrix_row()
{
    int * r;
    int jl, ju, inc;

    // Matrix A

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 0, &jl, &ju, &inc);
    gsl_test_int(jl, 0, "line %d", __LINE__);
    gsl_test_int(ju, 4, "line %d", __LINE__);
    gsl_test_int(r[0*inc], 11, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 12, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 13, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 14, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 1, &jl, &ju, &inc);
    gsl_test_int(jl, 0, "line %d", __LINE__);
    gsl_test_int(ju, 5, "line %d", __LINE__);
    gsl_test_int(r[0*inc], 21, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 22, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 23, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 24, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 25, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 2, &jl, &ju, &inc);
    gsl_test_int(jl, 0, "line %d", __LINE__);
    gsl_test_int(ju, 6, "line %d", __LINE__);
    gsl_test_int(r[0*inc], 31, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 32, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 33, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 34, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 35, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 36, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 3, &jl, &ju, &inc);
    gsl_test_int(jl, 1, "line %d", __LINE__);
    gsl_test_int(ju, 7, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 42, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 43, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 44, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 45, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 46, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 47, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 4, &jl, &ju, &inc);
    gsl_test_int(jl, 2, "line %d", __LINE__);
    gsl_test_int(ju, 8, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 53, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 54, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 55, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 56, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 57, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 58, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 5, &jl, &ju, &inc);
    gsl_test_int(jl, 3, "line %d", __LINE__);
    gsl_test_int(ju, 8, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 64, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 65, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 66, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 67, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 68, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 6, &jl, &ju, &inc);
    gsl_test_int(jl, 4, "line %d", __LINE__);
    gsl_test_int(ju, 8, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 75, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 76, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 77, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 78, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 7, &jl, &ju, &inc);
    gsl_test_int(jl, 5, "line %d", __LINE__);
    gsl_test_int(ju, 8, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 86, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 87, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 88, "line %d", __LINE__);

    r = row(am, an, al, au, (void *)a, la, sizeof(int), 8, &jl, &ju, &inc);
    gsl_test_int(jl, 6, "line %d", __LINE__);
    gsl_test_int(ju, 8, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 97, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 98, "line %d", __LINE__);


    // Matrix B

    r = row(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 0, &jl, &ju, &inc);
    gsl_test_int(jl, 0, "line %d", __LINE__);
    gsl_test_int(ju, 4, "line %d", __LINE__);
    gsl_test_int(r[0*inc], 11, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 12, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 13, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 14, "line %d", __LINE__);

    r = row(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 1, &jl, &ju, &inc);
    gsl_test_int(jl, 0, "line %d", __LINE__);
    gsl_test_int(ju, 5, "line %d", __LINE__);
    gsl_test_int(r[0*inc], 21, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 22, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 23, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 24, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 25, "line %d", __LINE__);

    r = row(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 2, &jl, &ju, &inc);
    gsl_test_int(jl, 0, "line %d", __LINE__);
    gsl_test_int(ju, 6, "line %d", __LINE__);
    gsl_test_int(r[0*inc], 31, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 32, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 33, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 34, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 35, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 36, "line %d", __LINE__);

    r = row(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 3, &jl, &ju, &inc);
    gsl_test_int(jl, 1, "line %d", __LINE__);
    gsl_test_int(ju, 7, "line %d", __LINE__);
    gsl_test_int(r[1*inc], 42, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 43, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 44, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 45, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 46, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 47, "line %d", __LINE__);

    r = row(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 4, &jl, &ju, &inc);
    gsl_test_int(jl, 2, "line %d", __LINE__);
    gsl_test_int(ju, 8, "line %d", __LINE__);
    gsl_test_int(r[2*inc], 53, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 54, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 55, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 56, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 57, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 58, "line %d", __LINE__);

    r = row(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 5, &jl, &ju, &inc);
    gsl_test_int(jl, 3, "line %d", __LINE__);
    gsl_test_int(ju, 9, "line %d", __LINE__);
    gsl_test_int(r[3*inc], 64, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 65, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 66, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 67, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 68, "line %d", __LINE__);
    gsl_test_int(r[8*inc], 69, "line %d", __LINE__);

    r = row(bm, bn, bl, bu, (void *)b, lb, sizeof(int), 6, &jl, &ju, &inc);
    gsl_test_int(jl, 4, "line %d", __LINE__);
    gsl_test_int(ju, 9, "line %d", __LINE__);
    gsl_test_int(r[4*inc], 75, "line %d", __LINE__);
    gsl_test_int(r[5*inc], 76, "line %d", __LINE__);
    gsl_test_int(r[6*inc], 77, "line %d", __LINE__);
    gsl_test_int(r[7*inc], 78, "line %d", __LINE__);
    gsl_test_int(r[8*inc], 79, "line %d", __LINE__);
}

int main()
{
    test_suzerain_gbmatrix_offset();
    test_suzerain_gbmatrix_in_band();
    test_suzerain_gbmatrix_col();
    test_suzerain_gbmatrix_row();

    exit(gsl_test_summary());
}
