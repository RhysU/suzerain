#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/gbmatrix.h>
#include <gsl/gsl_test.h>

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

int main()
{
    test_suzerain_gbmatrix_offset();

    exit(gsl_test_summary());
}
