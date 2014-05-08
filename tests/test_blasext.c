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

#include <suzerain/blas_et_al/blasext.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#ifdef SUZERAIN_HAVE_MKL
#include <mkl.h>
#include <mkl_service.h>
#endif

static
void
test_dgb_diag_scale_acc1(const double beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const double alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    double b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const double d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const double expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(0,10,0.10), R(0,20,0.20),
          -2.2,  -1.1, R(1, 1,0.01), R(1,11,0.11), R(1,21,0.21), R(1,31,0.31),
          -2.2,  -1.1, R(2,12,0.12), R(2,22,0.22), R(2,32,0.32), R(2,42,0.42),
          -2.2,  -1.1, R(3,23,0.23), R(3,33,0.33), R(3,43,0.43), R(3,53,0.53),
          -2.2,  -1.1, R(4,34,0.34), R(4,44,0.44), R(4,54,0.54), R(4,64,0.64),
          -2.2,  -1.1, R(5,45,0.45), R(5,55,0.55), R(5,65,0.65), R(5,75,0.75),
          -2.2,  -1.1, R(6,56,0.56), R(6,66,0.66), R(6,76,0.76), R(6,86,0.86),
          -2.2,  -1.1, R(7,67,0.67), R(7,77,0.77), R(7,87,0.87), R(7,97,0.97),
          -2.2,  -1.1, R(8,78,0.78), R(8,88,0.88), R(8,98,0.98),           -1,
          -2.2,  -1.1, R(9,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_dgb_diag_scale_acc('R', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_dgb_diag_scale_acc2(const double beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const double alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const double a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    double b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const double d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const double expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(1,10,0.10), R(2,20,0.20),
          -2.2,  -1.1, R(0, 1,0.01), R(1,11,0.11), R(2,21,0.21), R(3,31,0.31),
          -2.2,  -1.1, R(1,12,0.12), R(2,22,0.22), R(3,32,0.32), R(4,42,0.42),
          -2.2,  -1.1, R(2,23,0.23), R(3,33,0.33), R(4,43,0.43), R(5,53,0.53),
          -2.2,  -1.1, R(3,34,0.34), R(4,44,0.44), R(5,54,0.54), R(6,64,0.64),
          -2.2,  -1.1, R(4,45,0.45), R(5,55,0.55), R(6,65,0.65), R(7,75,0.75),
          -2.2,  -1.1, R(5,56,0.56), R(6,66,0.66), R(7,76,0.76), R(8,86,0.86),
          -2.2,  -1.1, R(6,67,0.67), R(7,77,0.77), R(8,87,0.87), R(9,97,0.97),
          -2.2,  -1.1, R(7,78,0.78), R(8,88,0.88), R(9,98,0.98),           -1,
          -2.2,  -1.1, R(8,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const double *a = a_data + lda-(ku+1+kl);
    double       *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_dgb_diag_scale_acc('L', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_DBL_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_sgb_diag_scale_acc1(const float  beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const float  alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float  a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    float  b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const float  d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const float  expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(0,10,0.10), R(0,20,0.20),
          -2.2,  -1.1, R(1, 1,0.01), R(1,11,0.11), R(1,21,0.21), R(1,31,0.31),
          -2.2,  -1.1, R(2,12,0.12), R(2,22,0.22), R(2,32,0.32), R(2,42,0.42),
          -2.2,  -1.1, R(3,23,0.23), R(3,33,0.33), R(3,43,0.43), R(3,53,0.53),
          -2.2,  -1.1, R(4,34,0.34), R(4,44,0.44), R(4,54,0.54), R(4,64,0.64),
          -2.2,  -1.1, R(5,45,0.45), R(5,55,0.55), R(5,65,0.65), R(5,75,0.75),
          -2.2,  -1.1, R(6,56,0.56), R(6,66,0.66), R(6,76,0.76), R(6,86,0.86),
          -2.2,  -1.1, R(7,67,0.67), R(7,77,0.77), R(7,87,0.87), R(7,97,0.97),
          -2.2,  -1.1, R(8,78,0.78), R(8,88,0.88), R(8,98,0.98),           -1,
          -2.2,  -1.1, R(9,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const float  *a = a_data + lda-(ku+1+kl);
    float        *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_sgb_diag_scale_acc('R', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_sgb_diag_scale_acc2(const float  beta)
{
    int i;

    const int    m     = 10;
    const int    n     = 10;
    const int    ku    = 1;
    const int    kl    = 2;
    const int    lda   = 5;
    const int    ldb   = 6;
    const float  alpha = 2.0;

    /* Negative one represents values outside of the band */
    /* Decimal parts flag regions outside the matrix entirely */
    const float  a_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
                 -1.1,     -1,       0,    10,      20,
                 -1.1,      1,      11,    21,      31,
                 -1.1,     12,      22,    32,      42,
                 -1.1,     23,      33,    43,      53,
                 -1.1,     34,      44,    54,      64,
                 -1.1,     45,      55,    65,      75,
                 -1.1,     56,      66,    76,      86,
                 -1.1,     67,      77,    87,      97,
                 -1.1,     78,      88,    98,      -1,
                 -1.1,     89,      99,    -1,      -1
    };

    float  b_data[] = {
        /*lda buffer*/ /*ku1*/ /*diag*/ /*kl1*/ /*kl2*/
          -2.2,  -1.1,     -1,    0.00,  0.10,    0.20,
          -2.2,  -1.1,   0.01,    0.11,  0.21,    0.31,
          -2.2,  -1.1,   0.12,    0.22,  0.32,    0.42,
          -2.2,  -1.1,   0.23,    0.33,  0.43,    0.53,
          -2.2,  -1.1,   0.34,    0.44,  0.54,    0.64,
          -2.2,  -1.1,   0.45,    0.55,  0.65,    0.75,
          -2.2,  -1.1,   0.56,    0.66,  0.76,    0.86,
          -2.2,  -1.1,   0.67,    0.77,  0.87,    0.97,
          -2.2,  -1.1,   0.78,    0.88,  0.98,      -1,
          -2.2,  -1.1,   0.89,    0.99,    -1,      -1
    };

    const float  d_data[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

#define R(d,a,b) (alpha*d_data[d]*a + beta*b)
    const float  expected_data[] = {
        /*lda buffer*/       /*ku1*/      /*diag*/       /*kl1*/       /*kl2*/
          -2.2,  -1.1,           -1, R(0, 0,0.00), R(1,10,0.10), R(2,20,0.20),
          -2.2,  -1.1, R(0, 1,0.01), R(1,11,0.11), R(2,21,0.21), R(3,31,0.31),
          -2.2,  -1.1, R(1,12,0.12), R(2,22,0.22), R(3,32,0.32), R(4,42,0.42),
          -2.2,  -1.1, R(2,23,0.23), R(3,33,0.33), R(4,43,0.43), R(5,53,0.53),
          -2.2,  -1.1, R(3,34,0.34), R(4,44,0.44), R(5,54,0.54), R(6,64,0.64),
          -2.2,  -1.1, R(4,45,0.45), R(5,55,0.55), R(6,65,0.65), R(7,75,0.75),
          -2.2,  -1.1, R(5,56,0.56), R(6,66,0.66), R(7,76,0.76), R(8,86,0.86),
          -2.2,  -1.1, R(6,67,0.67), R(7,77,0.77), R(8,87,0.87), R(9,97,0.97),
          -2.2,  -1.1, R(7,78,0.78), R(8,88,0.88), R(9,98,0.98),           -1,
          -2.2,  -1.1, R(8,89,0.89), R(9,99,0.99),         -1,             -1
    };
#undef R

    const float  *a = a_data + lda-(ku+1+kl);
    float        *b = b_data + ldb-(ku+1+kl);

    const int nb        = sizeof(b_data)/sizeof(b_data[0]);
    const int nexpected = sizeof(expected_data)/sizeof(expected_data[0]);
    gsl_test_int(nb, nexpected, "Expected results' length");

    suzerain_blasext_sgb_diag_scale_acc('L', m, n, kl, ku,
                                        alpha, d_data, 1, a, 1, lda,
                                        beta,             b, 1, ldb);
    for (i = 0; i < nexpected; ++i) {
        gsl_test_abs(b_data[i], expected_data[i], GSL_FLT_EPSILON,
                "%s:%d index %d", __func__, __LINE__, i);
    }
}

static
void
test_blasext_i2s_zaxpby2()
{
    { /* Test pure copy for single, stride one vectors */

        const int N = 4;
        double z[N*2];
        double z_re[N];
        double z_im[N];

        for (int i = 0; i < N; ++i) {
            z[2*i]   = i+1;   /* Real part */
            z[2*i+1] = i+1.5; /* Imag part */
        }

        suzerain_blasext_i2s_zaxpby2(
                1, N, NULL, z, 1, N,
                NULL, z_re, 1, N, z_im, 1, N);

        for (int i = 0; i < N; ++i) {
            gsl_test_abs(z_re[i], i+1, GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
            gsl_test_abs(z_im[i], i+1.5, GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
        }
    }

    { /* Test pure copy for multiple, strided vectors */
        const int M       = 3;
        const int N       = 4;
        const int incz    = 2;
        const int incz_re = 3;
        const int incz_im = 4;
        const int ldz     = N*incz + 2;
        const int ldz_re  = N*incz_re + 3;
        const int ldz_im  = N*incz_im + 5;
        double z[2*(N*incz)*(M*ldz)];
        double z_re[(N*incz_re)*(M*ldz_re)];
        double z_im[(N*incz_im)*(M*ldz_im)];

        memset(z_re, 'r', sizeof(z_re)/sizeof(z_re[0])); // Avoid NaNs
        memset(z_im, 'r', sizeof(z_im)/sizeof(z_im[0])); // Avoid NaNs

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                z[2*(i*incz + j*ldz)]   = i+1;   /* Real part */
                z[2*(i*incz + j*ldz)+1] = i+1.5; /* Imag part */
            }
        }

        suzerain_blasext_i2s_zaxpby2(
                M, N, NULL,
                z, incz, ldz,
                NULL,
                z_re, incz_re, ldz_re,
                z_im, incz_im, ldz_im);

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                gsl_test_abs(z_re[i*incz_re+j*ldz_re], i+1, GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im], i+1.5, GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
            }
        }
    }

    { /* Test real scaling for single, stride one vectors */
        const int N = 4;
        const double alpha[2] = {2.0, 0.0};
        double z[2*N];
        double z_re[N];
        double z_im[N];

        for (int i = 0; i < N; ++i) {
            z[2*i]   = i+1;   /* Real part */
            z[2*i+1] = i+1.5; /* Imag part */
        }

        suzerain_blasext_i2s_zaxpby2(
                1, N, alpha, z, 1, N,
                NULL,
                z_re, 1, N, z_im, 1, N);

        for (int i = 0; i < N; ++i) {
            gsl_test_abs(z_re[i], 2*(i+1), GSL_DBL_EPSILON,
                        "%s:%d real part of index",
                        __func__, __LINE__, i);
            gsl_test_abs(z_im[i], 2*(i+1.5), GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
        }
    }

    { /* Test real scaling for multiple, strided vectors */
        const int M       = 3;
        const int N       = 4;
        const double alpha[2] = {2.0, 0.0};
        const int incz    = 2;
        const int incz_re = 3;
        const int incz_im = 4;
        const int ldz     = N*incz + 2;
        const int ldz_re  = N*incz_re + 3;
        const int ldz_im  = N*incz_im + 5;
        double z[2*(N*incz)*(M*ldz)];
        double z_re[(N*incz_re)*(M*ldz_re)];
        double z_im[(N*incz_im)*(M*ldz_im)];

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                z[2*(i*incz + j*ldz)]   = i+1;   /* Real part */
                z[2*(i*incz + j*ldz)+1] = i+1.5; /* Imag part */
            }
        }

        suzerain_blasext_i2s_zaxpby2(
                M, N, alpha,
                z, incz, ldz,
                NULL,
                z_re, incz_re, ldz_re,
                z_im, incz_im, ldz_im);

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                gsl_test_abs(z_re[i*incz_re+j*ldz_re],
                        2*(i+1), GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im],
                        2*(i+1.5), GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
            }
        }
    }

    { /* Test complex scaling for contiguous vectors */
        const int N = 4;
        const double alpha[2] = {3.0, 2.0};
        const double beta[2]  = {4.0, 5.0};
        double z[2*N];
        double z_re[N];
        double z_im[N];

        for (int i = 0; i < N; ++i) {
            z[2*(i)] = (i+1)*(i+1);
            z[2*(i)+1] = i+1.5;
            z_re[i] = i+7;
            z_im[i] = i+11;
        }

        suzerain_blasext_i2s_zaxpby2(
                1, N, alpha, z, 1, N,
                beta, z_re, 1, N, z_im, 1, N);

        for (int i = 0; i < N; ++i) {
            gsl_test_abs(z_re[i], 3*i*i + 3*i - 27,
                        GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
            gsl_test_abs(z_im[i], 2*i*i + 16*i + 85.5,
                        GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
        }
    }

    { /* Test complex scaling for multiple, strided vectors */
        const int M       = 3;
        const int N       = 4;
        const double alpha[2] = {3.0, 2.0};
        const double beta[2]  = {4.0, 5.0};
        const int incz    = 2;
        const int incz_re = 3;
        const int incz_im = 4;
        const int ldz     = N*incz + 2;
        const int ldz_re  = N*incz_re + 3;
        const int ldz_im  = N*incz_im + 5;
        double z[2*(N*incz)*(M*ldz)];
        double z_re[(N*incz_re)*(M*ldz_re)];
        double z_im[(N*incz_im)*(M*ldz_im)];

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                z[2*(i*incz + j*ldz)] = (i+1)*(i+1);
                z[2*(i*incz + j*ldz)+1] = i+1.5;
                z_re[i*incz_re + j*ldz_re] = i+7;
                z_im[i*incz_im + j*ldz_im] = i+11;
            }
        }

        suzerain_blasext_i2s_zaxpby2(
                M, N, alpha,
                z, incz, ldz,
                beta,
                z_re, incz_re, ldz_re,
                z_im, incz_im, ldz_im);

        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                gsl_test_abs(z_re[i*incz_re+j*ldz_re], 3*i*i + 3*i - 27,
                        GSL_DBL_EPSILON,
                        "%s:%d real part of index %d",
                        __func__, __LINE__, i);
                gsl_test_abs(z_im[i*incz_im+j*ldz_im], 2*i*i + 16*i + 85.5,
                        GSL_DBL_EPSILON,
                        "%s:%d imag part of index %d",
                        __func__, __LINE__, i);
            }
        }
    }
}

static
void
test_blasext_sgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const float a[] = { 555, 555, 555, 555, 6,   -1,  -2,  -4,
                        555, 555, 555, -4,  10,  1,   -5,  555,
                        555, 555, 10,  10,  -7,  5,   555, 555,
                        555, -10, -1,  -9,  -2,  555, 555, 555,
                        -7,  3,   -5,  3,   555, 555, 555, 555 };

    float norm1;
    const int status = suzerain_blasext_sgbnorm1(m, n, kl, ku, a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, 32.0, GSL_FLT_EPSILON, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
void
test_blasext_dgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const double a[] = { 555, 555, 555, 555, 6,   -1,  -2,  -4,
                         555, 555, 555, -4,  10,  1,   -5,  555,
                         555, 555, 10,  10,  -7,  5,   555, 555,
                         555, -10, -1,  -9,  -2,  555, 555, 555,
                         -7,  3,   -5,  3,   555, 555, 555, 555 };

    double norm1;
    const int status = suzerain_blasext_dgbnorm1(m, n, kl, ku, a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, 32.0, GSL_DBL_EPSILON, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
void
test_blasext_cgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const float a[][2] = {
        {555,0},{555,0}, {555,0},{555,0}, {6,3}, {-1,5},  {-2,5},  {-4,-5},
        {555,0},{555,0}, {555,0},{-4,-4},{10,-7},{1,-2},  {-5,-7}, {555,0},
        {555,0},{555,0}, {10,-2},{10,-4},{-7,7}, {5,6},   {555,0}, {555,0},
        {555,0},{-10,-4},{-1,2}, {-9,0}, {-2,-7},{555,0}, {555,0}, {555,0},
        {-7,7}, {3,6},   {-5,-7},{3,-3}, {555,0},{555,0}, {555,0}, {555,0}
    };
    const float expected = 7.*sqrt(2.)+2.*sqrt(26.)+2.*sqrt(29.)+sqrt(61.);

    float norm1;
    const int status = suzerain_blasext_cgbnorm1(
            m, n, kl, ku, (complex_float *) a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, expected, GSL_FLT_EPSILON, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
void
test_blasext_zgbnorm1()
{
    const int m = 4, n = 5, kl = 3, ku = 4, lda = ku + 1 + kl;
    const double a[][2] = {
        {555,0},{555,0}, {555,0},{555,0}, {6,3}, {-1,5},  {-2,5},  {-4,-5},
        {555,0},{555,0}, {555,0},{-4,-4},{10,-7},{1,-2},  {-5,-7}, {555,0},
        {555,0},{555,0}, {10,-2},{10,-4},{-7,7}, {5,6},   {555,0}, {555,0},
        {555,0},{-10,-4},{-1,2}, {-9,0}, {-2,-7},{555,0}, {555,0}, {555,0},
        {-7,7}, {3,6},   {-5,-7},{3,-3}, {555,0},{555,0}, {555,0}, {555,0}
    };
    const double expected = 7.*sqrt(2.)+2.*sqrt(26.)+2.*sqrt(29.)+sqrt(61.);

    double norm1;
    const int status = suzerain_blasext_zgbnorm1(
            m, n, kl, ku, (complex_double *) a, lda, &norm1);
    gsl_test_int(status, 0, "%s:%d call success", __func__, __LINE__);
    gsl_test_abs(norm1, expected, GSL_DBL_EPSILON*100, "%s:%d norm result %d",
                 __func__, __LINE__, norm1);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_ddemote()
{
    int info;
    double x[]   = {1, 2, 3, 4, 5, 6, 7};
    float  e[]   = {1, 2, 3, 4, 5, 6, 7};
    const int n = sizeof(x)/sizeof(x[0]);

    info = suzerain_blasext_ddemote(n, x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_dpromote()
{
    int info;
    double  e[]   = {1, 2, 3, 4, 5, 6, 7};
    float   x[]   = {1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 0, 0, 0};
    const int n = sizeof(e)/sizeof(e[0]);
    gsl_test_int(2*n, sizeof(x)/sizeof(x[0]),
            "%s.%d consistency", __func__, __LINE__);

    info = suzerain_blasext_dpromote(n, (void *) x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_zdemote()
{
    int info;
    double x[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7}};
    float  e[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7}};
    const int n = sizeof(x)/sizeof(x[0]);

    info = suzerain_blasext_zdemote(n, (void *) x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

static
SUZERAIN_DONTINLINE
void
test_blasext_zpromote()
{
    int info;
    float  x[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7},
                       {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
    double e[][2]   = {{1,-1}, {2,-2}, {3,-3}, {4,-4}, {5,-5}, {6,-6}, {7,-7}};
    const int n = sizeof(e)/sizeof(e[0]);
    gsl_test_int(2*n, sizeof(x)/sizeof(x[0]),
            "%s.%d consistency", __func__, __LINE__);

    info = suzerain_blasext_zpromote(n, (void *) x);
    gsl_test_int(info, 0, "%s.%d report success", __func__, __LINE__);

    info = memcmp(x, e, sizeof(e));
    gsl_test_int(info, 0, "%s.%d matches expected", __func__, __LINE__);
}

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    gsl_ieee_env_setup();

    test_dgb_diag_scale_acc1(/* beta = */ 6.0);
    test_dgb_diag_scale_acc1(/* beta = */ 1.0);
    test_dgb_diag_scale_acc2(/* beta = */ 6.0);
    test_dgb_diag_scale_acc2(/* beta = */ 1.0);

    test_sgb_diag_scale_acc1(/* beta = */ 6.0);
    test_sgb_diag_scale_acc1(/* beta = */ 1.0);
    test_sgb_diag_scale_acc2(/* beta = */ 6.0);
    test_sgb_diag_scale_acc2(/* beta = */ 1.0);

    /* TODO Add suzerain_blasext_d{s,g}bmzv test cases                  */
    /* suzerain_blasext_d{s,g}bmzv exercised somewhat in test_bsplineop */
    /* via accumulate_complex and apply_complex                         */

    test_blasext_i2s_zaxpby2();
    test_blasext_sgbnorm1();
    test_blasext_dgbnorm1();
    test_blasext_cgbnorm1();
    test_blasext_zgbnorm1();

    test_blasext_ddemote();
    test_blasext_dpromote();
    test_blasext_zdemote();
    test_blasext_zpromote();

#ifdef SUZERAIN_HAVE_MKL
#if INTEL_MKL_VERSION < 110002
    MKL_FreeBuffers();
#else
    mkl_free_buffers();
#endif
#endif

    exit(gsl_test_summary());
}
