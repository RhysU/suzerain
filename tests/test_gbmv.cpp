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

#include <suzerain/blas_et_al/gbmv.h>

#include <boost/test/parameterized_test.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al/blas_et_al.h>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

static gsl_rng * rng;  // Managed within main()

using boost::algorithm::replace_first_copy;
using boost::unit_test::framework::master_test_suite;
using boost::unit_test::make_test_case;
using std::numeric_limits;
using std::size_t;

// I said button yer lip Mugsy.
#pragma warning(disable:1418 1572)

// Test suzerain_gbmv_?  against suzerain_blas_*gbmv_external
// Test suzerain_gbmv_?? against suzerain_blasext_?gbmzv_external
// Idea behind testing is that matching the BLAS is goodness

// For unary function-based test case registration
struct gbmv_tc_type
{
    char trans;
    int m, n, kl, ku;
    double alpha;
    int lda, incx;
    double beta;
    int incy;
};

// For unary function-based test case registration
struct gbmzv_tc_type
{
    char trans;
    int m, n, kl, ku;
    double alpha[2];
    int lda, incx;
    double beta[2];
    int incy;

    gbmzv_tc_type(const gbmv_tc_type& o)
        : trans(o.trans), m(o.m), n(o.n), kl(o.kl), ku(o.ku),
          lda(o.lda), incx(o.incx), incy(o.incy)
    {
        alpha[0] = o.alpha; alpha[1] = 0;
        beta[0]  = o.beta;  beta[1]  = 0;
    }
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbmv_tc_type& t)
{
    os << "{trans="  << t.trans
       << ", m="     << t.m
       << ", n="     << t.n
       << ", kl="    << t.kl
       << ", ku="    << t.ku
       << ", alpha=" << t.alpha
       << ", lda="   << t.lda
       << ", incx="  << t.incx
       << ", beta="  << t.beta
       << ", incy="  << t.incy
       << '}';
    return os;
}

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbmzv_tc_type& t)
{
    os << "{trans="  << t.trans
       << ", m="     << t.m
       << ", n="     << t.n
       << ", kl="    << t.kl
       << ", ku="    << t.ku
       << ", alpha=" << std::complex<double>(t.alpha[0], t.alpha[1])
       << ", lda="   << t.lda
       << ", incx="  << t.incx
       << ", beta="  << std::complex<double>(t.beta[0], t.beta[1])
       << ", incy="  << t.incy
       << '}';
    return os;
}

static void test_gbmv_s(const gbmv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.m*t.n*15;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const float alpha = (float) t.alpha;
    const float beta  = (float) t.beta;

    // Compute expected result using external BLAS
    suzerain_blas_sgbmv_external(t.trans, t.m, t.n, t.kl, t.ku,
                                   alpha, a.get(), t.lda, x.get(), t.incx,
                                   beta,                  e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbmv_s(
                t.trans, t.m, t.n, t.kl, t.ku,
                  alpha, a.get(), t.lda, x.get(), t.incx,
                  beta,                  y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbmv_d(const gbmv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.m*t.n*15;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Compute expected result using external BLAS
    suzerain_blas_dgbmv_external(t.trans, t.m, t.n, t.kl, t.ku,
                                 t.alpha, a.get(), t.lda, x.get(), t.incx,
                                 t.beta,                  e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbmv_d(
                t.trans, t.m, t.n, t.kl, t.ku,
                t.alpha, a.get(), t.lda, x.get(), t.incx,
                t.beta,                  y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbmv_scc(const gbmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.m*t.n*300;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = 2 * abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_float alpha(t.alpha[0], t.alpha[1]);
    const complex_float beta( t.beta[0],  t.beta[1]);

    // Compute expected result using external BLAS
    suzerain_blas_cgbmv_s_c_external(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (complex_float *) x.get(), t.incx,
            beta,                  (complex_float *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbmv_scc(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (complex_float *) x.get(), t.incx,
            beta,                  (complex_float *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbmv_dzz(const gbmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.m*t.n*250;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = 2 * abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_double alpha(t.alpha[0], t.alpha[1]);
    const complex_double beta( t.beta[0],  t.beta[1]);

    // Compute expected result using external BLAS
    suzerain_blas_zgbmv_d_z_external(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (complex_double *) x.get(), t.incx,
            beta,                  (complex_double *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbmv_dzz(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (complex_double *) x.get(), t.incx,
            beta,                  (complex_double *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbmv_ssc(const gbmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*2500;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = 2 * abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_float alpha(t.alpha[0], t.alpha[1]);
    const complex_float beta( t.beta[0],  t.beta[1]);

    // Compute expected result using scc implementation with Im(x) = 0
    for (int i = 0; i < lenx/2/abs(t.incx); i++) {
        x[2*i*abs(t.incx)+1] = 0;
    }
    suzerain_blas_cgbmv_s_c(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (complex_float *) x.get(), t.incx,
            beta,                  (complex_float *) e.get(), t.incy);

    // Compute observed result using ssc with a poisoned Im(x) = NaN
    for (int i = 0; i < lenx/2/abs(t.incx); i++) {
        x[2*i*abs(t.incx)+1] = std::numeric_limits<float>::quiet_NaN();
    }
    BOOST_REQUIRE_EQUAL(0, suzerain_blas_cgbmv_s_s(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda,                   x.get(), 2*t.incx,
            beta,                  (complex_float *) y.get(),   t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbmv_ddz(const gbmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*1000;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = 2 * abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_double alpha(t.alpha[0], t.alpha[1]);
    const complex_double beta( t.beta[0],  t.beta[1]);

    // Compute expected result using dzz implementation with Im(x) = 0
    for (int i = 0; i < lenx/2/abs(t.incx); i++) {
        x[2*i*abs(t.incx)+1] = 0;
    }
    suzerain_blas_zgbmv_d_z(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (complex_double *) x.get(), t.incx,
            beta,                  (complex_double *) e.get(), t.incy);

    // Compute observed result using ddz with a poisoned Im(x) = NaN
    for (int i = 0; i < lenx/2/abs(t.incx); i++) {
        x[2*i*abs(t.incx)+1] = std::numeric_limits<double>::quiet_NaN();
    }
    BOOST_REQUIRE_EQUAL(0, suzerain_blas_zgbmv_d_d(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda,                    x.get(), 2*t.incx,
            beta,                  (complex_double *) y.get(),   t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
bool init_unit_test_suite() {
#else
::boost::unit_test::test_suite* init_unit_test_suite( int, char* [] )   {
#endif

    master_test_suite().p_name.value = __FILE__;

    // -------------------------------------------------------
    // Register test cases designed for generalized bandwidths
    // -------------------------------------------------------

    const gbmv_tc_type gbmv_tc[] = {
        // trans,  m,  n, kl, ku, alpha, lda, incx, beta, incy
         {   'N', 17, 19,  3,  2,  -5.0,   7,    1,  0.0,    1} // Regular
        ,{   'N', 17, 19,  3,  2,   5.0,   7,   -1, 11.0,    1}
        ,{   'N', 17, 19,  3,  2,  -5.0,   6,    1,  0.0,    1}
        ,{   'N', 17, 19,  3,  2,   5.0,   6,    1, 11.0,   -1}
        ,{   'N', 17, 19,  3,  2,  -5.0,   7,    1,  0.0,    2}
        ,{   'N', 17, 19,  3,  2,   5.0,   7,   -1, 11.0,    2}
        ,{   'N', 17, 19,  3,  2,  -5.0,   6,    1,  0.0,    2}
        ,{   'N', 17, 19,  3,  2,   5.0,   6,    1, 11.0,   -2}
        ,{   'N', 17, 19,  3,  2,  -5.0,   7,    3,  0.0,    1}
        ,{   'N', 17, 19,  3,  2,   5.0,   7,   -3, 11.0,    1}
        ,{   'N', 17, 19,  3,  2,  -5.0,   6,    3,  0.0,    1}
        ,{   'N', 17, 19,  3,  2,   5.0,   6,    3, 11.0,   -1}
        ,{   'N', 17, 19,  3,  2,  -5.0,   7,    3,  0.0,    2}
        ,{   'N', 17, 19,  3,  2,   5.0,   7,   -3, 11.0,    2}
        ,{   'N', 17, 19,  3,  2,  -5.0,   6,    3,  0.0,    2}
        ,{   'N', 17, 19,  3,  2,   5.0,   6,    3, 11.0,   -2}
        ,{   'T', 17, 19,  3,  2,  -5.0,   7,    1,  0.0,    1}
        ,{   'T', 17, 19,  3,  2,   5.0,   7,   -1, 11.0,    1}
        ,{   'T', 17, 19,  3,  2,  -5.0,   6,    1,  0.0,    1}
        ,{   'T', 17, 19,  3,  2,   5.0,   6,    1, 11.0,   -1}
        ,{   'T', 17, 19,  3,  2,  -5.0,   7,    1,  0.0,    2}
        ,{   'T', 17, 19,  3,  2,   5.0,   7,   -1, 11.0,    2}
        ,{   'T', 17, 19,  3,  2,  -5.0,   6,    1,  0.0,    2}
        ,{   'T', 17, 19,  3,  2,   5.0,   6,    1, 11.0,   -2}
        ,{   'T', 17, 19,  3,  2,  -5.0,   7,    3,  0.0,    1}
        ,{   'T', 17, 19,  3,  2,   5.0,   7,   -3, 11.0,    1}
        ,{   'T', 17, 19,  3,  2,  -5.0,   6,    3,  0.0,    1}
        ,{   'T', 17, 19,  3,  2,   5.0,   6,    3, 11.0,   -1}
        ,{   'T', 17, 19,  3,  2,  -5.0,   7,    3,  0.0,    2}
        ,{   'T', 17, 19,  3,  2,   5.0,   7,   -3, 11.0,    2}
        ,{   'T', 17, 19,  3,  2,  -5.0,   6,    3,  0.0,    2}
        ,{   'T', 17, 19,  3,  2,   5.0,   6,    3, 11.0,   -2}
        ,{   'N', 19, 17,  0,  2,  -5.0,   4,    1,  0.0,    1} // kl == 0
        ,{   'N', 19, 17,  0,  2,   5.0,   4,   -1, 11.0,    1}
        ,{   'N', 19, 17,  0,  2,  -5.0,   3,    1,  0.0,    1}
        ,{   'N', 19, 17,  0,  2,   5.0,   3,    1, 11.0,   -1}
        ,{   'N', 19, 17,  0,  2,  -5.0,   4,    1,  0.0,    2}
        ,{   'N', 19, 17,  0,  2,   5.0,   4,   -1, 11.0,    2}
        ,{   'N', 19, 17,  0,  2,  -5.0,   3,    1,  0.0,    2}
        ,{   'N', 19, 17,  0,  2,   5.0,   3,    1, 11.0,   -2}
        ,{   'N', 19, 17,  0,  2,  -5.0,   4,    3,  0.0,    1}
        ,{   'N', 19, 17,  0,  2,   5.0,   4,   -3, 11.0,    1}
        ,{   'N', 19, 17,  0,  2,  -5.0,   3,    3,  0.0,    1}
        ,{   'N', 19, 17,  0,  2,   5.0,   3,    3, 11.0,   -1}
        ,{   'N', 19, 17,  0,  2,  -5.0,   4,    3,  0.0,    2}
        ,{   'N', 19, 17,  0,  2,   5.0,   4,   -3, 11.0,    2}
        ,{   'N', 19, 17,  0,  2,  -5.0,   3,    3,  0.0,    2}
        ,{   'N', 19, 17,  0,  2,   5.0,   3,    3, 11.0,   -2}
        ,{   'T', 19, 17,  0,  2,  -5.0,   4,    1,  0.0,    1}
        ,{   'T', 19, 17,  0,  2,   5.0,   4,   -1, 11.0,    1}
        ,{   'T', 19, 17,  0,  2,  -5.0,   3,    1,  0.0,    1}
        ,{   'T', 19, 17,  0,  2,   5.0,   3,    1, 11.0,   -1}
        ,{   'T', 19, 17,  0,  2,  -5.0,   4,    1,  0.0,    2}
        ,{   'T', 19, 17,  0,  2,   5.0,   4,   -1, 11.0,    2}
        ,{   'T', 19, 17,  0,  2,  -5.0,   3,    1,  0.0,    2}
        ,{   'T', 19, 17,  0,  2,   5.0,   3,    1, 11.0,   -2}
        ,{   'T', 19, 17,  0,  2,  -5.0,   4,    3,  0.0,    1}
        ,{   'T', 19, 17,  0,  2,   5.0,   4,   -3, 11.0,    1}
        ,{   'T', 19, 17,  0,  2,  -5.0,   3,    3,  0.0,    1}
        ,{   'T', 19, 17,  0,  2,   5.0,   3,    3, 11.0,   -1}
        ,{   'T', 19, 17,  0,  2,  -5.0,   4,    3,  0.0,    2}
        ,{   'T', 19, 17,  0,  2,   5.0,   4,   -3, 11.0,    2}
        ,{   'T', 19, 17,  0,  2,  -5.0,   3,    3,  0.0,    2}
        ,{   'T', 19, 17,  0,  2,   5.0,   4,    3, 11.0,   -2}
        ,{   'N', 19, 17,  3,  0,  -5.0,   5,    1,  0.0,    1} // ku == 0
        ,{   'N', 19, 17,  3,  0,   5.0,   5,   -1, 11.0,    1}
        ,{   'N', 19, 17,  3,  0,  -5.0,   4,    1,  0.0,    1}
        ,{   'N', 19, 17,  3,  0,   5.0,   4,    1, 11.0,   -1}
        ,{   'N', 19, 17,  3,  0,  -5.0,   5,    1,  0.0,    2}
        ,{   'N', 19, 17,  3,  0,   5.0,   5,   -1, 11.0,    2}
        ,{   'N', 19, 17,  3,  0,  -5.0,   4,    1,  0.0,    2}
        ,{   'N', 19, 17,  3,  0,   5.0,   4,    1, 11.0,   -2}
        ,{   'N', 19, 17,  3,  0,  -5.0,   5,    3,  0.0,    1}
        ,{   'N', 19, 17,  3,  0,   5.0,   5,   -3, 11.0,    1}
        ,{   'N', 19, 17,  3,  0,  -5.0,   4,    3,  0.0,    1}
        ,{   'N', 19, 17,  3,  0,   5.0,   4,    3, 11.0,   -1}
        ,{   'N', 19, 17,  3,  0,  -5.0,   5,    3,  0.0,    2}
        ,{   'N', 19, 17,  3,  0,   5.0,   5,   -3, 11.0,    2}
        ,{   'N', 19, 17,  3,  0,  -5.0,   4,    3,  0.0,    2}
        ,{   'N', 19, 17,  3,  0,   5.0,   4,    3, 11.0,   -2}
        ,{   'T', 19, 17,  3,  0,  -5.0,   5,    1,  0.0,    1}
        ,{   'T', 19, 17,  3,  0,   5.0,   5,   -1, 11.0,    1}
        ,{   'T', 19, 17,  3,  0,  -5.0,   4,    1,  0.0,    1}
        ,{   'T', 19, 17,  3,  0,   5.0,   4,    1, 11.0,   -1}
        ,{   'T', 19, 17,  3,  0,  -5.0,   5,    1,  0.0,    2}
        ,{   'T', 19, 17,  3,  0,   5.0,   5,   -1, 11.0,    2}
        ,{   'T', 19, 17,  3,  0,  -5.0,   4,    1,  0.0,    2}
        ,{   'T', 19, 17,  3,  0,   5.0,   4,    1, 11.0,   -2}
        ,{   'T', 19, 17,  3,  0,  -5.0,   5,    3,  0.0,    1}
        ,{   'T', 19, 17,  3,  0,   5.0,   5,   -3, 11.0,    1}
        ,{   'T', 19, 17,  3,  0,  -5.0,   4,    3,  0.0,    1}
        ,{   'T', 19, 17,  3,  0,   5.0,   4,    3, 11.0,   -1}
        ,{   'T', 19, 17,  3,  0,  -5.0,   5,    3,  0.0,    2}
        ,{   'T', 19, 17,  3,  0,   5.0,   5,   -3, 11.0,    2}
        ,{   'T', 19, 17,  3,  0,  -5.0,   4,    3,  0.0,    2}
        ,{   'T', 19, 17,  3,  0,   5.0,   4,    3, 11.0,   -2}
        ,{   'N',  4,  5,  3,  4,  -5.0,   8,    1,  0.0,    1} // Degenerate A
        ,{   'N',  4,  5,  3,  4,   5.0,   8,   -1, 11.0,    1}
        ,{   'N',  4,  5,  3,  4,  -5.0,   9,    1,  0.0,    1}
        ,{   'N',  4,  5,  3,  4,   5.0,   9,    1, 11.0,   -1}
        ,{   'N',  4,  5,  3,  4,  -5.0,   8,    1,  0.0,    2}
        ,{   'N',  4,  5,  3,  4,   5.0,   8,   -1, 11.0,    2}
        ,{   'N',  4,  5,  3,  4,  -5.0,   9,    1,  0.0,    2}
        ,{   'N',  4,  5,  3,  4,   5.0,   9,    1, 11.0,   -2}
        ,{   'N',  4,  5,  3,  4,  -5.0,   8,    3,  0.0,    1}
        ,{   'N',  4,  5,  3,  4,   5.0,   8,   -3, 11.0,    1}
        ,{   'N',  4,  5,  3,  4,  -5.0,   9,    3,  0.0,    1}
        ,{   'N',  4,  5,  3,  4,   5.0,   9,    3, 11.0,   -1}
        ,{   'N',  4,  5,  3,  4,  -5.0,   8,    3,  0.0,    2}
        ,{   'N',  4,  5,  3,  4,   5.0,   8,   -3, 11.0,    2}
        ,{   'N',  4,  5,  3,  4,  -5.0,   9,    3,  0.0,    2}
        ,{   'N',  4,  5,  3,  4,   5.0,   9,    3, 11.0,   -2}
        ,{   'T',  4,  5,  3,  4,  -5.0,   8,    1,  0.0,    1}
        ,{   'T',  4,  5,  3,  4,   5.0,   8,   -1, 11.0,    1}
        ,{   'T',  4,  5,  3,  4,  -5.0,   9,    1,  0.0,    1}
        ,{   'T',  4,  5,  3,  4,   5.0,   9,    1, 11.0,   -1}
        ,{   'T',  4,  5,  3,  4,  -5.0,   8,    1,  0.0,    2}
        ,{   'T',  4,  5,  3,  4,   5.0,   8,   -1, 11.0,    2}
        ,{   'T',  4,  5,  3,  4,  -5.0,   9,    1,  0.0,    2}
        ,{   'T',  4,  5,  3,  4,   5.0,   9,    1, 11.0,   -2}
        ,{   'T',  4,  5,  3,  4,  -5.0,   8,    3,  0.0,    1}
        ,{   'T',  4,  5,  3,  4,   5.0,   8,   -3, 11.0,    1}
        ,{   'T',  4,  5,  3,  4,  -5.0,   9,    3,  0.0,    1}
        ,{   'T',  4,  5,  3,  4,   5.0,   9,    3, 11.0,   -1}
        ,{   'T',  4,  5,  3,  4,  -5.0,   8,    3,  0.0,    2}
        ,{   'T',  4,  5,  3,  4,   5.0,   8,   -3, 11.0,    2}
        ,{   'T',  4,  5,  3,  4,  -5.0,   9,    3,  0.0,    2}
        ,{   'T',  4,  5,  3,  4,   5.0,   9,    3, 11.0,   -2}
        ,{   'N',  5,  4,  3,  2,  -5.0,   8,    1,  0.0,    1} // Degenerate B
        ,{   'N',  5,  4,  3,  2,   5.0,   8,   -1, 11.0,    1}
        ,{   'N',  5,  4,  3,  2,  -5.0,   9,    1,  0.0,    1}
        ,{   'N',  5,  4,  3,  2,   5.0,   9,    1, 11.0,   -1}
        ,{   'N',  5,  4,  3,  2,  -5.0,   8,    1,  0.0,    2}
        ,{   'N',  5,  4,  3,  2,   5.0,   8,   -1, 11.0,    2}
        ,{   'N',  5,  4,  3,  2,  -5.0,   9,    1,  0.0,    2}
        ,{   'N',  5,  4,  3,  2,   5.0,   9,    1, 11.0,   -2}
        ,{   'N',  5,  4,  3,  2,  -5.0,   8,    3,  0.0,    1}
        ,{   'N',  5,  4,  3,  2,   5.0,   8,   -3, 11.0,    1}
        ,{   'N',  5,  4,  3,  2,  -5.0,   9,    3,  0.0,    1}
        ,{   'N',  5,  4,  3,  2,   5.0,   9,    3, 11.0,   -1}
        ,{   'N',  5,  4,  3,  2,  -5.0,   8,    3,  0.0,    2}
        ,{   'N',  5,  4,  3,  2,   5.0,   8,   -3, 11.0,    2}
        ,{   'N',  5,  4,  3,  2,  -5.0,   9,    3,  0.0,    2}
        ,{   'N',  5,  4,  3,  2,   5.0,   9,    3, 11.0,   -2}
        ,{   'T',  5,  4,  3,  2,  -5.0,   8,    1,  0.0,    1}
        ,{   'T',  5,  4,  3,  2,   5.0,   8,   -1, 11.0,    1}
        ,{   'T',  5,  4,  3,  2,  -5.0,   9,    1,  0.0,    1}
        ,{   'T',  5,  4,  3,  2,   5.0,   9,    1, 11.0,   -1}
        ,{   'T',  5,  4,  3,  2,  -5.0,   8,    1,  0.0,    2}
        ,{   'T',  5,  4,  3,  2,   5.0,   8,   -1, 11.0,    2}
        ,{   'T',  5,  4,  3,  2,  -5.0,   9,    1,  0.0,    2}
        ,{   'T',  5,  4,  3,  2,   5.0,   9,    1, 11.0,   -2}
        ,{   'T',  5,  4,  3,  2,  -5.0,   8,    3,  0.0,    1}
        ,{   'T',  5,  4,  3,  2,   5.0,   8,   -3, 11.0,    1}
        ,{   'T',  5,  4,  3,  2,  -5.0,   9,    3,  0.0,    1}
        ,{   'T',  5,  4,  3,  2,   5.0,   9,    3, 11.0,   -1}
        ,{   'T',  5,  4,  3,  2,  -5.0,   8,    3,  0.0,    2}
        ,{   'T',  5,  4,  3,  2,   5.0,   8,   -3, 11.0,    2}
        ,{   'T',  5,  4,  3,  2,  -5.0,   9,    3,  0.0,    2}
        ,{   'T',  5,  4,  3,  2,   5.0,   9,    3, 11.0,   -2}
        ,{   'T', 19, 17,  3,  4,   0.0,   9,    3,  1.0,    1} // Quick
    };
    const size_t gcases = sizeof(gbmv_tc)/sizeof(gbmv_tc[0]);

    // Register test_gbmv_s cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbmv_s) << ' ' << gbmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbmv_s, name.str(), gbmv_tc + i, gbmv_tc + i + 1));
    }

    // Register test_gbmv_d cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbmv_d) << ' ' << gbmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbmv_d, name.str(), gbmv_tc + i, gbmv_tc + i + 1));
    }

    // Register test_gbmv_scc and test_gbmv_ssc cases
    for (size_t i = 0; i < gcases; ++i) {

        gbmzv_tc_type c(gbmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_scc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_scc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_scc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }
    }

    // Register test_gbmv_dzz and test_gbmv_ddz cases
    for (size_t i = 0; i < gcases; ++i) {

        gbmzv_tc_type c(gbmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_dzz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_dzz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_dzz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const gbmv_tc_type fixed_tc[] = {
        // kl, ku to be set while lda is a delta
        // trans,  m,  n, kl, ku, alpha, lda, incx, beta, incy
         {   'N', 17, 19,  0,  0,  -5.0,   1,    1,  0.0,    1} // Regular A
        ,{   'N', 17, 19,  0,  0,   5.0,   1,   -1, 11.0,    1}
        ,{   'N', 17, 19,  0,  0,  -5.0,   0,    1,  0.0,    1}
        ,{   'N', 17, 19,  0,  0,   5.0,   0,    1, 11.0,   -1}
        ,{   'N', 17, 19,  0,  0,  -5.0,   1,    1,  0.0,    2}
        ,{   'N', 17, 19,  0,  0,   5.0,   1,   -1, 11.0,    2}
        ,{   'N', 17, 19,  0,  0,  -5.0,   0,    1,  0.0,    2}
        ,{   'N', 17, 19,  0,  0,   5.0,   0,    1, 11.0,   -2}
        ,{   'N', 17, 19,  0,  0,  -5.0,   1,    3,  0.0,    1}
        ,{   'N', 17, 19,  0,  0,   5.0,   1,   -3, 11.0,    1}
        ,{   'N', 17, 19,  0,  0,  -5.0,   0,    3,  0.0,    1}
        ,{   'N', 17, 19,  0,  0,   5.0,   0,    3, 11.0,   -1}
        ,{   'N', 17, 19,  0,  0,  -5.0,   1,    3,  0.0,    2}
        ,{   'N', 17, 19,  0,  0,   5.0,   1,   -3, 11.0,    2}
        ,{   'N', 17, 19,  0,  0,  -5.0,   0,    3,  0.0,    2}
        ,{   'N', 17, 19,  0,  0,   5.0,   0,    3, 11.0,   -2}
        ,{   'T', 17, 19,  0,  0,  -5.0,   1,    1,  0.0,    1}
        ,{   'T', 17, 19,  0,  0,   5.0,   1,   -1, 11.0,    1}
        ,{   'T', 17, 19,  0,  0,  -5.0,   0,    1,  0.0,    1}
        ,{   'T', 17, 19,  0,  0,   5.0,   0,    1, 11.0,   -1}
        ,{   'T', 17, 19,  0,  0,  -5.0,   1,    1,  0.0,    2}
        ,{   'T', 17, 19,  0,  0,   5.0,   1,   -1, 11.0,    2}
        ,{   'T', 17, 19,  0,  0,  -5.0,   0,    1,  0.0,    2}
        ,{   'T', 17, 19,  0,  0,   5.0,   0,    1, 11.0,   -2}
        ,{   'T', 17, 19,  0,  0,  -5.0,   1,    3,  0.0,    1}
        ,{   'T', 17, 19,  0,  0,   5.0,   1,   -3, 11.0,    1}
        ,{   'T', 17, 19,  0,  0,  -5.0,   0,    3,  0.0,    1}
        ,{   'T', 17, 19,  0,  0,   5.0,   0,    3, 11.0,   -1}
        ,{   'T', 17, 19,  0,  0,  -5.0,   1,    3,  0.0,    2}
        ,{   'T', 17, 19,  0,  0,   5.0,   1,   -3, 11.0,    2}
        ,{   'T', 17, 19,  0,  0,  -5.0,   0,    3,  0.0,    2}
        ,{   'T', 17, 19,  0,  0,   5.0,   0,    3, 11.0,   -2}
        ,{   'N', 19, 17,  0,  0,  -5.0,   1,    1,  0.0,    1} // Regular B
        ,{   'N', 19, 17,  0,  0,   5.0,   1,   -1, 11.0,    1}
        ,{   'N', 19, 17,  0,  0,  -5.0,   0,    1,  0.0,    1}
        ,{   'N', 19, 17,  0,  0,   5.0,   0,    1, 11.0,   -1}
        ,{   'N', 19, 17,  0,  0,  -5.0,   1,    1,  0.0,    2}
        ,{   'N', 19, 17,  0,  0,   5.0,   1,   -1, 11.0,    2}
        ,{   'N', 19, 17,  0,  0,  -5.0,   0,    1,  0.0,    2}
        ,{   'N', 19, 17,  0,  0,   5.0,   0,    1, 11.0,   -2}
        ,{   'N', 19, 17,  0,  0,  -5.0,   1,    3,  0.0,    1}
        ,{   'N', 19, 17,  0,  0,   5.0,   1,   -3, 11.0,    1}
        ,{   'N', 19, 17,  0,  0,  -5.0,   0,    3,  0.0,    1}
        ,{   'N', 19, 17,  0,  0,   5.0,   0,    3, 11.0,   -1}
        ,{   'N', 19, 17,  0,  0,  -5.0,   1,    3,  0.0,    2}
        ,{   'N', 19, 17,  0,  0,   5.0,   1,   -3, 11.0,    2}
        ,{   'N', 19, 17,  0,  0,  -5.0,   0,    3,  0.0,    2}
        ,{   'N', 19, 17,  0,  0,   5.0,   0,    3, 11.0,   -2}
        ,{   'T', 19, 17,  0,  0,  -5.0,   1,    1,  0.0,    1}
        ,{   'T', 19, 17,  0,  0,   5.0,   1,   -1, 11.0,    1}
        ,{   'T', 19, 17,  0,  0,  -5.0,   0,    1,  0.0,    1}
        ,{   'T', 19, 17,  0,  0,   5.0,   0,    1, 11.0,   -1}
        ,{   'T', 19, 17,  0,  0,  -5.0,   1,    1,  0.0,    2}
        ,{   'T', 19, 17,  0,  0,   5.0,   1,   -1, 11.0,    2}
        ,{   'T', 19, 17,  0,  0,  -5.0,   0,    1,  0.0,    2}
        ,{   'T', 19, 17,  0,  0,   5.0,   0,    1, 11.0,   -2}
        ,{   'T', 19, 17,  0,  0,  -5.0,   1,    3,  0.0,    1}
        ,{   'T', 19, 17,  0,  0,   5.0,   1,   -3, 11.0,    1}
        ,{   'T', 19, 17,  0,  0,  -5.0,   0,    3,  0.0,    1}
        ,{   'T', 19, 17,  0,  0,   5.0,   0,    3, 11.0,   -1}
        ,{   'T', 19, 17,  0,  0,  -5.0,   1,    3,  0.0,    2}
        ,{   'T', 19, 17,  0,  0,   5.0,   1,   -3, 11.0,    2}
        ,{   'T', 19, 17,  0,  0,  -5.0,   0,    3,  0.0,    2}
        ,{   'T', 19, 17,  0,  0,   5.0,   0,    3, 11.0,   -2}
        ,{   'T', 19, 17,  0,  0,   0.0,   2,    3,  1.0,    1} // Quick
    };

    // Loop over fixed kl = ku = k bandwidths
    const size_t max_fixed_bandwidth = 20;
    for (size_t k = 0; k < max_fixed_bandwidth; ++k) {

        // Loop over fixed_tc and register both real and complex tests
        for (size_t i = 0; i < sizeof(fixed_tc)/sizeof(fixed_tc[0]); ++i) {

            // Prepare real-valued test details for bandwidth k
            gbmv_tc_type r(fixed_tc[i]);
            r.kl   = r.ku = k;
            r.lda += (r.kl + r.ku + 1);

            { // Register test_gbmv_s case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbmv_s) << ' ' << gbmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbmv_s, name.str(), &r, &r + 1));
            }

            { // Register test_gbmv_d case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbmv_d) << ' ' << gbmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbmv_d, name.str(), &r, &r + 1));
            }

            { // Register test_gbmv_scc and test_gbmv_ssc cases
                gbmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_scc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_scc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_scc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }
            }

            { // Register test_gbmv_dzz cases
                gbmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_dzz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbmv_ddz,
                            replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_dzz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbmv_ddz,
                            replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_dzz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbmv_ddz,
                            replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
                }
            }

        }
    }

#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
    return true;
}
#else
    return 0;
}
#endif

int main( int argc, char* argv[] )
{
    gsl_ieee_env_setup();
    if (!(rng = gsl_rng_alloc(gsl_rng_env_setup()))) return 1;
    const int retval = ::boost::unit_test::unit_test_main(
            &init_unit_test_suite, argc, argv);
    gsl_rng_free(rng);
    return retval;
}
