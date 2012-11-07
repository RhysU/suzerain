//--------------------------------------------------------------------------
//
// Copyright (C) 2011, 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
//--------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/gbdddmv.h>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>
#include <suzerain/blas_et_al.h>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

static gsl_rng * rng;  // Managed within main()

using boost::algorithm::replace_first_copy;
using boost::unit_test::framework::master_test_suite;
using boost::unit_test::make_test_case;
using std::numeric_limits;
using std::size_t;

#pragma warning(disable:1418 1572)

// Test suzerain_gbdddmv_?  against suzerain_blasext_*gbdddmv_external
// Test suzerain_gbdddmv_?? against suzerain_blasext_?gbdddmzv_external
// Idea behind testing is that matching the BLAS is goodness

// For unary function-based test case registration
struct gbdddmv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0;
    int ldd0;
    double alpha1;
    int ldd1;
    double alpha2;
    int ldd2, lda, incx;
    double beta;
    int incy;
};

// For unary function-based test case registration
struct gbdddmzv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0[2];
    int ldd0;
    double alpha1[2];
    int ldd1;
    double alpha2[2];
    int ldd2, lda, incx;
    double beta[2];
    int incy;

    gbdddmzv_tc_type(const gbdddmv_tc_type& o)
        : trans(o.trans), n(o.n), kl(o.kl), ku(o.ku),
          ldd0(o.ldd0), ldd1(o.ldd1), ldd2(o.ldd2),
          lda(o.lda), incx(o.incx), incy(o.incy)
    {
        alpha0[0] = o.alpha0; alpha0[1] = 0;
        alpha1[0] = o.alpha1; alpha1[1] = 0;
        alpha2[0] = o.alpha2; alpha2[1] = 0;
        beta[0]   = o.beta;   beta[1]   = 0;
    }
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbdddmv_tc_type& t)
{
    os << "{trans="   << t.trans
       << ", n="      << t.n
       << ", kl="     << t.kl
       << ", ku="     << t.ku
       << ", alpha0=" << t.alpha0
       << ", ldd0="   << t.ldd0
       << ", alpha1=" << t.alpha1
       << ", ldd1="   << t.ldd1
       << ", alpha2=" << t.alpha2
       << ", ldd2="   << t.ldd2
       << ", lda="    << t.lda
       << ", incx="   << t.incx
       << ", beta="   << t.beta
       << ", incy="   << t.incy
       << '}';
    return os;
}

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbdddmzv_tc_type& t)
{
    os << "{trans="   << t.trans
       << ", n="      << t.n
       << ", kl="     << t.kl
       << ", ku="     << t.ku
       << ", alpha0=" << std::complex<double>(t.alpha0[0], t.alpha0[1])
       << ", ldd0="   << t.ldd0
       << ", alpha1=" << std::complex<double>(t.alpha1[0], t.alpha1[1])
       << ", ldd1="   << t.ldd1
       << ", alpha2=" << std::complex<double>(t.alpha2[0], t.alpha2[1])
       << ", ldd2="   << t.ldd2
       << ", lda="    << t.lda
       << ", incx="   << t.incx
       << ", beta="   << std::complex<double>(t.beta[0], t.beta[1])
       << ", incy="   << t.incy
       << '}';
    return os;
}

static void test_gbdddmv_s(const gbdddmv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*3000;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = abs(t.incx) * t.n;
    const int leny  = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> d0(new float[lend0]);
    suzerain::scoped_array<float> d1(new float[lend1]);
    suzerain::scoped_array<float> d2(new float[lend2]);
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend1; ++i) d1[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend2; ++i) d2[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena;  ++i) a[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx;  ++i) x[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const float alpha0 = (float) t.alpha0;
    const float alpha1 = (float) t.alpha1;
    const float alpha2 = (float) t.alpha2;
    const float beta   = (float) t.beta;

    // Compute expected result using external BLAS
    suzerain_blasext_sgbdddmv_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, x.get(), t.incx,
            beta,           e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_s(
                  t.trans, t.n, t.kl, t.ku,
                  alpha0, d0.get(), t.ldd0,
                  alpha1, d1.get(), t.ldd1,
                  alpha2, d2.get(), t.ldd2,
                  a.get(), t.lda, x.get(), t.incx,
                  beta,           y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_d(const gbdddmv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*3500;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = abs(t.incx) * t.n;
    const int leny  = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> d0(new double[lend0]);
    suzerain::scoped_array<double> d1(new double[lend1]);
    suzerain::scoped_array<double> d2(new double[lend2]);
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend1; ++i) d1[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend2; ++i) d2[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena;  ++i) a[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx;  ++i) x[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = gsl_rng_uniform_pos(rng);

    // Compute expected result using external BLAS
    suzerain_blasext_dgbdddmv_external(
            t.trans, t.n, t.kl, t.ku,
            t.alpha0, d0.get(), t.ldd0,
            t.alpha1, d1.get(), t.ldd1,
            t.alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, x.get(), t.incx,
            t.beta,         e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_d(
                t.trans, t.n, t.kl, t.ku,
                t.alpha0, d0.get(), t.ldd0,
                t.alpha1, d1.get(), t.ldd1,
                t.alpha2, d2.get(), t.ldd2,
                a.get(), t.lda, x.get(), t.incx,
                t.beta,         y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_scc(const gbdddmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*5000;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> d0(new float[lend0]);
    suzerain::scoped_array<float> d1(new float[lend1]);
    suzerain::scoped_array<float> d2(new float[lend2]);
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend1; ++i) d1[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend2; ++i) d2[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena;  ++i) a[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx;  ++i) x[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_float alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_float alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_float alpha2( t.alpha2[0], t.alpha2[1] );
    const complex_float beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_cgbdddmv_s_c_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,           (      complex_float *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_scc(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,           (      complex_float *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_dzz(const gbdddmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*100000;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> d0(new double[lend0]);
    suzerain::scoped_array<double> d1(new double[lend1]);
    suzerain::scoped_array<double> d2(new double[lend2]);
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend1; ++i) d1[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend2; ++i) d2[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena;  ++i) a[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx;  ++i) x[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_double alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_double alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_double alpha2( t.alpha2[0], t.alpha2[1] );
    const complex_double beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_zgbdddmv_d_z_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,           (      complex_double *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_dzz(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,           (      complex_double *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_ssc(const gbdddmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*1000;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> d0(new float[lend0]);
    suzerain::scoped_array<float> d1(new float[lend1]);
    suzerain::scoped_array<float> d2(new float[lend2]);
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend1; ++i) d1[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend2; ++i) d2[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena;  ++i) a[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx;  ++i) x[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_float alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_float alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_float alpha2( t.alpha2[0], t.alpha2[1] );
    const complex_float beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using scc implementation with Im(x) = 0
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = 0;
    }
    suzerain_blasext_cgbdddmv_s_c(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,           (      complex_float *) e.get(), t.incy);

    // Compute observed result using ssc with a poisoned Im(x) = NaN
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = std::numeric_limits<float>::quiet_NaN();
    }
    BOOST_REQUIRE_EQUAL(0, suzerain_blasext_cgbdddmv_s_s(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda,                         x.get(), 2*t.incx,
            beta,           (      complex_float *) y.get(),   t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_ddz(const gbdddmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*1000;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> d0(new double[lend0]);
    suzerain::scoped_array<double> d1(new double[lend1]);
    suzerain::scoped_array<double> d2(new double[lend2]);
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend1; ++i) d1[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lend2; ++i) d2[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena;  ++i) a[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx;  ++i) x[i]  = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_double alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_double alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_double alpha2( t.alpha2[0], t.alpha2[1] );
    const complex_double beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using dzz implementation with Im(x) = 0
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = 0;
    }
    suzerain_blasext_zgbdddmv_d_z(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,           (      complex_double *) e.get(), t.incy);

    // Compute observed result using ddz with a poisoned Im(x) = NaN
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = std::numeric_limits<double>::quiet_NaN();
    }
    BOOST_REQUIRE_EQUAL(0, suzerain_blasext_zgbdddmv_d_d(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            a.get(), t.lda,                          x.get(), 2*t.incx,
            beta,           (      complex_double *) y.get(),   t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
bool init_unit_test_suite() {
#else
::boost::unit_test::test_suite* init_unit_test_suite( int, char* [] ) {
#endif

    master_test_suite().p_name.value = __FILE__;

    // -------------------------------------------------------
    // Register test cases designed for generalized bandwidths
    // -------------------------------------------------------

    const gbdddmv_tc_type gbdddmv_tc[] = {
        // trans,  n, kl, ku, alpha0, ld0, alpha1, ld1, alpha2, ld2, lda, incx, beta, incy
         {   'N', 19,  3,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,   7,    1,  0.0,    1} // Regular
        ,{   'N', 19,  3,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   7,   -1, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,   6,    1,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   6,    1, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,   7,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   7,   -1, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,   6,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   6,    1, 11.0,   -2}
        ,{   'N', 19,  3,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,   7,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   7,   -3, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   1,   -6.0,   1,    0.6,   2,   6,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   1,    6.0,   2,   -0.6,   3,   6,    3, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -7.0,   2,    0.7,   1,   7,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   2,    7.0,   3,   -0.7,   2,   7,   -3, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -8.0,   3,    0.8,   3,   6,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   3,    8.0,   1,   -0.8,   1,   6,    3, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -9.0,   2,    0.9,   2,   7,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   3,    9.0,   3,   -0.9,   3,   7,   -1, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,   6,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   6,    1, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,   7,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   7,   -1, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,   6,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   6,    1, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,   7,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   7,   -3, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,   6,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   6,    3, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -6.0,   1,    0.6,   2,   7,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   1,    6.0,   2,   -0.6,   3,   7,   -3, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   2,   -7.0,   2,    0.7,   1,   6,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   2,    7.0,   3,   -0.7,   2,   6,    3, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   2,   -8.0,   3,    0.8,   3,   4,    1,  0.0,    1} // kl == 0
        ,{   'N', 17,  0,  2,    5.0,   3,    8.0,   1,   -0.8,   1,   4,   -1, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -9.0,   2,    0.9,   2,   3,    1,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   3,    9.0,   3,   -0.9,   3,   3,    1, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,   4,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   4,   -1, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,   3,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   3,    1, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,   4,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   4,   -3, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,   3,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   3,    3, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,   4,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   4,   -3, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -6.0,   1,    0.6,   2,   3,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   1,    6.0,   2,   -0.6,   3,   3,    3, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -7.0,   2,    0.7,   1,   4,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    7.0,   3,   -0.7,   2,   4,   -1, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -8.0,   3,    0.8,   3,   3,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   3,    8.0,   1,   -0.8,   1,   3,    1, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -9.0,   2,    0.9,   2,   4,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   3,    9.0,   3,   -0.9,   3,   4,   -1, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,   3,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   3,    1, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,   4,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   4,   -3, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,   3,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   3,    3, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,   4,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   4,   -3, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,   3,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   4,    3, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -6.0,   1,    0.6,   2,   5,    1,  0.0,    1} // ku == 0
        ,{   'N', 17,  3,  0,    5.0,   1,    6.0,   2,   -0.6,   3,   5,   -1, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -7.0,   2,    0.7,   1,   4,    1,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   2,    7.0,   3,   -0.7,   2,   4,    1, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -8.0,   3,    0.8,   3,   5,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    8.0,   1,   -0.8,   1,   5,   -1, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   3,   -9.0,   2,    0.9,   2,   4,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    9.0,   3,   -0.9,   3,   4,    1, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -1.0,   1,    0.1,   1,   5,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   1,    1.0,   1,   -0.1,   2,   5,   -3, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -2.0,   2,    0.2,   3,   4,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   2,    2.0,   2,   -0.2,   1,   4,    3, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -3.0,   3,    0.3,   2,   5,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   2,    3.0,   3,   -0.3,   3,   5,   -3, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   3,   -4.0,   1,    0.4,   1,   4,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    4.0,   2,   -0.4,   2,   4,    3, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   3,   -5.0,   3,    0.5,   3,   5,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   1,    5.0,   1,   -0.5,   1,   5,   -1, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -6.0,   1,    0.6,   2,   4,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   1,    6.0,   2,   -0.6,   3,   4,    1, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -7.0,   2,    0.7,   1,   5,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    7.0,   3,   -0.7,   2,   5,   -1, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -8.0,   3,    0.8,   3,   4,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   3,    8.0,   1,   -0.8,   1,   4,    1, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   3,   -9.0,   2,    0.9,   2,   5,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   3,    9.0,   3,   -0.9,   3,   5,   -3, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -1.0,   1,    0.1,   1,   4,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   1,    1.0,   1,   -0.1,   2,   4,    3, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -2.0,   2,    0.2,   3,   5,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    2.0,   2,   -0.2,   1,   5,   -3, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -3.0,   3,    0.3,   2,   4,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    3.0,   3,   -0.3,   3,   4,    3, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -4.0,   1,    0.4,   1,   8,    1,  0.0,    1} // Degenerate
        ,{   'N',  5,  3,  4,    5.0,   3,    4.0,   2,   -0.4,   2,   8,   -1, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -5.0,   3,    0.5,   3,   9,    1,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   1,    5.0,   1,   -0.5,   1,   9,    1, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -6.0,   1,    0.6,   2,   8,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   1,    6.0,   2,   -0.6,   3,   8,   -1, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   2,   -7.0,   2,    0.7,   1,   9,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   2,    7.0,   3,   -0.7,   2,   9,    1, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   2,   -8.0,   3,    0.8,   3,   8,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   3,    8.0,   1,   -0.8,   1,   8,   -3, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -9.0,   2,    0.9,   2,   9,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   3,    9.0,   3,   -0.9,   3,   9,    3, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -1.0,   1,    0.1,   1,   8,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   1,    1.0,   1,   -0.1,   2,   8,   -3, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -2.0,   2,    0.2,   3,   9,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   2,    2.0,   2,   -0.2,   1,   9,    3, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -3.0,   3,    0.3,   2,   8,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   2,    3.0,   3,   -0.3,   3,   8,   -1, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -4.0,   1,    0.4,   1,   9,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   3,    4.0,   2,   -0.4,   2,   9,    1, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -5.0,   3,    0.5,   3,   8,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   1,    5.0,   1,   -0.5,   1,   8,   -1, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   1,   -6.0,   1,    0.6,   2,   9,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   1,    6.0,   2,   -0.6,   3,   9,    1, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -7.0,   2,    0.7,   1,   8,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   2,    7.0,   3,   -0.7,   2,   8,   -3, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -8.0,   3,    0.8,   3,   9,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   3,    8.0,   1,   -0.8,   1,   9,    3, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -9.0,   2,    0.9,   2,   8,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   3,    9.0,   3,   -0.9,   3,   8,   -3, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   1,   -1.0,   1,    0.1,   1,   9,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   1,    1.0,   2,   -0.1,   2,   9,    3, 11.0,   -2}
        ,{   'T', 17,  3,  4,    0.0,   1,    0.0,   3,    0.0,   3,   9,    3,  1.0,    1} // Quick
    };
    const size_t gcases = sizeof(gbdddmv_tc)/sizeof(gbdddmv_tc[0]);

    // Register test_gbdddmv_s cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbdddmv_s) << ' ' << gbdddmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbdddmv_s, name.str(), gbdddmv_tc + i, gbdddmv_tc + i + 1));
    }

    // Register test_gbdddmv_d cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbdddmv_d) << ' ' << gbdddmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbdddmv_d, name.str(), gbdddmv_tc + i, gbdddmv_tc + i + 1));
    }

    // Register test_gbdddmv_scc and test_gbdddmv_ssc cases
    for (size_t i = 0; i < gcases; ++i) {

        gbdddmzv_tc_type c(gbdddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_scc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdddmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.alpha2[0], c.alpha2[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_scc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdddmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.alpha2[0] += 3.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_scc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdddmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }
    }

    // Register test_gbdddmv_dzz and test_gbdddmv cases
    for (size_t i = 0; i < gcases; ++i) {

        gbdddmzv_tc_type c(gbdddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_dzz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdddmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.alpha2[0], c.alpha2[1]);
            std::swap(c.beta[0],   c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_dzz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdddmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.alpha2[0] += 3.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_dzz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdddmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const gbdddmv_tc_type fixed_tc[] = {
        // kl, ku to be set while lda is a delta
        // trans,  n, kl, ku, alpha0, ld0, alpha1, ld1, alpha2, ld2, lda, incx, beta, incy
         {   'N', 19,  0,  0,   -5.0,   1,   -5.1,   3,    1.5,   1,   1,    1,  0.0,    1} // Regular
        ,{   'N', 19,  0,  0,    5.0,   2,    5.1,   3,   -1.5,   2,   1,   -1, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.2,   2,    2.5,   1,   0,    1,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.2,   2,   -2.5,   2,   0,    1, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   2,   -5.3,   1,    3.5,   1,   1,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   3,    5.3,   1,   -3.5,   2,   1,   -1, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   1,   -5.4,   3,    4.5,   1,   0,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   2,    5.4,   2,   -4.5,   2,   0,    1, 11.0,   -2}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.5,   1,    5.5,   1,   1,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.5,   3,   -5.5,   2,   1,   -3, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   2,   -5.6,   2,    6.5,   1,   0,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   3,    5.6,   1,   -6.5,   2,   0,    3, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   1,   -5.7,   3,    7.5,   1,   1,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   2,    5.7,   3,   -7.5,   2,   1,   -3, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.8,   2,    8.5,   1,   0,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.8,   2,   -8.5,   2,   0,    3, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.9,   1,    9.5,   1,   1,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.9,   1,   -9.5,   2,   1,   -1, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.1,   3,    1.5,   1,   0,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.1,   2,   -1.5,   2,   0,    1, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   3,   -5.2,   1,    2.5,   1,   1,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   1,    5.2,   3,   -2.5,   2,   1,   -1, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.3,   2,    3.5,   1,   0,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.3,   1,   -3.5,   2,   0,    1, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.4,   3,    4.5,   1,   1,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.4,   3,   -4.5,   2,   1,   -3, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   3,   -5.5,   3,    5.5,   1,   0,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   1,    5.5,   2,   -5.5,   2,   0,    3, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.6,   2,    6.5,   1,   1,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.6,   2,   -6.5,   2,   1,   -3, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.7,   1,    7.5,   3,   0,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.7,   1,   -7.5,   3,   0,    3, 11.0,   -2}
        ,{   'T', 17,  0,  0,    0.0,   3,    0.0,   1,    0.0,   3,   2,    3,  1.0,    1} // Quick
    };

    // Loop over fixed kl = ku = k bandwidths
    const size_t max_fixed_bandwidth = 20;
    for (size_t k = 0; k < max_fixed_bandwidth; ++k) {

        // Loop over fixed_tc and register both real and complex tests
        for (size_t i = 0; i < sizeof(fixed_tc)/sizeof(fixed_tc[0]); ++i) {

            // Prepare real-valued test details for bandwidth k
            gbdddmv_tc_type r(fixed_tc[i]);
            r.kl   = r.ku = k;
            r.lda += (r.kl + r.ku + 1);

            { // Register test_gbdddmv_s case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbdddmv_s) << ' ' << gbdddmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbdddmv_s, name.str(), &r, &r + 1));
            }

            { // Register test_gbdddmv_d case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbdddmv_d) << ' ' << gbdddmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbdddmv_d, name.str(), &r, &r + 1));
            }

            { // Register test_gbdddmv_scc and test_gbdddmv_ssc cases
                gbdddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_scc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdddmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.alpha2[0], c.alpha2[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_scc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdddmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.alpha2[0] += 3.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_scc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdddmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }
            }

            { // Register test_gbdddmv_dzz and test_gbdddmv_ddz cases
                gbdddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_dzz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdddmv_ddz,
                            replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.alpha2[0], c.alpha2[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_dzz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdddmv_ddz,
                            replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.alpha2[0] += 3.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_dzz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdddmv_ddz,
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
