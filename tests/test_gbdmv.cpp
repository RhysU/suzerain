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

#include <suzerain/blas_et_al/gbdmv.h>

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

// That'll teach ya to get any ideers Mugsy.
#pragma warning(disable:1418 1572)

// Test suzerain_gbdmv_?  against suzerain_blasext_*gbdmv_external
// Test suzerain_gbdmv_?? against suzerain_blasext_?gbdmzv_external
// Idea behind testing is that matching the BLAS is goodness

// For unary function-based test case registration
struct gbdmv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha;
    int ldd, lda, incx;
    double beta;
    int incy;
};

// For unary function-based test case registration
struct gbdmzv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha[2];
    int ldd, lda, incx;
    double beta[2];
    int incy;

    gbdmzv_tc_type(const gbdmv_tc_type& o)
        : trans(o.trans), n(o.n), kl(o.kl), ku(o.ku),
          ldd(o.ldd), lda(o.lda), incx(o.incx), incy(o.incy)
    {
        alpha[0] = o.alpha; alpha[1] = 0;
        beta[0]  = o.beta;  beta[1]  = 0;
    }
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbdmv_tc_type& t)
{
    os << "{trans="  << t.trans
       << ", n="     << t.n
       << ", kl="    << t.kl
       << ", ku="    << t.ku
       << ", alpha=" << t.alpha
       << ", ldd="   << t.ldd
       << ", lda="   << t.lda
       << ", incx="  << t.incx
       << ", beta="  << t.beta
       << ", incy="  << t.incy
       << '}';
    return os;
}

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbdmzv_tc_type& t)
{
    os << "{trans="  << t.trans
       << ", n="     << t.n
       << ", kl="    << t.kl
       << ", ku="    << t.ku
       << ", alpha=" << std::complex<double>(t.alpha[0], t.alpha[1])
       << ", ldd="   << t.ldd
       << ", lda="   << t.lda
       << ", incx="  << t.incx
       << ", beta="  << std::complex<double>(t.beta[0], t.beta[1])
       << ", incy="  << t.incy
       << '}';
    return os;
}

static void test_gbdmv_s(const gbdmv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*15;
    const int lend = t.ldd * t.n;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * t.n;
    const int leny = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> d(new float[lend]);
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend; ++i) d[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const float alpha = (float) t.alpha;
    const float beta  = (float) t.beta;

    // Compute expected result using external BLAS
    suzerain_blasext_sgbdmv_external(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd, a.get(), t.lda, x.get(), t.incx,
            beta,                                  e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdmv_s(
                  t.trans, t.n, t.kl, t.ku,
                  alpha, d.get(), t.ldd,
                  a.get(), t.lda, x.get(), t.incx,
                  beta,           y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdmv_d(const gbdmv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*15;
    const int lend = t.ldd * t.n;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * t.n;
    const int leny = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> d(new double[lend]);
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend; ++i) d[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Compute expected result using external BLAS
    suzerain_blasext_dgbdmv_external(
            t.trans, t.n, t.kl, t.ku,
            t.alpha, d.get(), t.ldd, a.get(), t.lda, x.get(), t.incx,
            t.beta,                                  e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdmv_d(
                t.trans, t.n, t.kl, t.ku,
                t.alpha, d.get(), t.ldd,
                a.get(), t.lda, x.get(), t.incx,
                t.beta,         y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdmv_scc(const gbdmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*2500;
    const int lend = t.ldd * t.n;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> d(new float[lend]);
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend; ++i) d[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_float alpha( t.alpha[0], t.alpha[1] );
    const complex_float beta( t.beta[0],  t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_cgbdmv_s_c_external(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda, (complex_float *) x.get(), t.incx,
            beta,           (complex_float *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdmv_scc(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda, (complex_float *) x.get(), t.incx,
            beta,           (complex_float *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdmv_dzz(const gbdmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*10000;
    const int lend = t.ldd * t.n;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> d(new double[lend]);
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend; ++i) d[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_double alpha( t.alpha[0], t.alpha[1] );
    const complex_double beta( t.beta[0],  t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_zgbdmv_d_z_external(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda, (complex_double *) x.get(), t.incx,
            beta,           (complex_double *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdmv_dzz(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda, (complex_double *) x.get(), t.incx,
            beta,           (complex_double *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdmv_ssc(const gbdmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*1000;
    const int lend = t.ldd * t.n;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<float> d(new float[lend]);
    suzerain::scoped_array<float> a(new float[lena]);
    suzerain::scoped_array<float> x(new float[lenx]);
    suzerain::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend; ++i) d[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_float alpha( t.alpha[0], t.alpha[1] );
    const complex_float beta( t.beta[0],  t.beta[1]  );

    // Compute expected result using scc implementation with Im(x) = 0
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = 0;
    }
    suzerain_blasext_cgbdmv_s_c(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda, (complex_float *) x.get(), t.incx,
            beta,           (complex_float *) e.get(), t.incy);

    // Compute observed result using ssc with a poisoned Im(x) = NaN
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = std::numeric_limits<float>::quiet_NaN();
    }
    BOOST_REQUIRE_EQUAL(0, suzerain_blasext_cgbdmv_s_s(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda,                   x.get(), 2*t.incx,
            beta,           (complex_float *) y.get(),   t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdmv_ddz(const gbdmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*1000;
    const int lend = t.ldd * t.n;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    suzerain::scoped_array<double> d(new double[lend]);
    suzerain::scoped_array<double> a(new double[lena]);
    suzerain::scoped_array<double> x(new double[lenx]);
    suzerain::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend; ++i) d[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lena; ++i) a[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < lenx; ++i) x[i] = gsl_rng_uniform_pos(rng);
    for (int i = 0; i < leny; ++i) e[i] = y[i] = gsl_rng_uniform_pos(rng);

    // Get appropriately typed alpha and beta constants
    const complex_double alpha( t.alpha[0], t.alpha[1] );
    const complex_double beta( t.beta[0],  t.beta[1]  );

    // Compute expected result using dzz implementation with Im(x) = 0
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = 0;
    }
    suzerain_blasext_zgbdmv_d_z(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda, (complex_double *) x.get(), t.incx,
            beta,           (complex_double *) e.get(), t.incy);

    // Compute observed result using ddz with a poisoned Im(x) = NaN
    for (int i = 0; i < t.n; i++) {
        x[2*i*abs(t.incx)+1] = std::numeric_limits<double>::quiet_NaN();
    }
    BOOST_REQUIRE_EQUAL(0, suzerain_blasext_zgbdmv_d_d(
            t.trans, t.n, t.kl, t.ku,
            alpha, d.get(), t.ldd,
            a.get(), t.lda,                    x.get(), 2*t.incx,
            beta,           (complex_double *) y.get(),   t.incy));

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

    const gbdmv_tc_type gbdmv_tc[] = {
        // trans,  n, kl, ku, alpha, ldd, lda, incx, beta, incy
         {   'N', 19,  3,  2,  -5.0,   1, 7,    1,  0.0,    1} // Regular
        ,{   'N', 19,  3,  2,   5.0,   2, 7,   -1, 11.0,    1}
        ,{   'N', 19,  3,  2,  -5.0,   3, 6,    1,  0.0,    1}
        ,{   'N', 19,  3,  2,   5.0,   1, 6,    1, 11.0,   -1}
        ,{   'N', 19,  3,  2,  -5.0,   2, 7,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,   5.0,   3, 7,   -1, 11.0,    2}
        ,{   'N', 19,  3,  2,  -5.0,   1, 6,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,   5.0,   2, 6,    1, 11.0,   -2}
        ,{   'N', 19,  3,  2,  -5.0,   3, 7,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,   5.0,   1, 7,   -3, 11.0,    1}
        ,{   'N', 19,  3,  2,  -5.0,   2, 6,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,   5.0,   3, 6,    3, 11.0,   -1}
        ,{   'N', 19,  3,  2,  -5.0,   1, 7,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,   5.0,   2, 7,   -3, 11.0,    2}
        ,{   'N', 19,  3,  2,  -5.0,   3, 6,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,   5.0,   1, 6,    3, 11.0,   -2}
        ,{   'T', 19,  3,  2,  -5.0,   2, 7,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,   5.0,   3, 7,   -1, 11.0,    1}
        ,{   'T', 19,  3,  2,  -5.0,   1, 6,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,   5.0,   2, 6,    1, 11.0,   -1}
        ,{   'T', 19,  3,  2,  -5.0,   3, 7,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,   5.0,   1, 7,   -1, 11.0,    2}
        ,{   'T', 19,  3,  2,  -5.0,   2, 6,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,   5.0,   3, 6,    1, 11.0,   -2}
        ,{   'T', 19,  3,  2,  -5.0,   1, 7,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,   5.0,   2, 7,   -3, 11.0,    1}
        ,{   'T', 19,  3,  2,  -5.0,   3, 6,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,   5.0,   1, 6,    3, 11.0,   -1}
        ,{   'T', 19,  3,  2,  -5.0,   2, 7,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,   5.0,   3, 7,   -3, 11.0,    2}
        ,{   'T', 19,  3,  2,  -5.0,   1, 6,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,   5.0,   2, 6,    3, 11.0,   -2}
        ,{   'N', 17,  0,  2,  -5.0,   3, 4,    1,  0.0,    1} // kl == 0
        ,{   'N', 17,  0,  2,   5.0,   1, 4,   -1, 11.0,    1}
        ,{   'N', 17,  0,  2,  -5.0,   2, 3,    1,  0.0,    1}
        ,{   'N', 17,  0,  2,   5.0,   3, 3,    1, 11.0,   -1}
        ,{   'N', 17,  0,  2,  -5.0,   1, 4,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,   5.0,   2, 4,   -1, 11.0,    2}
        ,{   'N', 17,  0,  2,  -5.0,   3, 3,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,   5.0,   1, 3,    1, 11.0,   -2}
        ,{   'N', 17,  0,  2,  -5.0,   2, 4,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,   5.0,   3, 4,   -3, 11.0,    1}
        ,{   'N', 17,  0,  2,  -5.0,   1, 3,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,   5.0,   2, 3,    3, 11.0,   -1}
        ,{   'N', 17,  0,  2,  -5.0,   3, 4,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,   5.0,   1, 4,   -3, 11.0,    2}
        ,{   'N', 17,  0,  2,  -5.0,   2, 3,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,   5.0,   3, 3,    3, 11.0,   -2}
        ,{   'T', 17,  0,  2,  -5.0,   1, 4,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,   5.0,   2, 4,   -1, 11.0,    1}
        ,{   'T', 17,  0,  2,  -5.0,   3, 3,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,   5.0,   1, 3,    1, 11.0,   -1}
        ,{   'T', 17,  0,  2,  -5.0,   2, 4,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,   5.0,   3, 4,   -1, 11.0,    2}
        ,{   'T', 17,  0,  2,  -5.0,   1, 3,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,   5.0,   2, 3,    1, 11.0,   -2}
        ,{   'T', 17,  0,  2,  -5.0,   3, 4,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,   5.0,   1, 4,   -3, 11.0,    1}
        ,{   'T', 17,  0,  2,  -5.0,   2, 3,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,   5.0,   3, 3,    3, 11.0,   -1}
        ,{   'T', 17,  0,  2,  -5.0,   1, 4,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,   5.0,   2, 4,   -3, 11.0,    2}
        ,{   'T', 17,  0,  2,  -5.0,   3, 3,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,   5.0,   1, 4,    3, 11.0,   -2}
        ,{   'N', 17,  3,  0,  -5.0,   2, 5,    1,  0.0,    1} // ku == 0
        ,{   'N', 17,  3,  0,   5.0,   3, 5,   -1, 11.0,    1}
        ,{   'N', 17,  3,  0,  -5.0,   1, 4,    1,  0.0,    1}
        ,{   'N', 17,  3,  0,   5.0,   2, 4,    1, 11.0,   -1}
        ,{   'N', 17,  3,  0,  -5.0,   3, 5,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,   5.0,   1, 5,   -1, 11.0,    2}
        ,{   'N', 17,  3,  0,  -5.0,   2, 4,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,   5.0,   3, 4,    1, 11.0,   -2}
        ,{   'N', 17,  3,  0,  -5.0,   1, 5,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,   5.0,   2, 5,   -3, 11.0,    1}
        ,{   'N', 17,  3,  0,  -5.0,   3, 4,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,   5.0,   1, 4,    3, 11.0,   -1}
        ,{   'N', 17,  3,  0,  -5.0,   2, 5,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,   5.0,   3, 5,   -3, 11.0,    2}
        ,{   'N', 17,  3,  0,  -5.0,   1, 4,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,   5.0,   2, 4,    3, 11.0,   -2}
        ,{   'T', 17,  3,  0,  -5.0,   3, 5,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,   5.0,   1, 5,   -1, 11.0,    1}
        ,{   'T', 17,  3,  0,  -5.0,   2, 4,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,   5.0,   3, 4,    1, 11.0,   -1}
        ,{   'T', 17,  3,  0,  -5.0,   1, 5,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,   5.0,   2, 5,   -1, 11.0,    2}
        ,{   'T', 17,  3,  0,  -5.0,   3, 4,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,   5.0,   1, 4,    1, 11.0,   -2}
        ,{   'T', 17,  3,  0,  -5.0,   2, 5,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,   5.0,   3, 5,   -3, 11.0,    1}
        ,{   'T', 17,  3,  0,  -5.0,   1, 4,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,   5.0,   2, 4,    3, 11.0,   -1}
        ,{   'T', 17,  3,  0,  -5.0,   3, 5,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,   5.0,   1, 5,   -3, 11.0,    2}
        ,{   'T', 17,  3,  0,  -5.0,   2, 4,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,   5.0,   3, 4,    3, 11.0,   -2}
        ,{   'N',  5,  3,  4,  -5.0,   1, 8,    1,  0.0,    1} // Degenerate
        ,{   'N',  5,  3,  4,   5.0,   2, 8,   -1, 11.0,    1}
        ,{   'N',  5,  3,  4,  -5.0,   3, 9,    1,  0.0,    1}
        ,{   'N',  5,  3,  4,   5.0,   1, 9,    1, 11.0,   -1}
        ,{   'N',  5,  3,  4,  -5.0,   2, 8,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,   5.0,   3, 8,   -1, 11.0,    2}
        ,{   'N',  5,  3,  4,  -5.0,   1, 9,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,   5.0,   2, 9,    1, 11.0,   -2}
        ,{   'N',  5,  3,  4,  -5.0,   3, 8,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,   5.0,   1, 8,   -3, 11.0,    1}
        ,{   'N',  5,  3,  4,  -5.0,   2, 9,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,   5.0,   3, 9,    3, 11.0,   -1}
        ,{   'N',  5,  3,  4,  -5.0,   1, 8,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,   5.0,   2, 8,   -3, 11.0,    2}
        ,{   'N',  5,  3,  4,  -5.0,   3, 9,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,   5.0,   1, 9,    3, 11.0,   -2}
        ,{   'T',  5,  3,  4,  -5.0,   2, 8,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,   5.0,   3, 8,   -1, 11.0,    1}
        ,{   'T',  5,  3,  4,  -5.0,   1, 9,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,   5.0,   2, 9,    1, 11.0,   -1}
        ,{   'T',  5,  3,  4,  -5.0,   3, 8,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,   5.0,   1, 8,   -1, 11.0,    2}
        ,{   'T',  5,  3,  4,  -5.0,   2, 9,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,   5.0,   3, 9,    1, 11.0,   -2}
        ,{   'T',  5,  3,  4,  -5.0,   1, 8,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,   5.0,   2, 8,   -3, 11.0,    1}
        ,{   'T',  5,  3,  4,  -5.0,   3, 9,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,   5.0,   1, 9,    3, 11.0,   -1}
        ,{   'T',  5,  3,  4,  -5.0,   2, 8,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,   5.0,   3, 8,   -3, 11.0,    2}
        ,{   'T',  5,  3,  4,  -5.0,   1, 9,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,   5.0,   2, 9,    3, 11.0,   -2}
        ,{   'T', 17,  3,  4,   0.0,   3, 9,    3,  1.0,    1} // Quick
    };
    const size_t gcases = sizeof(gbdmv_tc)/sizeof(gbdmv_tc[0]);

    // Register test_gbdmv_s cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbdmv_s) << ' ' << gbdmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbdmv_s, name.str(), gbdmv_tc + i, gbdmv_tc + i + 1));
    }

    // Register test_gbdmv_d cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbdmv_d) << ' ' << gbdmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbdmv_d, name.str(), gbdmv_tc + i, gbdmv_tc + i + 1));
    }

    // Register test_gbdmv_scc test_gbdmv_ssc cases
    for (size_t i = 0; i < gcases; ++i) {

        gbdmzv_tc_type c(gbdmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdmv_scc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdmv_scc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdmv_scc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdmv_scc, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdmv_ssc,
                    replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
        }
    }

    // Register test_gbdmv_dzz and test_gbdmv_ddz cases
    for (size_t i = 0; i < gcases; ++i) {

        gbdmzv_tc_type c(gbdmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdmv_dzz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdmv_dzz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdmv_dzz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdmv_dzz, name.str(), &c, &c + 1));
            master_test_suite().add(make_test_case(&test_gbdmv_ddz,
                    replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const gbdmv_tc_type fixed_tc[] = {
        // kl, ku to be set while lda is a delta
        // trans,  n, kl, ku, alpha, ldd, lda, incx, beta, incy
         {   'N', 19,  0,  0,  -5.0,   3, 1,    1,  0.0,    1} // Regular
        ,{   'N', 19,  0,  0,   5.0,   2, 1,   -1, 11.0,    1}
        ,{   'N', 19,  0,  0,  -5.0,   1, 0,    1,  0.0,    1}
        ,{   'N', 19,  0,  0,   5.0,   3, 0,    1, 11.0,   -1}
        ,{   'N', 19,  0,  0,  -5.0,   2, 1,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,   5.0,   1, 1,   -1, 11.0,    2}
        ,{   'N', 19,  0,  0,  -5.0,   3, 0,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,   5.0,   2, 0,    1, 11.0,   -2}
        ,{   'N', 19,  0,  0,  -5.0,   1, 1,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,   5.0,   3, 1,   -3, 11.0,    1}
        ,{   'N', 19,  0,  0,  -5.0,   2, 0,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,   5.0,   1, 0,    3, 11.0,   -1}
        ,{   'N', 19,  0,  0,  -5.0,   3, 1,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,   5.0,   2, 1,   -3, 11.0,    2}
        ,{   'N', 19,  0,  0,  -5.0,   1, 0,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,   5.0,   3, 0,    3, 11.0,   -2}
        ,{   'T', 19,  0,  0,  -5.0,   2, 1,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,   5.0,   1, 1,   -1, 11.0,    1}
        ,{   'T', 19,  0,  0,  -5.0,   3, 0,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,   5.0,   2, 0,    1, 11.0,   -1}
        ,{   'T', 19,  0,  0,  -5.0,   1, 1,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,   5.0,   3, 1,   -1, 11.0,    2}
        ,{   'T', 19,  0,  0,  -5.0,   2, 0,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,   5.0,   1, 0,    1, 11.0,   -2}
        ,{   'T', 19,  0,  0,  -5.0,   3, 1,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,   5.0,   2, 1,   -3, 11.0,    1}
        ,{   'T', 19,  0,  0,  -5.0,   1, 0,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,   5.0,   3, 0,    3, 11.0,   -1}
        ,{   'T', 19,  0,  0,  -5.0,   2, 1,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,   5.0,   1, 1,   -3, 11.0,    2}
        ,{   'T', 19,  0,  0,  -5.0,   3, 0,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,   5.0,   2, 0,    3, 11.0,   -2}
        ,{   'T', 17,  0,  0,   0.0,   1, 2,    3,  1.0,    1} // Quick
    };

    // Loop over fixed kl = ku = k bandwidths
    const size_t max_fixed_bandwidth = 20;
    for (size_t k = 0; k < max_fixed_bandwidth; ++k) {

        // Loop over fixed_tc and register both real and complex tests
        for (size_t i = 0; i < sizeof(fixed_tc)/sizeof(fixed_tc[0]); ++i) {

            // Prepare real-valued test details for bandwidth k
            gbdmv_tc_type r(fixed_tc[i]);
            r.kl   = r.ku = k;
            r.lda += (r.kl + r.ku + 1);

            { // Register test_gbdmv_s case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbdmv_s) << ' ' << gbdmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbdmv_s, name.str(), &r, &r + 1));
            }

            { // Register test_gbdmv_d case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbdmv_d) << ' ' << gbdmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbdmv_d, name.str(), &r, &r + 1));
            }

            { // Register test_gbdmv_scc and test_gbdmv_ssc cases
                gbdmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdmv_scc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdmv_scc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdmv_scc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdmv_scc, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdmv_ssc,
                            replace_first_copy(name.str(), "scc", "ssc"), &c, &c + 1));
                }
            }

            { // Register test_gbdmv_dzz and test_gbdmv_ddz cases
                gbdmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdmv_dzz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdmv_ddz,
                            replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdmv_dzz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdmv_ddz,
                            replace_first_copy(name.str(), "dzz", "ddz"), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdmv_dzz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdmv_dzz, name.str(), &c, &c + 1));
                    master_test_suite().add(make_test_case(&test_gbdmv_ddz,
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
