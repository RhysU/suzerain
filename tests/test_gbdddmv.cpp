#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/blas_et_al.h>
#include <suzerain/gbdddmv.h>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

using boost::unit_test::framework::master_test_suite;
using boost::unit_test::make_test_case;
using std::numeric_limits;
using std::size_t;

// Test suzerain_gbdddmv_?  against suzerain_blasext_*gbdddmv_external
// Test suzerain_gbdddmv_?? against suzerain_blasext_?gbdddmzv_external
// Idea behind testing is that matching the BLAS is goodness

// For unary function-based test case registration
struct gbdddmv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0, alpha1, alpha2;
    int lda, incx;
    double beta;
    int incy;
};

// For unary function-based test case registration
struct gbdddmzv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0[2], alpha1[2], alpha2[2];
    int lda, incx;
    double beta[2];
    int incy;

    gbdddmzv_tc_type(const gbdddmv_tc_type& o)
        : trans(o.trans), n(o.n), kl(o.kl), ku(o.ku),
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
       << ", alpha1=" << t.alpha1
       << ", alpha2=" << t.alpha2
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
       << ", alpha1=" << std::complex<double>(t.alpha1[0], t.alpha1[1])
       << ", alpha2=" << std::complex<double>(t.alpha2[0], t.alpha2[1])
       << ", lda="    << t.lda
       << ", incx="   << t.incx
       << ", beta="   << std::complex<double>(t.beta[0], t.beta[1])
       << ", incy="   << t.incy
       << '}';
    return os;
}

static void test_gbdddmv_s(const gbdddmv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*2000;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lend = t.n;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * t.n;
    const int leny = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<float> d0(new float[lend]);
    boost::scoped_array<float> d1(new float[lend]);
    boost::scoped_array<float> d2(new float[lend]);
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lena; ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const float alpha0 = (float) t.alpha0;
    const float alpha1 = (float) t.alpha1;
    const float alpha2 = (float) t.alpha2;
    const float beta   = (float) t.beta;

    // Compute expected result using external BLAS
    suzerain_blasext_sgbdddmv_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), alpha1, d1.get(), alpha2, d2.get(),
            a.get(), t.lda, x.get(), t.incx,
            beta,           e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_s(
                  t.trans, t.n, t.kl, t.ku,
                  alpha0, d0.get(), alpha1, d1.get(), alpha2, d2.get(),
                  a.get(), t.lda, x.get(), t.incx,
                  beta,           y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_d(const gbdddmv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*1000;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lend = t.n;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * t.n;
    const int leny = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> d0(new double[lend]);
    boost::scoped_array<double> d1(new double[lend]);
    boost::scoped_array<double> d2(new double[lend]);
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lena; ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Compute expected result using external BLAS
    suzerain_blasext_dgbdddmv_external(
            t.trans, t.n, t.kl, t.ku,
            t.alpha0, d0.get(), t.alpha1, d1.get(), t.alpha2, d2.get(),
            a.get(), t.lda, x.get(), t.incx,
            t.beta,         e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_d(
                t.trans, t.n, t.kl, t.ku,
                t.alpha0, d0.get(), t.alpha1, d1.get(), t.alpha2, d2.get(),
                a.get(), t.lda, x.get(), t.incx,
                t.beta,         y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_sc(const gbdddmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*5000;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lend = t.n;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<float> d0(new float[lend]);
    boost::scoped_array<float> d1(new float[lend]);
    boost::scoped_array<float> d2(new float[lend]);
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lena; ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const float alpha0[2] = { t.alpha0[0], t.alpha0[1] };
    const float alpha1[2] = { t.alpha1[0], t.alpha1[1] };
    const float alpha2[2] = { t.alpha2[0], t.alpha2[1] };
    const float beta[2]   = { t.beta[0],   t.beta[1]  };

    // Compute expected result using external BLAS
    suzerain_blasext_sgbdddmzv_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), alpha1, d1.get(), alpha2, d2.get(),
            a.get(), t.lda, (const float(*)[2]) x.get(), t.incx,
            beta,           (      float(*)[2]) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_sc(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), alpha1, d1.get(), alpha2, d2.get(),
            a.get(), t.lda, (const float(*)[2]) x.get(), t.incx,
            beta,           (      float(*)[2]) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbdddmv_dz(const gbdddmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*100000;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lend = t.n;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> d0(new double[lend]);
    boost::scoped_array<double> d1(new double[lend]);
    boost::scoped_array<double> d2(new double[lend]);
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lena; ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Compute expected result using external BLAS
    suzerain_blasext_dgbdddmzv_external(
            t.trans, t.n, t.kl, t.ku,
            t.alpha0, d0.get(), t.alpha1, d1.get(), t.alpha2, d2.get(),
            a.get(), t.lda, (const double(*)[2]) x.get(), t.incx,
            t.beta,         (      double(*)[2]) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbdddmv_dz(
            t.trans, t.n, t.kl, t.ku,
            t.alpha0, d0.get(), t.alpha1, d1.get(), t.alpha2, d2.get(),
            a.get(), t.lda, (const double(*)[2]) x.get(), t.incx,
            t.beta,         (      double(*)[2]) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

boost::unit_test::test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    (void) argc; // Unused
    (void) argv; // Unused
    master_test_suite().p_name.value = __FILE__;

    // -------------------------------------------------------
    // Register test cases designed for generalized bandwidths
    // -------------------------------------------------------

    const gbdddmv_tc_type gbdddmv_tc[] = {
        // trans,  n, kl, ku, alpha0, alpha1, alpha2, lda, incx, beta, incy
         {   'N', 19,  3,  2,   -5.0,   -1.0,    0.1,   7,    1,  0.0,    1} // Regular
        ,{   'N', 19,  3,  2,    5.0,    1.0,   -0.1,   7,   -1, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   -2.0,    0.2,   6,    1,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,    2.0,   -0.2,   6,    1, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   -3.0,    0.3,   7,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,    3.0,   -0.3,   7,   -1, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   -4.0,    0.4,   6,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,    4.0,   -0.4,   6,    1, 11.0,   -2}
        ,{   'N', 19,  3,  2,   -5.0,   -5.0,    0.5,   7,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,    5.0,   -0.5,   7,   -3, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   -6.0,    0.6,   6,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,    6.0,   -0.6,   6,    3, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   -7.0,    0.7,   7,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,    7.0,   -0.7,   7,   -3, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   -8.0,    0.8,   6,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,    8.0,   -0.8,   6,    3, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   -9.0,    0.9,   7,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,    9.0,   -0.9,   7,   -1, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   -1.0,    0.1,   6,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,    1.0,   -0.1,   6,    1, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   -2.0,    0.2,   7,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,    2.0,   -0.2,   7,   -1, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   -3.0,    0.3,   6,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,    3.0,   -0.3,   6,    1, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   -4.0,    0.4,   7,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,    4.0,   -0.4,   7,   -3, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   -5.0,    0.5,   6,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,    5.0,   -0.5,   6,    3, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   -6.0,    0.6,   7,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,    6.0,   -0.6,   7,   -3, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   -7.0,    0.7,   6,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,    7.0,   -0.7,   6,    3, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   -8.0,    0.8,   4,    1,  0.0,    1} // kl == 0
        ,{   'N', 17,  0,  2,    5.0,    8.0,   -0.8,   4,   -1, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   -9.0,    0.9,   3,    1,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,    9.0,   -0.9,   3,    1, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   -1.0,    0.1,   4,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,    1.0,   -0.1,   4,   -1, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   -2.0,    0.2,   3,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,    2.0,   -0.2,   3,    1, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   -3.0,    0.3,   4,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,    3.0,   -0.3,   4,   -3, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   -4.0,    0.4,   3,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,    4.0,   -0.4,   3,    3, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   -5.0,    0.5,   4,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,    5.0,   -0.5,   4,   -3, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   -6.0,    0.6,   3,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,    6.0,   -0.6,   3,    3, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   -7.0,    0.7,   4,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,    7.0,   -0.7,   4,   -1, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   -8.0,    0.8,   3,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,    8.0,   -0.8,   3,    1, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   -9.0,    0.9,   4,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,    9.0,   -0.9,   4,   -1, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   -1.0,    0.1,   3,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,    1.0,   -0.1,   3,    1, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   -2.0,    0.2,   4,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,    2.0,   -0.2,   4,   -3, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   -3.0,    0.3,   3,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,    3.0,   -0.3,   3,    3, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   -4.0,    0.4,   4,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,    4.0,   -0.4,   4,   -3, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   -5.0,    0.5,   3,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,    5.0,   -0.5,   4,    3, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   -6.0,    0.6,   5,    1,  0.0,    1} // ku == 0
        ,{   'N', 17,  3,  0,    5.0,    6.0,   -0.6,   5,   -1, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   -7.0,    0.7,   4,    1,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,    7.0,   -0.7,   4,    1, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   -8.0,    0.8,   5,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,    8.0,   -0.8,   5,   -1, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   -9.0,    0.9,   4,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,    9.0,   -0.9,   4,    1, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   -1.0,    0.1,   5,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,    1.0,   -0.1,   5,   -3, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   -2.0,    0.2,   4,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,    2.0,   -0.2,   4,    3, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   -3.0,    0.3,   5,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,    3.0,   -0.3,   5,   -3, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   -4.0,    0.4,   4,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,    4.0,   -0.4,   4,    3, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   -5.0,    0.5,   5,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,    5.0,   -0.5,   5,   -1, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   -6.0,    0.6,   4,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,    6.0,   -0.6,   4,    1, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   -7.0,    0.7,   5,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,    7.0,   -0.7,   5,   -1, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   -8.0,    0.8,   4,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,    8.0,   -0.8,   4,    1, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   -9.0,    0.9,   5,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,    9.0,   -0.9,   5,   -3, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   -1.0,    0.1,   4,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,    1.0,   -0.1,   4,    3, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   -2.0,    0.2,   5,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,    2.0,   -0.2,   5,   -3, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   -3.0,    0.3,   4,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,    3.0,   -0.3,   4,    3, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   -4.0,    0.4,   8,    1,  0.0,    1} // Degenerate
        ,{   'N',  5,  3,  4,    5.0,    4.0,   -0.4,   8,   -1, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   -5.0,    0.5,   9,    1,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,    5.0,   -0.5,   9,    1, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   -6.0,    0.6,   8,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,    6.0,   -0.6,   8,   -1, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   -7.0,    0.7,   9,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,    7.0,   -0.7,   9,    1, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   -8.0,    0.8,   8,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,    8.0,   -0.8,   8,   -3, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   -9.0,    0.9,   9,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,    9.0,   -0.9,   9,    3, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   -1.0,    0.1,   8,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,    1.0,   -0.1,   8,   -3, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   -2.0,    0.2,   9,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,    2.0,   -0.2,   9,    3, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   -3.0,    0.3,   8,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,    3.0,   -0.3,   8,   -1, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   -4.0,    0.4,   9,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,    4.0,   -0.4,   9,    1, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   -5.0,    0.5,   8,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,    5.0,   -0.5,   8,   -1, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   -6.0,    0.6,   9,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,    6.0,   -0.6,   9,    1, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   -7.0,    0.7,   8,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,    7.0,   -0.7,   8,   -3, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   -8.0,    0.8,   9,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,    8.0,   -0.8,   9,    3, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   -9.0,    0.9,   8,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,    9.0,   -0.9,   8,   -3, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   -1.0,    0.1,   9,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,    1.0,   -0.1,   9,    3, 11.0,   -2}
        ,{   'T', 17,  3,  4,    0.0,    0.0,    0.0,   9,    3,  1.0,    1} // Quick
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

    // Register test_gbdddmv_sc cases
    for (size_t i = 0; i < gcases; ++i) {

        gbdddmzv_tc_type c(gbdddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_sc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_sc, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.alpha2[0], c.alpha2[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_sc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_sc, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.alpha2[0] += 3.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_sc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_sc, name.str(), &c, &c + 1));
        }
    }

    // Register test_gbdddmv_dz cases
    for (size_t i = 0; i < gcases; ++i) {

        gbdddmzv_tc_type c(gbdddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_dz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_dz, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.alpha2[0], c.alpha2[1]);
            std::swap(c.beta[0],   c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_dz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_dz, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.alpha2[0] += 3.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbdddmv_dz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbdddmv_dz, name.str(), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const gbdddmv_tc_type fixed_tc[] = {
        // kl, ku to be set while lda is a delta
        // trans,  n, kl, ku, alpha0, alpha1, alpha2, lda, incx, beta, incy
         {   'N', 19,  0,  0,   -5.0,   -5.1,    1.5,   1,    1,  0.0,    1} // Regular
        ,{   'N', 19,  0,  0,    5.0,    5.1,   -1.5,   1,   -1, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   -5.2,    2.5,   0,    1,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,    5.2,   -2.5,   0,    1, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   -5.3,    3.5,   1,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,    5.3,   -3.5,   1,   -1, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   -5.4,    4.5,   0,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,    5.4,   -4.5,   0,    1, 11.0,   -2}
        ,{   'N', 19,  0,  0,   -5.0,   -5.5,    5.5,   1,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,    5.5,   -5.5,   1,   -3, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   -5.6,    6.5,   0,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,    5.6,   -6.5,   0,    3, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   -5.7,    7.5,   1,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,    5.7,   -7.5,   1,   -3, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   -5.8,    8.5,   0,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,    5.8,   -8.5,   0,    3, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   -5.9,    9.5,   1,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,    5.9,   -9.5,   1,   -1, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   -5.1,    1.5,   0,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,    5.1,   -1.5,   0,    1, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   -5.2,    2.5,   1,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,    5.2,   -2.5,   1,   -1, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   -5.3,    3.5,   0,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,    5.3,   -3.5,   0,    1, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   -5.4,    4.5,   1,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,    5.4,   -4.5,   1,   -3, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   -5.5,    5.5,   0,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,    5.5,   -5.5,   0,    3, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   -5.6,    6.5,   1,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,    5.6,   -6.5,   1,   -3, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   -5.7,    7.5,   0,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,    5.7,   -7.5,   0,    3, 11.0,   -2}
        ,{   'T', 17,  0,  0,    0.0,    0.0,    0.0,   2,    3,  1.0,    1} // Quick
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

            { // Register test_gbdddmv_sc cases
                gbdddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_sc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_sc, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.alpha2[0], c.alpha2[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_sc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_sc, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.alpha2[0] += 3.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_sc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_sc, name.str(), &c, &c + 1));
                }
            }

            { // Register test_gbdddmv_dz cases
                gbdddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_dz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_dz, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.alpha2[0], c.alpha2[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_dz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_dz, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.alpha2[0] += 3.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbdddmv_dz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbdddmv_dz, name.str(), &c, &c + 1));
                }
            }

        }
    }

    return 0;
}
