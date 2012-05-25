#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/gbddmv.h>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/unit_test.hpp>
#include <suzerain/blas_et_al.h>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

using boost::unit_test::framework::master_test_suite;
using boost::unit_test::make_test_case;
using std::numeric_limits;
using std::size_t;

// Test suzerain_gbddmv_?  against suzerain_blasext_*gbddmv_external
// Test suzerain_gbddmv_?? against suzerain_blasext_?gbddmzv_external
// Idea behind testing is that matching the BLAS is goodness

// For unary function-based test case registration
struct gbddmv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0;
    int ldd0;
    double alpha1;
    int ldd1, lda, incx;
    double beta;
    int incy;
};

// For unary function-based test case registration
struct gbddmzv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0[2];
    int ldd0;
    double alpha1[2];
    int ldd1, lda, incx;
    double beta[2];
    int incy;

    gbddmzv_tc_type(const gbddmv_tc_type& o)
        : trans(o.trans), n(o.n), kl(o.kl), ku(o.ku),
          ldd0(o.ldd0), ldd1(o.ldd1), lda(o.lda), incx(o.incx), incy(o.incy)
    {
        alpha0[0] = o.alpha0; alpha0[1] = 0;
        alpha1[0] = o.alpha1; alpha1[1] = 0;
        beta[0]   = o.beta;   beta[1]   = 0;
    }
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbddmv_tc_type& t)
{
    os << "{trans="   << t.trans
       << ", n="      << t.n
       << ", kl="     << t.kl
       << ", ku="     << t.ku
       << ", alpha0=" << t.alpha0
       << ", ldd0="   << t.ldd0
       << ", alpha1=" << t.alpha1
       << ", ldd1="   << t.ldd1
       << ", lda="    << t.lda
       << ", incx="   << t.incx
       << ", beta="   << t.beta
       << ", incy="   << t.incy
       << '}';
    return os;
}

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbddmzv_tc_type& t)
{
    os << "{trans="   << t.trans
       << ", n="      << t.n
       << ", kl="     << t.kl
       << ", ku="     << t.ku
       << ", alpha0=" << std::complex<double>(t.alpha0[0], t.alpha0[1])
       << ", ldd0="   << t.ldd0
       << ", alpha1=" << std::complex<double>(t.alpha1[0], t.alpha1[1])
       << ", ldd1="   << t.ldd1
       << ", lda="    << t.lda
       << ", incx="   << t.incx
       << ", beta="   << std::complex<double>(t.beta[0], t.beta[1])
       << ", incy="   << t.incy
       << '}';
    return os;
}

static void test_gbddmv_s(const gbddmv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*15;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = abs(t.incx) * t.n;
    const int leny  = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<float> d0(new float[lend0]);
    boost::scoped_array<float> d1(new float[lend1]);
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const float alpha0 = (float) t.alpha0;
    const float alpha1 = (float) t.alpha1;
    const float beta   = (float) t.beta;

    // Compute expected result using external BLAS
    suzerain_blasext_sgbddmv_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0, alpha1, d1.get(), t.ldd1,
            a.get(), t.lda, x.get(), t.incx,
            beta,           e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddmv_s(
                  t.trans, t.n, t.kl, t.ku,
                  alpha0, d0.get(), t.ldd0, alpha1, d1.get(), t.ldd1,
                  a.get(), t.lda, x.get(), t.incx,
                  beta,           y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbddmv_d(const gbddmv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*15;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = abs(t.incx) * t.n;
    const int leny  = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> d0(new double[lend0]);
    boost::scoped_array<double> d1(new double[lend1]);
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Compute expected result using external BLAS
    suzerain_blasext_dgbddmv_external(
            t.trans, t.n, t.kl, t.ku,
            t.alpha0, d0.get(), t.ldd0, t.alpha1, d1.get(), t.ldd1,
            a.get(), t.lda, x.get(), t.incx,
            t.beta,         e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddmv_d(
                t.trans, t.n, t.kl, t.ku,
                t.alpha0, d0.get(), t.ldd0, t.alpha1, d1.get(), t.ldd1,
                a.get(), t.lda, x.get(), t.incx,
                t.beta,         y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbddmv_sc(const gbddmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*2500;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<float> d0(new float[lend0]);
    boost::scoped_array<float> d1(new float[lend1]);
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const complex_float alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_float alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_float beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_cgbddmv_s_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0, alpha1, d1.get(), t.ldd1,
            a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,           (      complex_float *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddmv_sc(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0, alpha1, d1.get(), t.ldd1,
            a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,           (      complex_float *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbddmv_dz(const gbddmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*5000;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lena  = t.lda * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> d0(new double[lend0]);
    boost::scoped_array<double> d1(new double[lend1]);
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const complex_double alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_double alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_double beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_zgbddmv_d_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0, alpha1, d1.get(), t.ldd1,
            a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,           (      complex_double *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddmv_dz(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0, alpha1, d1.get(), t.ldd1,
            a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,           (      complex_double *) y.get(), t.incy));

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

    const gbddmv_tc_type gbddmv_tc[] = {
        // trans,  n, kl, ku, alpha0, ld0, alpha1, ld1, lda, incx, beta, incy
         {   'N', 19,  3,  2,   -5.0,   1,   -1.0,   1,   7,    1,  0.0,    1} // Regular
        ,{   'N', 19,  3,  2,    5.0,   1,    1.0,   2,   7,   -1, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -2.0,   3,   6,    1,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   2,    2.0,   1,   6,    1, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   3,   -3.0,   2,   7,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   3,    3.0,   3,   7,   -1, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   1,   -4.0,   1,   6,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   1,    4.0,   2,   6,    1, 11.0,   -2}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -5.0,   3,   7,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   2,    5.0,   1,   7,   -3, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   3,   -6.0,   2,   6,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   3,    6.0,   3,   6,    3, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   1,   -7.0,   1,   7,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   1,    7.0,   2,   7,   -3, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -8.0,   3,   6,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   2,    8.0,   1,   6,    3, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -9.0,   2,   7,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   3,    9.0,   3,   7,   -1, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -1.0,   1,   6,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   1,    1.0,   2,   6,    1, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   2,   -2.0,   3,   7,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   2,    2.0,   1,   7,   -1, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -3.0,   2,   6,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   3,    3.0,   3,   6,    1, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -4.0,   1,   7,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   1,    4.0,   2,   7,   -3, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   2,   -5.0,   3,   6,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   2,    5.0,   1,   6,    3, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -6.0,   2,   7,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   3,    6.0,   3,   7,   -3, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -7.0,   1,   6,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   1,    7.0,   2,   6,    3, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   2,   -8.0,   3,   4,    1,  0.0,    1} // kl == 0
        ,{   'N', 17,  0,  2,    5.0,   2,    8.0,   1,   4,   -1, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -9.0,   2,   3,    1,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   3,    9.0,   3,   3,    1, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -1.0,   1,   4,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   1,    1.0,   2,   4,   -1, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   2,   -2.0,   3,   3,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   2,    2.0,   1,   3,    1, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -3.0,   2,   4,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   3,    3.0,   3,   4,   -3, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -4.0,   1,   3,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   1,    4.0,   2,   3,    3, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   2,   -5.0,   3,   4,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   2,    5.0,   1,   4,   -3, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -6.0,   2,   3,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   3,    6.0,   3,   3,    3, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   1,   -7.0,   1,   4,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   1,    7.0,   2,   4,   -1, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -8.0,   3,   3,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    8.0,   1,   3,    1, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -9.0,   2,   4,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   3,    9.0,   3,   4,   -1, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   1,   -1.0,   1,   3,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   1,    1.0,   2,   3,    1, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -2.0,   3,   4,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    2.0,   1,   4,   -3, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -3.0,   2,   3,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   3,    3.0,   3,   3,    3, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   1,   -4.0,   1,   4,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   1,    4.0,   2,   4,   -3, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -5.0,   3,   3,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   2,    5.0,   1,   4,    3, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   3,   -6.0,   2,   5,    1,  0.0,    1} // ku == 0
        ,{   'N', 17,  3,  0,    5.0,   3,    6.0,   3,   5,   -1, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -7.0,   1,   4,    1,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   1,    7.0,   2,   4,    1, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -8.0,   3,   5,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   2,    8.0,   1,   5,   -1, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   3,   -9.0,   2,   4,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    9.0,   3,   4,    1, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -1.0,   1,   5,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   1,    1.0,   2,   5,   -3, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -2.0,   3,   4,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   2,    2.0,   1,   4,    3, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   3,   -3.0,   2,   5,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    3.0,   3,   5,   -3, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -4.0,   1,   4,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   1,    4.0,   2,   4,    3, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -5.0,   3,   5,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   2,    5.0,   1,   5,   -1, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   3,   -6.0,   2,   4,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   3,    6.0,   3,   4,    1, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -7.0,   1,   5,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   1,    7.0,   2,   5,   -1, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -8.0,   3,   4,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    8.0,   1,   4,    1, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   3,   -9.0,   2,   5,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   3,    9.0,   3,   5,   -3, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -1.0,   1,   4,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   1,    1.0,   2,   4,    3, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -2.0,   3,   5,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    2.0,   1,   5,   -3, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   3,   -3.0,   2,   4,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   3,    3.0,   3,   4,    3, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -4.0,   1,   8,    1,  0.0,    1} // Degenerate
        ,{   'N',  5,  3,  4,    5.0,   1,    4.0,   2,   8,   -1, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   2,   -5.0,   3,   9,    1,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   2,    5.0,   1,   9,    1, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -6.0,   2,   8,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   3,    6.0,   3,   8,   -1, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -7.0,   1,   9,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   1,    7.0,   2,   9,    1, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   2,   -8.0,   3,   8,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   2,    8.0,   1,   8,   -3, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -9.0,   2,   9,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   3,    9.0,   3,   9,    3, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -1.0,   1,   8,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   1,    1.0,   2,   8,   -3, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   2,   -2.0,   3,   9,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   2,    2.0,   1,   9,    3, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -3.0,   2,   8,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   3,    3.0,   3,   8,   -1, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   1,   -4.0,   1,   9,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   1,    4.0,   2,   9,    1, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -5.0,   3,   8,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   2,    5.0,   1,   8,   -1, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -6.0,   2,   9,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   3,    6.0,   3,   9,    1, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   1,   -7.0,   1,   8,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   1,    7.0,   2,   8,   -3, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -8.0,   3,   9,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   2,    8.0,   1,   9,    3, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -9.0,   2,   8,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   3,    9.0,   3,   8,   -3, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   1,   -1.0,   1,   9,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   1,    1.0,   2,   9,    3, 11.0,   -2}
        ,{   'T', 17,  3,  4,    0.0,   2,    0.0,   3,   9,    3,  1.0,    1} // Quick
    };
    const size_t gcases = sizeof(gbddmv_tc)/sizeof(gbddmv_tc[0]);

    // Register test_gbddmv_s cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbddmv_s) << ' ' << gbddmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbddmv_s, name.str(), gbddmv_tc + i, gbddmv_tc + i + 1));
    }

    // Register test_gbddmv_d cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbddmv_d) << ' ' << gbddmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbddmv_d, name.str(), gbddmv_tc + i, gbddmv_tc + i + 1));
    }

    // Register test_gbddmv_sc cases
    for (size_t i = 0; i < gcases; ++i) {

        gbddmzv_tc_type c(gbddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddmv_sc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddmv_sc, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddmv_sc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddmv_sc, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddmv_sc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddmv_sc, name.str(), &c, &c + 1));
        }
    }

    // Register test_gbddmv_dz cases
    for (size_t i = 0; i < gcases; ++i) {

        gbddmzv_tc_type c(gbddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddmv_dz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddmv_dz, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.beta[0],   c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddmv_dz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddmv_dz, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddmv_dz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddmv_dz, name.str(), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const gbddmv_tc_type fixed_tc[] = {
        // kl, ku to be set while lda is a delta
        // trans,  n, kl, ku, alpha0, ld0, alpha1, ld1, lda, incx, beta, incy
         {   'N', 19,  0,  0,   -5.0,   1,   -5.1,   1,   1,    1,  0.0,    1} // Regular
        ,{   'N', 19,  0,  0,    5.0,   2,    5.1,   1,   1,   -1, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.2,   2,   0,    1,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.2,   2,   0,    1, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   2,   -5.3,   3,   1,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   3,    5.3,   3,   1,   -1, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   1,   -5.4,   1,   0,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   2,    5.4,   1,   0,    1, 11.0,   -2}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.5,   2,   1,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.5,   2,   1,   -3, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   2,   -5.6,   3,   0,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   3,    5.6,   3,   0,    3, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   1,   -5.7,   1,   1,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   2,    5.7,   1,   1,   -3, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.8,   2,   0,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.8,   2,   0,    3, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.9,   3,   1,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.9,   3,   1,   -1, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.1,   1,   0,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.1,   1,   0,    1, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   3,   -5.2,   2,   1,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   1,    5.2,   2,   1,   -1, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.3,   3,   0,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.3,   3,   0,    1, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.4,   1,   1,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.4,   1,   1,   -3, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   3,   -5.5,   2,   0,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   1,    5.5,   2,   0,    3, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.6,   3,   1,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.6,   3,   1,   -3, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.7,   1,   0,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.7,   1,   0,    3, 11.0,   -2}
        ,{   'T', 17,  0,  0,    0.0,   3,    0.0,   2,   2,    3,  1.0,    1} // Quick
    };

    // Loop over fixed kl = ku = k bandwidths
    const size_t max_fixed_bandwidth = 20;
    for (size_t k = 0; k < max_fixed_bandwidth; ++k) {

        // Loop over fixed_tc and register both real and complex tests
        for (size_t i = 0; i < sizeof(fixed_tc)/sizeof(fixed_tc[0]); ++i) {

            // Prepare real-valued test details for bandwidth k
            gbddmv_tc_type r(fixed_tc[i]);
            r.kl   = r.ku = k;
            r.lda += (r.kl + r.ku + 1);

            { // Register test_gbddmv_s case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbddmv_s) << ' ' << gbddmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbddmv_s, name.str(), &r, &r + 1));
            }

            { // Register test_gbddmv_d case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbddmv_d) << ' ' << gbddmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbddmv_d, name.str(), &r, &r + 1));
            }

            { // Register test_gbddmv_sc cases
                gbddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddmv_sc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddmv_sc, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddmv_sc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddmv_sc, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddmv_sc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddmv_sc, name.str(), &c, &c + 1));
                }
            }

            { // Register test_gbddmv_dz cases
                gbddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddmv_dz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddmv_dz, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddmv_dz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddmv_dz, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddmv_dz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddmv_dz, name.str(), &c, &c + 1));
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
    return ::boost::unit_test::unit_test_main( &init_unit_test_suite, argc, argv );
}
