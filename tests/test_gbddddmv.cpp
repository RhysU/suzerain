#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/gbddddmv.h>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/unit_test.hpp>
#include <suzerain/blas_et_al.h>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

using boost::unit_test::framework::master_test_suite;
using boost::unit_test::make_test_case;
using std::numeric_limits;
using std::size_t;

// Test suzerain_gbddddmv_?  against suzerain_blasext_*gbddddmv_external
// Test suzerain_gbddddmv_?? against suzerain_blasext_?gbddddmzv_external
// Idea behind testing is that matching the BLAS is goodness

// For unary function-based test case registration
struct gbddddmv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0;
    int ldd0;
    double alpha1;
    int ldd1;
    double alpha2;
    int ldd2;
    double alpha3;
    int ldd3, lda, incx;
    double beta;
    int incy;
};

// For unary function-based test case registration
struct gbddddmzv_tc_type
{
    char trans;
    int n, kl, ku;
    double alpha0[2];
    int ldd0;
    double alpha1[2];
    int ldd1;
    double alpha2[2];
    int ldd2;
    double alpha3[2];
    int ldd3, lda, incx;
    double beta[2];
    int incy;

    gbddddmzv_tc_type(const gbddddmv_tc_type& o)
        : trans(o.trans), n(o.n), kl(o.kl), ku(o.ku),
          ldd0(o.ldd0), ldd1(o.ldd1), ldd2(o.ldd2), ldd3(o.ldd3),
          lda(o.lda), incx(o.incx), incy(o.incy)
    {
        alpha0[0] = o.alpha0; alpha0[1] = 0;
        alpha1[0] = o.alpha1; alpha1[1] = 0;
        alpha2[0] = o.alpha2; alpha2[1] = 0;
        alpha3[0] = o.alpha3; alpha3[1] = 0;
        beta[0]   = o.beta;   beta[1]   = 0;
    }
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbddddmv_tc_type& t)
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
       << ", alpha3=" << t.alpha3
       << ", ldd3="   << t.ldd3
       << ", lda="    << t.lda
       << ", incx="   << t.incx
       << ", beta="   << t.beta
       << ", incy="   << t.incy
       << '}';
    return os;
}

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const gbddddmzv_tc_type& t)
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
       << ", alpha3=" << std::complex<double>(t.alpha3[0], t.alpha3[1])
       << ", ldd3="   << t.ldd3
       << ", lda="    << t.lda
       << ", incx="   << t.incx
       << ", beta="   << std::complex<double>(t.beta[0], t.beta[1])
       << ", incy="   << t.incy
       << '}';
    return os;
}

static void test_gbddddmv_s(const gbddddmv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*3000;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lend3 = t.ldd3 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = abs(t.incx) * t.n;
    const int leny  = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<float> d0(new float[lend0]);
    boost::scoped_array<float> d1(new float[lend1]);
    boost::scoped_array<float> d2(new float[lend2]);
    boost::scoped_array<float> d3(new float[lend3]);
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend2; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lend3; ++i) d3[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const float alpha0 = (float) t.alpha0;
    const float alpha1 = (float) t.alpha1;
    const float alpha2 = (float) t.alpha2;
    const float alpha3 = (float) t.alpha3;
    const float beta   = (float) t.beta;

    // Compute expected result using external BLAS
    suzerain_blasext_sgbddddmv_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            alpha3, d3.get(), t.ldd3,
            a.get(), t.lda, x.get(), t.incx,
            beta,           e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddddmv_s(
                  t.trans, t.n, t.kl, t.ku,
                  alpha0, d0.get(), t.ldd0,
                  alpha1, d1.get(), t.ldd1,
                  alpha2, d2.get(), t.ldd2,
                  alpha3, d3.get(), t.ldd3,
                  a.get(), t.lda, x.get(), t.incx,
                  beta,           y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbddddmv_d(const gbddddmv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*4000;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lend3 = t.ldd3 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = abs(t.incx) * t.n;
    const int leny  = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> d0(new double[lend0]);
    boost::scoped_array<double> d1(new double[lend1]);
    boost::scoped_array<double> d2(new double[lend2]);
    boost::scoped_array<double> d3(new double[lend3]);
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend2; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lend3; ++i) d3[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Compute expected result using external BLAS
    suzerain_blasext_dgbddddmv_external(
            t.trans, t.n, t.kl, t.ku,
            t.alpha0, d0.get(), t.ldd0,
            t.alpha1, d1.get(), t.ldd1,
            t.alpha2, d2.get(), t.ldd2,
            t.alpha3, d3.get(), t.ldd3,
            a.get(), t.lda, x.get(), t.incx,
            t.beta,         e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddddmv_d(
                t.trans, t.n, t.kl, t.ku,
                t.alpha0, d0.get(), t.ldd0,
                t.alpha1, d1.get(), t.ldd1,
                t.alpha2, d2.get(), t.ldd2,
                t.alpha3, d3.get(), t.ldd3,
                a.get(), t.lda, x.get(), t.incx,
                t.beta,         y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbddddmv_scc(const gbddddmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*6000;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lend3 = t.ldd3 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<float> d0(new float[lend0]);
    boost::scoped_array<float> d1(new float[lend1]);
    boost::scoped_array<float> d2(new float[lend2]);
    boost::scoped_array<float> d3(new float[lend3]);
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend2; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lend3; ++i) d3[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const complex_float alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_float alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_float alpha2( t.alpha2[0], t.alpha2[1] );
    const complex_float alpha3( t.alpha3[0], t.alpha3[1] );
    const complex_float beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_cgbddddmv_s_c_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            alpha3, d3.get(), t.ldd3,
            a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,           (      complex_float *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddddmv_scc(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            alpha3, d3.get(), t.ldd3,
            a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,           (      complex_float *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbddddmv_dzz(const gbddddmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*100000;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lend0 = t.ldd0 * t.n;
    const int lend1 = t.ldd1 * t.n;
    const int lend2 = t.ldd2 * t.n;
    const int lend3 = t.ldd3 * t.n;
    const int lena  = t.lda  * t.n;
    const int lenx  = 2 * abs(t.incx) * t.n;
    const int leny  = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> d0(new double[lend0]);
    boost::scoped_array<double> d1(new double[lend1]);
    boost::scoped_array<double> d2(new double[lend2]);
    boost::scoped_array<double> d3(new double[lend3]);
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lend0; ++i) d0[i] = random() * inv_rand_max;
    for (int i = 0; i < lend1; ++i) d1[i] = random() * inv_rand_max;
    for (int i = 0; i < lend2; ++i) d2[i] = random() * inv_rand_max;
    for (int i = 0; i < lend3; ++i) d3[i] = random() * inv_rand_max;
    for (int i = 0; i < lena;  ++i) a[i]  = random() * inv_rand_max;
    for (int i = 0; i < lenx;  ++i) x[i]  = random() * inv_rand_max;
    for (int i = 0; i < leny;  ++i) e[i]  = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const complex_double alpha0( t.alpha0[0], t.alpha0[1] );
    const complex_double alpha1( t.alpha1[0], t.alpha1[1] );
    const complex_double alpha2( t.alpha2[0], t.alpha2[1] );
    const complex_double alpha3( t.alpha3[0], t.alpha3[1] );
    const complex_double beta  ( t.beta[0],   t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blasext_zgbddddmv_d_z_external(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            alpha3, d3.get(), t.ldd3,
            a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,           (      complex_double *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbddddmv_dzz(
            t.trans, t.n, t.kl, t.ku,
            alpha0, d0.get(), t.ldd0,
            alpha1, d1.get(), t.ldd1,
            alpha2, d2.get(), t.ldd2,
            alpha3, d3.get(), t.ldd3,
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

    const gbddddmv_tc_type gbddddmv_tc[] = {
        // trans,  n, kl, ku, alpha0, ld0, alpha1, ld1, alpha2, ld2, alpha3, ld3, lda, incx, beta, incy
         {   'N', 19,  3,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   4,   7,    1,  0.0,    1} // Regular
        ,{   'N', 19,  3,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   -7.1,   3,   7,   -1, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,    7.2,   2,   6,    1,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   -7.2,   1,   6,    1, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,    7.3,   4,   7,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   -7.3,   1,   7,   -1, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,    7.4,   2,   6,    1,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   -7.4,   3,   6,    1, 11.0,   -2}
        ,{   'N', 19,  3,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,    7.5,   3,   7,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   -7.5,   1,   7,   -3, 11.0,    1}
        ,{   'N', 19,  3,  2,   -5.0,   1,   -6.0,   1,    0.6,   2,    7.6,   2,   6,    3,  0.0,    1}
        ,{   'N', 19,  3,  2,    5.0,   1,    6.0,   2,   -0.6,   3,   -7.6,   4,   6,    3, 11.0,   -1}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -7.0,   2,    0.7,   1,    7.7,   4,   7,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   2,    7.0,   3,   -0.7,   2,   -7.7,   3,   7,   -3, 11.0,    2}
        ,{   'N', 19,  3,  2,   -5.0,   2,   -8.0,   3,    0.8,   3,    7.8,   2,   6,    3,  0.0,    2}
        ,{   'N', 19,  3,  2,    5.0,   3,    8.0,   1,   -0.8,   1,   -7.8,   1,   6,    3, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -9.0,   2,    0.9,   2,    7.9,   4,   7,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   3,    9.0,   3,   -0.9,   3,   -7.9,   1,   7,   -1, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   2,   6,    1,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   -7.1,   3,   6,    1, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,    7.2,   3,   7,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   -7.2,   1,   7,   -1, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,    7.3,   2,   6,    1,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   -7.3,   4,   6,    1, 11.0,   -2}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,    7.4,   4,   7,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   -7.4,   3,   7,   -3, 11.0,    1}
        ,{   'T', 19,  3,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,    7.5,   2,   6,    3,  0.0,    1}
        ,{   'T', 19,  3,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   -7.5,   1,   6,    3, 11.0,   -1}
        ,{   'T', 19,  3,  2,   -5.0,   1,   -6.0,   1,    0.6,   2,    7.6,   4,   7,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   1,    6.0,   2,   -0.6,   3,   -7.6,   1,   7,   -3, 11.0,    2}
        ,{   'T', 19,  3,  2,   -5.0,   2,   -7.0,   2,    0.7,   1,    7.7,   2,   6,    3,  0.0,    2}
        ,{   'T', 19,  3,  2,    5.0,   2,    7.0,   3,   -0.7,   2,   -7.7,   3,   6,    3, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   2,   -8.0,   3,    0.8,   3,    7.8,   3,   4,    1,  0.0,    1} // kl == 0
        ,{   'N', 17,  0,  2,    5.0,   3,    8.0,   1,   -0.8,   1,   -7.8,   1,   4,   -1, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -9.0,   2,    0.9,   2,    7.9,   2,   3,    1,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   3,    9.0,   3,   -0.9,   3,   -7.9,   4,   3,    1, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   4,   4,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   -7.1,   3,   4,   -1, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,    7.2,   2,   3,    1,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   -7.2,   1,   3,    1, 11.0,   -2}
        ,{   'N', 17,  0,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,    7.3,   4,   4,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   -7.3,   1,   4,   -3, 11.0,    1}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,    7.4,   2,   3,    3,  0.0,    1}
        ,{   'N', 17,  0,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   -7.4,   3,   3,    3, 11.0,   -1}
        ,{   'N', 17,  0,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,    7.5,   3,   4,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   -7.5,   1,   4,   -3, 11.0,    2}
        ,{   'N', 17,  0,  2,   -5.0,   1,   -6.0,   1,    0.6,   2,    7.6,   2,   3,    3,  0.0,    2}
        ,{   'N', 17,  0,  2,    5.0,   1,    6.0,   2,   -0.6,   3,   -7.6,   4,   3,    3, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -7.0,   2,    0.7,   1,    7.7,   4,   4,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    7.0,   3,   -0.7,   2,   -7.7,   3,   4,   -1, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -8.0,   3,    0.8,   3,    7.8,   2,   3,    1,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   3,    8.0,   1,   -0.8,   1,   -7.8,   1,   3,    1, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -9.0,   2,    0.9,   2,    7.9,   4,   4,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   3,    9.0,   3,   -0.9,   3,   -7.9,   1,   4,   -1, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   2,   3,    1,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   1,    1.0,   1,   -0.1,   2,   -7.1,   3,   3,    1, 11.0,   -2}
        ,{   'T', 17,  0,  2,   -5.0,   1,   -2.0,   2,    0.2,   3,    7.2,   3,   4,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    2.0,   2,   -0.2,   1,   -7.2,   1,   4,   -3, 11.0,    1}
        ,{   'T', 17,  0,  2,   -5.0,   2,   -3.0,   3,    0.3,   2,    7.3,   2,   3,    3,  0.0,    1}
        ,{   'T', 17,  0,  2,    5.0,   2,    3.0,   3,   -0.3,   3,   -7.3,   4,   3,    3, 11.0,   -1}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -4.0,   1,    0.4,   1,    7.4,   4,   4,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   3,    4.0,   2,   -0.4,   2,   -7.4,   3,   4,   -3, 11.0,    2}
        ,{   'T', 17,  0,  2,   -5.0,   3,   -5.0,   3,    0.5,   3,    7.5,   2,   3,    3,  0.0,    2}
        ,{   'T', 17,  0,  2,    5.0,   1,    5.0,   1,   -0.5,   1,   -7.5,   1,   4,    3, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -6.0,   1,    0.6,   2,    7.6,   4,   5,    1,  0.0,    1} // ku == 0
        ,{   'N', 17,  3,  0,    5.0,   1,    6.0,   2,   -0.6,   3,   -7.6,   1,   5,   -1, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -7.0,   2,    0.7,   1,    7.7,   2,   4,    1,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   2,    7.0,   3,   -0.7,   2,   -7.7,   3,   4,    1, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -8.0,   3,    0.8,   3,    7.8,   3,   5,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    8.0,   1,   -0.8,   1,   -7.8,   1,   5,   -1, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   3,   -9.0,   2,    0.9,   2,    7.9,   2,   4,    1,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    9.0,   3,   -0.9,   3,   -7.9,   4,   4,    1, 11.0,   -2}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   4,   5,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   1,    1.0,   1,   -0.1,   2,   -7.1,   3,   5,   -3, 11.0,    1}
        ,{   'N', 17,  3,  0,   -5.0,   1,   -2.0,   2,    0.2,   3,    7.2,   2,   4,    3,  0.0,    1}
        ,{   'N', 17,  3,  0,    5.0,   2,    2.0,   2,   -0.2,   1,   -7.2,   1,   4,    3, 11.0,   -1}
        ,{   'N', 17,  3,  0,   -5.0,   2,   -3.0,   3,    0.3,   2,    7.3,   4,   5,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   2,    3.0,   3,   -0.3,   3,   -7.3,   1,   5,   -3, 11.0,    2}
        ,{   'N', 17,  3,  0,   -5.0,   3,   -4.0,   1,    0.4,   1,    7.4,   2,   4,    3,  0.0,    2}
        ,{   'N', 17,  3,  0,    5.0,   3,    4.0,   2,   -0.4,   2,   -7.4,   3,   4,    3, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   3,   -5.0,   3,    0.5,   3,    7.5,   3,   5,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   1,    5.0,   1,   -0.5,   1,   -7.5,   1,   5,   -1, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -6.0,   1,    0.6,   2,    7.6,   2,   4,    1,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   1,    6.0,   2,   -0.6,   3,   -7.6,   4,   4,    1, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -7.0,   2,    0.7,   1,    7.7,   4,   5,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    7.0,   3,   -0.7,   2,   -7.7,   3,   5,   -1, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -8.0,   3,    0.8,   3,    7.8,   2,   4,    1,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   3,    8.0,   1,   -0.8,   1,   -7.8,   1,   4,    1, 11.0,   -2}
        ,{   'T', 17,  3,  0,   -5.0,   3,   -9.0,   2,    0.9,   2,    7.9,   4,   5,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   3,    9.0,   3,   -0.9,   3,   -7.9,   1,   5,   -3, 11.0,    1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   2,   4,    3,  0.0,    1}
        ,{   'T', 17,  3,  0,    5.0,   1,    1.0,   1,   -0.1,   2,   -7.1,   3,   4,    3, 11.0,   -1}
        ,{   'T', 17,  3,  0,   -5.0,   1,   -2.0,   2,    0.2,   3,    7.2,   3,   5,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    2.0,   2,   -0.2,   1,   -7.2,   1,   5,   -3, 11.0,    2}
        ,{   'T', 17,  3,  0,   -5.0,   2,   -3.0,   3,    0.3,   2,    7.3,   2,   4,    3,  0.0,    2}
        ,{   'T', 17,  3,  0,    5.0,   2,    3.0,   3,   -0.3,   3,   -7.3,   4,   4,    3, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -4.0,   1,    0.4,   1,    7.4,   4,   8,    1,  0.0,    1} // Degenerate
        ,{   'N',  5,  3,  4,    5.0,   3,    4.0,   2,   -0.4,   2,   -7.4,   3,   8,   -1, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -5.0,   3,    0.5,   3,    7.5,   2,   9,    1,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   1,    5.0,   1,   -0.5,   1,   -7.5,   1,   9,    1, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -6.0,   1,    0.6,   2,    7.6,   4,   8,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   1,    6.0,   2,   -0.6,   3,   -7.6,   1,   8,   -1, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   2,   -7.0,   2,    0.7,   1,    7.7,   2,   9,    1,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   2,    7.0,   3,   -0.7,   2,   -7.7,   3,   9,    1, 11.0,   -2}
        ,{   'N',  5,  3,  4,   -5.0,   2,   -8.0,   3,    0.8,   3,    7.8,   3,   8,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   3,    8.0,   1,   -0.8,   1,   -7.8,   1,   8,   -3, 11.0,    1}
        ,{   'N',  5,  3,  4,   -5.0,   3,   -9.0,   2,    0.9,   2,    7.9,   2,   9,    3,  0.0,    1}
        ,{   'N',  5,  3,  4,    5.0,   3,    9.0,   3,   -0.9,   3,   -7.9,   4,   9,    3, 11.0,   -1}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   4,   8,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   1,    1.0,   1,   -0.1,   2,   -7.1,   3,   8,   -3, 11.0,    2}
        ,{   'N',  5,  3,  4,   -5.0,   1,   -2.0,   2,    0.2,   3,    7.2,   2,   9,    3,  0.0,    2}
        ,{   'N',  5,  3,  4,    5.0,   2,    2.0,   2,   -0.2,   1,   -7.2,   1,   9,    3, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -3.0,   3,    0.3,   2,    7.3,   4,   8,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   2,    3.0,   3,   -0.3,   3,   -7.3,   1,   8,   -1, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -4.0,   1,    0.4,   1,    7.4,   2,   9,    1,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   3,    4.0,   2,   -0.4,   2,   -7.4,   3,   9,    1, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -5.0,   3,    0.5,   3,    7.5,   3,   8,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   1,    5.0,   1,   -0.5,   1,   -7.5,   1,   8,   -1, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   1,   -6.0,   1,    0.6,   2,    7.6,   2,   9,    1,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   1,    6.0,   2,   -0.6,   3,   -7.6,   4,   9,    1, 11.0,   -2}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -7.0,   2,    0.7,   1,    7.7,   1,   8,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   2,    7.0,   3,   -0.7,   2,   -7.7,   2,   8,   -3, 11.0,    1}
        ,{   'T',  5,  3,  4,   -5.0,   2,   -8.0,   3,    0.8,   3,    7.8,   3,   9,    3,  0.0,    1}
        ,{   'T',  5,  3,  4,    5.0,   3,    8.0,   1,   -0.8,   1,   -7.8,   1,   9,    3, 11.0,   -1}
        ,{   'T',  5,  3,  4,   -5.0,   3,   -9.0,   2,    0.9,   2,    7.9,   2,   8,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   3,    9.0,   3,   -0.9,   3,   -7.9,   3,   8,   -3, 11.0,    2}
        ,{   'T',  5,  3,  4,   -5.0,   1,   -1.0,   1,    0.1,   1,    7.1,   1,   9,    3,  0.0,    2}
        ,{   'T',  5,  3,  4,    5.0,   1,    1.0,   2,   -0.1,   2,   -7.1,   2,   9,    3, 11.0,   -2}
        ,{   'T', 17,  3,  4,    0.0,   1,    0.0,   3,    0.0,   3,    0.0,   3,   9,    3,  1.0,    1} // Quick
    };
    const size_t gcases = sizeof(gbddddmv_tc)/sizeof(gbddddmv_tc[0]);

    // Register test_gbddddmv_s cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbddddmv_s) << ' ' << gbddddmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbddddmv_s, name.str(), gbddddmv_tc + i, gbddddmv_tc + i + 1));
    }

    // Register test_gbddddmv_d cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_gbddddmv_d) << ' ' << gbddddmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_gbddddmv_d, name.str(), gbddddmv_tc + i, gbddddmv_tc + i + 1));
    }

    // Register test_gbddddmv_scc cases
    for (size_t i = 0; i < gcases; ++i) {

        gbddddmzv_tc_type c(gbddddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddddmv_scc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddddmv_scc, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.alpha2[0], c.alpha2[1]);
            std::swap(c.alpha3[0], c.alpha3[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddddmv_scc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddddmv_scc, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.alpha2[0] += 3.5;
            c.alpha3[0] += 4.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddddmv_scc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddddmv_scc, name.str(), &c, &c + 1));
        }
    }

    // Register test_gbddddmv_dzz cases
    for (size_t i = 0; i < gcases; ++i) {

        gbddddmzv_tc_type c(gbddddmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddddmv_dzz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddddmv_dzz, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha0[0], c.alpha0[1]);
            std::swap(c.alpha1[0], c.alpha1[1]);
            std::swap(c.alpha2[0], c.alpha2[1]);
            std::swap(c.alpha3[0], c.alpha3[1]);
            std::swap(c.beta[0],   c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddddmv_dzz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddddmv_dzz, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha0[0] += 1.5;
            c.alpha1[0] += 2.5;
            c.alpha2[0] += 3.5;
            c.alpha3[0] += 4.5;
            c.beta[0]   -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbddddmv_dzz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbddddmv_dzz, name.str(), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const gbddddmv_tc_type fixed_tc[] = {
        // kl, ku to be set while lda is a delta
        // trans,  n, kl, ku, alpha0, ld0, alpha1, ld1, alpha2, ld2, alpha3, ld3, lda, incx, beta, incy
         {   'N', 19,  0,  0,   -5.0,   1,   -5.1,   3,    1.5,   1,    1.5,   4,   1,    1,  0.0,    1} // Regular
        ,{   'N', 19,  0,  0,    5.0,   2,    5.1,   3,   -1.5,   2,   -1.5,   4,   1,   -1, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.2,   2,    2.5,   1,    2.5,   3,   0,    1,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.2,   2,   -2.5,   2,   -2.5,   3,   0,    1, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   2,   -5.3,   1,    3.5,   1,    3.5,   2,   1,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   3,    5.3,   1,   -3.5,   2,   -3.5,   2,   1,   -1, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   1,   -5.4,   3,    4.5,   1,    4.5,   1,   0,    1,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   2,    5.4,   2,   -4.5,   2,   -4.5,   1,   0,    1, 11.0,   -2}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.5,   1,    5.5,   1,    5.5,   4,   1,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.5,   3,   -5.5,   2,   -5.5,   4,   1,   -3, 11.0,    1}
        ,{   'N', 19,  0,  0,   -5.0,   2,   -5.6,   2,    6.5,   1,    6.5,   3,   0,    3,  0.0,    1}
        ,{   'N', 19,  0,  0,    5.0,   3,    5.6,   1,   -6.5,   2,   -6.5,   3,   0,    3, 11.0,   -1}
        ,{   'N', 19,  0,  0,   -5.0,   1,   -5.7,   3,    7.5,   1,    7.5,   2,   1,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   2,    5.7,   3,   -7.5,   2,   -7.5,   2,   1,   -3, 11.0,    2}
        ,{   'N', 19,  0,  0,   -5.0,   3,   -5.8,   2,    8.5,   1,    8.5,   1,   0,    3,  0.0,    2}
        ,{   'N', 19,  0,  0,    5.0,   1,    5.8,   2,   -8.5,   2,   -8.5,   1,   0,    3, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.9,   1,    9.5,   1,    9.5,   4,   1,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.9,   1,   -9.5,   2,   -9.5,   4,   1,   -1, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.1,   3,    1.5,   1,    1.5,   3,   0,    1,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.1,   2,   -1.5,   2,   -1.5,   3,   0,    1, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   3,   -5.2,   1,    2.5,   1,    2.5,   2,   1,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   1,    5.2,   3,   -2.5,   2,   -2.5,   2,   1,   -1, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.3,   2,    3.5,   1,    3.5,   1,   0,    1,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.3,   1,   -3.5,   2,   -3.5,   1,   0,    1, 11.0,   -2}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.4,   3,    4.5,   1,    4.5,   4,   1,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.4,   3,   -4.5,   2,   -4.5,   4,   1,   -3, 11.0,    1}
        ,{   'T', 19,  0,  0,   -5.0,   3,   -5.5,   3,    5.5,   1,    5.5,   3,   0,    3,  0.0,    1}
        ,{   'T', 19,  0,  0,    5.0,   1,    5.5,   2,   -5.5,   2,   -5.5,   3,   0,    3, 11.0,   -1}
        ,{   'T', 19,  0,  0,   -5.0,   2,   -5.6,   2,    6.5,   1,    6.5,   2,   1,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   3,    5.6,   2,   -6.5,   2,   -6.5,   2,   1,   -3, 11.0,    2}
        ,{   'T', 19,  0,  0,   -5.0,   1,   -5.7,   1,    7.5,   3,    7.5,   1,   0,    3,  0.0,    2}
        ,{   'T', 19,  0,  0,    5.0,   2,    5.7,   1,   -7.5,   3,   -7.5,   1,   0,    3, 11.0,   -2}
        ,{   'T', 17,  0,  0,    0.0,   3,    0.0,   1,    0.0,   3,    0.0,   3,   2,    3,  1.0,    1} // Quick
    };

    // Loop over fixed kl = ku = k bandwidths
    const size_t max_fixed_bandwidth = 20;
    for (size_t k = 0; k < max_fixed_bandwidth; ++k) {

        // Loop over fixed_tc and register both real and complex tests
        for (size_t i = 0; i < sizeof(fixed_tc)/sizeof(fixed_tc[0]); ++i) {

            // Prepare real-valued test details for bandwidth k
            gbddddmv_tc_type r(fixed_tc[i]);
            r.kl   = r.ku = k;
            r.lda += (r.kl + r.ku + 1);

            { // Register test_gbddddmv_s case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbddddmv_s) << ' ' << gbddddmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbddddmv_s, name.str(), &r, &r + 1));
            }

            { // Register test_gbddddmv_d case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_gbddddmv_d) << ' ' << gbddddmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_gbddddmv_d, name.str(), &r, &r + 1));
            }

            { // Register test_gbddddmv_scc cases
                gbddddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddddmv_scc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddddmv_scc, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.alpha2[0], c.alpha2[1]);
                    std::swap(c.alpha3[0], c.alpha3[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddddmv_scc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddddmv_scc, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.alpha2[0] += 3.5;
                    c.alpha3[0] += 4.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddddmv_scc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddddmv_scc, name.str(), &c, &c + 1));
                }
            }

            { // Register test_gbddddmv_dzz cases
                gbddddmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddddmv_dzz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddddmv_dzz, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha0[0], c.alpha0[1]);
                    std::swap(c.alpha1[0], c.alpha1[1]);
                    std::swap(c.alpha2[0], c.alpha2[1]);
                    std::swap(c.alpha3[0], c.alpha3[1]);
                    std::swap(c.beta[0],   c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddddmv_dzz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddddmv_dzz, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha0[0] += 1.5;
                    c.alpha1[0] += 2.5;
                    c.alpha2[0] += 3.5;
                    c.alpha3[0] += 4.5;
                    c.beta[0]   -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbddddmv_dzz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbddddmv_dzz, name.str(), &c, &c + 1));
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
