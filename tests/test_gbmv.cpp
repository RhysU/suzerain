#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/blas_et_al.h>
#include <suzerain/gbmv.h>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

using boost::unit_test::framework::master_test_suite;
using boost::unit_test::make_test_case;
using std::numeric_limits;
using std::size_t;

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
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lena; ++i) a[i] = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i] = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i] = y[i] = random() * inv_rand_max;

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
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lena; ++i) a[i] = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i] = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i] = y[i] = random() * inv_rand_max;

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

static void test_gbmv_sc(const gbmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.m*t.n*20;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = 2 * abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lena; ++i) a[i] = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i] = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i] = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const float alpha[2] = { t.alpha[0], t.alpha[1] };
    const float beta[2]  = { t.beta[0],  t.beta[1]  };

    // Compute expected result using external BLAS
    suzerain_blasext_sgbmzv_external(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (const float(*)[2]) x.get(), t.incx,
            beta,                  (      float(*)[2]) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbmv_sc(
            t.trans, t.m, t.n, t.kl, t.ku,
            alpha, a.get(), t.lda, (const float(*)[2]) x.get(), t.incx,
            beta,                  (      float(*)[2]) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbmv_dz(const gbmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.m*t.n*25;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * (toupper(t.trans) == 'N' ? t.n : t.m);
    const int leny = 2 * abs(t.incy) * (toupper(t.trans) == 'N' ? t.m : t.n);

    // Allocate random data for testing purposes
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lena; ++i) a[i] = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i] = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i] = y[i] = random() * inv_rand_max;

    // Compute expected result using external BLAS
    suzerain_blasext_dgbmzv_external(
            t.trans, t.m, t.n, t.kl, t.ku,
            t.alpha, a.get(), t.lda, (const double(*)[2]) x.get(), t.incx,
            t.beta,                  (      double(*)[2]) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_gbmv_dz(
            t.trans, t.m, t.n, t.kl, t.ku,
            t.alpha, a.get(), t.lda, (const double(*)[2]) x.get(), t.incx,
            t.beta,                  (      double(*)[2]) y.get(), t.incy));

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
        ,{   'N',  4,  5,  3,  4,  -5.0,   8,    1,  0.0,    1} // Degenerate
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

    // Register test_gbmv_sc cases
    for (size_t i = 0; i < gcases; ++i) {

        gbmzv_tc_type c(gbmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_sc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_sc, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_sc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_sc, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_sc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_sc, name.str(), &c, &c + 1));
        }
    }

    // Register test_gbmv_dz cases
    for (size_t i = 0; i < gcases; ++i) {

        gbmzv_tc_type c(gbmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_dz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_dz, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_dz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_dz, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_gbmv_dz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_gbmv_dz, name.str(), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const gbmv_tc_type fixed_tc[] = {
        // kl, ku to be set while lda is a delta
        // trans,  m,  n, kl, ku, alpha, lda, incx, beta, incy
         {   'N', 17, 19,  0,  0,  -5.0,   1,    1,  0.0,    1} // Regular
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
        ,{   'T', 19, 17,  0,  0,   0.0,   2,    3,  1.0,    1} // Quick
    };

    // Loop over fixed kl = ku = k bandwidths
    const size_t max_fixed_bandwidth = 9;
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

            { // Register test_gbmv_sc cases
                gbmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_sc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_sc, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_sc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_sc, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_sc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_sc, name.str(), &c, &c + 1));
                }
            }

            { // Register test_gbmv_dz cases
                gbmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_dz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_dz, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_dz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_dz, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_gbmv_dz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_gbmv_dz, name.str(), &c, &c + 1));
                }
            }

        }
    }

    return 0;
}
