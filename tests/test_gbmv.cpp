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

using std::numeric_limits;

// Test suzerain_gbmv_?  against suzerain_blas_*gbmv_external
// Test suzerain_gbmb_?? against suzerain_blasext_?gbmzv_external
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
    BOOST_TEST_MESSAGE("\tscenario " << t);

    const float close_enough = numeric_limits<float>::epsilon()*t.m*t.n*10;
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
    BOOST_REQUIRE_EQUAL(0,
    suzerain_gbmv_s(t.trans, t.m, t.n, t.kl, t.ku,
                      alpha, a.get(), t.lda, x.get(), t.incx,
                      beta,                  y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_gbmv_d(const gbmv_tc_type& t)
{
    BOOST_TEST_MESSAGE("\tscenario " << t);

    const double close_enough = numeric_limits<double>::epsilon()*t.m*t.n;
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
    suzerain_gbmv_d(t.trans, t.m, t.n, t.kl, t.ku,
                    t.alpha, a.get(), t.lda, x.get(), t.incx,
                    t.beta,                  y.get(), t.incy);

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

boost::unit_test::test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    (void) argc; // Unused
    (void) argv; // Unused
    boost::unit_test::framework::master_test_suite().p_name.value = __FILE__;

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
    };

    {
        const std::size_t ncases = sizeof(gbmv_tc)/sizeof(gbmv_tc[0]);

        // Register test_gbmv_s cases
        boost::unit_test::framework::master_test_suite().add(
                BOOST_PARAM_TEST_CASE(&test_gbmv_s, gbmv_tc, gbmv_tc + ncases),
                /* timeout in seconds */ 10);

        // Register test_gbmv_d cases
        boost::unit_test::framework::master_test_suite().add(
                BOOST_PARAM_TEST_CASE(&test_gbmv_d, gbmv_tc, gbmv_tc + ncases),
                /* timeout in seconds */ 10);
    }

    return 0;
}
