#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/blas_et_al.h>
#include <suzerain/sbmv.h>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

using boost::unit_test::framework::master_test_suite;
using boost::unit_test::make_test_case;
using std::numeric_limits;
using std::size_t;

// Most basic of upper storage test cases to help weed on glaring errors
// If we cannot nail baby_steps, we've got no right to do everything else.
BOOST_AUTO_TEST_CASE( baby_steps )
{
    const float aU[] = { 555, 555,  11,   // Upper triangular storage
                         555,  12,  22,
                          13,  23,  33,
                          24,  34,  44,
                          35,  45,  55 };
    const float aL[] = {  11,  12,  13,   // Equivalent lower storage
                          22,  23,  24,
                          33,  34,  35,
                          44,  45, 555,
                          55, 555, 555 };
    const int n         = 5;
    const int k         = 2;
    const float alpha   = 3;
    const int lda       = 3;
    const float x[]     = { 1, 2, 3, 4, 5 };
    const int incx      = 1;
    const float beta    = 1;
    const float yinit[] = { 10, 11, 12, 13, 14 };
    const int incy      = 1;
    const float e[]     = { 232, 674, 1419, 1666, 1694 };
    const float close_enough = numeric_limits<float>::epsilon()*n*n*15;

    // Scratch storage
    float y[sizeof(yinit)/sizeof(yinit[0])];

    // Check that BLAS kicks back what we expect for y := alpha*aU*x + beta*y
    suzerain_blas_scopy(n, yinit, incy, y, incy);
    suzerain_blas_ssbmv_external(
            'U', n, k, alpha, aU, lda, x, incx, beta, y, incy);
    check_close_collections(e, e + n, y, y + n, close_enough);

    // Check that BLAS kicks back what we expect for y := alpha*aL*x + beta*y
    suzerain_blas_scopy(n, yinit, incy, y, incy);
    suzerain_blas_ssbmv_external(
            'L', n, k, alpha, aL, lda, x, incx, beta, y, incy);
    check_close_collections(e, e + n, y, y + n, close_enough);

    // Does our implementation return the expected values for 'U'?
    BOOST_TEST_MESSAGE("Baby steps: 'U'");
    suzerain_blas_scopy(n, yinit, incy, y, incy);
    suzerain_sbmv_s('U', n, k, alpha, aU, lda, x, incx, beta, y, incy);
    check_close_collections(e, e + n, y, y + n, close_enough);


    // Does our implementation return the expected values for 'L'?
    BOOST_TEST_MESSAGE("Baby steps: 'L'");
    suzerain_blas_scopy(n, yinit, incy, y, incy);
    suzerain_sbmv_s('L', n, k, alpha, aL, lda, x, incx, beta, y, incy);
    check_close_collections(e, e + n, y, y + n, close_enough);
}


// Test suzerain_sbmv_?  against suzerain_blas_*sbmv_external
// Test suzerain_sbmv_?? against suzerain_blasext_?sbmzv_external
// Idea behind testing is that matching the BLAS is goodness

// For unary function-based test case registration
struct sbmv_tc_type
{
    char uplo;
    int n, k;
    double alpha;
    int lda, incx;
    double beta;
    int incy;
};

// For unary function-based test case registration
struct sbmzv_tc_type
{
    char uplo;
    int n, k;
    double alpha[2];
    int lda, incx;
    double beta[2];
    int incy;

    sbmzv_tc_type(const sbmv_tc_type& o)
        : uplo(o.uplo), n(o.n), k(o.k),
          lda(o.lda), incx(o.incx), incy(o.incy)
    {
        alpha[0] = o.alpha; alpha[1] = 0;
        beta[0]  = o.beta;  beta[1]  = 0;
    }
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const sbmv_tc_type& t)
{
    os << "{uplo="   << t.uplo
       << ", n="     << t.n
       << ", k="     << t.k
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
        std::basic_ostream<charT,traits> &os, const sbmzv_tc_type& t)
{
    os << "{uplo="   << t.uplo
       << ", n="     << t.n
       << ", k="     << t.k
       << ", alpha=" << std::complex<double>(t.alpha[0], t.alpha[1])
       << ", lda="   << t.lda
       << ", incx="  << t.incx
       << ", beta="  << std::complex<double>(t.beta[0], t.beta[1])
       << ", incy="  << t.incy
       << '}';
    return os;
}

static void test_sbmv_s(const sbmv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*15;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * t.n;
    const int leny = abs(t.incy) * t.n;

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
    suzerain_blas_ssbmv_external(t.uplo, t.n, t.k,
                                 alpha, a.get(), t.lda, x.get(), t.incx,
                                 beta,                  e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_sbmv_s(
                t.uplo, t.n, t.k,
                alpha, a.get(), t.lda, x.get(), t.incx,
                beta,                  y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_sbmv_d(const sbmv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*15;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = abs(t.incx) * t.n;
    const int leny = abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lena; ++i) a[i] = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i] = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i] = y[i] = random() * inv_rand_max;

    // Compute expected result using external BLAS
    suzerain_blas_dsbmv_external(t.uplo, t.n, t.k,
                                 t.alpha, a.get(), t.lda, x.get(), t.incx,
                                 t.beta,                  e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_sbmv_d(
                t.uplo, t.n, t.k,
                t.alpha, a.get(), t.lda, x.get(), t.incx,
                t.beta,                  y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_sbmv_sc(const sbmzv_tc_type& t)
{
    const float close_enough = numeric_limits<float>::epsilon()*t.n*t.n*50;
    const float inv_rand_max = float(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<float> a(new float[lena]);
    boost::scoped_array<float> x(new float[lenx]);
    boost::scoped_array<float> y(new float[leny]), e(new float[leny]);
    for (int i = 0; i < lena; ++i) a[i] = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i] = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i] = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const complex_float alpha( t.alpha[0], t.alpha[1] );
    const complex_float beta ( t.beta[0],  t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blas_csbmv_s_external(
            t.uplo, t.n, t.k,
            alpha, a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,                  (      complex_float *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_sbmv_sc(
            t.uplo, t.n, t.k,
            alpha, a.get(), t.lda, (const complex_float *) x.get(), t.incx,
            beta,                  (      complex_float *) y.get(), t.incy));

    check_close_collections(e.get(), e.get() + leny,
                            y.get(), y.get() + leny,
                            close_enough);
}

static void test_sbmv_dz(const sbmzv_tc_type& t)
{
    const double close_enough = numeric_limits<double>::epsilon()*t.n*t.n*50;
    const double inv_rand_max = double(1) / RAND_MAX;
    const int lena = t.lda * t.n;
    const int lenx = 2 * abs(t.incx) * t.n;
    const int leny = 2 * abs(t.incy) * t.n;

    // Allocate random data for testing purposes
    boost::scoped_array<double> a(new double[lena]);
    boost::scoped_array<double> x(new double[lenx]);
    boost::scoped_array<double> y(new double[leny]), e(new double[leny]);
    for (int i = 0; i < lena; ++i) a[i] = random() * inv_rand_max;
    for (int i = 0; i < lenx; ++i) x[i] = random() * inv_rand_max;
    for (int i = 0; i < leny; ++i) e[i] = y[i] = random() * inv_rand_max;

    // Get appropriately typed alpha and beta constants
    const complex_double alpha( t.alpha[0], t.alpha[1] );
    const complex_double beta ( t.beta[0],  t.beta[1]  );

    // Compute expected result using external BLAS
    suzerain_blas_zsbmv_d_external(
            t.uplo, t.n, t.k,
            alpha, a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,                  (      complex_double *) e.get(), t.incy);

    // Compute observed result using our implementation
    BOOST_REQUIRE_EQUAL(0, suzerain_sbmv_dz(
            t.uplo, t.n, t.k,
            alpha, a.get(), t.lda, (const complex_double *) x.get(), t.incx,
            beta,                  (      complex_double *) y.get(), t.incy));

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

    const sbmv_tc_type sbmv_tc[] = {
        //  uplo,  n,  k, alpha, lda, incx, beta, incy
         {   'U', 19,  3,  -5.0,   7,    1,  0.0,    1} // Regular
        ,{   'U', 19,  3,   5.0,   7,   -1, 11.0,    1}
        ,{   'U', 19,  3,  -5.0,   6,    1,  0.0,    1}
        ,{   'U', 19,  3,   5.0,   6,    1, 11.0,   -1}
        ,{   'U', 19,  3,  -5.0,   7,    1,  0.0,    2}
        ,{   'U', 19,  3,   5.0,   7,   -1, 11.0,    2}
        ,{   'U', 19,  3,  -5.0,   6,    1,  0.0,    2}
        ,{   'U', 19,  3,   5.0,   6,    1, 11.0,   -2}
        ,{   'U', 19,  3,  -5.0,   7,    3,  0.0,    1}
        ,{   'U', 19,  3,   5.0,   7,   -3, 11.0,    1}
        ,{   'U', 19,  3,  -5.0,   6,    3,  0.0,    1}
        ,{   'U', 19,  3,   5.0,   6,    3, 11.0,   -1}
        ,{   'U', 19,  3,  -5.0,   7,    3,  0.0,    2}
        ,{   'U', 19,  3,   5.0,   7,   -3, 11.0,    2}
        ,{   'U', 19,  3,  -5.0,   6,    3,  0.0,    2}
        ,{   'U', 19,  3,   5.0,   6,    3, 11.0,   -2}
        ,{   'L', 19,  3,  -5.0,   7,    1,  0.0,    1}
        ,{   'L', 19,  3,   5.0,   7,   -1, 11.0,    1}
        ,{   'L', 19,  3,  -5.0,   6,    1,  0.0,    1}
        ,{   'L', 19,  3,   5.0,   6,    1, 11.0,   -1}
        ,{   'L', 19,  3,  -5.0,   7,    1,  0.0,    2}
        ,{   'L', 19,  3,   5.0,   7,   -1, 11.0,    2}
        ,{   'L', 19,  3,  -5.0,   6,    1,  0.0,    2}
        ,{   'L', 19,  3,   5.0,   6,    1, 11.0,   -2}
        ,{   'L', 19,  3,  -5.0,   7,    3,  0.0,    1}
        ,{   'L', 19,  3,   5.0,   7,   -3, 11.0,    1}
        ,{   'L', 19,  3,  -5.0,   6,    3,  0.0,    1}
        ,{   'L', 19,  3,   5.0,   6,    3, 11.0,   -1}
        ,{   'L', 19,  3,  -5.0,   7,    3,  0.0,    2}
        ,{   'L', 19,  3,   5.0,   7,   -3, 11.0,    2}
        ,{   'L', 19,  3,  -5.0,   6,    3,  0.0,    2}
        ,{   'L', 19,  3,   5.0,   6,    3, 11.0,   -2}
        ,{   'U', 17,  0,  -5.0,   4,    1,  0.0,    1} // k == 0
        ,{   'U', 17,  0,   5.0,   4,   -1, 11.0,    1}
        ,{   'U', 17,  0,  -5.0,   3,    1,  0.0,    1}
        ,{   'U', 17,  0,   5.0,   3,    1, 11.0,   -1}
        ,{   'U', 17,  0,  -5.0,   4,    1,  0.0,    2}
        ,{   'U', 17,  0,   5.0,   4,   -1, 11.0,    2}
        ,{   'U', 17,  0,  -5.0,   3,    1,  0.0,    2}
        ,{   'U', 17,  0,   5.0,   3,    1, 11.0,   -2}
        ,{   'U', 17,  0,  -5.0,   4,    3,  0.0,    1}
        ,{   'U', 17,  0,   5.0,   4,   -3, 11.0,    1}
        ,{   'U', 17,  0,  -5.0,   3,    3,  0.0,    1}
        ,{   'U', 17,  0,   5.0,   3,    3, 11.0,   -1}
        ,{   'U', 17,  0,  -5.0,   4,    3,  0.0,    2}
        ,{   'U', 17,  0,   5.0,   4,   -3, 11.0,    2}
        ,{   'U', 17,  0,  -5.0,   3,    3,  0.0,    2}
        ,{   'U', 17,  0,   5.0,   3,    3, 11.0,   -2}
        ,{   'L', 17,  0,  -5.0,   4,    1,  0.0,    1}
        ,{   'L', 17,  0,   5.0,   4,   -1, 11.0,    1}
        ,{   'L', 17,  0,  -5.0,   3,    1,  0.0,    1}
        ,{   'L', 17,  0,   5.0,   3,    1, 11.0,   -1}
        ,{   'L', 17,  0,  -5.0,   4,    1,  0.0,    2}
        ,{   'L', 17,  0,   5.0,   4,   -1, 11.0,    2}
        ,{   'L', 17,  0,  -5.0,   3,    1,  0.0,    2}
        ,{   'L', 17,  0,   5.0,   3,    1, 11.0,   -2}
        ,{   'L', 17,  0,  -5.0,   4,    3,  0.0,    1}
        ,{   'L', 17,  0,   5.0,   4,   -3, 11.0,    1}
        ,{   'L', 17,  0,  -5.0,   3,    3,  0.0,    1}
        ,{   'L', 17,  0,   5.0,   3,    3, 11.0,   -1}
        ,{   'L', 17,  0,  -5.0,   4,    3,  0.0,    2}
        ,{   'L', 17,  0,   5.0,   4,   -3, 11.0,    2}
        ,{   'L', 17,  0,  -5.0,   3,    3,  0.0,    2}
        ,{   'L', 17,  0,   5.0,   4,    3, 11.0,   -2}
        ,{   'U',  5,  3,  -5.0,   8,    1,  0.0,    1} // Degenerate
        ,{   'U',  5,  3,   5.0,   8,   -1, 11.0,    1}
        ,{   'U',  5,  3,  -5.0,   9,    1,  0.0,    1}
        ,{   'U',  5,  3,   5.0,   9,    1, 11.0,   -1}
        ,{   'U',  5,  3,  -5.0,   8,    1,  0.0,    2}
        ,{   'U',  5,  3,   5.0,   8,   -1, 11.0,    2}
        ,{   'U',  5,  3,  -5.0,   9,    1,  0.0,    2}
        ,{   'U',  5,  3,   5.0,   9,    1, 11.0,   -2}
        ,{   'U',  5,  3,  -5.0,   8,    3,  0.0,    1}
        ,{   'U',  5,  3,   5.0,   8,   -3, 11.0,    1}
        ,{   'U',  5,  3,  -5.0,   9,    3,  0.0,    1}
        ,{   'U',  5,  3,   5.0,   9,    3, 11.0,   -1}
        ,{   'U',  5,  3,  -5.0,   8,    3,  0.0,    2}
        ,{   'U',  5,  3,   5.0,   8,   -3, 11.0,    2}
        ,{   'U',  5,  3,  -5.0,   9,    3,  0.0,    2}
        ,{   'U',  5,  3,   5.0,   9,    3, 11.0,   -2}
        ,{   'L',  5,  3,  -5.0,   8,    1,  0.0,    1}
        ,{   'L',  5,  3,   5.0,   8,   -1, 11.0,    1}
        ,{   'L',  5,  3,  -5.0,   9,    1,  0.0,    1}
        ,{   'L',  5,  3,   5.0,   9,    1, 11.0,   -1}
        ,{   'L',  5,  3,  -5.0,   8,    1,  0.0,    2}
        ,{   'L',  5,  3,   5.0,   8,   -1, 11.0,    2}
        ,{   'L',  5,  3,  -5.0,   9,    1,  0.0,    2}
        ,{   'L',  5,  3,   5.0,   9,    1, 11.0,   -2}
        ,{   'L',  5,  3,  -5.0,   8,    3,  0.0,    1}
        ,{   'L',  5,  3,   5.0,   8,   -3, 11.0,    1}
        ,{   'L',  5,  3,  -5.0,   9,    3,  0.0,    1}
        ,{   'L',  5,  3,   5.0,   9,    3, 11.0,   -1}
        ,{   'L',  5,  3,  -5.0,   8,    3,  0.0,    2}
        ,{   'L',  5,  3,   5.0,   8,   -3, 11.0,    2}
        ,{   'L',  5,  3,  -5.0,   9,    3,  0.0,    2}
        ,{   'L',  5,  3,   5.0,   9,    3, 11.0,   -2}
        ,{   'L', 17,  3,   0.0,   9,    3,  1.0,    1} // Quick
    };
    const size_t gcases = sizeof(sbmv_tc)/sizeof(sbmv_tc[0]);

    // Register test_sbmv_s cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_sbmv_s) << ' ' << sbmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_sbmv_s, name.str(), sbmv_tc + i, sbmv_tc + i + 1));
    }

    // Register test_sbmv_d cases
    for (size_t i = 0; i < gcases; ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(test_sbmv_d) << ' ' << sbmv_tc[i];
        master_test_suite().add(make_test_case(
                &test_sbmv_d, name.str(), sbmv_tc + i, sbmv_tc + i + 1));
    }

    // Register test_sbmv_sc cases
    for (size_t i = 0; i < gcases; ++i) {

        sbmzv_tc_type c(sbmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_sbmv_sc) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_sbmv_sc, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_sbmv_sc) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_sbmv_sc, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_sbmv_sc) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_sbmv_sc, name.str(), &c, &c + 1));
        }
    }

    // Register test_sbmv_dz cases
    for (size_t i = 0; i < gcases; ++i) {

        sbmzv_tc_type c(sbmv_tc[i]);

        { // Real-valued alpha, beta
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_sbmv_dz) << " real " << c;
            master_test_suite().add(make_test_case(
                    &test_sbmv_dz, name.str(), &c, &c + 1));
        }

        { // Imaginary-valued alpha, beta
            std::swap(c.alpha[0], c.alpha[1]);
            std::swap(c.beta[0],  c.beta[1]);
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_sbmv_dz) << " imag " << c;
            master_test_suite().add(make_test_case(
                    &test_sbmv_dz, name.str(), &c, &c + 1));
        }

        { // Truly complex alpha, beta
            c.alpha[0] += 1.5;
            c.beta[0]  -= 1.5;
            std::ostringstream name;
            name << BOOST_TEST_STRINGIZE(test_sbmv_dz) << " complex " << c;
            master_test_suite().add(make_test_case(
                    &test_sbmv_dz, name.str(), &c, &c + 1));
        }
    }

    // -------------------------------------------------------
    // Register test cases designed for fixed bandwidths
    // -------------------------------------------------------

    const sbmv_tc_type fixed_tc[] = {
        // k to be set while lda is a delta
        //  uplo,  n,  k, alpha, lda, incx, beta, incy
         {   'U', 19,  0,  -5.0,   1,    1,  0.0,    1} // Regular
        ,{   'U', 19,  0,   5.0,   1,   -1, 11.0,    1}
        ,{   'U', 19,  0,  -5.0,   0,    1,  0.0,    1}
        ,{   'U', 19,  0,   5.0,   0,    1, 11.0,   -1}
        ,{   'U', 19,  0,  -5.0,   1,    1,  0.0,    2}
        ,{   'U', 19,  0,   5.0,   1,   -1, 11.0,    2}
        ,{   'U', 19,  0,  -5.0,   0,    1,  0.0,    2}
        ,{   'U', 19,  0,   5.0,   0,    1, 11.0,   -2}
        ,{   'U', 19,  0,  -5.0,   1,    3,  0.0,    1}
        ,{   'U', 19,  0,   5.0,   1,   -3, 11.0,    1}
        ,{   'U', 19,  0,  -5.0,   0,    3,  0.0,    1}
        ,{   'U', 19,  0,   5.0,   0,    3, 11.0,   -1}
        ,{   'U', 19,  0,  -5.0,   1,    3,  0.0,    2}
        ,{   'U', 19,  0,   5.0,   1,   -3, 11.0,    2}
        ,{   'U', 19,  0,  -5.0,   0,    3,  0.0,    2}
        ,{   'U', 19,  0,   5.0,   0,    3, 11.0,   -2}
        ,{   'L', 19,  0,  -5.0,   1,    1,  0.0,    1}
        ,{   'L', 19,  0,   5.0,   1,   -1, 11.0,    1}
        ,{   'L', 19,  0,  -5.0,   0,    1,  0.0,    1}
        ,{   'L', 19,  0,   5.0,   0,    1, 11.0,   -1}
        ,{   'L', 19,  0,  -5.0,   1,    1,  0.0,    2}
        ,{   'L', 19,  0,   5.0,   1,   -1, 11.0,    2}
        ,{   'L', 19,  0,  -5.0,   0,    1,  0.0,    2}
        ,{   'L', 19,  0,   5.0,   0,    1, 11.0,   -2}
        ,{   'L', 19,  0,  -5.0,   1,    3,  0.0,    1}
        ,{   'L', 19,  0,   5.0,   1,   -3, 11.0,    1}
        ,{   'L', 19,  0,  -5.0,   0,    3,  0.0,    1}
        ,{   'L', 19,  0,   5.0,   0,    3, 11.0,   -1}
        ,{   'L', 19,  0,  -5.0,   1,    3,  0.0,    2}
        ,{   'L', 19,  0,   5.0,   1,   -3, 11.0,    2}
        ,{   'L', 19,  0,  -5.0,   0,    3,  0.0,    2}
        ,{   'L', 19,  0,   5.0,   0,    3, 11.0,   -2}
        ,{   'L', 17,  0,   0.0,   2,    3,  1.0,    1} // Quick
    };

    // Loop over fixed bandwidths
    const size_t max_fixed_bandwidth = 20;
    for (size_t k = 0; k < max_fixed_bandwidth; ++k) {

        // Loop over fixed_tc and register both real and complex tests
        for (size_t i = 0; i < sizeof(fixed_tc)/sizeof(fixed_tc[0]); ++i) {

            // Prepare real-valued test details for bandwidth k
            sbmv_tc_type r(fixed_tc[i]);
            r.k = k;
            r.lda += (r.k + 1);

            { // Register test_sbmv_s case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_sbmv_s) << ' ' << sbmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_sbmv_s, name.str(), &r, &r + 1));
            }

            { // Register test_sbmv_d case
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(test_sbmv_d) << ' ' << sbmv_tc[i];
                master_test_suite().add(make_test_case(
                        &test_sbmv_d, name.str(), &r, &r + 1));
            }

            { // Register test_sbmv_sc cases
                sbmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_sbmv_sc) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_sbmv_sc, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_sbmv_sc) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_sbmv_sc, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_sbmv_sc) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_sbmv_sc, name.str(), &c, &c + 1));
                }
            }

            { // Register test_sbmv_dz cases
                sbmzv_tc_type c(r);

                { // Real-valued alpha, beta
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_sbmv_dz) << " real " << c;
                    master_test_suite().add(make_test_case(
                            &test_sbmv_dz, name.str(), &c, &c + 1));
                }

                { // Imaginary-valued alpha, beta
                    std::swap(c.alpha[0], c.alpha[1]);
                    std::swap(c.beta[0],  c.beta[1]);
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_sbmv_dz) << " imag " << c;
                    master_test_suite().add(make_test_case(
                            &test_sbmv_dz, name.str(), &c, &c + 1));
                }

                { // Truly complex alpha, beta
                    c.alpha[0] += 1.5;
                    c.beta[0]  -= 1.5;
                    std::ostringstream name;
                    name << BOOST_TEST_STRINGIZE(test_sbmv_dz) << " complex " << c;
                    master_test_suite().add(make_test_case(
                            &test_sbmv_dz, name.str(), &c, &c + 1));
                }
            }

        }
    }

    return 0;
}
