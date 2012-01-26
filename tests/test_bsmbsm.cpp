#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <suzerain/bsmbsm.h>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/complex.hpp>
#include <boost/test/included/unit_test.hpp>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

BOOST_AUTO_TEST_SUITE( permutation )

BOOST_AUTO_TEST_CASE( q_S5n9 )
{
    // Compare against Octave-based results for
    //      q    = mod(0:44,5)*9 + floor((0:44)/5)
    // for case when S = 5, n = 9
    const int S = 5, n = 9;
    BOOST_CHECK_EQUAL( 0, suzerain_bsmbsm_q(S, n,  0));
    BOOST_CHECK_EQUAL( 9, suzerain_bsmbsm_q(S, n,  1));
    BOOST_CHECK_EQUAL(18, suzerain_bsmbsm_q(S, n,  2));
    BOOST_CHECK_EQUAL(27, suzerain_bsmbsm_q(S, n,  3));
    BOOST_CHECK_EQUAL(36, suzerain_bsmbsm_q(S, n,  4));
    BOOST_CHECK_EQUAL( 1, suzerain_bsmbsm_q(S, n,  5));
    BOOST_CHECK_EQUAL(10, suzerain_bsmbsm_q(S, n,  6));
    BOOST_CHECK_EQUAL(19, suzerain_bsmbsm_q(S, n,  7));
    BOOST_CHECK_EQUAL(28, suzerain_bsmbsm_q(S, n,  8));
    BOOST_CHECK_EQUAL(37, suzerain_bsmbsm_q(S, n,  9));
    BOOST_CHECK_EQUAL( 2, suzerain_bsmbsm_q(S, n, 10));
    BOOST_CHECK_EQUAL(11, suzerain_bsmbsm_q(S, n, 11));
    BOOST_CHECK_EQUAL(20, suzerain_bsmbsm_q(S, n, 12));
    BOOST_CHECK_EQUAL(29, suzerain_bsmbsm_q(S, n, 13));
    BOOST_CHECK_EQUAL(38, suzerain_bsmbsm_q(S, n, 14));
    BOOST_CHECK_EQUAL( 3, suzerain_bsmbsm_q(S, n, 15));
    BOOST_CHECK_EQUAL(12, suzerain_bsmbsm_q(S, n, 16));
    BOOST_CHECK_EQUAL(21, suzerain_bsmbsm_q(S, n, 17));
    BOOST_CHECK_EQUAL(30, suzerain_bsmbsm_q(S, n, 18));
    BOOST_CHECK_EQUAL(39, suzerain_bsmbsm_q(S, n, 19));
    BOOST_CHECK_EQUAL( 4, suzerain_bsmbsm_q(S, n, 20));
    BOOST_CHECK_EQUAL(13, suzerain_bsmbsm_q(S, n, 21));
    BOOST_CHECK_EQUAL(22, suzerain_bsmbsm_q(S, n, 22));
    BOOST_CHECK_EQUAL(31, suzerain_bsmbsm_q(S, n, 23));
    BOOST_CHECK_EQUAL(40, suzerain_bsmbsm_q(S, n, 24));
    BOOST_CHECK_EQUAL( 5, suzerain_bsmbsm_q(S, n, 25));
    BOOST_CHECK_EQUAL(14, suzerain_bsmbsm_q(S, n, 26));
    BOOST_CHECK_EQUAL(23, suzerain_bsmbsm_q(S, n, 27));
    BOOST_CHECK_EQUAL(32, suzerain_bsmbsm_q(S, n, 28));
    BOOST_CHECK_EQUAL(41, suzerain_bsmbsm_q(S, n, 29));
    BOOST_CHECK_EQUAL( 6, suzerain_bsmbsm_q(S, n, 30));
    BOOST_CHECK_EQUAL(15, suzerain_bsmbsm_q(S, n, 31));
    BOOST_CHECK_EQUAL(24, suzerain_bsmbsm_q(S, n, 32));
    BOOST_CHECK_EQUAL(33, suzerain_bsmbsm_q(S, n, 33));
    BOOST_CHECK_EQUAL(42, suzerain_bsmbsm_q(S, n, 34));
    BOOST_CHECK_EQUAL( 7, suzerain_bsmbsm_q(S, n, 35));
    BOOST_CHECK_EQUAL(16, suzerain_bsmbsm_q(S, n, 36));
    BOOST_CHECK_EQUAL(25, suzerain_bsmbsm_q(S, n, 37));
    BOOST_CHECK_EQUAL(34, suzerain_bsmbsm_q(S, n, 38));
    BOOST_CHECK_EQUAL(43, suzerain_bsmbsm_q(S, n, 39));
    BOOST_CHECK_EQUAL( 8, suzerain_bsmbsm_q(S, n, 40));
    BOOST_CHECK_EQUAL(17, suzerain_bsmbsm_q(S, n, 41));
    BOOST_CHECK_EQUAL(26, suzerain_bsmbsm_q(S, n, 42));
    BOOST_CHECK_EQUAL(35, suzerain_bsmbsm_q(S, n, 43));
    BOOST_CHECK_EQUAL(44, suzerain_bsmbsm_q(S, n, 44));
}

BOOST_AUTO_TEST_CASE( qinv_S5n9 )
{
    // Compare against Octave-based results for
    //      qinv = mod(0:44,9)*5 + floor((0:44)/9);
    // for case when S = 5, n = 9
    const int S = 5, n = 9;
    BOOST_CHECK_EQUAL( 0, suzerain_bsmbsm_qinv(S, n,  0));
    BOOST_CHECK_EQUAL( 5, suzerain_bsmbsm_qinv(S, n,  1));
    BOOST_CHECK_EQUAL(10, suzerain_bsmbsm_qinv(S, n,  2));
    BOOST_CHECK_EQUAL(15, suzerain_bsmbsm_qinv(S, n,  3));
    BOOST_CHECK_EQUAL(20, suzerain_bsmbsm_qinv(S, n,  4));
    BOOST_CHECK_EQUAL(25, suzerain_bsmbsm_qinv(S, n,  5));
    BOOST_CHECK_EQUAL(30, suzerain_bsmbsm_qinv(S, n,  6));
    BOOST_CHECK_EQUAL(35, suzerain_bsmbsm_qinv(S, n,  7));
    BOOST_CHECK_EQUAL(40, suzerain_bsmbsm_qinv(S, n,  8));
    BOOST_CHECK_EQUAL( 1, suzerain_bsmbsm_qinv(S, n,  9));
    BOOST_CHECK_EQUAL( 6, suzerain_bsmbsm_qinv(S, n, 10));
    BOOST_CHECK_EQUAL(11, suzerain_bsmbsm_qinv(S, n, 11));
    BOOST_CHECK_EQUAL(16, suzerain_bsmbsm_qinv(S, n, 12));
    BOOST_CHECK_EQUAL(21, suzerain_bsmbsm_qinv(S, n, 13));
    BOOST_CHECK_EQUAL(26, suzerain_bsmbsm_qinv(S, n, 14));
    BOOST_CHECK_EQUAL(31, suzerain_bsmbsm_qinv(S, n, 15));
    BOOST_CHECK_EQUAL(36, suzerain_bsmbsm_qinv(S, n, 16));
    BOOST_CHECK_EQUAL(41, suzerain_bsmbsm_qinv(S, n, 17));
    BOOST_CHECK_EQUAL( 2, suzerain_bsmbsm_qinv(S, n, 18));
    BOOST_CHECK_EQUAL( 7, suzerain_bsmbsm_qinv(S, n, 19));
    BOOST_CHECK_EQUAL(12, suzerain_bsmbsm_qinv(S, n, 20));
    BOOST_CHECK_EQUAL(17, suzerain_bsmbsm_qinv(S, n, 21));
    BOOST_CHECK_EQUAL(22, suzerain_bsmbsm_qinv(S, n, 22));
    BOOST_CHECK_EQUAL(27, suzerain_bsmbsm_qinv(S, n, 23));
    BOOST_CHECK_EQUAL(32, suzerain_bsmbsm_qinv(S, n, 24));
    BOOST_CHECK_EQUAL(37, suzerain_bsmbsm_qinv(S, n, 25));
    BOOST_CHECK_EQUAL(42, suzerain_bsmbsm_qinv(S, n, 26));
    BOOST_CHECK_EQUAL( 3, suzerain_bsmbsm_qinv(S, n, 27));
    BOOST_CHECK_EQUAL( 8, suzerain_bsmbsm_qinv(S, n, 28));
    BOOST_CHECK_EQUAL(13, suzerain_bsmbsm_qinv(S, n, 29));
    BOOST_CHECK_EQUAL(18, suzerain_bsmbsm_qinv(S, n, 30));
    BOOST_CHECK_EQUAL(23, suzerain_bsmbsm_qinv(S, n, 31));
    BOOST_CHECK_EQUAL(28, suzerain_bsmbsm_qinv(S, n, 32));
    BOOST_CHECK_EQUAL(33, suzerain_bsmbsm_qinv(S, n, 33));
    BOOST_CHECK_EQUAL(38, suzerain_bsmbsm_qinv(S, n, 34));
    BOOST_CHECK_EQUAL(43, suzerain_bsmbsm_qinv(S, n, 35));
    BOOST_CHECK_EQUAL( 4, suzerain_bsmbsm_qinv(S, n, 36));
    BOOST_CHECK_EQUAL( 9, suzerain_bsmbsm_qinv(S, n, 37));
    BOOST_CHECK_EQUAL(14, suzerain_bsmbsm_qinv(S, n, 38));
    BOOST_CHECK_EQUAL(19, suzerain_bsmbsm_qinv(S, n, 39));
    BOOST_CHECK_EQUAL(24, suzerain_bsmbsm_qinv(S, n, 40));
    BOOST_CHECK_EQUAL(29, suzerain_bsmbsm_qinv(S, n, 41));
    BOOST_CHECK_EQUAL(34, suzerain_bsmbsm_qinv(S, n, 42));
    BOOST_CHECK_EQUAL(39, suzerain_bsmbsm_qinv(S, n, 43));
    BOOST_CHECK_EQUAL(44, suzerain_bsmbsm_qinv(S, n, 44));
}

BOOST_AUTO_TEST_CASE( gsl_permutation_equivalence )
{
   const int S = 5, n = 9, N = S*n;
   int data[N];

   boost::scoped_ptr<gsl_permutation> p(suzerain_bsmbsm_permutation(S,n));
   BOOST_REQUIRE(p);

   for (int i = 0; i < N; ++i) data[i] = i;
   gsl_permute_int(p->data, data, 1, N);
   for (int i = 0; i < N; ++i) {
      BOOST_CHECK_EQUAL(data[i], suzerain_bsmbsm_q(S, n, i));
   }

   gsl_permute_int_inverse(p->data, data, 1, N);
   for (int i = 0; i < N; ++i) {
      BOOST_CHECK_EQUAL(data[i], i);
   }

   boost::scoped_ptr<gsl_permutation> pinv(gsl_permutation_alloc(N));
   BOOST_REQUIRE(pinv);
   gsl_permutation_inverse(pinv.get(), p.get());

   for (int i = 0; i < N; ++i) data[i] = i;
   gsl_permute_int(pinv->data, data, 1, N);
   for (int i = 0; i < N; ++i) {
      BOOST_CHECK_EQUAL(data[i], suzerain_bsmbsm_qinv(S, n, i));
   }
}

BOOST_AUTO_TEST_CASE( identity_relation )
{
   // Checks that qinv inverts q for a large operating S, n
   for (int S = 1; S < 5; ++S) {
      for (int n = 1; n < 10; ++n) {
         for (int i = 0; i < S*n; ++i) {
            BOOST_CHECK_EQUAL(i,
               suzerain_bsmbsm_qinv(S, n,
                  suzerain_bsmbsm_q(S, n, i)));
         }
      }
   }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( aPxbpy )

// Templated type encapsulating a particular problem and precision
template<typename Scalar>
struct Problem {
    char trans;
    int S;
    int n;
    Scalar alpha;
    int incx;
    Scalar beta;
    int incy;
};

// Helper outputting Problem<Scalar>
template< typename charT, typename traits, typename Scalar >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os,
        const Problem<Scalar> &p)
{
    os << '{'
       << "trans = "  << p.trans
       << ",S = "     << p.S
       << ",n = "     << p.n
       << ",alpha = " << p.alpha
       << ",incx = "  << p.incx
       << ",beta = "  << p.beta
       << ",incy = "  << p.incy
       << '}';
    return os;
}

// Precision-specific dispatch for floats
void aPxpby(const Problem<float> &p, const float *x, float *y)
{
   suzerain_bsmbsm_saPxpby(
         p.trans, p.S, p.n, p.alpha, x, p.incx, p.beta, y, p.incy);
}

// Precision-specific dispatch for doubles
void aPxpby(const Problem<double> &p, const double *x, double *y)
{
   suzerain_bsmbsm_daPxpby(
         p.trans, p.S, p.n, p.alpha, x, p.incx, p.beta, y, p.incy);
}

// Precision-specific dispatch for complex floats
void aPxpby(const Problem<std::complex<float> > &p,
            const std::complex<float> *x,
                  std::complex<float> *y)
{
   float alpha[2], beta[2];
   memcpy(alpha, &p.alpha, sizeof(alpha));
   memcpy(beta,  &p.beta,  sizeof(beta));
   suzerain_bsmbsm_caPxpby(p.trans, p.S, p.n,
                           alpha, (const float (*)[2]) x, p.incx,
                           beta,  (      float (*)[2]) y, p.incy);
}

// Precision-specific dispatch for complex doubles
void aPxpby(const Problem<std::complex<double> > &p,
            const std::complex<double> *x,
                  std::complex<double> *y)
{
   double alpha[2], beta[2];
   memcpy(alpha, &p.alpha, sizeof(alpha));
   memcpy(beta,  &p.beta,  sizeof(beta));
   suzerain_bsmbsm_zaPxpby(p.trans, p.S, p.n,
                           alpha, (const double (*)[2]) x, p.incx,
                           beta,  (      double (*)[2]) y, p.incy);
}

// Precision-specific dispatch for floats
void permute(const Problem<float> &p, float *x)
{
   boost::scoped_ptr<gsl_permutation> g(suzerain_bsmbsm_permutation(p.S,p.n));
   switch (toupper(p.trans)) {
      case 'N': gsl_permute_float        (g->data, x, p.incx, p.S*p.n);
                break;
      case 'T': gsl_permute_float_inverse(g->data, x, p.incx, p.S*p.n);
                break;
      default:  BOOST_FAIL("Unknown p.trans");
   }
}

// Precision-specific dispatch for doubles
void permute(const Problem<double> &p, double *x)
{
   boost::scoped_ptr<gsl_permutation> g(suzerain_bsmbsm_permutation(p.S,p.n));
   switch (toupper(p.trans)) {
      case 'N': gsl_permute        (g->data, x, p.incx, p.S*p.n);
                break;
      case 'T': gsl_permute_inverse(g->data, x, p.incx, p.S*p.n);
                break;
      default:  BOOST_FAIL("Unknown p.trans");
   }
}

// Precision-specific dispatch for complex floats
void permute(const Problem<std::complex<float> > &p,
             std::complex<float> *x)
{
   boost::scoped_ptr<gsl_permutation> g(suzerain_bsmbsm_permutation(p.S,p.n));
   switch (toupper(p.trans)) {
      case 'N': gsl_permute_complex_float        (g->data, (float *)x, p.incx, p.S*p.n);
                break;
      case 'T': gsl_permute_complex_float_inverse(g->data, (float *)x, p.incx, p.S*p.n);
                break;
      default:  BOOST_FAIL("Unknown p.trans");
   }
}

// Precision-specific dispatch for complex doubles
void permute(const Problem<std::complex<double> > &p,
             std::complex<double> *x)
{
   boost::scoped_ptr<gsl_permutation> g(suzerain_bsmbsm_permutation(p.S,p.n));
   switch (toupper(p.trans)) {
      case 'N': gsl_permute_complex        (g->data, (double *)x, p.incx, p.S*p.n);
                break;
      case 'T': gsl_permute_complex_inverse(g->data, (double *)x, p.incx, p.S*p.n);
                break;
      default:  BOOST_FAIL("Unknown p.trans");
   }
}

template< typename Scalar >
void test(const Problem<Scalar> &p)
{
   const int N = p.S*p.n;

   // Allocate working storage
   boost::scoped_array<Scalar> x(new Scalar[N*p.incx]);
   boost::scoped_array<Scalar> y(new Scalar[N*p.incy]);
   boost::scoped_array<Scalar> r(new Scalar[N*p.incy]);

   // Create test data
   std::fill(x.get(), x.get() + N*p.incx, 2*p.alpha+1);
   std::fill(y.get(), y.get() + N*p.incy,   p.alpha-7);
   std::accumulate(x.get(), x.get() + N*p.incx, 0);
   std::accumulate(y.get(), y.get() + N*p.incy, 0);
   std::fill(y.get(), y.get() + N*p.incy, r.get());

   // Compute the single-shot result using BSMBSM
   aPxpby(p, x.get(), r.get());

   // Compute same result by permuting followed by axpby
   permute(p, x.get());
   suzerain::blas::axpby(N, p.alpha, x.get(), p.incx,
                            p.beta(), y.get(), y.incy);

   // Do r (from aPxpby) and y (from permute/axpby) agree?
   namespace traits = suzerain::complex::traits;
   typedef typename traits::real<Scalar>::type real_type;
   real_type tol = std::numeric_limits<real_type>::epsilon();
   if (traits::is_complex<Scalar>::value) tol *= 4;

   check_close_collections(r.get(), r.get() + N*p.incy,
                           y.get(), y.get() + N*p.incy,
                           tol);
}

BOOST_AUTO_TEST_SUITE_END()
