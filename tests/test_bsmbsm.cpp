//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
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

#include <suzerain/bsmbsm.h>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_permute_vector_complex_double.h>
#include <gsl/gsl_permute_vector_complex_float.h>
#include <gsl/gsl_permute_vector_float.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_vector.h>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al/blas_et_al.hpp>
#include <suzerain/countof.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/traits.hpp>

#include "test_tools.hpp"

#pragma warning(disable:1418 1572)

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

    suzerain::shared_ptr<gsl_permutation> p(suzerain_bsmbsm_permutation(S,n),
                                            &gsl_permutation_free);
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

    suzerain::shared_ptr<gsl_permutation> pinv(gsl_permutation_alloc(N),
                                               &gsl_permutation_free);
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
class Problem {
public:

    Problem(char trans, int S, int n,
            const Scalar &alpha, int incx,
            const Scalar &beta, int incy)
       : trans(trans), S(S), n(n),
         alpha(alpha), incx(incx),
         beta(beta), incy(incy)
    {}

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
       << ", S = "     << p.S
       << ", n = "     << p.n
       << ", alpha = " << p.alpha
       << ", incx = "  << p.incx
       << ", beta = "  << p.beta
       << ", incy = "  << p.incy
       << '}';
    return os;
}

// Precision-specific dispatch for floats
static void aPxpby(const Problem<float> &p, const float *x, float *y)
{
    suzerain_bsmbsm_saPxpby(
        p.trans, p.S, p.n, p.alpha, x, p.incx, p.beta, y, p.incy);
}

// Precision-specific dispatch for doubles
static void aPxpby(const Problem<double> &p, const double *x, double *y)
{
    suzerain_bsmbsm_daPxpby(
        p.trans, p.S, p.n, p.alpha, x, p.incx, p.beta, y, p.incy);
}

// Precision-specific dispatch for complex floats
static void aPxpby(const Problem<std::complex<float> > &p,
                   const std::complex<float> *x,
                         std::complex<float> *y)
{
    suzerain_bsmbsm_caPxpby(p.trans, p.S, p.n, p.alpha, x, p.incx,
                                               p.beta,  y, p.incy);
}

// Precision-specific dispatch for complex doubles
static void aPxpby(const Problem<std::complex<double> > &p,
                   const std::complex<double> *x,
                         std::complex<double> *y)
{
    suzerain_bsmbsm_zaPxpby(p.trans, p.S, p.n, p.alpha, x, p.incx,
                                               p.beta,  y, p.incy);
}

// Precision-specific dispatch for floats
// GSL cannot permute a vector with negative strides.
static void permute(const Problem<float> &p, gsl_permutation *g, float *x)
{
    gsl_vector_float_view v
        = gsl_vector_float_view_array_with_stride(x, abs(p.incx), p.S*p.n);
    if (p.incx < 0) gsl_vector_float_reverse(&v.vector);
    switch (toupper(p.trans)) {
    case 'N': gsl_permute_vector_float        (g, &v.vector);
              break;
    case 'T': gsl_permute_vector_float_inverse(g, &v.vector);
              break;
    default:  BOOST_FAIL("Unknown p.trans");
    }
    if (p.incx < 0) gsl_vector_float_reverse(&v.vector);
}

// Precision-specific dispatch for floats
// GSL cannot permute a vector with negative strides
static void permute(const Problem<double> &p, gsl_permutation *g, double *x)
{
    gsl_vector_view v
        = gsl_vector_view_array_with_stride(x, abs(p.incx), p.S*p.n);
    if (p.incx < 0) gsl_vector_reverse(&v.vector);
    switch (toupper(p.trans)) {
    case 'N': gsl_permute_vector        (g, &v.vector);
              break;
    case 'T': gsl_permute_vector_inverse(g, &v.vector);
              break;
    default:  BOOST_FAIL("Unknown p.trans");
    }
    if (p.incx < 0) gsl_vector_reverse(&v.vector);
}

// Precision-specific dispatch for complex floats
// GSL cannot permute a vector with negative strides
static void permute(const Problem<std::complex<float> > &p,
                    gsl_permutation *g,
                    std::complex<float> *x)
{
    gsl_vector_complex_float_view v
        = gsl_vector_complex_float_view_array_with_stride(
                (float *) x, abs(p.incx), p.S*p.n);
    if (p.incx < 0) gsl_vector_complex_float_reverse(&v.vector);
    switch (toupper(p.trans)) {
    case 'N': gsl_permute_vector_complex_float        (g, &v.vector);
              break;
    case 'T': gsl_permute_vector_complex_float_inverse(g, &v.vector);
              break;
    default:  BOOST_FAIL("Unknown p.trans");
    }
    if (p.incx < 0) gsl_vector_complex_float_reverse(&v.vector);
}

// Precision-specific dispatch for complex doubles
// GSL cannot permute a vector with negative strides
static void permute(const Problem<std::complex<double> > &p,
                    gsl_permutation *g,
                    std::complex<double> *x)
{
    gsl_vector_complex_view v
        = gsl_vector_complex_view_array_with_stride(
                (double *) x, abs(p.incx), p.S*p.n);
    if (p.incx < 0) gsl_vector_complex_reverse(&v.vector);
    switch (toupper(p.trans)) {
    case 'N': gsl_permute_vector_complex        (g, &v.vector);
              break;
    case 'T': gsl_permute_vector_complex_inverse(g, &v.vector);
              break;
    default:  BOOST_FAIL("Unknown p.trans");
    }
    if (p.incx < 0) gsl_vector_complex_reverse(&v.vector);
}

template< typename Scalar >
bool test(const Problem<Scalar> &p)
{
    BOOST_TEST_MESSAGE("Testing problem " << p);

    typedef typename suzerain::traits::component<Scalar>::type component_type;
    using suzerain::complex::traits::is_complex;
    using std::abs;
    using std::copy;
    using std::fill;
    using std::partial_sum;
    using std::sqrt;

    const int N = p.S*p.n;

    // Allocate working storage
    suzerain::scoped_array<Scalar> x(new Scalar[N*abs(p.incx)]);
    suzerain::scoped_array<Scalar> y(new Scalar[N*abs(p.incy)]);
    suzerain::scoped_array<Scalar> r(new Scalar[N*abs(p.incy)]);

    // Synthesize test data
    fill(x.get(), x.get() + N*abs(p.incx), p.alpha+p.alpha+Scalar(1));
    if (is_complex<Scalar>::value) {
         fill(y.get(), y.get() + N*abs(p.incy),
              p.alpha-Scalar(7) + sqrt(Scalar(-1)));
    } else {
         fill(y.get(), y.get() + N*abs(p.incy),
              p.alpha-Scalar(7));
    }
    partial_sum(x.get(), x.get() + N*abs(p.incx), x.get());
    partial_sum(y.get(), y.get() + N*abs(p.incy), y.get());
    copy(y.get(), y.get() + N*abs(p.incy), r.get());

    // Compute the single-shot result using BSMBSM routines
    aPxpby(p, x.get(), r.get());

    // Compute same result by permuting followed by axpby
    suzerain::shared_ptr<gsl_permutation> g(suzerain_bsmbsm_permutation(p.S,p.n),
                                            &gsl_permutation_free);
    permute(p, g.get(), x.get());
    suzerain::blas::axpby(N, p.alpha, x.get(), p.incx,
                             p.beta,  y.get(), p.incy);

    // Cook a reasonable agreement tolerance
    component_type tol = std::numeric_limits<component_type>::epsilon();
    if (is_complex<Scalar>::value) tol *= 4;

    // Do r (from aPxpby) and y (from permute/axpby) agree?
    return check_close_collections(r.get(), r.get() + N*abs(p.incy),
                                   y.get(), y.get() + N*abs(p.incy),
                                   tol);
}

BOOST_AUTO_TEST_CASE( spot_checks )
{
    // Cases which have proven problematic to draft implementations
    BOOST_REQUIRE((test(Problem<float>('N', 2, 7, 1, -1, 0, 1))));
}

typedef boost::mpl::list<float, double> component_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( real_valued, Component, component_types )
{
    // Test conditions
    const char      trans[] = { 'N', 'T' };
    const int       S[]     = { 1, 2, 5 };
    const int       n[]     = { 1, 7, 10 };
    const Component alpha[] = { 1, 3, -2 };
    const int       incx[]  = { 1, 2, -1, -2 };
    const Component beta[]  = { 0, 1, -3 };
    const int       incy[]  = { 1, 3, -1, -3 };

    // Outer product of all test conditions
    for (size_t ii = 0; ii < SUZERAIN_COUNTOF(trans); ++ii)
    for (size_t ij = 0; ij < SUZERAIN_COUNTOF(S);     ++ij)
    for (size_t ik = 0; ik < SUZERAIN_COUNTOF(n);     ++ik)
    for (size_t il = 0; il < SUZERAIN_COUNTOF(alpha); ++il)
    for (size_t im = 0; im < SUZERAIN_COUNTOF(incx);  ++im)
    for (size_t in = 0; in < SUZERAIN_COUNTOF(beta);  ++in)
    for (size_t io = 0; io < SUZERAIN_COUNTOF(incy);  ++io)
    {
       Problem<Component> P(trans[ii],
                            S    [ij],
                            n    [ik],
                            alpha[il],
                            incx [im],
                            beta [in],
                            incy [io]);
       BOOST_REQUIRE_MESSAGE(test(P), "Stopping due to failure in " << P);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( complex_valued, Component, component_types )
{
    typedef typename std::complex<Component> C;

    // Test conditions
    const char trans[] = { 'N', 'T' };
    const int  S[]     = { 1, 2, 5 };
    const int  n[]     = { 1, 7, 10 };
    const C    alpha[] = { 1, 3, C(-2,2) };
    const int  incx[]  = { 1, 2, -1, -2 };
    const C    beta[]  = { 0, 1, -3, C(4,-2) };
    const int  incy[]  = { 1, 3, -1, -3 };

    // Outer product of all test conditions
    for (size_t ii = 0; ii < SUZERAIN_COUNTOF(trans); ++ii)
    for (size_t ij = 0; ij < SUZERAIN_COUNTOF(S);     ++ij)
    for (size_t ik = 0; ik < SUZERAIN_COUNTOF(n);     ++ik)
    for (size_t il = 0; il < SUZERAIN_COUNTOF(alpha); ++il)
    for (size_t im = 0; im < SUZERAIN_COUNTOF(incx);  ++im)
    for (size_t in = 0; in < SUZERAIN_COUNTOF(beta);  ++in)
    for (size_t io = 0; io < SUZERAIN_COUNTOF(incy);  ++io)
    {
       Problem<C> P(trans[ii],
                    S    [ij],
                    n    [ik],
                    alpha[il],
                    incx [im],
                    beta [in],
                    incy [io]);
       BOOST_REQUIRE_MESSAGE(test(P), "Stopping due to failure in " << P);
    }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( submatrix_packing )

// Precision-specific dispatch for floats
template<typename Scalar>
int packc(const suzerain_bsmbsm *A, int ihat, int jhat, Scalar alpha,
          const float *b, float *papt)
{
    return suzerain_bsmbsm_spackc(A, ihat, jhat, alpha, b, papt);
}

// Precision-specific dispatch for doubles
template<typename Scalar>
int packc(const suzerain_bsmbsm *A, int ihat, int jhat, Scalar alpha,
          const double *b, double *papt)
{
    return suzerain_bsmbsm_dpackc(A, ihat, jhat, alpha, b, papt);
}

// Precision-specific dispatch for complex floats
template<typename Scalar>
int packc(const suzerain_bsmbsm *A, int ihat, int jhat, Scalar alpha,
          const std::complex<float> *b, std::complex<float> *papt)
{
    return suzerain_bsmbsm_cpackc(A, ihat, jhat, alpha, b, papt);
}

// Precision-specific dispatch for complex doubles
template<typename Scalar>
int packc(const suzerain_bsmbsm *A, int ihat, int jhat, Scalar alpha,
          const std::complex<double> *b, std::complex<double> *papt)
{
    return suzerain_bsmbsm_zpackc(A, ihat, jhat, alpha, b, papt);
}

typedef boost::mpl::list<
        float,
        double,
        std::complex<float>,
        std::complex<double>
    > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( degenerate, Scalar, test_types )
{
    static const Scalar zero = 0;
    static const Scalar one  = 1;
    static const Scalar two  = 2;
    typedef typename suzerain::traits::component<Scalar>::type component_type;
    using suzerain::complex::traits::is_complex;

    const suzerain_bsmbsm A = suzerain_bsmbsm_construct(1, 10, 2, 3);

    // Allocate source and target storage
    suzerain::scoped_array<Scalar> b(new Scalar[A.n*A.ld]);
    suzerain::scoped_array<Scalar> papt(new Scalar[A.N*A.LD]);

    // Generate well-defined source data
    if (suzerain::complex::traits::is_complex<Scalar>::value) {
        std::fill(b.get(), b.get() + A.n*A.ld, one - std::sqrt(-one));
    } else {
        std::fill(b.get(), b.get() + A.n*A.ld, one);
    }
    std::partial_sum(b.get(), b.get() + A.n*A.ld, b.get());

    // Fill target with NaNs and perform the degenerate pack operation
    // Check that the operation was indeed nothing but a copy
    BOOST_TEST_MESSAGE("Case alpha == one");
    std::fill(((component_type*) papt.get()),
              ((component_type*) papt.get())
                + A.N*A.LD*sizeof(Scalar)/sizeof(component_type),
              std::numeric_limits<component_type>::quiet_NaN());
    packc(&A, 0, 0, one, b.get(), papt.get());
    CHECK_GBMATRIX_CLOSE(A.n, A.n, A.kl, A.ku, b.get(),    A.ld,
                         A.N, A.N, A.KL, A.KU, papt.get(), A.LD,
                         std::numeric_limits<component_type>::epsilon());

    // Scale the source data before and after pack
    // Fill target with NaNs and perform the degenerate pack operation
    // Check that the operation was indeed nothing but a scaled copy
    BOOST_TEST_MESSAGE("Case alpha == two");
    std::fill(((component_type*) papt.get()),
              ((component_type*) papt.get())
                + A.N*A.LD*sizeof(Scalar)/sizeof(component_type),
              std::numeric_limits<component_type>::quiet_NaN());
    for (int i = 0; i < A.n*A.ld; ++i) b[i] /= two;
    packc(&A, 0, 0, two, b.get(), papt.get());
    for (int i = 0; i < A.n*A.ld; ++i) b[i] *= two;
    CHECK_GBMATRIX_CLOSE(A.n, A.n, A.kl, A.ku, b.get(),    A.ld,
                         A.N, A.N, A.KL, A.KU, papt.get(), A.LD,
                         std::numeric_limits<component_type>::epsilon());

    // Fill target with NaNs and perform the degenerate pack operation
    // Check that the operation was indeed nothing but a zero fill
    BOOST_TEST_MESSAGE("Case alpha == zero");
    std::fill(((component_type*) papt.get()),
              ((component_type*) papt.get())
                + A.N*A.LD*sizeof(Scalar)/sizeof(component_type),
              std::numeric_limits<component_type>::quiet_NaN());
    packc(&A, 0, 0, zero, b.get(), papt.get());
    for (int i = 0; i < A.n*A.ld; ++i) b[i] = zero;
    CHECK_GBMATRIX_CLOSE(A.n, A.n, A.kl, A.ku, b.get(),    A.ld,
                         A.N, A.N, A.KL, A.KU, papt.get(), A.LD,
                         std::numeric_limits<component_type>::epsilon());

}

BOOST_AUTO_TEST_CASE_TEMPLATE( minimal_two_by_two, Scalar, test_types )
{
    static const Scalar one = 1;
    typedef typename suzerain::traits::component<Scalar>::type component_type;
    using suzerain::complex::traits::is_complex;

    // Prepare a 2x2 problem with 5x5 submatrices
    const suzerain_bsmbsm A = suzerain_bsmbsm_construct(2, 2, 0, 0);
    suzerain::scoped_array<Scalar> b(new Scalar[A.n*A.ld]);
    suzerain::scoped_array<Scalar> papt(new Scalar[A.N*A.LD]);

    // Generate well-defined source data
    if (suzerain::complex::traits::is_complex<Scalar>::value) {
        std::fill(b.get(), b.get() + A.n*A.ld, one - std::sqrt(-one));
    } else {
        std::fill(b.get(), b.get() + A.n*A.ld, one);
    }
    std::partial_sum(b.get(), b.get() + A.n*A.ld, b.get());

    // Fill target with NaNs for detecting untouched values
    std::fill(((component_type*) papt.get()),
              ((component_type*) papt.get())
                + A.N*A.LD*sizeof(Scalar)/sizeof(component_type),
              std::numeric_limits<component_type>::quiet_NaN());

    // Pack the 4 submatrices with distinct scaling coefficients
    packc(&A, 0, 0, 2, b.get(), papt.get());
    packc(&A, 0, 1, 3, b.get(), papt.get());
    packc(&A, 1, 0, 5, b.get(), papt.get());
    packc(&A, 1, 1, 7, b.get(), papt.get());

    // Test that no NaN values remain inside papt's bandwidth
    for (int i = 0; i < A.N; ++i) {
        const int qi = suzerain_bsmbsm_q(A.S, A.n, i);
        for (int j = 0; j < A.N; ++j) {
            const int qj = suzerain_bsmbsm_q(A.S, A.n, j);
            if (suzerain_gbmatrix_in_band(A.KL,A.KU,i,j)) {
                int o = suzerain_gbmatrix_offset(A.LD,A.KL,A.KU,i,j);
                BOOST_CHECK_MESSAGE(papt[o] == papt[o],
                    "NaN PAP^T_{"<<i<<","<<j<<"} from "
                    <<"submatrix ("<<qi/A.n<<","<<qj/A.n<<") "
                    <<"element ("<<qi%A.n<<","<<qj%A.n<<")");
            }
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( three_by_three, Scalar, test_types )
{
    static const Scalar one = 1;
    typedef typename suzerain::traits::component<Scalar>::type component_type;
    using suzerain::complex::traits::is_complex;

    // Prepare a 3x3 problem with 17x17 submatrices
    const suzerain_bsmbsm A = suzerain_bsmbsm_construct(3, 17, 4, 3);
    suzerain::scoped_array<Scalar> b(new Scalar[A.n*A.ld]);
    suzerain::scoped_array<Scalar> papt(new Scalar[A.N*A.LD]);

    // Generate well-defined source data
    if (suzerain::complex::traits::is_complex<Scalar>::value) {
        std::fill(b.get(), b.get() + A.n*A.ld, one - std::sqrt(-one));
    } else {
        std::fill(b.get(), b.get() + A.n*A.ld, one);
    }
    std::partial_sum(b.get(), b.get() + A.n*A.ld, b.get());

    // Fill target with NaNs for detecting untouched values
    std::fill(((component_type*) papt.get()),
              ((component_type*) papt.get())
                + A.N*A.LD*sizeof(Scalar)/sizeof(component_type),
              std::numeric_limits<component_type>::quiet_NaN());

    // Pack the 9 submatrices with distinct scaling coefficients
    packc(&A, 0, 0,  2, b.get(), papt.get());
    packc(&A, 0, 1,  3, b.get(), papt.get());
    packc(&A, 0, 2,  5, b.get(), papt.get());
    packc(&A, 1, 0,  7, b.get(), papt.get());
    packc(&A, 1, 1, 11, b.get(), papt.get());
    packc(&A, 1, 2, 13, b.get(), papt.get());
    packc(&A, 2, 0, 17, b.get(), papt.get());
    packc(&A, 2, 1, 19, b.get(), papt.get());
    packc(&A, 2, 2, 23, b.get(), papt.get());

    // Test that no NaN values remain inside papt's bandwidth
    for (int i = 0; i < A.N; ++i) {
        const int qi = suzerain_bsmbsm_q(A.S, A.n, i);
        for (int j = 0; j < A.N; ++j) {
            const int qj = suzerain_bsmbsm_q(A.S, A.n, j);
            if (suzerain_gbmatrix_in_band(A.KL,A.KU,i,j)) {
                int o = suzerain_gbmatrix_offset(A.LD,A.KL,A.KU,i,j);
                BOOST_CHECK_MESSAGE(papt[o] == papt[o],
                    "NaN PAP^T_{"<<i<<","<<j<<"} from "
                    <<"submatrix ("<<qi/A.n<<","<<qj/A.n<<") "
                    <<"element ("<<qi%A.n<<","<<qj%A.n<<")");
            }
        }
    }
}

// Employing helper functions from http://snipt.net/RhysU/tag/mathematica,
// this test problem comes from Mathematica like
//     b = {0, 1/4, 3/4, 1, 5/4, 2};
//     M  = CollocationMatrix[0, 6, b];
//     D1 = CollocationMatrix[1, 6, b];
//     D2 = CollocationMatrix[2, 6, b];
//     Z = ConstantArray[0, {10, 10}];
//     op = ArrayFlatten[{{M,D1/5,Z}, {D2/7,2M,D2/14}, {Z,D1/10,4M}}];
// where op has 2-norm of roughly 100.
//
// Submatrix band expressions are available in C-like storage using
//     Flatten[Transpose[GeneralBandStorage[M,  4, 4]]] // InputForm
//     Flatten[Transpose[GeneralBandStorage[D1, 4, 4]]] // InputForm
//     Flatten[Transpose[GeneralBandStorage[D2, 4, 4]]] // InputForm

static const double M[] = {
        0.0, 0.0, 0.0, 0.0, 1.0, 1024.0/3125.0, 1.0/3125.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 144496.0/253125.0, 80404.0/253125.0,
        16807.0/506250.0, 16.0/253125.0, 0.0, 0.0, 0.0, 0.0,
        310429.0/3240000.0, 24064.0/50625.0, 71362.0/253125.0,
        449693.0/16200000.0, 0.0, 0.0, 0.0, 0.0, 747409.0/135000000.0,
        384544.0/2109375.0, 7545889.0/16875000.0, 34903837.0/135000000.0,
        128.0/78125.0, 0.0, 0.0, 0.0, 2771.0/22500000.0, 16672.0/703125.0,
        1064753.0/4921875.0, 83241521.0/157500000.0, 4050373.0/17500000.0,
        2592.0/109375.0, 2401.0/1500000.0, 81.0/3500000.0, 1.0/1500000.0,
        32.0/46875.0, 48961.0/2296875.0, 12950477.0/73500000.0,
        62041527.0/122500000.0, 749088.0/3828125.0, 726131.0/22500000.0,
        165159.0/122500000.0, 0.0, 0.0, 81.0/1225000.0, 4096.0/459375.0,
        306864937.0/1225000000.0, 8953344.0/19140625.0,
        141590743.0/675000000.0, 32713209.0/1225000000.0, 0.0,
        0.0, 0.0, 0.0, 994703.0/100000000.0, 106584.0/390625.0,
        3415772899.0/8100000000.0, 18014271.0/100000000.0, 0.0, 0.0,
        0.0, 0.0, 1.0/800000.0, 124.0/3125.0, 56669767.0/194400000.0,
        371281.0/800000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/3125.0,
        32768.0/759375.0, 1024.0/3125.0, 1.0, 0.0, 0.0, 0.0, 0.0
    };

static const double D1[] = {
        0.0, 0.0, 0.0, 0.0, -20.0, -1024.0/125.0, -4.0/125.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 20.0, 47584.0/10125.0, -28796.0/10125.0,
        -4802.0/10125.0, -32.0/10125.0, 0.0, 0.0, 0.0, 0.0,
        102769.0/32400.0, 1264.0/2025.0, -15784.0/10125.0,
        -62779.0/162000.0, 0.0, 0.0, 0.0, 0.0, 419509.0/1350000.0,
        152944.0/84375.0, 60683.0/168750.0, -1809851.0/1350000.0,
        -128.0/3125.0, 0.0, 0.0, 0.0, 2171.0/225000.0, 11872.0/28125.0,
        278764.0/196875.0, 956717.0/1575000.0, -198087.0/175000.0,
        -864.0/4375.0, -343.0/15000.0, -27.0/35000.0, 1.0/15000.0,
        32.0/1875.0, 23468.0/91875.0, 248043.0/245000.0,
        -240533.0/1225000.0, -148896.0/153125.0, -74333.0/225000.0,
        -42453.0/1225000.0, 0.0, 0.0, 27.0/12250.0, 2048.0/18375.0,
        14747177.0/12250000.0, -372048.0/765625.0, -8305549.0/6750000.0,
        -5841303.0/12250000.0, 0.0, 0.0, 0.0, 0.0, 165263.0/1000000.0,
        18372.0/15625.0, -29829457.0/81000000.0, -1735857.0/1000000.0,
        0.0, 0.0, 0.0, 0.0, 1.0/8000.0, 176.0/375.0, 2746019.0/1944000.0,
        -11581.0/24000.0, -20.0/3.0, 0.0, 0.0, 0.0, 0.0, 4.0/375.0,
        16384.0/30375.0, 1024.0/375.0, 20.0/3.0, 0.0, 0.0, 0.0, 0.0
    };

static const double D2[] = {
        0.0, 0.0, 0.0, 0.0, 320.0, 4096.0/25.0, 64.0/25.0, 0.0, 0.0,
        0.0, 0.0, 0.0, -1280.0/3.0, -409856.0/2025.0, 34816.0/2025.0,
        10976.0/2025.0, 256.0/2025.0, 0.0, 0.0, 0.0, 320.0/3.0,
        11029.0/405.0, -11264.0/405.0, 2752.0/2025.0, 8237.0/2025.0,
        0.0, 0.0, 0.0, 0.0, 181609.0/16875.0, 45376.0/16875.0,
        -196792.0/16875.0, 3373.0/16875.0, 512.0/625.0, 0.0,
        0.0, 0.0, 3142.0/5625.0, 28288.0/5625.0, 100928.0/39375.0,
        -327182.0/39375.0, 13306.0/4375.0, 1152.0/875.0, 98.0/375.0,
        18.0/875.0, 2.0/375.0, 128.0/375.0, 41536.0/18375.0,
        51466.0/18375.0, -223586.0/30625.0, 64128.0/30625.0,
        12838.0/5625.0, 19902.0/30625.0, 0.0, 0.0, 72.0/1225.0,
        4096.0/3675.0, 194217.0/153125.0, -1027136.0/153125.0,
        113407.0/84375.0, 749401.0/153125.0, 0.0, 0.0, 0.0, 0.0,
        27023.0/12500.0, -9088.0/9375.0, -9796949.0/1012500.0,
        -3043.0/37500.0, 80.0/3.0, 0.0, 0.0, 0.0, 1.0/100.0, 896.0/225.0,
        9583.0/24300.0, -21319.0/900.0, -560.0/9.0, 0.0, 0.0, 0.0, 0.0,
        64.0/225.0, 32768.0/6075.0, 4096.0/225.0, 320.0/9.0, 0.0, 0.0,
        0.0, 0.0
    };

static const double BR[] = {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0,
        23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0
    };

static const double XR[] = { // Found with N[Solve[op, BR], 30]
        -17.1205296638145575157187077973,
        -18.2963281236408675358367770456,
        -17.7937495693148766965188650215,
        - 5.24873030297202382738973786214,
          7.29869863887962759089597440773,
         10.9351171600530258544716323088,
          6.83821683636469177866241781694,
          6.84399935549535219401212203209,
          7.68282008181871798408896236238,
          8.39869754941780737393079084220,
        -25.2918333858194511840992191907,
        -20.7617009698658118051695422414,
        - 4.80320627915998609921722468469,
         14.0898063526294883933158530571,
         18.2079079742776334930218292431,
         10.1670604181122878920891517858,
          5.77253066966409009521644168134,
          7.46243669806219462513064720412,
          8.47491284977710917644942981328,
          9.67588968771375364600133668163,
          2.98493379202318031053516152534,
          3.01485895262305588445744594234,
          3.12754283059044674270610469246,
          4.87104392757836263762851591752,
          6.55548510997155995215443480015,
          7.11292221620541167955200620204,
          6.70305515823495426680783148237,
          6.85543754777172484660916203578,
          7.06268704531070322713866827606,
          7.29983719367722592174134885528
    };

BOOST_AUTO_TEST_CASE( solve_real )
{
    namespace blas = suzerain::blas;
    const suzerain_bsmbsm A = suzerain_bsmbsm_construct(3, 10, 4, 4);

    // Sanity checks on compile-time information
    BOOST_REQUIRE_EQUAL(A.S,       3);
    BOOST_REQUIRE_EQUAL(A.n,      10);
    BOOST_REQUIRE_EQUAL(A.kl,      4);
    BOOST_REQUIRE_EQUAL(A.ku,      4);
    BOOST_REQUIRE_EQUAL(A.ld,  4+1+4);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(M ), (unsigned) A.ld*A.n);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(D1), (unsigned) A.ld*A.n);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(D2), (unsigned) A.ld*A.n);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(BR), (unsigned) A.N);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(XR), (unsigned) A.N);

    // Allocate working buffers for accumulating submatrices
    suzerain::scoped_array<double> b   (new double[A.n*A.ld]);
    suzerain::scoped_array<double> papt(new double[A.N*(A.LD+A.KL)]);
    std::fill(papt.get(), papt.get() + A.N*(A.LD+A.KL),
              std::numeric_limits<double>::quiet_NaN());

    // Build the packed operator "op" per Mathematica above
    // B^{0,0} is M
    blas::copy(SUZERAIN_COUNTOF(M), M, 1, b.get(), 1);
    suzerain_bsmbsm_dpackf(&A, 0, 0, 1.0, b.get(), papt.get());
    // B^{1,1} is 2*M
    suzerain_bsmbsm_dpackf(&A, 1, 1, 2.0, b.get(), papt.get());
    // B^{2,2} is 4*M
    suzerain_bsmbsm_dpackf(&A, 2, 2, 4.0, b.get(), papt.get());

    // B^{0,1} is (1/5)*D1
    blas::copy(SUZERAIN_COUNTOF(D1), D1, 1, b.get(), 1);
    blas::scal(SUZERAIN_COUNTOF(D1), 1.0/5.0, b.get(), 1);
    suzerain_bsmbsm_dpackf(&A, 0, 1, 1.0, b.get(), papt.get());
    // B^{2,1} is (1/10)*D1
    blas::scal(SUZERAIN_COUNTOF(D1), 1.0/2.0, b.get(), 1);
    suzerain_bsmbsm_dpackf(&A, 2, 1, 1.0, b.get(), papt.get());

    // B^{1,0} is (1/7)*D2
    blas::copy(SUZERAIN_COUNTOF(D2), D2, 1, b.get(), 1);
    blas::scal(SUZERAIN_COUNTOF(D2), 1.0/7.0, b.get(), 1);
    suzerain_bsmbsm_dpackf(&A, 1, 0, 1.0, b.get(), papt.get());
    // B^{1,2} is (1/14)*D2
    blas::scal(SUZERAIN_COUNTOF(D2), 1.0/2.0, b.get(), 1);
    suzerain_bsmbsm_dpackf(&A, 1, 2, 1.0, b.get(), papt.get());

    // Zero out submatrices B^{0,2}, B^{2,0} using pack
    // Strictly necessary since we filled papt with NaNs above
    suzerain_bsmbsm_dpackf(&A, 0, 2, 0.0, NULL, papt.get());
    suzerain_bsmbsm_dpackf(&A, 2, 0, 0.0, NULL, papt.get());

    // Reuse working buffer to permute right hand side for solve
    b.reset(new double[2*A.N]);
    suzerain_bsmbsm_daPxpby(
            'N', A.S, A.n, 1.0, BR, 1, 0, b.get(), 1);

    // Solve in place
    suzerain::scoped_array<int> ipiv(new int[A.N]);
    BOOST_REQUIRE_EQUAL(0, suzerain::lapack::gbsv(
        A.N, A.KL, A.KU, 1, papt.get(), A.LD + A.KL, ipiv.get(),
        b.get(), A.N));

    // Unpermute the right hand side out-of-place
    suzerain_bsmbsm_daPxpby(
            'T', A.S, A.n, 1.0, b.get(), 1, 0, b.get() + A.N, 1);

    // Do we match the expected solution?
    check_close_collections(XR, XR + A.N, b.get() + A.N, b.get() + 2*A.N,
                            std::numeric_limits<double>::epsilon()*1e4);
}

// Complex operator used for this test is i*op used for solve_real.
// This decouples the real- and imaginary equations but is useful
// for ensuring the data movement (what we care about) remains correct.
BOOST_AUTO_TEST_CASE( solve_complex )
{
    typedef double               real_t;
    typedef std::complex<real_t> complex_t;
    namespace blas = suzerain::blas;
    const suzerain_bsmbsm A = suzerain_bsmbsm_construct(3, 10, 4, 4);

    // Sanity checks on compile-time information
    BOOST_REQUIRE_EQUAL(A.S,       3);
    BOOST_REQUIRE_EQUAL(A.n,      10);
    BOOST_REQUIRE_EQUAL(A.kl,      4);
    BOOST_REQUIRE_EQUAL(A.ku,      4);
    BOOST_REQUIRE_EQUAL(A.ld,  4+1+4);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(M ), (unsigned) A.ld*A.n);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(D1), (unsigned) A.ld*A.n);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(D2), (unsigned) A.ld*A.n);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(BR), (unsigned) A.N);
    BOOST_REQUIRE_EQUAL(SUZERAIN_COUNTOF(XR), (unsigned) A.N);

    // Allocate working buffers for accumulating submatrices
    suzerain::scoped_array<real_t> b   (new real_t[2*A.n*A.ld]);
    suzerain::scoped_array<real_t> papt(new real_t[2*A.N*(A.LD+A.KL)]);
    std::fill(b.get(),    b.get()    + 2*A.n*A.ld,        0);
    std::fill(papt.get(), papt.get() + 2*A.N*(A.LD+A.KL),
              std::numeric_limits<real_t>::quiet_NaN());

    // Build the packed, complex operator "i*op" per Mathematica above.
    // Notice that we only manipulate the imaginary portions of b.get().
    // B^{0,0} is i*M
    blas::copy(SUZERAIN_COUNTOF(M), M, 1, b.get()+1, 2);
    suzerain_bsmbsm_zdpackf(&A, 0, 0, std::complex<double>(0,1), M,
                                                (complex_t *) papt.get());
    // B^{1,1} is i*2*M
    suzerain_bsmbsm_zpackf(&A, 1, 1, 2.0, (const complex_t *) b.get(),
                                          (      complex_t *) papt.get());
    // B^{2,2} is i*4*M
    suzerain_bsmbsm_zpackf(&A, 2, 2, 4.0, (const complex_t *) b.get(),
                                          (      complex_t *) papt.get());

    // B^{0,1} is i*(1/5)*D1
    blas::copy(SUZERAIN_COUNTOF(D1), D1, 1, b.get()+1, 2);
    blas::scal(SUZERAIN_COUNTOF(D1), 1.0/5.0, b.get()+1, 2);
    suzerain_bsmbsm_zpackf(&A, 0, 1, 1.0, (const complex_t *) b.get(),
                                          (      complex_t *) papt.get());
    // B^{2,1} is i*(1/10)*D1
    blas::scal(SUZERAIN_COUNTOF(D1), 1.0/2.0, b.get()+1, 2);
    suzerain_bsmbsm_zpackf(&A, 2, 1, 1.0, (const complex_t *) b.get(),
                                          (      complex_t *) papt.get());

    // B^{1,0} is i*(1/7)*D2
    blas::copy(SUZERAIN_COUNTOF(D2), D2, 1,   b.get()+1, 2);
    blas::scal(SUZERAIN_COUNTOF(D2), 1.0/7.0, b.get()+1, 2);
    suzerain_bsmbsm_zpackf(&A, 1, 0, 1.0, (const complex_t *) b.get(),
                                          (      complex_t *) papt.get());
    // B^{1,2} is i*(1/14)*D2
    blas::scal(SUZERAIN_COUNTOF(D2), 1.0/2.0, b.get()+1, 2);
    suzerain_bsmbsm_zpackf(&A, 1, 2, 1.0, (const complex_t *) b.get(),
                                          (      complex_t *) papt.get());

    // B^{2,0} is identically zero so b argument is not referenced
    suzerain_bsmbsm_zpackf(&A, 2, 0, 0.0, NULL, (complex_t *) papt.get());

    // B^{0,2} is identically zero so b argument is not referenced
    suzerain_bsmbsm_zpackf(&A, 0, 2, 0.0, NULL, (complex_t *) papt.get());

    // Reuse working buffer to permute right hand side for solve
    b.reset(new real_t[2*(2*A.N)]);
    suzerain_bsmbsm_daPxpby(
            'N', A.S, A.n, -1.0, BR, 1, 0, b.get(),   2); // Re(-BR)
    suzerain_bsmbsm_daPxpby(
            'N', A.S, A.n,  1.0, BR, 1, 0, b.get()+1, 2); // Im( BR)

    // Solve in place
    suzerain::scoped_array<int> ipiv(new int[A.N]);
    BOOST_REQUIRE_EQUAL(0, suzerain::lapack::gbsv(
        A.N, A.KL, A.KU, 1, (complex_t *) papt.get(), A.LD + A.KL, ipiv.get(),
        (complex_t *) b.get(), A.N));

    // Unpermute the right hand side out-of-place
    suzerain_bsmbsm_daPxpby(
            'T', A.S, A.n, 1.0, b.get(),   2, 0, b.get() + 2*A.N, 1); // Re(X)
    suzerain_bsmbsm_daPxpby(
            'T', A.S, A.n, 1.0, b.get()+1, 2, 0, b.get() + 3*A.N, 1); // Im(X)

    // Do we match the expected solution?
    check_close_collections(XR, XR + A.N, b.get() + 2*A.N, b.get() + 3*A.N,
                            std::numeric_limits<double>::epsilon()*1e4);
    check_close_collections(XR, XR + A.N, b.get() + 3*A.N, b.get() + 4*A.N,
                            std::numeric_limits<double>::epsilon()*1e4);
}

BOOST_AUTO_TEST_SUITE_END()
