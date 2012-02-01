#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/bspline.hpp>
#include <suzerain/rholut_imexop.h>
#include <suzerain/countof.h>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

// Check if the apply and pack operations are consistent in that
//     (M+\varphi{}L)^(-1) * (M+\varphi{}L) * I == I
// holds to within reasonable floating point error.
BOOST_AUTO_TEST_CASE( operator_consistency )
{
    typedef double real_t;
    typedef std::complex<double> complex_t;

    // Initialize discrete B-spline operators
    using suzerain::bspline;
    using suzerain::bsplineop;
    const real_t breakpts[] = { 0.0, 0.5, 1.5, 2.0, 2.5, 2.75, 3.0 };
    suzerain::bspline b(8, suzerain::bspline::from_breakpoints(),
                        SUZERAIN_COUNTOF(breakpts), breakpts);
    suzerain::bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    const int S = 5;      // Number of state scalars
    const int n = op.n(); // DOF per submatrix
    const int N = S*n;    // DOF for global operator

    // Initialize scenario parameters
    const double phi = M_SQRT2;
    const double km  = M_E;
    const double kn  = M_PI;

    suzerain_rholut_imexop_scenario s;
    s.Re    = 3000;
    s.Pr    = 0.7;
    s.Ma    = 1.5;
    s.alpha = 3;
    s.gamma = 1.4;

    // Initialize reference quantities with random, nonzero values
    suzerain_rholut_imexop_ref r;
    const std::size_t nrefs = sizeof(r)/sizeof(real_t *);
    boost::scoped_array<real_t> refs(new real_t[n*nrefs]);
    for (std::size_t i = 0; i < n*nrefs; ++i) {
        refs[i] = 1 + random() / RAND_MAX;
    }
    suzerain_rholut_imexop_refld ld;
    std::fill((int *)&ld, (int *)(&ld + 1), 0);

    // Allocate state storage and initialize B1 to eye(N)
    boost::scoped_array<complex_t> B1(new complex_t[N*N]);
    boost::scoped_array<complex_t> B2(new complex_t[N*N]);
    std::fill(B1.get(), B1.get() + N*N, 0);
    std::fill(B2.get(), B2.get() + N*N, 0);
    for (int i = 0; i < N; ++i) B1[i + i*N] = 1;

    // Accumulate (M + \varphi{} L) B1 into B2
    for (int j = 0; j < N; ++j) {
        const int jN = j*N;
        suzerain_rholut_imexop_apply(phi, km, kn, &s, &r, &ld, op.get(),
               &B1[0*n+jN], &B1[1*n+jN], &B1[2*n+jN], &B1[3*n+jN], &B1[4*n+jN],
            0, &B2[0*n+jN], &B2[1*n+jN], &B2[2*n+jN], &B2[3*n+jN], &B2[4*n+jN]);
    }

    // Compute B1 = P*B2
    for (int j = 0; j < N; ++j) {
        suzerain_bsmbsm_zaPxpby('N', S, n, 1.0, &B2[j*N], 1, 0.0, &B1[j*N], 1);
    }

    // Form PAP^T
    suzerain_bsmbsm A;
    boost::scoped_array<complex_t> buf(
        new complex_t[n*(op.max_kl() + 1 + op.max_ku())]);
    boost::scoped_array<complex_t> papt(
        new complex_t[n*(S*(op.max_kl() + op.max_ku() + 2) - 1)]);
    suzerain_rholut_imexop_packc(phi, km, kn, &s, &r, &ld, op.get(),
                                 0, 1, 2, 3, 4, buf.get(), &A, papt.get());

    // "Invert"
    // TODO

    // Compute B2 = P^T*B1
    for (int j = 0; j < N; ++j) {
        suzerain_bsmbsm_zaPxpby('T', S, n, 1.0, &B1[j*N], 1, 0.0, &B2[j*N], 1);
    }

    // TODO Check result
}
