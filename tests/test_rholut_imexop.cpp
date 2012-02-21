#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/bspline.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/countof.h>
#include <suzerain/gbmatrix.h>
#include <suzerain/rholut_imexop.h>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

typedef double real_t;
typedef std::complex<real_t> complex_t;

struct parameters
{
    real_t km;
    real_t kn;
    bool   zerorefs;
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const parameters& p)
{
    return os << "{km="        << p.km
              << ", kn="       << p.kn
              << ", zerorefs=" << p.zerorefs << '}';
}

// Free function checking if the apply and pack operations are
// consistent in that
//     (M+\varphi{}L)^(-1) * (M+\varphi{}L) * I == I
// holds to within reasonable floating point error for some parameters.
static void operator_consistency(const parameters& p)
{
    using suzerain::blas::iamax;
    using suzerain::blas::iamin;

    // Initialize discrete B-spline operators
    using suzerain::bspline;
    using suzerain::bsplineop;
    const real_t breakpts[] = { 0.0, 0.5, 1.5, 2.25, 2.5, 3.0, 3.35 };
    suzerain::bspline b(8, suzerain::bspline::from_breakpoints(),
                        SUZERAIN_COUNTOF(breakpts), breakpts);
    suzerain::bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    const int S = 5;      // Number of state scalars
    const int n = op.n(); // DOF per submatrix
    const int N = S*n;    // DOF for global operator

    // Initialize scenario parameters
    const complex_t phi(M_SQRT2/5, M_LOG2E);
    const real_t&   km = p.km;
    const real_t&   kn = p.kn;

    suzerain_rholut_imexop_scenario s;
    s.Re    = 3000;
    s.Pr    = 0.7;
    s.Ma    = 1.5;
    s.alpha = 3;
    s.gamma = 1.4;

    // Initialize reference quantities
    // Abuses what we know about the structure of ...imexop_ref{,ld}.
    suzerain_rholut_imexop_ref r;
    const std::size_t nrefs = sizeof(r)/sizeof(real_t *);
    boost::scoped_array<real_t> refs(new real_t[n*nrefs]);
    if (p.zerorefs) {  // Either use zero reference values...
        for (std::size_t i = 0; i < n*nrefs; ++i) {
            refs[i] = 0;
        }
    } else {           // ...or create interesting garbage
        for (std::size_t i = 0; i < n*nrefs; ++i) {
            refs[i] = 1 + (i / 100.0) + (random() / RAND_MAX);
        }
    }
    for (std::size_t i = 0; i < nrefs; ++i) {   // Establish refs
        ((real_t **)&r)[i] = &refs[i*n];
    }
    suzerain_rholut_imexop_refld ld;
    std::fill((int *)&ld, (int *)(&ld + 1), 1); // Establish lds

    // Allocate state storage and initialize B1 to eye(N)
    boost::scoped_array<complex_t> B1(new complex_t[N*N]);
    boost::scoped_array<complex_t> B2(new complex_t[N*N]);
    std::fill(B1.get(), B1.get() + N*N, 0);
    std::fill(B2.get(), B2.get() + N*N, 0);
    for (int i = 0; i < N; ++i) {
        B1[i + i*N] = 1;
    }

    // Accumulate (M + \varphi{} L) B1 into B2
    for (int j = 0; j < N; ++j) {
        const int jN = j*N;
        suzerain_rholut_imexop_accumulate(phi, km, kn, &s, &r, &ld, op.get(),
               &B1[0*n+jN], &B1[1*n+jN], &B1[2*n+jN], &B1[3*n+jN], &B1[4*n+jN],
            0, &B2[0*n+jN], &B2[1*n+jN], &B2[2*n+jN], &B2[3*n+jN], &B2[4*n+jN]);
    }

    // Compute B1 = P*B2
    for (int j = 0; j < N; ++j) {
        suzerain_bsmbsm_zaPxpby('N', S, n, 1.0, &B2[j*N], 1, 0.0, &B1[j*N], 1);
    }

    // Form PAP^T (with size computations to check documentation is right)
    const int bufsize  = n*(op.max_kl() + 1 + op.max_ku());
    const int paptsize = (S*n)*(S*(op.max_kl() + op.max_ku() + 2) - 1);
    const int lusize   = (S*n)*(S*(2*op.max_kl() + op.max_ku() + 3) - 2);
    suzerain_bsmbsm A = suzerain_bsmbsm_construct(
            S, n, op.max_kl(), op.max_ku());
    BOOST_REQUIRE_EQUAL(bufsize,  A.ld * A.n);
    BOOST_REQUIRE_EQUAL(paptsize, A.LD * A.N);
    BOOST_REQUIRE_EQUAL(lusize,   (A.LD + A.KL) * A.N);
    boost::scoped_array<complex_t> buf(new complex_t[bufsize]);
    boost::scoped_array<complex_t> papt(new complex_t[paptsize]);

    // Fill all working storage with NaNs, invoke suzerain_rholut_imexop_packc,
    // and be sure we get a matrix lacking NaNs on the band as a result.
    using suzerain::complex::NaN;
    std::fill(buf.get(),  buf.get()  + bufsize,  NaN<real_t>());
    std::fill(papt.get(), papt.get() + paptsize, NaN<real_t>());
    suzerain_rholut_imexop_packc(phi, km, kn, &s, &r, &ld, op.get(),
                                 0, 1, 2, 3, 4, buf.get(), &A, papt.get());
    for (int i = 0; i < A.N; ++i) {
        const int qi = suzerain_bsmbsm_q(A.S, A.n, i);
        for (int j = 0; j < A.N; ++j) {
            const int qj = suzerain_bsmbsm_q(A.S, A.n, j);
            if (suzerain_gbmatrix_in_band(A.LD,A.KL,A.KU,i,j)) {
                int o = suzerain_gbmatrix_offset(A.LD,A.KL,A.KU,i,j);
                BOOST_CHECK_MESSAGE(papt[o] == papt[o],
                    "NaN PAP^T_{"<<i<<","<<j<<"} from "
                    <<"submatrix ("<<qi/A.n<<","<<qj/A.n<<") "
                    <<"element ("<<qi%A.n<<","<<qj%A.n<<")");
            }
        }
    }

    // Factor LU = PAP^T and solve (LU)^{-1} B2 = B1
    boost::scoped_array<complex_t> lu(new complex_t[(A.LD + A.KL)*A.N]);
    std::fill(lu.get(), lu.get() + lusize, NaN<real_t>());
    boost::scoped_array<int>       ipiv(new int[A.N]);
    char equed;
    boost::scoped_array<real_t>    scale_r(new real_t[A.N]);
    boost::scoped_array<real_t>    scale_c(new real_t[A.N]);
    real_t rcond;
    boost::scoped_array<real_t>    ferr(new real_t[N]);
    boost::scoped_array<real_t>    berr(new real_t[N]);
    boost::scoped_array<complex_t> work(new complex_t[2*A.N]);
    boost::scoped_array<real_t>    rwork(new real_t[A.N]);

    BOOST_REQUIRE_EQUAL(0, suzerain_lapack_zgbsvx(/* fact  */ 'E',
                                                  /* trans */ 'N',
                                                  /* n     */ A.N,
                                                  /* kl    */ A.KL,
                                                  /* ku    */ A.KU,
                                                  /* nrhs  */ N,
                                                  /* ab    */ papt.get(),
                                                  /* ldab  */ A.LD,
                                                  /* afb   */ lu.get(),
                                                  /* ldafb */ A.LD + A.KL,
                                                  /* ipiv  */ ipiv.get(),
                                                  /* equed */ &equed,
                                                  /* r     */ scale_r.get(),
                                                  /* c     */ scale_c.get(),
                                                  /* b     */ B1.get(),
                                                  /* ldb   */ N,
                                                  /* x     */ B2.get(),
                                                  /* ldx   */ N,
                                                  /* rcond */ &rcond,
                                                  /* ferr  */ ferr.get(),
                                                  /* berr  */ berr.get(),
                                                  /* work  */ work.get(),
                                                  /* rwork */ rwork.get()));
    BOOST_TEST_MESSAGE("suzerain_lapack_zgbsvx returned equed = " <<   equed);
    BOOST_TEST_MESSAGE("suzerain_lapack_zgbsvx returned cond  = " << 1/rcond);

    // Compute B1 = P^T*B2
    for (int j = 0; j < N; ++j) {
        suzerain_bsmbsm_zaPxpby('T', S, n, 1.0, &B2[j*N], 1, 0.0, &B1[j*N], 1);
    }

    // Expected result is the identity so subtract identity matrix from B1
    for (int i = 0; i < N; ++i) {
        B1[i + i*N] -= 1;
    }

    // Expected result is now a matrix containing only zeros.
    // Check that the maximum absolute deviation from zero is small.
    const int imaxabs = iamax(N*N, B1.get(), 1);
    BOOST_TEST_MESSAGE("maximum absolute error is " << B1[imaxabs]);
    BOOST_CHECK_LT(std::abs(B1[imaxabs]),
                   N*N*N*std::numeric_limits<real_t>::epsilon());
}

boost::unit_test::test_suite*
init_unit_test_suite( int argc, char* argv[] )
{

    (void) argc; // Unused
    (void) argv; // Unused

    using boost::unit_test::framework::master_test_suite;
    master_test_suite().p_name.value = __FILE__;

    // Use a mixture of non-zero and zero wavenumbers in each of X, Z
    // Done as zero-wavenumbers hit degenerate portions of BLAS-like calls.
    //
    // Use both zero and non-zero reference values to tickle
    // degenerate-but-expected behavior.
    const parameters p[] = { {7*M_E, 3*M_PI, true },
                             {7*M_E,      0, true },
                             {    0, 3*M_PI, true },
                             {    0,      0, true },
                             {7*M_E, 3*M_PI, false},
                             {7*M_E,      0, false},
                             {    0, 3*M_PI, false},
                             {    0,      0, false} };
    for (size_t i = 0; i < sizeof(p)/sizeof(p[0]); ++i) {
        std::ostringstream name;
        name << BOOST_TEST_STRINGIZE(operator_consistency) << ' ' << p[i];
        master_test_suite().add(boost::unit_test::make_test_case(
                &operator_consistency, name.str(), p + i, p + i + 1));
    }

    return 0;
}
