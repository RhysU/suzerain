/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012, 2013 Rhys Ulerich
 * Copyright (C) 2012, 2013 The PECOS Development Team
 * Please see http://pecos.ices.utexas.edu for more information on PECOS.
 *
 * This file is part of Suzerain.
 *
 * Suzerain is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Suzerain is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 */

#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif

#include <suzerain/reacting_imexop.h>

#include <boost/test/parameterized_test.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/bspline.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/countof.h>
#include <suzerain/gbmatrix.h>

#include "test_tools.hpp"

// FIXME Testing loops over unused wavenumbers km, kn.  Remove.

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

typedef double real_t;
typedef std::complex<real_t> complex_t;

// Machine precision constant
static const real_t macheps = std::numeric_limits<real_t>::epsilon();

// Abuse what we know about suzerain_reacting_imexop_ref{,ld}.
static const size_t NREFS = sizeof(suzerain_reacting_imexop_ref)/sizeof(real_t*);
BOOST_STATIC_ASSERT(NREFS == sizeof(suzerain_reacting_imexop_refld)/sizeof(int));

struct parameters
{
    real_t km;
    real_t kn;
    bool   imagzero;
    int    refndx;
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const parameters& p)
{
    return os << "{km="        << p.km
              << ", kn="       << p.kn
              << ", imagzero=" << p.imagzero
              << ", refndx="   << p.refndx << '}';
}

// Free function checking if the apply and pack operations are
// consistent in that
//     (M+\varphi{}L)^(-1) * (M+\varphi{}L) * I == I
// holds to within reasonable floating point error for some parameters.
static void operator_consistency(const parameters& p)
{
    using std::abs;
    using std::fill;
    using std::sqrt;
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
    const complex_t phi(M_SQRT2/11, M_LOG2E/13);
    const real_t&   km = p.km;  SUZERAIN_UNUSED(km);
    const real_t&   kn = p.kn;  SUZERAIN_UNUSED(kn);

    suzerain_reacting_imexop_scenario s;
    s.alpha = 3;

    // Initialize nonzero reference quantities for p.refndx
    // Abuses what we know about the structure of ...imexop_ref{,ld}.
    suzerain_reacting_imexop_ref r;
    suzerain::scoped_array<real_t> refs(new real_t[n*NREFS]);
    for (int i = 0; i < (int) NREFS; ++i) {
        if (i == p.refndx) {
            for (int j = 0; j < n; ++j) {
                refs[i*n + j] = 1 + (random() / RAND_MAX);
            }
        } else {
            for (int j = 0; j < n; ++j) {
                refs[i*n + j] = 0;
            }
        }
    }
    for (size_t i = 0; i < NREFS; ++i) {        // Establish refs
        ((real_t **)&r)[i] = &refs[i*n];
    }
    suzerain_reacting_imexop_refld ld;
    fill((int *)&ld, (int *)(&ld + 1), 1); // Establish lds

    // Allocate state storage and initialize B1 to eye(N)
    // along with a poison NaN entry whenever p.imagzero is true.
    suzerain::scoped_array<complex_t> B1(new complex_t[N*N]);
    suzerain::scoped_array<complex_t> B2(new complex_t[N*N]);
    fill(B1.get(), B1.get() + N*N, 0);
    fill(B2.get(), B2.get() + N*N, 0);
    for (int i = 0; i < N; ++i) {
        B1[i + i*N] = complex_t(1,   p.imagzero
                                   ? std::numeric_limits<real_t>::quiet_NaN()
                                   : 0);
    }

    // Accumulate (M + \varphi{} L) B1 into B2
    for (int j = 0; j < N; ++j) {
        const int jN = j*N;
        suzerain_reacting_flow_imexop_accumulate(
            phi, /* km, kn, */ &s, &r, &ld, op.get(), p.imagzero,
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
    suzerain::scoped_array<complex_t> buf(new complex_t[bufsize]);
    suzerain::scoped_array<complex_t> papt(new complex_t[paptsize]);

    // Fill all working storage with NaNs, invoke suzerain_reacting_imexop_packc,
    // and be sure we get a matrix lacking NaNs on the band as a result.
    using suzerain::complex::NaN;
    fill(buf.get(),  buf.get()  + bufsize,  NaN<real_t>());
    fill(papt.get(), papt.get() + paptsize, NaN<real_t>());
    suzerain_reacting_flow_imexop_packc(phi, /* km, kn, */ &s, &r, &ld, op.get(),
                                        0, 1, 2, 3, 4, buf.get(), &A, papt.get());
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

    // Factor LU = PAP^T and solve (LU)^{-1} B2 = B1
    suzerain::scoped_array<complex_t> lu(new complex_t[(A.LD + A.KL)*A.N]);
    fill(lu.get(), lu.get() + lusize, NaN<real_t>());
    suzerain::scoped_array<int>       ipiv(new int[A.N]);
    char equed;
    suzerain::scoped_array<real_t>    scale_r(new real_t[A.N]);
    suzerain::scoped_array<real_t>    scale_c(new real_t[A.N]);
    real_t rcond;
    suzerain::scoped_array<real_t>    ferr(new real_t[N]);
    suzerain::scoped_array<real_t>    berr(new real_t[N]);
    suzerain::scoped_array<complex_t> work(new complex_t[2*A.N]);
    suzerain::scoped_array<real_t>    rwork(new real_t[A.N]);

    BOOST_REQUIRE_EQUAL(0, suzerain_lapack_zgbsvx(/* fact  */ 'E',
                                                  /* trans */ 'T',
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
    // "Small" defined as machine epsilon modified by conditioning estimate.
    const int imaxabs = iamax(N*N, B1.get(), 1);
    BOOST_TEST_MESSAGE("maximum absolute error is " << B1[imaxabs]);
    BOOST_WARN_LT(1/rcond, 1/sqrt(macheps));         // Squack on bad rcond
    BOOST_CHECK_LT(abs(B1[imaxabs]), macheps/rcond); // Check on estimate
}


// Free function checking if the apply and pack operations are
// consistent in that
//     (M+\varphi{}L)^(-1) * (M+\varphi{}L) * I == I
// holds to within reasonable floating point error for some parameters.
static void species_operator_consistency(const parameters& p)
{
    using std::abs;
    using std::fill;
    using std::sqrt;
    using suzerain::blas::iamax;
    using suzerain::blas::iamin;

    // Initialize discrete B-spline operators
    using suzerain::bspline;
    using suzerain::bsplineop;
    const real_t breakpts[] = { 0.0, 0.5, 1.5, 2.25, 2.5, 3.0, 3.35 };
    suzerain::bspline b(8, suzerain::bspline::from_breakpoints(),
                        SUZERAIN_COUNTOF(breakpts), breakpts);
    suzerain::bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    const int S = 1;      // Number of state scalars
    const int n = op.n(); // DOF per submatrix
    const int N = S*n;    // DOF for global operator

    // Initialize scenario parameters
    const complex_t phi(M_SQRT2/11, M_LOG2E/13);
    const real_t&   km = p.km;  SUZERAIN_UNUSED(km);
    const real_t&   kn = p.kn;  SUZERAIN_UNUSED(kn);

    suzerain_reacting_imexop_scenario s;
    s.alpha = 3;

    // Initialize nonzero reference quantities for p.refndx
    // Abuses what we know about the structure of ...imexop_ref{,ld}.
    suzerain_reacting_imexop_ref r;
    suzerain::scoped_array<real_t> refs(new real_t[n*NREFS]);
    for (int i = 0; i < (int) NREFS; ++i) {
        if (i == p.refndx) {
            for (int j = 0; j < n; ++j) {
                refs[i*n + j] = 1 + (random() / RAND_MAX);
            }
        } else {
            for (int j = 0; j < n; ++j) {
                refs[i*n + j] = 0;
            }
        }
    }
    for (size_t i = 0; i < NREFS; ++i) {        // Establish refs
        ((real_t **)&r)[i] = &refs[i*n];
    }
    suzerain_reacting_imexop_refld ld;
    fill((int *)&ld, (int *)(&ld + 1), 1); // Establish lds

    // Allocate state storage and initialize B1 to eye(N)
    // along with a poison NaN entry whenever p.imagzero is true.
    suzerain::scoped_array<complex_t> B1(new complex_t[N*N]);
    suzerain::scoped_array<complex_t> B2(new complex_t[N*N]);
    fill(B1.get(), B1.get() + N*N, 0);
    fill(B2.get(), B2.get() + N*N, 0);
    for (int i = 0; i < N; ++i) {
        B1[i + i*N] = complex_t(1,   p.imagzero
                                   ? std::numeric_limits<real_t>::quiet_NaN()
                                   : 0);
    }

    // Accumulate (M + \varphi{} L) B1 into B2
    for (int j = 0; j < N; ++j) {
        const int jN = j*N;
        suzerain_reacting_species_imexop_accumulate(
            phi, /* km, kn, */ &s, &r, &ld, op.get(), p.imagzero,
               &B1[0*n+jN],
            0, &B2[0*n+jN]);
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
    suzerain::scoped_array<complex_t> buf(new complex_t[bufsize]);
    suzerain::scoped_array<complex_t> papt(new complex_t[paptsize]);

    // Fill all working storage with NaNs, invoke suzerain_reacting_imexop_packc,
    // and be sure we get a matrix lacking NaNs on the band as a result.
    using suzerain::complex::NaN;
    fill(buf.get(),  buf.get()  + bufsize,  NaN<real_t>());
    fill(papt.get(), papt.get() + paptsize, NaN<real_t>());
    suzerain_reacting_species_imexop_packc(phi, /* km, kn, */ &s, &r, &ld, op.get(),
                                           buf.get(), &A, papt.get());
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

    // Factor LU = PAP^T and solve (LU)^{-1} B2 = B1
    suzerain::scoped_array<complex_t> lu(new complex_t[(A.LD + A.KL)*A.N]);
    fill(lu.get(), lu.get() + lusize, NaN<real_t>());
    suzerain::scoped_array<int>       ipiv(new int[A.N]);
    char equed;
    suzerain::scoped_array<real_t>    scale_r(new real_t[A.N]);
    suzerain::scoped_array<real_t>    scale_c(new real_t[A.N]);
    real_t rcond;
    suzerain::scoped_array<real_t>    ferr(new real_t[N]);
    suzerain::scoped_array<real_t>    berr(new real_t[N]);
    suzerain::scoped_array<complex_t> work(new complex_t[2*A.N]);
    suzerain::scoped_array<real_t>    rwork(new real_t[A.N]);

    BOOST_REQUIRE_EQUAL(0, suzerain_lapack_zgbsvx(/* fact  */ 'E',
                                                  /* trans */ 'T',
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
    // "Small" defined as machine epsilon modified by conditioning estimate.
    const int imaxabs = iamax(N*N, B1.get(), 1);
    BOOST_TEST_MESSAGE("maximum absolute error is " << B1[imaxabs]);
    BOOST_WARN_LT(1/rcond, 1/sqrt(macheps));         // Squack on bad rcond
    BOOST_CHECK_LT(abs(B1[imaxabs]), macheps/rcond); // Check on estimate
}

#ifdef BOOST_TEST_ALTERNATIVE_INIT_API
bool init_unit_test_suite() {
#else
::boost::unit_test::test_suite* init_unit_test_suite( int, char* [] )   {
#endif

    using boost::unit_test::framework::master_test_suite;
    master_test_suite().p_name.value = __FILE__;

    // Use a mixture of non-zero and zero wavenumbers in each of X, Z
    // Done as zero-wavenumbers hit degenerate portions of BLAS-like calls.
    //
    // First register all-zero reference values to tickle degenerate cases.
    // Then register nonzero references one-by-one to ensure consistency.
    parameters p[] = { {    0,      0, /*imagzero*/0, /*refndx*/-1},
                       {    0, 3*M_PI, /*imagzero*/0, /*refndx*/-1},
                       {7*M_E,      0, /*imagzero*/0, /*refndx*/-1},
                       {    1, 3*M_PI, /*imagzero*/0, /*refndx*/-1},
                       {7*M_E,      1, /*imagzero*/0, /*refndx*/-1},
                       {7*M_E, 3*M_PI, /*imagzero*/0, /*refndx*/-1} };
    for (int i = -1; i < (int) NREFS; ++i) {
        for (size_t j = 0; j < sizeof(p)/sizeof(p[0]); ++j) {
            p[j].refndx = i;

            p[j].imagzero = 0;
            {
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(operator_consistency)
                     << ' ' << p[j];
                master_test_suite().add(boost::unit_test::make_test_case(
                        &operator_consistency, name.str(), p + j, p + j + 1));
            }

            {
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(species_operator_consistency)
                     << ' ' << p[j];
                master_test_suite().add(boost::unit_test::make_test_case(
                        &species_operator_consistency, name.str(), p + j, p + j + 1));
            }

            p[j].imagzero = 1;
            {
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(operator_consistency)
                     << ' ' << p[j];
                master_test_suite().add(boost::unit_test::make_test_case(
                        &operator_consistency, name.str(), p + j, p + j + 1));
            }

            {
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(species_operator_consistency)
                     << ' ' << p[j];
                master_test_suite().add(boost::unit_test::make_test_case(
                        &species_operator_consistency, name.str(), p + j, p + j + 1));
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
