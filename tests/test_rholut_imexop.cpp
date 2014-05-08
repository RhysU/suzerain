/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012-2014 Rhys Ulerich
 * Copyright (C) 2012-2014 The PECOS Development Team
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

#include <suzerain/rholut_imexop.h>

#include <boost/test/parameterized_test.hpp>
#include <boost/test/unit_test.hpp>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al/blas_et_al.hpp>
#include <suzerain/bsmbsm.h>
#include <suzerain/bspline.hpp>
#include <suzerain/complex.hpp>
#include <suzerain/countof.h>
#include <suzerain/gbmatrix.h>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

typedef double real_t;
typedef std::complex<real_t> complex_t;

// Machine precision constant
static const real_t macheps = std::numeric_limits<real_t>::epsilon();

// Abuse what we know about suzerain_rholut_imexop_ref{,ld}.
static const size_t NREFS = sizeof(suzerain_rholut_imexop_ref)/sizeof(real_t*);
BOOST_STATIC_ASSERT(NREFS == sizeof(suzerain_rholut_imexop_refld)/sizeof(int));

// Abuse what we know about NRBC matrices supplied to suzerain_rholut_XXX
static const size_t ASIZE = 5 * 5;
static const size_t BSIZE = 5 * 5;
static const size_t CSIZE = 5 * 5;

struct parameters
{
    real_t km;
    real_t kn;
    int    refndx;
    int    andx;   // Negative -1 to disable NRBC "a" matrix
                   // Zero to provide a "a" matrix full of zeros
                   // Positive to set nonzero (andx - 1) entry
    int    bndx;   // As andx but for NRBC "b" matrix
    int    cndx;   // As andx but for NRBC "c" matrix
};

template< typename charT, typename traits >
std::basic_ostream<charT,traits>& operator<<(
        std::basic_ostream<charT,traits> &os, const parameters& p)
{
    return os << "{km="      << p.km
              << ", kn="     << p.kn
              << ", refndx=" << p.refndx
              << ", andx="   << p.andx
              << ", bndx="   << p.bndx
              << ", cndx="   << p.cndx
              << '}';
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
    suzerain::bspline b(6, suzerain::bspline::from_breakpoints(),
                        SUZERAIN_COUNTOF(breakpts), breakpts);
    suzerain::bsplineop op(b, 2, SUZERAIN_BSPLINEOP_COLLOCATION_GREVILLE);
    const int S = 5;      // Number of state scalars
    const int n = op.n(); // DOF per submatrix
    const int N = S*n;    // DOF for global operator

    // Initialize scenario parameters
    const complex_t phi(M_SQRT2/11, M_LOG2E/13);
    const real_t&   km = p.km;
    const real_t&   kn = p.kn;

    suzerain_rholut_imexop_scenario s;
    s.Re    = 3000;
    s.Pr    = 0.7;
    s.Ma    = 1.5;
    s.alpha = 3;
    s.gamma = 1.4;

    // Initialize nonzero reference quantities for p.refndx
    // Abuses what we know about the structure of ...imexop_ref{,ld}.
    suzerain_rholut_imexop_ref r;
    suzerain::scoped_array<real_t> refs(new real_t[n*NREFS]);
    for (int i = 0; i < (int) NREFS; ++i) {
        if (i == p.refndx) {
            for (int j = 0; j < n; ++j) {
                refs[i*n + j] = 1;
            }
        } else {
            for (int j = 0; j < n; ++j) {
                refs[i*n + j] = 0;
            }
        }
    }
    for (size_t i = 0; i < NREFS; ++i) {    // Establish refs
        ((real_t **)&r)[i] = &refs[i*n];
    }
    suzerain_rholut_imexop_refld ld;
    fill((int *)&ld, (int *)(&ld + 1), 1);  // Establish lds

    suzerain::scoped_array<real_t> a55;     // Establish NRBC "a55" matrix
    if (p.andx != -1) {
        a55.reset(new real_t[ASIZE]);
        fill(a55.get(), a55.get() + ASIZE, 0);
        if (p.andx > 0) a55[p.andx - 1] = 1;
    }

    suzerain::scoped_array<real_t> b55;     // Establish NRBC "b" matrix
    if (p.bndx != -1) {
        b55.reset(new real_t[BSIZE]);
        fill(b55.get(), b55.get() + BSIZE, 0);
        if (p.bndx > 0) b55[p.bndx - 1] = 1;
    }

    suzerain::scoped_array<real_t> c55;     // Establish NRBC "c55" matrix
    if (p.cndx != -1) {
        c55.reset(new real_t[CSIZE]);
        fill(c55.get(), c55.get() + CSIZE, 0);
        if (p.cndx > 0) c55[p.cndx - 1] = 1;
    }

    // Allocate state storage and initialize B1 to eye(N)
    suzerain::scoped_array<complex_t> B1(new complex_t[N*N]);
    suzerain::scoped_array<complex_t> B2(new complex_t[N*N]);
    fill(B1.get(), B1.get() + N*N, 0);
    fill(B2.get(), B2.get() + N*N, 0);
    for (int i = 0; i < N; ++i) B1[i + i*N] = 1;

    // Accumulate (M + \varphi{} L) B1 into B2
    for (int j = 0; j < N; ++j) {
        const int jN = j*N;
        suzerain_rholut_imexop_accumulate(
            phi, km, kn, &s, &r, &ld, op.get(),
               &B1[0*n+jN], &B1[1*n+jN], &B1[2*n+jN], &B1[3*n+jN], &B1[4*n+jN],
            0, &B2[0*n+jN], &B2[1*n+jN], &B2[2*n+jN], &B2[3*n+jN], &B2[4*n+jN],
            a55.get(), b55.get(), c55.get());
    }

    // Compute B1 = P*B2
    for (int j = 0; j < N; ++j) {
        suzerain_bsmbsm_zaPxpby('N', S, n, 1.0, &B2[j*N], 1, 0.0, &B1[j*N], 1);
    }

    // Form PAP^T (with size computations to check documentation is right)
    const int bufsize  = std::max(75, n*(op.max_kl() + 1 + op.max_ku()));
    const int paptsize = (S*n)*(S*(op.max_kl() + op.max_ku() + 2) - 1);
    const int lusize   = (S*n)*(S*(2*op.max_kl() + op.max_ku() + 3) - 2);
    suzerain_bsmbsm A = suzerain_bsmbsm_construct(
            S, n, op.max_kl(), op.max_ku());
    BOOST_REQUIRE_EQUAL(bufsize,  A.ld * A.n);
    BOOST_REQUIRE_EQUAL(paptsize, A.LD * A.N);
    BOOST_REQUIRE_EQUAL(lusize,   (A.LD + A.KL) * A.N);
    suzerain::scoped_array<complex_t> buf(new complex_t[bufsize]);
    suzerain::scoped_array<complex_t> papt(new complex_t[paptsize]);

    // Fill all working storage with NaNs, invoke suzerain_rholut_imexop_packc,
    // and be sure we get a matrix lacking NaNs on the band as a result.
    using suzerain::complex::NaN;
    fill(buf.get(),  buf.get()  + bufsize,  NaN<complex_t>());
    fill(papt.get(), papt.get() + paptsize, NaN<complex_t>());
    suzerain_rholut_imexop_packc(phi, km, kn, &s, &r, &ld, op.get(),
                                 0, 1, 2, 3, 4, buf.get(), &A, papt.get(),
                                 a55.get(), b55.get(), c55.get());
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

    // Check that packf assembly produces the same result as packc assembly.
    // Else, switching from --solver=zgbsvx to --solver=zgbsv breaks badly.
    // Added due to Redmine #3073 where --solver=zgbsv with NRBC was NaNing.
    {
        suzerain_bsmbsm AF;
        const int packfsize = (S*n)*(S*(2*op.max_kl() + op.max_ku() + 3) - 2);
        suzerain::scoped_array<complex_t> packf(new complex_t[packfsize]);
        fill(packf.get(), packf.get() + packfsize, NaN<complex_t>());
        suzerain_rholut_imexop_packf(phi, km, kn, &s, &r, &ld, op.get(),
                                     0, 1, 2, 3, 4, buf.get(), &AF, packf.get(),
                                     a55.get(), b55.get(), c55.get());
        CHECK_GBMATRIX_CLOSE(
            S*n, S*n, A .KL, A .KU, &papt [   0],  A .LD,
            S*n, S*n, AF.KL, AF.KU, &packf[AF.KL], AF.LD + AF.KL,
            std::sqrt(std::numeric_limits<real_t>::epsilon()));
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
    BOOST_WARN_LT(1/rcond, 1/sqrt(macheps));         // Squawk on bad rcond
    BOOST_CHECK_LT(abs(B1[imaxabs]), macheps/rcond); // Check on estimate
    BOOST_TEST_MESSAGE("maximum absolute error is from value " << B1[imaxabs]
                       << " at index (" << (imaxabs % N)
                       << ", " << (imaxabs / N) << ")");

////// DEBUG
////std::cout << "B1=[";
////std::cout.setf(std::ios::showpos);
////for (int i = 0; i < N; ++i) {
////    for (int j = 0; j < N; ++j) {
////        std::cout << B1[i + j*N].real() << B1[i + j*N].imag() << "i";
////        if (j != N-1) std::cout << ",";
////    }
////    if (i != N-1) std::cout << ";\n";
////}
////std::cout << "]\n";
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
    // Due to excessive test runtime, the wavenumber sweep is avoided
    // unless we are testing the operator sans NRBC matrices a, b, and c.
    //
    // First register all-zero reference values to tickle degenerate cases.
    // Then register nonzero references one-by-one to ensure consistency.
    const real_t wavenumbers[][2] = {
         /*kx*/   /*kz*/
#ifdef TEST_RHOLUT_IMEXOP
        {     0,       0},
        {     0,  3*M_PI},
        { 7*M_E,       0},
        {     1,  3*M_PI},
        { 7*M_E,       1},
#endif
        { 7*M_E,  3*M_PI},
    };

    // TODO Add andx, bndx, cndx
    parameters p[] = {
          /*kx*/   /*kz*/  /*refndx*/  /*andx*/ /*bndx*/ /*cndx*/
#ifdef TEST_RHOLUT_IMEXOP
        {    -1,      -1,         -1,        -1,     -1,      -1},
#endif
#ifdef TEST_RHOLUT_IMEXOPA
        {    -1,      -1,         -1,         0,     -1,      -1},
        {    -1,      -1,         -1,         1,     -1,      -1},
        {    -1,      -1,         -1,         2,     -1,      -1},
        {    -1,      -1,         -1,         3,     -1,      -1},
        {    -1,      -1,         -1,         4,     -1,      -1},
        {    -1,      -1,         -1,         5,     -1,      -1},
        {    -1,      -1,         -1,         6,     -1,      -1},
        {    -1,      -1,         -1,         7,     -1,      -1},
        {    -1,      -1,         -1,         8,     -1,      -1},
        {    -1,      -1,         -1,         9,     -1,      -1},
        {    -1,      -1,         -1,        10,     -1,      -1},
        {    -1,      -1,         -1,        11,     -1,      -1},
        {    -1,      -1,         -1,        12,     -1,      -1},
        {    -1,      -1,         -1,        13,     -1,      -1},
        {    -1,      -1,         -1,        14,     -1,      -1},
        {    -1,      -1,         -1,        15,     -1,      -1},
        {    -1,      -1,         -1,        16,     -1,      -1},
        {    -1,      -1,         -1,        17,     -1,      -1},
        {    -1,      -1,         -1,        18,     -1,      -1},
        {    -1,      -1,         -1,        19,     -1,      -1},
        {    -1,      -1,         -1,        20,     -1,      -1},
        {    -1,      -1,         -1,        21,     -1,      -1},
        {    -1,      -1,         -1,        22,     -1,      -1},
        {    -1,      -1,         -1,        23,     -1,      -1},
        {    -1,      -1,         -1,        24,     -1,      -1},
        {    -1,      -1,         -1,        25,     -1,      -1},
#endif
#ifdef TEST_RHOLUT_IMEXOPB
        {    -1,      -1,         -1,        -1,      0,      -1},
        {    -1,      -1,         -1,        -1,      1,      -1},
        {    -1,      -1,         -1,        -1,      2,      -1},
        {    -1,      -1,         -1,        -1,      3,      -1},
        {    -1,      -1,         -1,        -1,      4,      -1},
        {    -1,      -1,         -1,        -1,      5,      -1},
        {    -1,      -1,         -1,        -1,      6,      -1},
        {    -1,      -1,         -1,        -1,      7,      -1},
        {    -1,      -1,         -1,        -1,      8,      -1},
        {    -1,      -1,         -1,        -1,      9,      -1},
        {    -1,      -1,         -1,        -1,     10,      -1},
        {    -1,      -1,         -1,        -1,     11,      -1},
        {    -1,      -1,         -1,        -1,     12,      -1},
        {    -1,      -1,         -1,        -1,     13,      -1},
        {    -1,      -1,         -1,        -1,     14,      -1},
        {    -1,      -1,         -1,        -1,     15,      -1},
        {    -1,      -1,         -1,        -1,     16,      -1},
        {    -1,      -1,         -1,        -1,     17,      -1},
        {    -1,      -1,         -1,        -1,     18,      -1},
        {    -1,      -1,         -1,        -1,     19,      -1},
        {    -1,      -1,         -1,        -1,     20,      -1},
        {    -1,      -1,         -1,        -1,     21,      -1},
        {    -1,      -1,         -1,        -1,     22,      -1},
        {    -1,      -1,         -1,        -1,     23,      -1},
        {    -1,      -1,         -1,        -1,     24,      -1},
        {    -1,      -1,         -1,        -1,     25,      -1},
#endif
#ifdef TEST_RHOLUT_IMEXOPC
        {    -1,      -1,         -1,        -1,     -1,       0},
        {    -1,      -1,         -1,        -1,     -1,       1},
        {    -1,      -1,         -1,        -1,     -1,       2},
        {    -1,      -1,         -1,        -1,     -1,       3},
        {    -1,      -1,         -1,        -1,     -1,       4},
        {    -1,      -1,         -1,        -1,     -1,       5},
        {    -1,      -1,         -1,        -1,     -1,       6},
        {    -1,      -1,         -1,        -1,     -1,       7},
        {    -1,      -1,         -1,        -1,     -1,       8},
        {    -1,      -1,         -1,        -1,     -1,       9},
        {    -1,      -1,         -1,        -1,     -1,      10},
        {    -1,      -1,         -1,        -1,     -1,      11},
        {    -1,      -1,         -1,        -1,     -1,      12},
        {    -1,      -1,         -1,        -1,     -1,      13},
        {    -1,      -1,         -1,        -1,     -1,      14},
        {    -1,      -1,         -1,        -1,     -1,      15},
        {    -1,      -1,         -1,        -1,     -1,      16},
        {    -1,      -1,         -1,        -1,     -1,      17},
        {    -1,      -1,         -1,        -1,     -1,      18},
        {    -1,      -1,         -1,        -1,     -1,      19},
        {    -1,      -1,         -1,        -1,     -1,      20},
        {    -1,      -1,         -1,        -1,     -1,      21},
        {    -1,      -1,         -1,        -1,     -1,      22},
        {    -1,      -1,         -1,        -1,     -1,      23},
        {    -1,      -1,         -1,        -1,     -1,      24},
        {    -1,      -1,         -1,        -1,     -1,      25},
#endif
    };
    for (size_t i = 0; i < sizeof(wavenumbers)/sizeof(wavenumbers[0]); ++i) {
        for (int j = -1; j < (int) NREFS; ++j) {
            for (size_t k = 0; k < sizeof(p)/sizeof(p[0]); ++k) {
                p[k].km = wavenumbers[i][0];
                p[k].kn = wavenumbers[i][1];
                p[k].refndx = j;
                std::ostringstream name;
                name << BOOST_TEST_STRINGIZE(operator_consistency)
                    << ' ' << p[k];
                master_test_suite().add(boost::unit_test::make_test_case(
                        &operator_consistency, name.str(), p + k, p + k + 1));
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

int main(int argc, char* argv[])
{
    return boost::unit_test::unit_test_main(&init_unit_test_suite, argc, argv);
}
