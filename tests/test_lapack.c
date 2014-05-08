//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#include <suzerain/blas_et_al/lapack.h>

#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_test.h>

#ifdef SUZERAIN_HAVE_MKL
#include <mkl.h>
#include <mkl_service.h>
#endif

static
void
test_lapack_dgbcon()
{
    /*
     * Test matrix is
     *      5  -3   0  0
     *      0   5  -3  0
     *      0   0   5 -3
     *      0   0   3 -1
     * which has a 1-norm of 11.  It has one subdiagonal but we'll use kl = 2
     * to test storage-padding details.  Extra kl = 2 rows of padding present
     * to store factorization multipliers.
     */

    const int kl = 2;
    const int ku = 1;
    const int n  = 4;
    double ab[] = { 0, 0,  0,  5, 0, 0,
                    0, 0, -3,  5, 0, 0,
                    0, 0, -3,  5, 3, 0,
                    0, 0, -3, -1, 0, 0 };
    int ipiv[n];
    const double norm1 = 11;

    // Prepare factorization
    const int f = suzerain_lapack_dgbtrf(n, n, kl, ku, ab, 2*kl + ku + 1, ipiv);
    gsl_test_int(f, 0, "%s:%d factorization success", __func__, __LINE__);

    // Compute reciprocal of condition number
    double rcond = -555;
    double work[3*n];
    int    iwork[n];
    const int g = suzerain_lapack_dgbcon('1', n, kl, ku, ab, 2*kl + ku + 1,
                                         ipiv, norm1, &rcond, work, iwork);
    gsl_test_int(g, 0, "%s:%d condition number estimation success",
                 __func__, __LINE__);

    // Check result against expected
    gsl_test_abs(rcond, 25.0/748.0, GSL_DBL_EPSILON,
            "%s:%d condition number estimation result %d",
            __func__, __LINE__, rcond);
}

int
main(int argc, char **argv)
{
    SUZERAIN_UNUSED(argc);
    SUZERAIN_UNUSED(argv);

    gsl_ieee_env_setup();

    /* TODO Add test_lapack_{c,z}gbtr{f,s} */
    /* Already exercised to some extent in test_bsplineop */

    test_lapack_dgbcon();

    /* TODO Add test_lapack_{s,c,z}gbcon */
    /* zgbcon already exercised to some extent in test_bsplineop */

#ifdef SUZERAIN_HAVE_MKL
#if INTEL_MKL_VERSION < 110002
    MKL_FreeBuffers();
#else
    mkl_free_buffers();
#endif
#endif

    exit(gsl_test_summary());
}
