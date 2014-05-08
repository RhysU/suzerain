/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2014 Rhys Ulerich
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

#include <suzerain/diffwave.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_pow_int.h>

#include <suzerain/common.hpp>
#include <suzerain/blas_et_al.h>
#include <suzerain/inorder.h>

#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);


BOOST_AUTO_TEST_SUITE( prereqs )

BOOST_AUTO_TEST_CASE( gsl_sf_pow_int_zero_to_zero )
{
    // We rely on 0.0^0 == 1.0 according to gsl_sf_pow_int
    // If it changes, then some rework is required in diffwave.c
    BOOST_CHECK_EQUAL(gsl_sf_pow_int(0.0, 0), 1.0);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( accumulate )

// Tests mainly that the implementation accesses memory in the way that we
// anticipate.  The implementation walks memory linearly with little branching
// and utilizes the BLAS.  Here we inefficiently check that it did what we
// expect.  Functional correctness in terms of representing known functions is
// handled in test_diffwave_p3dfft.
static void test_accumulate_helper(const int dxcnt, const int dzcnt,
                                   const double Lx, const double Lz,
                                   const int Ny,
                                   const int Nx, const int dNx,
                                   const int Nz, const int dNz,
                                   const double small_enough)
{
    // Allocate test arrays
    const int nelem = Ny*(dNx/2+1)*dNz;
    double (* const x)[2] = (double (*)[2]) malloc(nelem*sizeof(x[0]));
    double (* const y)[2] = (double (*)[2]) malloc(nelem*sizeof(y[0]));

    // Load up a synthetic field
    {
        double (*p)[2] = x, (*q)[2] = y;
        for (int n = 0; n < dNz; ++n) {
            for (int m = 0; m < (dNx/2+1); ++m) {
                for (int l = 0; l < Ny; ++l) {
                    (*p)[0] =  (l+1+ 2)*(m+1+ 3)*(n+1+ 5); // Fill x
                    (*p)[1] = -(l+1+ 7)*(m+1+11)*(n+1+13);
                    p++;
                    (*q)[0] =  (l+1+17)*(m+1+19)*(n+1+23); // Fill y
                    (*q)[1] = -(l+1+29)*(m+1+31)*(n+1+37);
                    q++;
                }
            }
        }
        BOOST_REQUIRE_EQUAL(p - x, nelem);  // Covered all of x?
        BOOST_REQUIRE_EQUAL(q - y, nelem);  // Covered all of y?
    }

    // Call the function under test
    const gsl_complex alpha = gsl_complex_rect(2, 0);
    const gsl_complex beta  = gsl_complex_rect(3, 0);
    suzerain::diffwave::accumulate(dxcnt, dzcnt, alpha.dat, x, beta.dat, y,
            Lx, Lz, Ny, Nx, dNx, 0, (dNx/2+1), Nz, dNz, 0, dNz);

    // Ensure the results match the synthetic fields
    // Horribly inefficient test code
    const gsl_complex itwopioverLx = gsl_complex_rect(0, 2*M_PI/Lx);
    const gsl_complex itwopioverLz = gsl_complex_rect(0, 2*M_PI/Lz);

    double (*q)[2] = y;
    for (int n = 0; n < dNz; ++n) {
        // Z-wavenumber-dependent scaling
        const int nfreqidx = suzerain_inorder_wavenumber_diff(Nz, dNz, n);
        const gsl_complex zfactor = gsl_complex_mul_real(
                itwopioverLz, nfreqidx);
        const gsl_complex zscale = (dzcnt == 0)
                ? gsl_complex_rect(1, 0)
                : gsl_complex_pow_real(zfactor, dzcnt);

        for (int m = 0; m < (dNx/2+1); ++m) {
            // X-wavenumber-dependent scaling
            const int mfreqidx = suzerain_inorder_wavenumber_diff(Nx, dNx, m);
            const gsl_complex xfactor = gsl_complex_mul_real(
                    itwopioverLx, mfreqidx);
            const gsl_complex xscale = (dxcnt == 0)
                    ? gsl_complex_rect(1, 0)
                    : gsl_complex_pow_real(xfactor, dxcnt);

            // Combine Z- and X-wavenumber-dependent scaling
            const gsl_complex xzscale = gsl_complex_mul(xscale, zscale);

            for (int l = 0; l < Ny; ++l) {
                const gsl_complex observed = gsl_complex_rect((*q)[0],(*q)[1]);

                gsl_complex xsrc = gsl_complex_rect(
                     (l+1+ 2)*(m+1+ 3)*(n+1+ 5), -(l+1+ 7)*(m+1+11)*(n+1+13));
                const gsl_complex ysrc = gsl_complex_rect(
                     (l+1+17)*(m+1+19)*(n+1+23), -(l+1+29)*(m+1+31)*(n+1+37));

                const gsl_complex alpha_D_x
                    = gsl_complex_mul(gsl_complex_mul(xzscale, alpha), xsrc);
                const gsl_complex beta_y
                    = gsl_complex_mul(beta, ysrc);

                const bool n_keeper = suzerain_inorder_wavenumber_abs(dNz, n)
                                   <= suzerain_inorder_wavenumber_absmin(Nz);
                const bool m_keeper = suzerain_inorder_wavenumber_abs(dNx, m)
                                   <= suzerain_inorder_wavenumber_absmin(Nx);

                gsl_complex expected;
                if (n_keeper && m_keeper) {
                    expected = gsl_complex_add(alpha_D_x, beta_y);
                } else {
                    expected = beta_y;
                }

                const double diff = gsl_complex_abs(
                        gsl_complex_sub(expected, observed));
                BOOST_CHECK_SMALL(diff, small_enough);

                ++q;
            }
        }
    }

    // Deallocate test arrays
    free(x);
    free(y);
}

BOOST_AUTO_TEST_CASE( accumulate )
{
    const int MAX_DXCNT_INCLUSIVE = 4;
    const int MAX_DZCNT_INCLUSIVE = 4;

    suzerain::array<int,7> c[] = {
        /* Lx, Lz, Ny, Nx, dNx, Nz, dNz */
        // Beat on the Z direction in quasi-1D cases
        {{   5,  7,  1,  1,   1,  6,   9  }}
       ,{{   5,  7,  1,  1,   1,  6,   8  }}
       ,{{   5,  7,  1,  1,   1,  5,   8  }}
       ,{{   5,  7,  1,  1,   1,  5,   9  }}
        // Beat on the X direction in quasi-1D cases
       ,{{   5,  7,  1,  6,   9,  1,   1  }}
       ,{{   5,  7,  1,  6,   8,  1,   1  }}
       ,{{   5,  7,  1,  5,   8,  1,   1  }}
       ,{{   5,  7,  1,  5,   9,  1,   1  }}
        // Beat on the X and Z directions in quasi-2D cases
       ,{{   5,  7,  1,  6,   9,  6,   9  }}
       ,{{   5,  7,  1,  6,   8,  6,   8  }}
       ,{{   5,  7,  1,  5,   8,  5,   8  }}
       ,{{   5,  7,  1,  5,   9,  5,   9  }}
        // Beat on everything in full 3D cases
       ,{{   5,  7,  3,  4,   4,  4,   4  }}
       ,{{   5,  7,  3,  8,   8,  8,   8  }}
       ,{{   5,  7,  3,  7,   7,  7,   7  }}
       ,{{   5,  7,  3,  8,  12, 16,  24  }}
    };

    for (int dxcnt = 0; dxcnt <= MAX_DXCNT_INCLUSIVE; ++dxcnt) {
        for (int dzcnt = 0; dzcnt <= MAX_DZCNT_INCLUSIVE; ++dzcnt) {
            for (int k = 0; k < (int) (sizeof(c)/sizeof(c[0])); ++k) {

                // Empirical tolerance choice: maybe too small, maybe not.
                const double small = 7*std::pow(10, -10 + (dxcnt+dzcnt)/2.5);
                BOOST_TEST_MESSAGE("Testing dxcnt = " << dxcnt
                                                    << ", dzcnt = " << dzcnt
                                                    << " for params " << c[k]
                                                    << " using tol " << small);

                test_accumulate_helper(dxcnt, dzcnt, c[k][0], c[k][1],
                        c[k][2], c[k][3], c[k][4], c[k][5], c[k][6], small);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( apply )

// Tests mainly that the implementation accesses memory in the way that we
// anticipate.  The implementation walks memory linearly with little branching
// and utilizes the BLAS.  Here we inefficiently check that it did what we
// expect.  Functional correctness in terms of representing known functions is
// handled in test_diffwave_p3dfft.
static void test_apply_helper(const int dxcnt, const int dzcnt,
                              const double alpha_re, const double alpha_im,
                              const double Lx, const double Lz,
                              const int Ny,
                              const int Nx, const int dNx,
                              const int Nz, const int dNz,
                              const double small_enough)
{
    // Allocate test array
    const int nelem = Ny*(dNx/2+1)*dNz;
    double (* const x)[2] = (double (*)[2]) malloc(nelem*sizeof(x[0]));

    // Load up a synthetic field
    {
        double (*p)[2] = x;
        for (int n = 0; n < dNz; ++n) {
            for (int m = 0; m < (dNx/2+1); ++m) {
                for (int l = 0; l < Ny; ++l) {
                    (*p)[0] =  (l+1+ 2)*(m+1+ 3)*(n+1+ 5); // Fill x
                    (*p)[1] = -(l+1+ 7)*(m+1+11)*(n+1+13);
                    p++;
                }
            }
        }
        BOOST_REQUIRE_EQUAL(p - x, nelem);  // Covered all of x?
    }

    // Call the function under test
    const gsl_complex alpha = gsl_complex_rect(alpha_re, alpha_im);
    suzerain::diffwave::apply(dxcnt, dzcnt, alpha.dat, x,
            Lx, Lz, Ny, Nx, dNx, 0, (dNx/2+1), Nz, dNz, 0, dNz);

    // Ensure the results match the synthetic fields
    // Horribly inefficient test code
    const gsl_complex itwopioverLx = gsl_complex_rect(0, 2*M_PI/Lx);
    const gsl_complex itwopioverLz = gsl_complex_rect(0, 2*M_PI/Lz);

    double (*p)[2] = x;
    for (int n = 0; n < dNz; ++n) {
        // Z-wavenumber-dependent scaling
        const int nfreqidx = suzerain_inorder_wavenumber_diff(Nz, dNz, n);
        const gsl_complex zfactor = gsl_complex_mul_real(
                itwopioverLz, nfreqidx);
        const gsl_complex zscale = (dzcnt == 0)
                ? gsl_complex_rect(1, 0)
                : gsl_complex_pow_real(zfactor, dzcnt);

        for (int m = 0; m < (dNx/2+1); ++m) {
            // X-wavenumber-dependent scaling
            const int mfreqidx = suzerain_inorder_wavenumber_diff(Nx, dNx, m);
            const gsl_complex xfactor = gsl_complex_mul_real(
                    itwopioverLx, mfreqidx);
            const gsl_complex xscale = (dxcnt == 0)
                    ? gsl_complex_rect(1, 0)
                    : gsl_complex_pow_real(xfactor, dxcnt);

            // Combine Z- and X-wavenumber-dependent scaling
            const gsl_complex xzscale = gsl_complex_mul(xscale, zscale);

            for (int l = 0; l < Ny; ++l) {
                const gsl_complex observed = gsl_complex_rect((*p)[0],(*p)[1]);

                gsl_complex xsrc = gsl_complex_rect(
                     (l+1+ 2)*(m+1+ 3)*(n+1+ 5), -(l+1+ 7)*(m+1+11)*(n+1+13));

                const bool n_keeper = suzerain_inorder_wavenumber_abs(dNz, n)
                                   <= suzerain_inorder_wavenumber_absmin(Nz);
                const bool m_keeper = suzerain_inorder_wavenumber_abs(dNx, m)
                                   <= suzerain_inorder_wavenumber_absmin(Nx);

                const gsl_complex alpha_D_x
                    = gsl_complex_mul(gsl_complex_mul(xzscale, alpha), xsrc);

                gsl_complex expected;
                if (n_keeper && m_keeper) {
                    expected = alpha_D_x;
                } else {
                    expected = gsl_complex_rect(0, 0);
                }

                const double diff = gsl_complex_abs(
                        gsl_complex_sub(expected, observed));
                BOOST_CHECK_SMALL(diff, small_enough);

                ++p;
            }
        }
    }

    // Check that apply is equivalent to accumulate to within small
    // tolerances (Prepare identical fields in x and y, apply to y,
    // accumulate x to negative y, ensure the difference is small).
    {
        double (*p)[2] = x;
        for (int n = 0; n < dNz; ++n) {
            for (int m = 0; m < (dNx/2+1); ++m) {
                for (int l = 0; l < Ny; ++l) {
                    (*p)[0] =  (l+1+ 2)*(m+1+ 3)*(n+1+ 5); // Fill x
                    (*p)[1] = -(l+1+ 7)*(m+1+11)*(n+1+13);
                    p++;
                }
            }
        }
        BOOST_REQUIRE_EQUAL(p - x, nelem);  // Covered all of x?
    }
    double (* const y)[2] = (double (*)[2]) malloc(nelem*sizeof(y[0]));
    memcpy(y, x, nelem*sizeof(y[0]));
    suzerain::diffwave::apply(dxcnt, dzcnt, alpha.dat, y,
            Lx, Lz, Ny, Nx, dNx, 0, (dNx/2+1), Nz, dNz, 0, dNz);
    suzerain::diffwave::accumulate(dxcnt, dzcnt, alpha.dat, x, -1.0, y,
            Lx, Lz, Ny, Nx, dNx, 0, (dNx/2+1), Nz, dNz, 0, dNz);
    const int idamax[2] = { suzerain_blas_idamax(nelem, &y[0][0], 2),
                            suzerain_blas_idamax(nelem, &y[0][1], 2) };
    const double damax[2] = { y[idamax[0]][0], y[idamax[1]][1] };
    BOOST_CHECK_SMALL(damax[0], small_enough);
    BOOST_CHECK_SMALL(damax[1], small_enough);

    // Deallocate test arrays
    free(x);
    free(y);
}

BOOST_AUTO_TEST_CASE( apply )
{
    const int MAX_DXCNT_INCLUSIVE = 4;
    const int MAX_DZCNT_INCLUSIVE = 4;

    suzerain::array<int,7> c[] = {
        /* Lx, Lz, Ny, Nx, dNx, Nz, dNz */
        // Beat on the Z direction in quasi-1D cases
        {{   5,  7,  1,  1,   1,  6,   9  }}
       ,{{   5,  7,  1,  1,   1,  6,   8  }}
       ,{{   5,  7,  1,  1,   1,  5,   8  }}
       ,{{   5,  7,  1,  1,   1,  5,   9  }}
        // Beat on the X direction in quasi-1D cases
       ,{{   5,  7,  1,  6,   9,  1,   1  }}
       ,{{   5,  7,  1,  6,   8,  1,   1  }}
       ,{{   5,  7,  1,  5,   8,  1,   1  }}
       ,{{   5,  7,  1,  5,   9,  1,   1  }}
        // Beat on the X and Z directions in quasi-2D cases
       ,{{   5,  7,  1,  6,   9,  6,   9  }}
       ,{{   5,  7,  1,  6,   8,  6,   8  }}
       ,{{   5,  7,  1,  5,   8,  5,   8  }}
       ,{{   5,  7,  1,  5,   9,  5,   9  }}
        // Beat on everything in full 3D cases
       ,{{   5,  7,  3,  4,   4,  4,   4  }}
       ,{{   5,  7,  3,  8,   8,  8,   8  }}
       ,{{   5,  7,  3,  7,   7,  7,   7  }}
       ,{{   5,  7,  3,  8,  12, 16,  24  }}
    };

    for (int dxcnt = 0; dxcnt <= MAX_DXCNT_INCLUSIVE; ++dxcnt) {
        for (int dzcnt = 0; dzcnt <= MAX_DZCNT_INCLUSIVE; ++dzcnt) {
            for (int k = 0; k < (int) (sizeof(c)/sizeof(c[0])); ++k) {

                // Empirical tolerance choice: maybe too small, maybe not.
                double small = 7*std::pow(10, -10 + (dxcnt+dzcnt)/2.5);

                // Tickles the general logic using "arbitrary" alpha
                BOOST_TEST_MESSAGE("Testing dxcnt = " << dxcnt
                                                    << ", dzcnt = " << dzcnt
                                                    << ", alpha = 2+1i"
                                                    << " for params " << c[k]
                                                    << " using tol " << small);
                test_apply_helper(dxcnt, dzcnt, 2, 1, c[k][0], c[k][1],
                        c[k][2], c[k][3], c[k][4], c[k][5], c[k][6], small);

                // Tickles possibly specialized logic for alpha = 1
                if (dxcnt == 0 && dzcnt == 0) {
                    // No derivatives, scaling => no precision loss allowed
                    small = std::numeric_limits<double>::epsilon();
                }
                BOOST_TEST_MESSAGE("Testing dxcnt = " << dxcnt
                                                    << ", dzcnt = " << dzcnt
                                                    << ", alpha = 1+0i"
                                                    << " for params " << c[k]
                                                    << " using tol " << small);
                test_apply_helper(dxcnt, dzcnt, 1, 0, c[k][0], c[k][1],
                        c[k][2], c[k][3], c[k][4], c[k][5], c[k][6], small);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
