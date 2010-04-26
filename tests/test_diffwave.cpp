#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <suzerain/diffwave.h>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

BOOST_AUTO_TEST_CASE( freqindex_nodealiasing )
{
    const int expected[][8] = { {0},
                                {0,               0},
                                {0,1,            -1},
                                {0,1,    0,      -1},
                                {0,1,2,       -2,-1},
                                {0,1,2,  0,   -2,-1},
                                {0,1,2,3,  -3,-2,-1},
                                {0,1,2,3,0,-3,-2,-1} };

    for (int i = 0; i < sizeof(expected)/sizeof(expected[0]); ++i) {
        const int N = i + 1, dN = N;
        for (int j = 0; j < dN; ++j) {
            BOOST_CHECK_EQUAL(expected[i][j],
                              suzerain_diffwave_freqindex(N, dN, j));
        }
    }
}

BOOST_AUTO_TEST_CASE( freqindex_dealiasing )
{
    {
        const int expected[][12] = { {0,0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,1,0,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,0,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,2,0,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,0,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,3,0,0,0,0,0,-3,-2,-1},
                                     {0,1,2,3,0,0,0,0,0,-3,-2,-1} };

        for (int i = 0; i < sizeof(expected)/sizeof(expected[0]); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
                BOOST_CHECK_EQUAL(expected[i][j],
                                suzerain_diffwave_freqindex(N, dN, j));
            }
        }
    }

    {
        const int expected[][11] = { {0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,0,0,0,0,0,0,0, 0, 0, 0},
                                     {0,1,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,0,0,0,0,0,0, 0, 0,-1},
                                     {0,1,2,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,0,0,0,0,0, 0,-2,-1},
                                     {0,1,2,3,0,0,0,0,-3,-2,-1},
                                     {0,1,2,3,0,0,0,0,-3,-2,-1} };

        for (int i = 0; i < sizeof(expected)/sizeof(expected[0]); ++i) {
            const int N  = i + 1;
            const int dN = sizeof(expected[0])/sizeof(expected[0][0]);
            for (int j = 0; j < dN; ++j) {
                BOOST_CHECK_EQUAL(expected[i][j],
                                suzerain_diffwave_freqindex(N, dN, j));
            }
        }
    }
}

static void test_accumulate_helper(const int dxcnt, const int dzcnt,
                                   const double Lx, const double Lz,
                                   const int Ny,
                                   const int Nx, const int dNx,
                                   const int Nz, const int dNz)
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
    }

    // Call the function under test
    const gsl_complex alpha = gsl_complex_rect(2, 0);
    const gsl_complex beta  = gsl_complex_rect(3, 0);
    suzerain_diffwave_accumulate(dxcnt, dzcnt, alpha.dat, x, beta.dat, y,
            Lx, Lz, Ny, Nx, dNx, 0, (dNx/2+1), Nz, dNz, 0, dNz);

    // Ensure the results match the synthetic fields
    // Horribly inefficient test code
    const gsl_complex itwopioverLx = gsl_complex_rect(0, 2*M_PI/Lx);
    const gsl_complex itwopioverLz = gsl_complex_rect(0, 2*M_PI/Lz);

    // Empirical tolerance choice: maybe too small, maybe not.
    const double small_enough = 7*std::pow(10, -10 + (dxcnt+dzcnt)/2.5);
    BOOST_TEST_MESSAGE("Small enough is " << small_enough);

    double (*q)[2] = y;
    for (int n = 0; n < dNz; ++n) {
        const int nfreqidx = suzerain_diffwave_freqindex(Nz, dNz, n);
        for (int m = 0; m < (dNx/2+1); ++m) {
            const int mfreqidx = suzerain_diffwave_freqindex(Nx, dNx, m);
            for (int l = 0; l < Ny; ++l) {
                const gsl_complex observed = gsl_complex_rect((*q)[0],(*q)[1]);

                const gsl_complex xsrc = gsl_complex_rect(
                     (l+1+2)*(m+1+ 3)*(n+1+ 5),  -(l+1+ 7)*(m+1+11)*(n+1+13));
                const gsl_complex ysrc = gsl_complex_rect(
                     (l+1+17)*(m+1+19)*(n+1+23), -(l+1+29)*(m+1+31)*(n+1+37));

                gsl_complex expected;

                if (dxcnt != 0 || dzcnt != 0) {
                    if (   (n <= (Nz+1)/2 || n >= (dNz - (Nz-1)/2))
                        && (m <= (Nx+1)/2 || m >= (dNx - (Nx-1)/2))) {

                        const gsl_complex xfactor
                            = gsl_complex_mul_real(itwopioverLx, mfreqidx);
                        const gsl_complex xscale
                            = gsl_complex_pow_real(xfactor, dxcnt);

                        const gsl_complex zfactor
                            = gsl_complex_mul_real(itwopioverLz, nfreqidx);
                        const gsl_complex zscale
                            = gsl_complex_pow_real(zfactor, dzcnt);

                        const gsl_complex scale
                            = gsl_complex_mul(xscale, zscale);
                        const gsl_complex alpha_D_x
                            = gsl_complex_mul(gsl_complex_mul(scale, alpha),
                                              xsrc);

                        expected = gsl_complex_add(gsl_complex_mul(beta, ysrc),
                                                   alpha_D_x);

                    } else {
                        expected = gsl_complex_mul(beta, ysrc);
                    }
                } else {
                    // Special handling for the no derivative case
                    // since we do not nuke the zeroth and Nyquist frequencies.
                    if (   (n <= (Nz+1)/2 || n >= (dNz - (Nz-1)/2))
                        && (m <= (Nx+1)/2 || m >= (dNx - (Nx-1)/2))) {

                        expected = gsl_complex_add(
                                gsl_complex_mul(alpha, xsrc),
                                gsl_complex_mul(beta,  ysrc));
                    } else {
                        expected = gsl_complex_mul(beta, ysrc);
                    }
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

    boost::array<int,7> c[] = {
        /* Lx, Lz, Ny, Nx, dNx, Nz, dNz */
        {   5,  7,  3,  4,   4,  4,   4  },
        {   5,  7,  3,  8,   8,  8,   8  },
        {   5,  7,  3,  7,   7,  7,   7  },
        {   5,  7,  3,  8,  12, 16,  24  }
    };

    for (int dxcnt = 0; dxcnt <= MAX_DXCNT_INCLUSIVE; ++dxcnt) {
        for (int dzcnt = 0; dzcnt <= MAX_DZCNT_INCLUSIVE; ++dzcnt) {
            for (int k = 0; k < sizeof(c)/sizeof(c[0]); ++k) {
                BOOST_TEST_MESSAGE("Testing dxcnt = " << dxcnt
                                                    << ", dzcnt = " << dzcnt
                                                    << " for params " << c[k]);
                test_accumulate_helper(dxcnt, dzcnt, c[k][0], c[k][1],
                        c[k][2], c[k][3], c[k][4], c[k][5], c[k][6]);
            }
        }
    }
}
