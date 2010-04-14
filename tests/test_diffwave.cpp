#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <suzerain/diffwave.h>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

#pragma warning(disable:383)

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
            const int N = i + 1, dN = sizeof(expected[0])/sizeof(expected[0][0]);
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
            const int N = i + 1, dN = sizeof(expected[0])/sizeof(expected[0][0]);
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
    const int nelem = Ny*(dNx/2+1)*dNz;
    double (* const x)[2] = (double (*)[2]) malloc(nelem*sizeof(x[0]));
    double (* const y)[2] = (double (*)[2]) malloc(nelem*sizeof(y[0]));
    {
        double (*p)[2] = x, (*q)[2] = y;
        for (int n = 0; n < dNz; ++n) {
            for (int m = 0; m < (dNx/2+1); ++m) {
                for (int l = 0; l < Ny; ++l) {
                    (*p)[0] =  (l+2)*(m+3)*(n+5);   // Fill x
                    (*p)[1] = -(l+2)*(m+3)*(n+5);
                    p++;
                    (*q)[0] =  (l+7)*(m+11)*(n+13); // Fill y
                    (*q)[1] = -(l+7)*(m+11)*(n+13);
                    q++;
                }
            }
        }
    }

    const gsl_complex alpha = gsl_complex_rect(2, 0);
    const gsl_complex beta  = gsl_complex_rect(3, 0);
    suzerain_diffwave_accumulate(dxcnt, dzcnt, alpha.dat, x, beta.dat, y,
            Lx, Lz, Ny, Nx, dNx, 0, (dNx/2+1), Nz, dNz, 0, dNz);

    const double close_enough = std::numeric_limits<double>::epsilon()*100;
    double (*q)[2] = y;
    for (int n = 0; n < dNz; ++n) {
        for (int m = 0; m < (dNx/2+1); ++m) {
            for (int l = 0; l < Ny; ++l) {
                const gsl_complex observed = gsl_complex_rect((*q)[0],(*q)[1]);

                const gsl_complex xsrc = gsl_complex_rect(
                            (l+2)*(m+3)*(n+5),   -(l+2)*(m+3)*(n+5));
                const gsl_complex ysrc = gsl_complex_rect(
                            (l+7)*(m+11)*(n+13), -(l+7)*(m+11)*(n+13));

                gsl_complex expected;

                if (   (n <= (Nz+1)/2 || n >= (dNz - (Nz-1)/2))
                    && (m <= (Nx+1)/2 || m >= (dNx - (Nx-1)/2))) {

                    expected = gsl_complex_add(gsl_complex_mul(alpha, xsrc),
                                               gsl_complex_mul(beta,  ysrc));
                } else {
                    expected = gsl_complex_mul(beta, ysrc);
                }

                const double diff = gsl_complex_abs(
                        gsl_complex_sub(expected, observed));
                BOOST_CHECK_SMALL(diff, close_enough);

                ++q;
            }
        }
    }

    free(x);
    free(y);
}

BOOST_AUTO_TEST_CASE( accumulate )
{
    boost::array<int,9> c[] = {
        /* dxcnt, dzcnt,  Lx,  Lz, Ny, Nx, dNx, Nz, dNz */
        {      0,     0, 555, 777,  3,  4,   4,  4,   4  },
        {      0,     0, 555, 777,  3,  8,   8,  8,   8  },
        {      0,     0, 555, 777,  3,  7,   7,  7,   7  },
        {      0,     0, 555, 777,  3,  8,  12, 16,  24  }
    };

    for (int i = 0; i < sizeof(c)/sizeof(c[0]); ++i) {
        BOOST_TEST_MESSAGE("Testing " << c[i]);
        test_accumulate_helper(
                c[i][0], c[i][1], c[i][2], c[i][3], c[i][4],
                c[i][5], c[i][6], c[i][7], c[i][8]);
    }
}
