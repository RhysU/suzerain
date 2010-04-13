#ifdef HAVE_CONFIG_H
#include <suzerain/config.h>
#endif
#include <suzerain/common.hpp>
#pragma hdrstop
#define BOOST_TEST_MODULE $Id$
#include <boost/test/included/unit_test.hpp>
#include <suzerain/diffwave.h>

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

BOOST_AUTO_TEST_CASE( y0x0z0 )
{
    const int Ny = 2, Nx = 4, Nz = 4;
    double * const x = (double *) malloc(2*Ny*(Nx/2+1)*Nz*sizeof(x[0]));
    double * const y = (double *) malloc(2*Ny*(Nx/2+1)*Nz*sizeof(y[0]));
    {
        double *p = x, *q = y;
        for (int n = 0; n < Nz; ++n) {
            for (int m = 0; m < (Nx/2+1); ++m) {
                for (int l = 0; l < Ny; ++l) {
                    *p++ =  (l+1)*(m+1)*(n+1); // Fill x
                    *p++ = -(l+1)*(m+1)*(n+1);
                    *q++ =  (l+2)*(m+3)*(n+4); // Fill y
                    *q++ = -(l+2)*(m+3)*(n+4);
                }
            }
        }
    }

    suzerain_diffwave_accumulate_y0x0z0(
            2, x, 3, y, 555, 777, Ny, Nx, Nx, 0, (Nx/2+1), Nz, Nz, 0, Nz);

    {
        double *q = y;
        for (int n = 0; n < Nz; ++n) {
            for (int m = 0; m < (Nx/2+1); ++m) {
                for (int l = 0; l < Ny; ++l) {
                    BOOST_CHECK_EQUAL(
                            *q++,  2*(l+1)*(m+1)*(n+1) +  3*(l+2)*(m+3)*(n+4));
                    BOOST_CHECK_EQUAL(
                            *q++, -2*(l+1)*(m+1)*(n+1) + -3*(l+2)*(m+3)*(n+4));
                }
            }
        }
    }

    free(x);
    free(y);
}
