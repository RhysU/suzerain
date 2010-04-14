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

static void test_y0x0z0_helper(const int Ny,
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
                    (*p)[0] =  (l+1)*(m+1)*(n+1); // Fill x
                    (*p)[1] = -(l+1)*(m+1)*(n+1);
                    p++;
                    (*q)[0] =  (l+2)*(m+3)*(n+4); // Fill y
                    (*q)[1] = -(l+2)*(m+3)*(n+4);
                    q++;
                }
            }
        }
    }

    suzerain_diffwave_accumulate_y0x0z0(
            2, x, 3, y, 555, 777, Ny, Nx, dNx, 0, (dNx/2+1), Nz, dNz, 0, dNz);

    {
        double (*q)[2] = y;
        for (int n = 0; n < dNz; ++n) {
            for (int m = 0; m < (dNx/2+1); ++m) {
                for (int l = 0; l < Ny; ++l) {
                    if (   (n <= (Nz+1)/2 || n >= (dNz - (Nz-1)/2))
                        && (m <= (Nx+1)/2 || m >= (dNx - (Nx-1)/2))) {
                        BOOST_CHECK_EQUAL(
                            (*q)[0], 2*(l+1)*(m+1)*(n+1)+ 3*(l+2)*(m+3)*(n+4));
                        BOOST_CHECK_EQUAL(
                            (*q)[1],-2*(l+1)*(m+1)*(n+1)+-3*(l+2)*(m+3)*(n+4));
                        q++;
                    } else {
                        BOOST_CHECK_EQUAL((*q)[0], 0);
                        BOOST_CHECK_EQUAL((*q)[1], 0);
                        q++;
                    }
                }
            }
        }
    }

    free(x);
    free(y);
}

BOOST_AUTO_TEST_CASE( y0x0z0 )
{
    boost::array<int,5> c[] = {
        /* Ny, Nx, dNx, Nz, dNz */
        {   3,  8,   8,  8,   8  },
        {   3,  7,   7,  7,   7  },
        {   3,  8,  12, 16,  24  }
    };

    for (int i = 0; i < sizeof(c)/sizeof(c[0]); ++i) {
        BOOST_TEST_MESSAGE("Testing " << c[i]);
        test_y0x0z0_helper(c[i][0], c[i][1], c[i][2], c[i][3], c[i][4]);
    }
}
