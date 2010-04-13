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
