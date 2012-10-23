/*--------------------------------------------------------------------------
 *
 * Copyright (C) 2012 The PECOS Development Team
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
#include <suzerain/common.hpp>
#pragma hdrstop
#include <suzerain/filterop.h>
#define BOOST_TEST_MAIN
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include "test_tools.hpp"

BOOST_GLOBAL_FIXTURE(BlasCleanupFixture);

struct FilteropFixture {
    suzerain_filterop_workspace *w;
    FilteropFixture()  : w(NULL) {}
    ~FilteropFixture() { suzerain_filterop_free(w); }
};

BOOST_FIXTURE_TEST_SUITE( basic_suite, FilteropFixture )

BOOST_AUTO_TEST_CASE( filterop_matrices )
{
    // 'suzerain_filterop_workspace *w' handled by FilteropFixture
    w = suzerain_filterop_alloc(9, SUZERAIN_FILTEROP_COOKCABOT2005,
                                /* default params */ NULL,
                                SUZERAIN_FILTEROP_BOUNDARY_IGNORE,
                                SUZERAIN_FILTEROP_BOUNDARY_IGNORE);

    // Ensure workspace looks like the allocation call we requested
    BOOST_REQUIRE( w );
    BOOST_CHECK_EQUAL( w->n,    9 );
    BOOST_CHECK_EQUAL( w->klat, 2 );
    BOOST_CHECK_EQUAL( w->kuat, 2 );
    BOOST_CHECK_EQUAL( w->ldat, 2*2 + 1 + 2 );
    BOOST_CHECK_EQUAL( w->klbt, 4 );
    BOOST_CHECK_EQUAL( w->kubt, 4 );
    BOOST_CHECK_EQUAL( w->ldbt, 4 + 1 + 4);
    BOOST_CHECK( w->A_T );
    BOOST_CHECK( w->ipiva );
    BOOST_CHECK( w->B_T );

    // Coefficients for CookCabot2005 filter 
    // with default value of alpha_1
    const double alpha_0 = 1.0;
    const double alpha_1 = 66624./100000.;
    const double alpha_2 = (  1.-     alpha_1)/  2.;
    const double a_0     = ( 58.+105.*alpha_1)/128.;
    const double a_1     = ( 14.+ 11.*alpha_1)/ 32.;
    const double a_2     = ( 18.- 11.*alpha_1)/ 64.;
    const double a_3     = (  2.-  3.*alpha_1)/ 32.;
    const double a_4     = (- 2.+  3.*alpha_1)/256.;


    // Ensure A_T looks as expected prior to factorization
    // Notice: -5555 marks storage outside of matrix which is ignored
    // Notice: Values look "transposed" because C reads right-to-left
    const double good_A_T[] = {
        //  ku2       ku1      diag       kl1       kl2
          -5555,    -5555,  alpha_0,  alpha_1,  alpha_2,
          -5555,  alpha_1,  alpha_0,  alpha_1,  alpha_2,
        alpha_2,  alpha_1,  alpha_0,  alpha_1,  alpha_2,
        alpha_2,  alpha_1,  alpha_0,  alpha_1,  alpha_2,
        alpha_2,  alpha_1,  alpha_0,  alpha_1,  alpha_2,
        alpha_2,  alpha_1,  alpha_0,  alpha_1,  alpha_2,
        alpha_2,  alpha_1,  alpha_0,  alpha_1,  alpha_2,
        alpha_2,  alpha_1,  alpha_0,  alpha_1,    -5555,
        alpha_2,  alpha_1,  alpha_0,    -5555,    -5555,
    };


    // See test_tools.hpp for the macro signature
    // Notice: Though A_T has two additional superdiagonals for factorization
    //         we do not require them for call to CHECK_GBMATRIX_CLOSE
    //         because "w->A_T + w->klat" was used instead of "w->A_T".
    CHECK_GBMATRIX_CLOSE(
           9,    9,       2,       2,         good_A_T,       5,
        w->n, w->n, w->klat, w->kuat, w->A_T + w->klat, w->ldat,
        std::numeric_limits<double>::epsilon());



    // Ensure B_T looks as expected prior to factorization
    // Notice: -5555 marks storage outside of matrix which is ignored
    // Notice: Values look "transposed" because C reads right-to-left
    const double good_B_T[] = {
      //   ku4     ku3     ku2     ku1    diag     kl1     kl2     kl3     kl4  
         -5555,  -5555,  -5555,  -5555,    a_0,    a_1,    a_2,    a_3,    a_4,
         -5555,  -5555,  -5555,    a_1,    a_0,    a_1,    a_2,    a_3,    a_4,
         -5555,  -5555,    a_2,    a_1,    a_0,    a_1,    a_2,    a_3,    a_4,
         -5555,    a_3,    a_2,    a_1,    a_0,    a_1,    a_2,    a_3,    a_4,
           a_4,    a_3,    a_2,    a_1,    a_0,    a_1,    a_2,    a_3,    a_4,
           a_4,    a_3,    a_2,    a_1,    a_0,    a_1,    a_2,    a_3,  -5555,
           a_4,    a_3,    a_2,    a_1,    a_0,    a_1,    a_2,  -5555,  -5555,
           a_4,    a_3,    a_2,    a_1,    a_0,    a_1,  -5555,  -5555,  -5555,
           a_4,    a_3,    a_2,    a_1,    a_0,  -5555,  -5555,  -5555,  -5555,
    };


    // See test_tools.hpp for the macro signature
    CHECK_GBMATRIX_CLOSE(
           9,    9,       4,       4, good_B_T,       9,
        w->n, w->n, w->klbt, w->kubt,   w->B_T, w->ldbt,
        std::numeric_limits<double>::epsilon());

    // 'suzerain_filterop_free(w)' handled by FilteropFixture
}



BOOST_AUTO_TEST_CASE( filterop_nofilter )
{
    // 'suzerain_filterop_workspace *w' handled by FilteropFixture
    w = suzerain_filterop_alloc(16, SUZERAIN_FILTEROP_COOKCABOT2005,
                                /* default params */ NULL,
                                SUZERAIN_FILTEROP_BOUNDARY_NOFILTER,
                                SUZERAIN_FILTEROP_BOUNDARY_NOFILTER);

    // Ensure workspace looks like the allocation call we requested
    BOOST_REQUIRE( w );
    BOOST_CHECK_EQUAL( w->n,    16 );
    BOOST_CHECK_EQUAL( w->klat,  2 );
    BOOST_CHECK_EQUAL( w->kuat,  2 );
    BOOST_CHECK_EQUAL( w->ldat,  2*2 + 1 + 2 );
    BOOST_CHECK_EQUAL( w->klbt,  4 );
    BOOST_CHECK_EQUAL( w->kubt,  4 );
    BOOST_CHECK_EQUAL( w->ldbt,  4 + 1 + 4);
    BOOST_CHECK( w->A_T );
    BOOST_CHECK( w->ipiva );
    BOOST_CHECK( w->B_T );

    // Factorize A_T
    suzerain_filterop_factorize(w);

    // Test with a 7th order Chebyshev T polynomial
    // Notice: The Cook and Cabot filter produces no effect on this polynomial
    // Declare test function 
    const double ChebyshevT7_test_function[] = { 
               -1.           ,  148730387./170859375.,   -85000619./170859375.,
           -76443./78125.    ,  -43434727./170859375.,        1511./     2187.,
            77111./78125.    ,   76924511./170859375.,   -76924511./170859375.,
           -77111./78125.    ,      -1511./     2187.,    43434727./170859375.,
            76443./78125.    ,   85000619./170859375.,  -148730387./170859375.,
                1.
    };


    // Reference right-hand side result
    // Notice: In this case it is the same as the test function
    const double ChebyshevT7_good_rhs[] = { 
           -1., 148730387./170859375., -85000619./170859375.,
           -76443./78125., -38882317124./106787109375., 
           116501768948./106787109375., 58011238244./35595703125.,
           404959172084./533935546875., -404959172084./533935546875.,
           -58011238244./35595703125., -116501768948./106787109375.,
           38882317124./106787109375., 76443./78125., 
           85000619./170859375., -148730387./170859375., 1.
    };


    // Reference result
    // Notice: In this case it is the same as the test function
    const double ChebyshevT7_good_filtered[] = { 
               -1.           ,  148730387./170859375.,   -85000619./170859375.,
           -76443./78125.    ,  -43434727./170859375.,        1511./     2187.,
            77111./78125.    ,   76924511./170859375.,   -76924511./170859375.,
           -77111./78125.    ,      -1511./     2187.,    43434727./170859375.,
            76443./78125.    ,   85000619./170859375.,  -148730387./170859375.,
                1.
    };


    // Declare and initialize results vector
    double ChebyshevT7_result[] = { 
           0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. 
    };

    // Compute right-hand side
    // Apply filter to test function
    suzerain_filterop_apply(
           1., ChebyshevT7_test_function, /* contiguous */ 1, 
           ChebyshevT7_result, 1, w);

    // Check result for rhs
    // Notice: Tolerance adjusted manually
    check_close_collections(
            ChebyshevT7_result  , ChebyshevT7_result   + w->n,
            ChebyshevT7_good_rhs, ChebyshevT7_good_rhs + w->n,
            std::numeric_limits<double>::epsilon()*100.);


    // Compute result using suzerain_filterop_solve
    suzerain_filterop_solve(/* one vector */ 1,
           ChebyshevT7_result, w->n, w);

    // Check result
    // Notice: Tolerance adjusted manually
    check_close_collections(
            ChebyshevT7_result,        ChebyshevT7_result + w->n,
            ChebyshevT7_good_filtered, ChebyshevT7_good_filtered + w->n,
            std::numeric_limits<double>::epsilon()*w->n*w->n*100.);


    // Test with an 8th order Chebyshev T polynomial
    // Notice: The Cook and Cabot filter produces no effect on this polynomial
    // Declare test function 
    const double ChebyshevT8_test_function[] = {
        1., -1304175487./2562890625., 2446513793./2562890625., 164833./390625., 
        -1888197247./2562890625., -5983./6561., -15647./390625., 
        2206433153./2562890625., 2206433153./2562890625., 
        -15647./390625., -5983./6561., -1888197247./2562890625., 
        164833./390625., 2446513793./2562890625., -1304175487./2562890625., 1.
    };

    // Reference result
    const double ChebyshevT8_good_filtered[] = {
        1., -1304175487./2562890625., 
        2446513793./2562890625., 
        164833./390625., 
        -14247886892786344169./19338889058996484375., 
        -17635186346974875113./19338889058996484375., 
        -86076619383575777./2148765450999609375., 
        16649146794802214167./19338889058996484375., 
        16649146794802214167./19338889058996484375., 
        -86076619383575777./2148765450999609375., 
        -17635186346974875113./19338889058996484375., 
        -14247886892786344169./19338889058996484375., 
        164833./390625.,
        2446513793./2562890625., 
        -1304175487./2562890625., 1.
    };

    // Declare and initialize results vector
    double ChebyshevT8_result[] = { 
           0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. 
    };

    // Apply filter to test function
    suzerain_filterop_filter(1., ChebyshevT8_test_function, /* contiguous */ 1, 
                             ChebyshevT8_result, w);

    // Check result
    // Notice: Tolerance adjusted manually
    check_close_collections(
            ChebyshevT8_result,        ChebyshevT8_result + w->n,
            ChebyshevT8_good_filtered, ChebyshevT8_good_filtered + w->n,
            std::numeric_limits<double>::epsilon()*w->n*w->n*1000.);

    // 'suzerain_filterop_free(w)' handled by FilteropFixture
}



BOOST_AUTO_TEST_SUITE_END()
