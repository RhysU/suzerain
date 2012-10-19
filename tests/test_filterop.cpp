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
    const double a_0     = ( 58.-105.*alpha_1)/128.;
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
                                SUZERAIN_FILTEROP_BOUNDARY_IGNORE);

    // FIXME: Apply operator and check result


    // 'suzerain_filterop_free(w)' handled by FilteropFixture
}



BOOST_AUTO_TEST_SUITE_END()
