#define BOOST_TEST_MODULE $Id$

#include <suzerain/config.h>

#include <boost/format.hpp>
#include <boost/test/included/unit_test.hpp>
#include <suzerain/bspline_eval.h>

#include "test_tools.hpp"

BOOST_AUTO_TEST_CASE( test_simple_knot_search_piecewise_quadratic )
{
    const int    k   = 3;
    const double t[] = { 0, 0, 0, 1, 2, 3, 3, 3 };
    const int nt  = sizeof(t)/sizeof(t[0]);

    /* Mid-interval search results */
    BOOST_CHECK_EQUAL( 0, suzerain_bspline_series_span_search(k, nt, t, 0.5));
    BOOST_CHECK_EQUAL( 1, suzerain_bspline_series_span_search(k, nt, t, 1.5));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 2.5));

    /* Left endpoint of interval search results */
    BOOST_CHECK_EQUAL( 0, suzerain_bspline_series_span_search(k, nt, t, 0));
    BOOST_CHECK_EQUAL( 1, suzerain_bspline_series_span_search(k, nt, t, 1));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 2));

    /* Right endpoint of interval search results */
    BOOST_CHECK_EQUAL( 1, suzerain_bspline_series_span_search(k, nt, t, 1));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 2));
    /* Check for left <= x <= right continuity in the rightmost interval */
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 3));
}

BOOST_AUTO_TEST_CASE( test_simple_knot_search_piecewise_cubic )
{
    const int    k   = 4;
    const double t[] = { 0, 0, 0, 0, 1, 2, 3, 3, 3, 3 };
    const int nt  = sizeof(t)/sizeof(t[0]);

    /* Mid-interval search results */
    BOOST_CHECK_EQUAL( 0, suzerain_bspline_series_span_search(k, nt, t, 0.5));
    BOOST_CHECK_EQUAL( 1, suzerain_bspline_series_span_search(k, nt, t, 1.5));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 2.5));

    /* Left endpoint of interval search results */
    BOOST_CHECK_EQUAL( 0, suzerain_bspline_series_span_search(k, nt, t, 0));
    BOOST_CHECK_EQUAL( 1, suzerain_bspline_series_span_search(k, nt, t, 1));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 2));

    /* Right endpoint of interval search results */
    BOOST_CHECK_EQUAL( 1, suzerain_bspline_series_span_search(k, nt, t, 1));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 2));
    /* Check for left <= x <= right continuity in the rightmost interval */
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 3));
}

BOOST_AUTO_TEST_CASE( test_repeated_knot_search_piecewise_quadratic )
{
    const int    k   = 3;
    const double t[] = { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 };
    const int    nt  = sizeof(t)/sizeof(t[0]);

    /* Mid-interval search results */
    BOOST_CHECK_EQUAL( 0, suzerain_bspline_series_span_search(k, nt, t, 0.5));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 1.5));
    BOOST_CHECK_EQUAL( 4, suzerain_bspline_series_span_search(k, nt, t, 2.5));

    /* Left endpoint of interval search results */
    BOOST_CHECK_EQUAL( 0, suzerain_bspline_series_span_search(k, nt, t, 0));
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 1));
    BOOST_CHECK_EQUAL( 4, suzerain_bspline_series_span_search(k, nt, t, 2));

    /* Right endpoint of interval search results */
    BOOST_CHECK_EQUAL( 2, suzerain_bspline_series_span_search(k, nt, t, 1));
    BOOST_CHECK_EQUAL( 4, suzerain_bspline_series_span_search(k, nt, t, 2));
    /* Check for left <= x <= right continuity in the rightmost interval */
    BOOST_CHECK_EQUAL( 4, suzerain_bspline_series_span_search(k, nt, t, 3));
}

BOOST_AUTO_TEST_CASE( test_piecewise_constants )
{
    const int    k    = 1;
    const double t[]  = {0, 1, 2, 3};
    const int    nt   = sizeof(t)/sizeof(t[0]);
    const int    n    = (nt - 2*(k-1) /* nbreakpoints */) + k - 2;
    const double C[n] = {1, 2, 3};

    double Q[k];

    /* Mid-interval evaluation results */
    suzerain_bspline_series_evaluate(n, k, t, 0, 0.5, C, Q);
    BOOST_CHECK_EQUAL(1, Q[0]);
    suzerain_bspline_series_evaluate(n, k, t, 0, 1.5, C, Q);
    BOOST_CHECK_EQUAL(2, Q[0]);
    suzerain_bspline_series_evaluate(n, k, t, 0, 2.5, C, Q);
    BOOST_CHECK_EQUAL(3, Q[0]);

    /* Left endpoint of interval evaluation results */
    suzerain_bspline_series_evaluate(n, k, t, 0, 0.0, C, Q);
    BOOST_CHECK_EQUAL(1, Q[0]);
    suzerain_bspline_series_evaluate(n, k, t, 0, 1.0, C, Q);
    BOOST_CHECK_EQUAL(2, Q[0]);
    suzerain_bspline_series_evaluate(n, k, t, 0, 2.0, C, Q);
    BOOST_CHECK_EQUAL(3, Q[0]);

    /* Right endpoint of interval evaluation results */
    suzerain_bspline_series_evaluate(n, k, t, 0, 1.0, C, Q);
    BOOST_CHECK_EQUAL(2, Q[0]);
    suzerain_bspline_series_evaluate(n, k, t, 0, 2.0, C, Q);
    BOOST_CHECK_EQUAL(3, Q[0]);
    suzerain_bspline_series_evaluate(n, k, t, 0, 3.0, C, Q);
    BOOST_CHECK_EQUAL(3, Q[0]);
}

BOOST_AUTO_TEST_CASE( spot_check_derivative_of_piecewise_linears )
{
    const int    k    = 2;
    const double t[]  = {0, 0, 2, 7, 7};
    const int    nt   = sizeof(t)/sizeof(t[0]);
    const int    n    = (nt - 2*(k-1) /* nbreakpoints */) + k - 2;
    double       C[n];
    double       Q[k];

    const double xloc[4] = { 0.0, 1.0, 6.0, 7.0 };
    const double deriv[4][3] =
    {
      { -1.0/2.0,  1.0/2.0, 0.0     },
      { -1.0/2.0,  1.0/2.0, 0.0     },
      {      0.0, -1.0/5.0, 1.0/5.0 },
      {      0.0, -1.0/5.0, 1.0/5.0 }
    };

    for (int i = 0; i < sizeof(xloc)/sizeof(xloc[0]); ++i) {
        for (int j = 0; j < n; ++j) {
            // Put garbage in Q
            for (int l = 0; l < k; ++l) Q[l] = 789e123;

            // Set only C[j] = 1.0 to obtain a single basis function
            for (int l = 0; l < n; ++l) C[l] = 0.0;
            C[j] = 1.0;

            // Check derivative of basis function j at point xloc[i]
            suzerain_bspline_series_evaluate(n, k, t, k-1, xloc[i], C, Q);
            BOOST_CHECK_EQUAL(Q[1], deriv[i][j]);
        }
    }
}

BOOST_AUTO_TEST_CASE( spot_check_evaluation_and_derivatives_piecewise_quadratics )
{
    const int    k    = 3;
    const double t[]  = {-3, -3, -3, 2, 9, 12, 21, 21, 21};
    const int    nt   = sizeof(t)/sizeof(t[0]);
    const int    n    = (nt - 2*(k-1) /* nbreakpoints */) + k - 2;
    double       C[n];
    double       Q[k];

    const double xloc[5]     =  { 0.0, 5.0, 9.0, 12.0, 15.0 };
    const double eval[5][6]  =
    {
      { 4./25.,  69./100.,   3./ 20. ,  0.    , 0.   , 0.    },
      { 0.     ,  4./21. , 143./210. ,  9./70., 0.   , 0.    },
      { 0.     ,  0.     ,   3./ 10. ,  7./10., 0.   , 0.    },
      { 0.     ,  0.     ,   0.      ,  3./4. , 1./4., 0.    },
      { 0.     ,  0.     ,   0.      ,  1./3. , 5./9., 1./9. }
    };
    const double deriv[5][6] =
    {
      { -4./25.,  3./50.,   1./ 10.,  0.    , 0.    , 0.      },
      {  0.    , -2./21.,   1./105.,  3./35., 0.    , 0.      },
      {  0.    ,  0.    ,  -1./5.  ,  1./ 5., 0.    , 0.      },
      {  0.    ,  0.    ,   0.     , -1./ 6., 1./6. , 0.      },
      {  0.    ,  0.    ,   0.     , -1./ 9., 1./27., 2./27. }
    };
    const double deriv2[5][6] =
    {
      { 2./25., -17./150.,   1.0/30.0 ,  0.0     ,  0.     , 0.     },
      { 0.    ,   1./ 42., -11.0/210.0,  1.0/35.0,  0.     , 0.     },
      { 0.    ,   0.     ,   1.0/15.0 ,-11.0/90.0,  1./18. , 0.     },
      { 0.    ,   0.     ,   0.0      ,  1.0/54.0, -7./162., 2./81. },
      { 0.    ,   0.     ,   0.0      ,  1.0/54.0, -7./162., 2./81. }
    };

    for (int i = 0; i < sizeof(xloc)/sizeof(xloc[0]); ++i) {
        for (int j = 0; j < n; ++j) {
            // Put garbage in Q
            for (int l = 0; l < k; ++l) Q[l] = 789e123;

            // Set only C[j] = 1.0 to obtain a single basis function
            for (int l = 0; l < n; ++l) C[l] = 0.0;
            C[j] = 1.0;

            // Check values and derivatives of basis j at point xloc[i]
            suzerain_bspline_series_evaluate(n, k, t, k-1, xloc[i], C, Q);
            BOOST_CHECK_CLOSE(Q[0], eval[i][j],   5.0e-12);
            BOOST_CHECK_CLOSE(Q[1], deriv[i][j],  5.0e-12);
            BOOST_CHECK_CLOSE(Q[2], deriv2[i][j], 5.0e-12);
        }
    }
}

BOOST_AUTO_TEST_CASE( partition_of_unity_for_piecewise_linears )
{
    const int    k    = 2;
    const double t[]  = {0, 0, 1, 2, 3, 3};
    const int    nt   = sizeof(t)/sizeof(t[0]);
    const int    n    = (nt - 2*(k-1) /* nbreakpoints */) + k - 2;
    const double C[n] = {1, 1, 1, 1};

    double Q[k];

    const int    nsamples = 25;
    const double incr_x   = (t[nt-1] - t[0])/nsamples;
    for (int i = 0; i <= nsamples; ++i) {
        const double x = t[0] + nsamples * incr_x;
        suzerain_bspline_series_evaluate(n, k, t, k-1, x, C, Q);

        BOOST_CHECK_EQUAL(1.0, Q[0]);
        BOOST_CHECK_EQUAL(0.0, Q[1]);
    }
}

BOOST_AUTO_TEST_CASE( partition_of_unity_for_piecewise_quadratics )
{
    const int    k    = 3;
    const double t[]  = {0, 0, 0, 1, 2, 3, 3, 3};
    const int    nt   = sizeof(t)/sizeof(t[0]);
    const int    n    = (nt - 2*(k-1) /* nbreakpoints */) + k - 2;
    const double C[n] = {1, 1, 1, 1, 1};

    double Q[k];

    const int    nsamples = 25;
    const double incr_x   = (t[nt-1] - t[0])/nsamples;
    for (int i = 0; i <= nsamples; ++i) {
        const double x = t[0] + nsamples * incr_x;
        suzerain_bspline_series_evaluate(n, k, t, k-1, x, C, Q);

        BOOST_CHECK_EQUAL(1.0, Q[0]);
        BOOST_CHECK_EQUAL(0.0, Q[1]);
        BOOST_CHECK_EQUAL(0.0, Q[2]);
    }
}

BOOST_AUTO_TEST_CASE( partition_of_unity_for_piecewise_cubics )
{
    const int    k    = 4;
    const double t[]  = {0, 0, 0, 0, 1, 2, 3, 3, 3, 3};
    const int    nt   = sizeof(t)/sizeof(t[0]);
    const int    n    = (nt - 2*(k-1) /* nbreakpoints */) + k - 2;
    const double C[n] = {1, 1, 1, 1, 1, 1};

    double Q[k];

    const int    nsamples = 25;
    const double incr_x   = (t[nt-1] - t[0])/nsamples;
    for (int i = 0; i <= nsamples; ++i) {
        const double x = t[0] + nsamples * incr_x;
        suzerain_bspline_series_evaluate(n, k, t, k-1, x, C, Q);

        BOOST_CHECK_EQUAL(1.0, Q[0]);
        BOOST_CHECK_EQUAL(0.0, Q[1]);
        BOOST_CHECK_EQUAL(0.0, Q[2]);
        BOOST_CHECK_EQUAL(0.0, Q[3]);
    }
}
