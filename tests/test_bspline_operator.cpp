#define BOOST_TEST_MODULE $Id$

#include <config.h>

#include <boost/format.hpp>
#include <boost/test/included/unit_test.hpp>
#include <log4cxx/logger.h>
#include <suzerain/bspline_operator.h>

using boost::format;
log4cxx::LoggerPtr logger = log4cxx::Logger::getRootLogger();

BOOST_AUTO_TEST_CASE( allocation_okay )
{
    const int order  = 4, nbreak = 10, nderiv = 2;

    suzerain_bspline_operator_workspace * w
        = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);
    suzerain_bspline_operator_free(w);
}


// Check a simple piecewise linear case's general banded storage
BOOST_AUTO_TEST_CASE( memory_layout )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 2;
    const int nderiv = 1;

    suzerain_bspline_operator_workspace * w
        = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    suzerain_bspline_operator_create(breakpoints, w);

    /* Check w->D[0], the mass matrix, against known good solution:
     *   1              0              0              0
     *   0              1              0              0
     *   0              0              1              0
     *   0              0              0              1
     * Known good is in general banded matrix column-major order.
     */
    const double good_D0[] = { /*DK*/   1,     0,
                                   0,   1,     0,
                                   0,   1,     0,
                                   0,   1  /*DK*/ };
    BOOST_CHECK_EQUAL_COLLECTIONS(
        good_D0, good_D0 + sizeof(good_D0)/sizeof(good_D0[0]),
        w->D[0] + 1, w->D[0] + w->storagesize - 1);

    /* Check w->D[1], the first derivative matrix, against known good:
     *       -1              1              0              0
     *        0             -1              1              0
     *        0              0             -1              1
     *        0              0             -1              1
     * Known good is in general banded matrix column-major order.
     */
    const double good_D1[] = { /*DK*/  -1,     0,
                                   1,  -1,     0,
                                   1,  -1,    -1,
                                   1,   1  /*DK*/ };
    BOOST_CHECK_EQUAL_COLLECTIONS(
        good_D1, good_D1 + sizeof(good_D1)/sizeof(good_D1[0]),
        w->D[1] + 1, w->D[1] + w->storagesize - 1);

    suzerain_bspline_operator_free(w);
}
