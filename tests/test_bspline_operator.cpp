#define BOOST_TEST_MODULE $Id$

#include <config.h>

#include <boost/test/included/unit_test.hpp>
#include <suzerain/bspline_operator.h>

BOOST_AUTO_TEST_CASE( allocation_okay )
{
    const int order  = 4;
    const int nbreak = 10;
    const int nderiv = 2;
    suzerain_bspline_operator_workspace * w
        = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);
    suzerain_bspline_operator_free(w);
}


BOOST_AUTO_TEST_CASE( memory_layout )
{
    const int order  = 3;
    const int nbreak = 5;
    const int nderiv = 2;

    suzerain_bspline_operator_workspace * w
        = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    suzerain_bspline_operator_create(NULL, w); // FIXME

    suzerain_bspline_operator_free(w);
}
