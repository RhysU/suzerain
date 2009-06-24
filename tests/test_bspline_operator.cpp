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


BOOST_AUTO_TEST_CASE( memory_layout )
{
    const double breakpoints[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    const int nbreak = sizeof(breakpoints)/sizeof(breakpoints[0]);
    const int order  = 3;
    const int nderiv = 2;

    LOG4CXX_DEBUG(logger, format("nbreak=%d, order=%d, nderiv=%d")
                          % nbreak % order % nderiv);

    suzerain_bspline_operator_workspace * w
        = suzerain_bspline_operator_alloc(order, nbreak, nderiv,
            SUZERAIN_BSPLINE_OPERATOR_COLLOCATION_GREVILLE);

    LOG4CXX_DEBUG(logger, format("w->nbreakpoints=%d") % w->nbreakpoints);
    LOG4CXX_DEBUG(logger, format("w->order       =%d") % w->order);
    LOG4CXX_DEBUG(logger, format("w->nderivatives=%d") % w->nderivatives);
    LOG4CXX_DEBUG(logger, format("w->n           =%d") % w->n);

    suzerain_bspline_operator_create(breakpoints, w);

    suzerain_bspline_operator_free(w);
}
