#include "config.h"

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <boost/foreach.hpp>
#include <log4cxx/logger.h>

using namespace log4cxx;

int
main()
{
    LoggerPtr logger = Logger::getRootLogger();

    const size_t nbreak = 5, order = 2;
    const double left_end = 0.0, right_end = 1.0;
    const double eval_points[] = { 1.0 / 8.0, 1.0 / 16.0, 1.0 / 32.0, 0 };

    gsl_bspline_workspace *bw = gsl_bspline_alloc(order, nbreak);

    gsl_bspline_knots_uniform(left_end, right_end, bw);

    if (logger->isDebugEnabled()) {
        for (size_t i = 0; i < bw->knots->size; i++) {
            LOG4CXX_DEBUG(logger,
                          "knot[" << i << "] at " << gsl_vector_get(bw->knots, i)
                         )
        }
    }

    gsl_vector *v = gsl_vector_alloc(gsl_bspline_ncoeffs(bw));

    BOOST_FOREACH(double eval_point, eval_points) {
        gsl_bspline_eval(eval_point, v, bw);

        if (logger->isDebugEnabled()) {
            for (size_t i = 0; i < v->size; i++) {
                LOG4CXX_DEBUG(logger,
                              "At x = " << eval_point
                              << " basis[" << i << "] = "
                              << gsl_vector_get(v, i)
                             )
            }
        }
    }

    gsl_vector_free(v);
    gsl_bspline_free(bw);

    return 0;
}
