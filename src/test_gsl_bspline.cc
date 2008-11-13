#include <stdio.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>

using namespace log4cxx;

int
main()
{
    BasicConfigurator::configure();
    LoggerPtr logger = Logger::getRootLogger();

    const size_t nbreak= 5, order = 2;
    const double left_end = 0.0, right_end = 1.0;
    boost::array<double,4> eval_points = {{ 1.0/8.0, 1.0/16.0, 1.0/32.0, 0 }};

    gsl_bspline_workspace *bw = gsl_bspline_alloc(order, nbreak);

    gsl_bspline_knots_uniform(left_end, right_end, bw);
    LOG4CXX_DEBUG(logger, "Knots:");
    gsl_vector_fprintf(stdout, bw->knots, "\t%0.5f");

    gsl_vector *v = gsl_vector_alloc(gsl_bspline_ncoeffs(bw));
    BOOST_FOREACH( double eval_point, eval_points ) {
        gsl_bspline_eval(eval_point, v, bw);
        LOG4CXX_DEBUG(logger, boost::format("Values at %0.5f:") % eval_point);
        gsl_vector_fprintf(stdout, v, "\t%0.5f");
    }

    gsl_vector_free(v);
    gsl_bspline_free(bw);

    return 0;
}
