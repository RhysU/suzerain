#define BOOST_TEST_MODULE $Id$

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <boost/format.hpp>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <log4cxx/logger.h>

using namespace log4cxx;

LoggerPtr logger = Logger::getRootLogger();

double f(double x, void *p) {
    double retval = 0.0;

    if (0.0 < x && x <= 1.0) {
        retval = x;
    } else if (1.0 < x && x <= 2.0) {
        retval = 1.0 - x;
    }
    retval *= retval;

    return retval;
}

int main(int argc, char **argv)
{
    double abserr = 0.0;
    double relerr = 1.0;
    double left = 0.0;
    double right = 2.0;

    const gsl_function gsl_f = { f, NULL };

    if (argc >= 2) abserr = atof(argv[1]);
    if (argc >= 3) relerr = atof(argv[2]);
    if (argc >= 4) left   = atof(argv[2]);
    if (argc >= 5) right  = atof(argv[2]);

    LOG4CXX_DEBUG(logger, boost::format("abserr:           %g") % abserr);
    LOG4CXX_DEBUG(logger, boost::format("relerr:           %g") % relerr);
    LOG4CXX_DEBUG(logger, boost::format("left:             %g") % left);
    LOG4CXX_DEBUG(logger, boost::format("right:            %g") % right);

    {
        double result1, result2;
        double esterr1, esterr2;
        size_t neval1, neval2;

        const double center = (left+right)/2.0;
        gsl_integration_qng(&gsl_f, left, center, abserr, relerr,
                            &result1, &esterr1, &neval1);
        gsl_integration_qng(&gsl_f, center, right, abserr, relerr,
                            &result2, &esterr2, &neval2);

        LOG4CXX_DEBUG(logger, "Two intervals...");
        LOG4CXX_DEBUG(logger, boost::format("Result:           %g")  % (result1 + result2));
        LOG4CXX_DEBUG(logger, boost::format("Estimated error:  %g")  % (esterr1 + esterr2));
        LOG4CXX_DEBUG(logger, boost::format("# of evaluations: %zu") % (neval1 + neval2));
    }

    {
        double result;
        double esterr;
        size_t neval;

        gsl_integration_qng(&gsl_f, left, right, abserr, relerr,
                            &result, &esterr, &neval);

        LOG4CXX_DEBUG(logger, "One interval...");
        LOG4CXX_DEBUG(logger, boost::format("Result:           %g")  % result);
        LOG4CXX_DEBUG(logger, boost::format("Estimated error:  %g")  % esterr);
        LOG4CXX_DEBUG(logger, boost::format("# of evaluations: %zu") % neval);
    }

    return 0;
}
